# Train a GAT attention graph from the MixFind celltype annotated object
rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

if (!requireNamespace("torch", quietly = TRUE)) {
  stop(
    "The R package 'torch' is required for GAT attention training. ",
    "Install/load torch on the server, then rerun this script."
  )
}
suppressPackageStartupMessages(library(torch))

# =========================
# 0. Parameters
# =========================
output_root <- "/media/zenglab/result/lingyuan/STEM/RNA-proteinRegulatoryModules/imputation/attention/magic_stagate_like_pipeline"
input_dir <- file.path(output_root, "02_mixfind_annotation")
out_dir <- file.path(output_root, "03_attention")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_rds <- file.path(input_dir, "02_mixfind_celltype_annotation.rds")
out_rds <- file.path(out_dir, "03_attention_graph.rds")
out_rdata <- file.path(out_dir, "03_attention_graph.RData")

wnn_dims <- 1:20
candidate_spatial_k <- 30
attention_k <- 20
magic_k <- 10
min_same_celltype_neighbors <- 5
celltype_prior_col <- "celltype_annotation"
wnn_graph_name <- "expression_wsnn"
spatial_sample_col <- "type"

gat_hidden_dim <- 32
gat_latent_dim <- 16
gat_epochs <- 400
gat_learning_rate <- 1e-3
gat_weight_decay <- 1e-4
gat_dropout <- 0.05
gat_device <- if (torch::cuda_is_available()) "cuda" else "cpu"
gat_seed <- 20260622

# Biological / structural priors are added as fixed edge bias before softmax.
# Celltype is intentionally dominant because downstream MAGIC should diffuse
# mainly inside annotated cell types.
prior_celltype_same_bonus <- 8
prior_celltype_cross_penalty <- 12
prior_wnn_similarity_weight <- 3
prior_spatial_distance_penalty <- 2
prior_same_expression_cluster_bonus <- 1
prior_same_wnn_mix_cluster_bonus <- 1.5
prior_same_z_bonus <- 0.4

# =========================
# 1. Helpers
# =========================
scale_matrix_no_na <- function(mat) {
  mat <- scale(mat)
  mat[is.na(mat)] <- 0
  mat
}

safe_same <- function(x, i, j) {
  same <- x[i] == x[j]
  same[is.na(same)] <- FALSE
  same
}

is_real_celltype <- function(x) {
  !is.na(x) & !x %in% c("Mix", "Unknown", "")
}

safe_mean_num <- function(x) {
  if (length(x) == 0) return(NA_real_)
  mean(x, na.rm = TRUE)
}

get_reduction_features <- function(seu, reductions = c("rbrna.harmony", "ntrna.harmony"), dims = 1:20) {
  missing_reductions <- setdiff(reductions, names(seu@reductions))
  if (length(missing_reductions) > 0) {
    stop("Missing reductions: ", paste(missing_reductions, collapse = ", "))
  }
  cell_names <- Cells(seu)
  feature_list <- lapply(reductions, function(reduction_name) {
    emb <- Embeddings(seu, reduction = reduction_name)
    use_dims <- dims[dims <= ncol(emb)]
    if (length(use_dims) == 0) stop("No usable dimensions in reduction: ", reduction_name)
    emb[cell_names, use_dims, drop = FALSE]
  })
  scale_matrix_no_na(do.call(cbind, feature_list))
}

get_sparse_graph_values <- function(seu, graph_name, from_idx, to_idx) {
  if (!graph_name %in% names(seu@graphs)) {
    warning("WNN graph '", graph_name, "' is missing. WNN similarity prior will be set to 0.")
    return(rep(0, length(from_idx)))
  }
  graph_mat <- seu@graphs[[graph_name]]
  values_ij <- graph_mat[cbind(from_idx, to_idx)]
  values_ji <- graph_mat[cbind(to_idx, from_idx)]
  pmax(as.numeric(values_ij), as.numeric(values_ji), na.rm = TRUE)
}

build_prior_edges <- function(
    seu,
    spatial_coords,
    celltype_prior_col,
    candidate_spatial_k,
    attention_k,
    min_same_celltype_neighbors,
    wnn_graph_name,
    spatial_sample_col = NULL) {
  cell_names <- Cells(seu)
  n_cells <- length(cell_names)
  spatial_coords <- spatial_coords[cell_names, , drop = FALSE]
  if (anyNA(spatial_coords)) {
    stop("Spatial coordinates contain NA after aligning to Seurat cells.")
  }
  if (!celltype_prior_col %in% colnames(seu@meta.data)) {
    stop(sprintf("Column '%s' is missing from seu@meta.data.", celltype_prior_col))
  }

  celltype_label <- as.character(seu@meta.data[[celltype_prior_col]])
  expression_cluster <- if ("expression_cluster" %in% colnames(seu@meta.data)) {
    as.character(seu@meta.data$expression_cluster)
  } else {
    rep(NA_character_, n_cells)
  }
  wnn_cluster_mix <- if ("wnn_cluster_mix" %in% colnames(seu@meta.data)) {
    as.character(seu@meta.data$wnn_cluster_mix)
  } else {
    rep(NA_character_, n_cells)
  }

  meta_prior_cols <- intersect(c("z"), colnames(seu@meta.data))
  meta_prior <- lapply(meta_prior_cols, function(col) as.character(seu@meta.data[[col]]))
  names(meta_prior) <- meta_prior_cols

  if (!is.null(spatial_sample_col) && spatial_sample_col %in% colnames(seu@meta.data)) {
    spatial_sample <- as.character(seu@meta.data[[spatial_sample_col]])
    spatial_sample[is.na(spatial_sample) | spatial_sample == ""] <- "UnknownSample"
  } else {
    if (!is.null(spatial_sample_col)) {
      warning("Column '", spatial_sample_col, "' is missing. Spatial KNN will be built on all cells together.")
    }
    spatial_sample <- rep("all_cells", n_cells)
  }

  edge_chunks <- vector("list", n_cells)
  for (sample_id in sort(unique(spatial_sample))) {
    sample_idx <- which(spatial_sample == sample_id)
    if (length(sample_idx) <= 1) next

    sample_k <- min(candidate_spatial_k + 1, length(sample_idx))
    spatial_nn <- FindNeighbors(
      object = spatial_coords[sample_idx, , drop = FALSE],
      k.param = sample_k,
      return.neighbor = TRUE,
      compute.SNN = FALSE,
      verbose = FALSE
    )
    nn_idx_raw <- spatial_nn@nn.idx
    nn_dist_raw <- spatial_nn@nn.dist

    for (local_i in seq_along(sample_idx)) {
      i <- sample_idx[local_i]
      keep <- nn_idx_raw[local_i, ] != local_i
      nbr_local_idx <- nn_idx_raw[local_i, keep]
      nbr_dist <- nn_dist_raw[local_i, keep]
      valid <- !is.na(nbr_local_idx) & nbr_local_idx > 0
      nbr_local_idx <- nbr_local_idx[valid]
      nbr_dist <- nbr_dist[valid]
      if (length(nbr_local_idx) == 0) next

      nbr_idx <- sample_idx[nbr_local_idx]
      same_raw_celltype <- safe_same(celltype_label, i, nbr_idx)
      real_query_celltype <- is_real_celltype(celltype_label[i])
      real_neighbor_celltype <- is_real_celltype(celltype_label[nbr_idx])
      same_celltype <- same_raw_celltype & real_query_celltype & real_neighbor_celltype
      same_available <- real_query_celltype && sum(same_celltype) >= min_same_celltype_neighbors
      if (same_available) {
        keep_nbr <- same_celltype
      } else {
        keep_nbr <- rep(TRUE, length(nbr_idx))
      }
      nbr_idx <- nbr_idx[keep_nbr]
      nbr_dist <- nbr_dist[keep_nbr]
      same_raw_celltype <- same_raw_celltype[keep_nbr]
      real_neighbor_celltype <- real_neighbor_celltype[keep_nbr]
      same_celltype <- same_celltype[keep_nbr]

      ord <- order(!same_celltype, nbr_dist)
      ord <- ord[seq_len(min(attention_k, length(ord)))]
      nbr_idx <- nbr_idx[ord]
      nbr_dist <- nbr_dist[ord]
      same_raw_celltype <- same_raw_celltype[ord]
      real_neighbor_celltype <- real_neighbor_celltype[ord]
      same_celltype <- same_celltype[ord]

      edge_chunks[[i]] <- data.frame(
        from_idx = i,
        to_idx = nbr_idx,
        Cell1 = cell_names[i],
        Cell2 = cell_names[nbr_idx],
        SpatialSample = sample_id,
        Celltype1 = celltype_label[i],
        Celltype2 = celltype_label[nbr_idx],
        SpatialDistance = nbr_dist,
        SameCelltype = same_celltype,
        SameRawCelltype = same_raw_celltype,
        CelltypePriorEligible = real_query_celltype & real_neighbor_celltype,
        SameExpressionCluster = safe_same(expression_cluster, i, nbr_idx),
        SameWNNMixCluster = safe_same(wnn_cluster_mix, i, nbr_idx),
        stringsAsFactors = FALSE
      )

      for (col in meta_prior_cols) {
        edge_chunks[[i]][[paste0("Same_", col)]] <- safe_same(meta_prior[[col]], i, nbr_idx)
      }
    }
  }

  edges <- bind_rows(edge_chunks)
  if (nrow(edges) == 0) {
    stop("No spatial prior edges were constructed.")
  }

  edges$WNNSimilarity <- get_sparse_graph_values(seu, wnn_graph_name, edges$from_idx, edges$to_idx)
  edges$SpatialDistanceScaled <- edges$SpatialDistance / max(edges$SpatialDistance, na.rm = TRUE)
  if (max(edges$WNNSimilarity, na.rm = TRUE) > 0) {
    edges$WNNSimilarityScaled <- edges$WNNSimilarity / max(edges$WNNSimilarity, na.rm = TRUE)
  } else {
    edges$WNNSimilarityScaled <- 0
  }

  cross_real_celltype <- edges$CelltypePriorEligible & !edges$SameCelltype

  edges$PriorScore <- 0
  edges$PriorScore <- edges$PriorScore +
    prior_celltype_same_bonus * as.numeric(edges$SameCelltype) -
    prior_celltype_cross_penalty * as.numeric(cross_real_celltype) +
    prior_wnn_similarity_weight * edges$WNNSimilarityScaled -
    prior_spatial_distance_penalty * edges$SpatialDistanceScaled +
    prior_same_expression_cluster_bonus * as.numeric(edges$SameExpressionCluster) +
    prior_same_wnn_mix_cluster_bonus * as.numeric(edges$SameWNNMixCluster)

  if ("Same_z" %in% colnames(edges)) edges$PriorScore <- edges$PriorScore + prior_same_z_bonus * as.numeric(edges$Same_z)

  edges
}

gat_autoencoder <- nn_module(
  "GATAutoencoder",
  initialize = function(input_dim, hidden_dim, latent_dim, dropout = 0.05) {
    self$feature_proj <- nn_linear(input_dim, hidden_dim, bias = FALSE)
    self$attn_vec <- nn_parameter(torch_randn(c(hidden_dim * 2, 1)) * 0.01)
    self$latent_proj <- nn_linear(hidden_dim, latent_dim)
    self$decoder <- nn_linear(latent_dim, input_dim)
    self$dropout <- dropout
  },
  forward = function(x, edge_src, edge_dst, edge_bias, n_nodes) {
    h <- self$feature_proj(x)
    h <- nnf_elu(h)
    h_src <- h$index_select(1, edge_src)
    h_dst <- h$index_select(1, edge_dst)
    edge_h <- torch_cat(list(h_src, h_dst), dim = 2)
    edge_score <- torch_matmul(edge_h, self$attn_vec)$squeeze(2)
    edge_score <- nnf_leaky_relu(edge_score, negative_slope = 0.2) + edge_bias
    edge_score <- torch_clamp(edge_score, min = -30, max = 30)

    exp_score <- torch_exp(edge_score)
    denom <- torch_zeros(c(n_nodes), device = x$device)
    denom$index_add_(1, edge_src, exp_score)
    alpha <- exp_score / (denom$index_select(1, edge_src) + 1e-12)
    alpha <- nnf_dropout(alpha, p = self$dropout, training = self$training)

    msg <- h_dst * alpha$unsqueeze(2)
    agg <- torch_zeros(c(n_nodes, h$size(2)), device = x$device)
    agg$index_add_(1, edge_src, msg)

    z <- self$latent_proj(agg)
    z <- nnf_elu(z)
    recon <- self$decoder(z)
    list(recon = recon, z = z, alpha = alpha)
  }
)

train_gat_attention <- function(feature_mat, edges, hidden_dim, latent_dim, epochs, lr, weight_decay, dropout, device, seed) {
  torch_manual_seed(seed)
  n_nodes <- nrow(feature_mat)
  input_dim <- ncol(feature_mat)

  x <- torch_tensor(feature_mat, dtype = torch_float(), device = device)
  edge_src <- torch_tensor(edges$from_idx, dtype = torch_long(), device = device)
  edge_dst <- torch_tensor(edges$to_idx, dtype = torch_long(), device = device)
  edge_bias <- torch_tensor(edges$PriorScore, dtype = torch_float(), device = device)

  model <- gat_autoencoder(input_dim, hidden_dim, latent_dim, dropout = dropout)$to(device = device)
  optimizer <- optim_adam(model$parameters, lr = lr, weight_decay = weight_decay)

  loss_history <- data.frame(epoch = integer(), reconstruction_loss = numeric())
  for (epoch in seq_len(epochs)) {
    model$train()
    optimizer$zero_grad()
    out <- model(x, edge_src, edge_dst, edge_bias, n_nodes)
    loss <- nnf_mse_loss(out$recon, x)
    loss$backward()
    optimizer$step()

    loss_value <- as.numeric(loss$item())
    loss_history <- rbind(loss_history, data.frame(epoch = epoch, reconstruction_loss = loss_value))
    if (epoch == 1 || epoch %% 10 == 0 || epoch == epochs) {
      message(sprintf("GAT epoch %03d/%03d | reconstruction MSE = %.6f", epoch, epochs, loss_value))
    }
  }

  model$eval()
  with_no_grad({
    out <- model(x, edge_src, edge_dst, edge_bias, n_nodes)
    alpha <- as.numeric(as_array(out$alpha$to(device = "cpu")))
    latent <- as_array(out$z$to(device = "cpu"))
  })
  rownames(latent) <- rownames(feature_mat)
  colnames(latent) <- paste0("GAT_", seq_len(ncol(latent)))

  list(alpha = alpha, latent = latent, loss_history = loss_history)
}

build_magic_transition_from_attention <- function(attention_edges, cell_names, k_magic = 10) {
  cell_id <- setNames(seq_along(cell_names), cell_names)
  attention_edges <- attention_edges %>%
    filter(Cell1 %in% cell_names, Cell2 %in% cell_names, Cell1 != Cell2, Weight > 0) %>%
    bind_rows(transmute(., Cell1 = Cell2, Cell2 = Cell1, Weight = Weight)) %>%
    group_by(Cell1, Cell2) %>%
    summarise(Weight = max(Weight), .groups = "drop") %>%
    group_by(Cell1) %>%
    slice_max(Weight, n = k_magic, with_ties = FALSE) %>%
    ungroup()

  split_edges <- split(attention_edges, attention_edges$Cell1)
  nn_idx <- matrix(NA_integer_, nrow = length(cell_names), ncol = k_magic, dimnames = list(cell_names, NULL))
  nn_weight <- matrix(0, nrow = length(cell_names), ncol = k_magic, dimnames = list(cell_names, NULL))
  for (cell in names(split_edges)) {
    row_i <- cell_id[[cell]]
    cell_edges <- split_edges[[cell]]
    ord <- order(cell_edges$Weight, decreasing = TRUE)
    cell_edges <- cell_edges[ord[seq_len(min(k_magic, nrow(cell_edges)))], ]
    nn_idx[row_i, seq_len(nrow(cell_edges))] <- unname(cell_id[cell_edges$Cell2])
    nn_weight[row_i, seq_len(nrow(cell_edges))] <- cell_edges$Weight
  }

  A <- sparseMatrix(
    i = unname(cell_id[attention_edges$Cell1]),
    j = unname(cell_id[attention_edges$Cell2]),
    x = attention_edges$Weight,
    dims = c(length(cell_names), length(cell_names)),
    dimnames = list(cell_names, cell_names)
  )
  A_sym <- A + t(A)
  if (length(A_sym@x) == 0) {
    stop("The GAT attention graph did not return any positive edges for MAGIC.")
  }
  diag(A_sym) <- max(A_sym@x)
  row_sums <- Matrix::rowSums(A_sym)
  row_sums[row_sums == 0] <- 1
  M <- Diagonal(x = 1 / row_sums) %*% A_sym

  list(
    nn = list(nn.idx = nn_idx, nn.weight = nn_weight, nn.dist = 1 / (nn_weight + 1e-8)),
    affinity = A_sym,
    transition = M
  )
}

# =========================
# 2. Train GAT attention graph
# =========================
states <- readRDS(input_rds)
if (is.null(states@misc$spatial_coords)) {
  stop("states@misc$spatial_coords is missing. Rerun 00_Harmony_preprocess.R and 01_WNN_clustering.R.")
}

feature_mat <- get_reduction_features(
  states,
  reductions = c("rbrna.harmony", "ntrna.harmony"),
  dims = wnn_dims
)

prior_edges <- build_prior_edges(
  seu = states,
  spatial_coords = states@misc$spatial_coords,
  celltype_prior_col = celltype_prior_col,
  candidate_spatial_k = candidate_spatial_k,
  attention_k = attention_k,
  min_same_celltype_neighbors = min_same_celltype_neighbors,
  wnn_graph_name = wnn_graph_name,
  spatial_sample_col = spatial_sample_col
)

write.csv(prior_edges, file.path(out_dir, "03_gat_attention_prior_edges.csv"), row.names = FALSE)

gat_result <- train_gat_attention(
  feature_mat = feature_mat,
  edges = prior_edges,
  hidden_dim = gat_hidden_dim,
  latent_dim = gat_latent_dim,
  epochs = gat_epochs,
  lr = gat_learning_rate,
  weight_decay = gat_weight_decay,
  dropout = gat_dropout,
  device = gat_device,
  seed = gat_seed
)

attention_edges <- prior_edges
attention_edges$Weight <- gat_result$alpha
attention_edges <- attention_edges %>%
  group_by(Cell1) %>%
  mutate(Weight = Weight / sum(Weight)) %>%
  ungroup()
attention_edges$CelltypeEdgeClass <- ifelse(
  !attention_edges$CelltypePriorEligible,
  "Mix/Unknown involved",
  ifelse(attention_edges$SameCelltype, "Same real celltype", "Cross real celltype")
)

magic_graph <- build_magic_transition_from_attention(
  attention_edges = attention_edges,
  cell_names = Cells(states),
  k_magic = magic_k
)

states[["gat.attention"]] <- CreateDimReducObject(
  embeddings = gat_result$latent,
  key = "GAT_",
  assay = DefaultAssay(states)
)
states@misc$gat_attention <- list(
  edges = attention_edges,
  prior_edges = prior_edges,
  loss_history = gat_result$loss_history,
  params = list(
    candidate_spatial_k = candidate_spatial_k,
    attention_k = attention_k,
    magic_k = magic_k,
    gat_hidden_dim = gat_hidden_dim,
    gat_latent_dim = gat_latent_dim,
    gat_epochs = gat_epochs,
    gat_learning_rate = gat_learning_rate,
    gat_weight_decay = gat_weight_decay,
    gat_dropout = gat_dropout,
    gat_device = gat_device,
    spatial_sample_col = spatial_sample_col,
    prior_celltype_same_bonus = prior_celltype_same_bonus,
    prior_celltype_cross_penalty = prior_celltype_cross_penalty,
    prior_wnn_similarity_weight = prior_wnn_similarity_weight,
    prior_spatial_distance_penalty = prior_spatial_distance_penalty
  )
)
states@misc$magic_attention_graph <- magic_graph
states@misc$magic_transition <- magic_graph$transition

write.csv(attention_edges, file.path(out_dir, "03_gat_attention_edges.csv"), row.names = FALSE)
write.csv(gat_result$loss_history, file.path(out_dir, "03_gat_training_loss.csv"), row.names = FALSE)

p_loss <- ggplot(gat_result$loss_history, aes(x = epoch, y = reconstruction_loss)) +
  geom_line(color = "#2c3e50", linewidth = 1) +
  theme_minimal() +
  labs(title = "GAT autoencoder training loss", x = "Epoch", y = "Reconstruction MSE")
ggsave(file.path(out_dir, "03_gat_training_loss.png"), p_loss, width = 7, height = 5, dpi = 300)

p_attention <- ggplot(attention_edges, aes(x = Weight, fill = CelltypeEdgeClass)) +
  geom_histogram(bins = 100, color = "white", position = "identity", alpha = 0.7) +
  theme_minimal() +
  labs(title = "GAT attention weight distribution", x = "Attention weight", y = "Edges", fill = "Celltype edge class")
ggsave(file.path(out_dir, "03_attention_weight_distribution.png"), p_attention, width = 7, height = 5, dpi = 300)

edge_summary <- attention_edges %>%
  summarise(
    n_edges = n(),
    celltype_prior_eligible_fraction = mean(CelltypePriorEligible),
    same_real_celltype_fraction = safe_mean_num(SameCelltype[CelltypePriorEligible]),
    mix_or_unknown_edge_fraction = mean(!CelltypePriorEligible),
    mean_weight_same_real_celltype = safe_mean_num(Weight[CelltypePriorEligible & SameCelltype]),
    mean_weight_cross_real_celltype = safe_mean_num(Weight[CelltypePriorEligible & !SameCelltype]),
    mean_weight_mix_or_unknown_involved = safe_mean_num(Weight[!CelltypePriorEligible]),
    median_spatial_distance = median(SpatialDistance),
    mean_wnn_similarity = mean(WNNSimilarity)
  )
write.csv(edge_summary, file.path(out_dir, "03_attention_edge_summary.csv"), row.names = FALSE)

save(states, attention_edges, prior_edges, gat_result, magic_graph, file = out_rdata)
saveRDS(states, file = out_rds)
message("Done. Saved GAT attention graph object to: ", out_rds)
