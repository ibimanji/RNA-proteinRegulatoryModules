# Train a GAT attention graph from the WNN object using h5ad metadata celltypes
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
output_root <- "/media/zenglab/result/lingyuan/STEM/imputation/magic_stagate_like_pipeline"
input_dir <- file.path(output_root, "01_wnn_clustering")
out_dir <- file.path(output_root, "03_attention")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
attention_plot_dir <- file.path(out_dir, "gat_umap_spatial_qc")
dir.create(attention_plot_dir, recursive = TRUE, showWarnings = FALSE)

input_rds <- file.path(input_dir, "01_wnn_clustering.rds")
out_rds <- file.path(out_dir, "03_attention_graph.rds")
out_rdata <- file.path(out_dir, "03_attention_graph.RData")

wnn_dims <- 1:20
candidate_spatial_k <- 30
attention_k <- 20
magic_k <- 10
min_same_celltype_neighbors <- 5
celltype_prior_col <- "rna_nn_alg1_label3"
wnn_graph_name <- "expression_wsnn"
batch_var <- "protocol-replicate"
spatial_sample_col <- make.names(batch_var)
gat_cluster_col <- "gat_attention_cluster"
gat_neighbor_graph_name <- "gat_attention_snn"
gat_cluster_resolution <- 0.5

celltype_colors <- c(
  "AC1" = "#ccba33",      
  "AC2" = "#ffbe85",   
  "AC3" = "#e3782b",  
  "CHOR" = "#7f52a9",     
  "EPEN" = "#c4b0d4",    
  "CHO/PEP" = "#97f4f7", 
  "INH_Sst" = "#96abeb",      
  "INH_Pvalb" = "#96cad4", 
  "INH_Cnr1_Vip" = "#a8e1eb",
  "MLG" = "#8597c6",     
  "OPC" = "#667872", 
  "OLG1" = "#e4f768",      
  "OLG2" = "#e6db17",   
  "VLMC" = "#1f76b3" ,      
  "VSMC" = "#00aeef" ,
  "Peri/VEC" = "#d3a59c",
  "DE/MEN" = "#b274e8", 
  "MSN" = "#D96DA1" ,
  "DGGRC" = "#a6e8a6",  
  "TEGLU CA1" = "#77ed8f",      
  "TEGLU CA2" = "#82ad2d",      
  "TEGLU CA3" = "#28330b",      
  "TEGLU L2/3" = "#cbfc60",
  "TEGLU L4" = "#96db00",
  "TEGLU L5" = "#04b361",
  "TEGLU L5/6" = "#40d102", 
  "TEGLU L6" = "#32a630",
  "TEGLU L6b" = "#406e27",
  "TEGLU Mix" = "#c5fcc5",      
  "Mix"="#F5F5F5"
)

# Estimate celltype-specific spatial constraints from the annotation itself.
# A larger same-celltype / all-cell kth-neighbor distance ratio means the
# celltype is spatially scattered, so it receives lower distance penalty and
# a wider candidate neighborhood. Compact celltypes receive stronger spatial
# restriction. Mix/Unknown are not treated as biological cell types.
adaptive_celltype_spatial <- TRUE
adaptive_aggregation_k <- 10
adaptive_min_cells_per_sample <- 50
adaptive_min_samples_per_celltype <- 1
adaptive_penalty_range <- c(1.2, 3.8)
adaptive_candidate_spatial_k_range <- c(20, 60)
adaptive_attention_k_range <- c(10, 30)
adaptive_exclude_celltypes <- c("Mix", "Unknown", "")

# Optional manual overrides. Leave empty to use fully adaptive values.
manual_celltype_candidate_spatial_k <- c()
manual_celltype_attention_k <- c()
manual_celltype_spatial_distance_penalty <- c()

gat_hidden_dim <- 16
gat_latent_dim <- 16
gat_epochs <- 400
gat_learning_rate <- 1e-3
gat_weight_decay <- 1e-4
gat_dropout <- 0.05
# CPU is slower but much safer for this large graph. CUDA can be enabled
# manually if the server has enough GPU memory.
gat_device <- "cpu"
gat_seed <- 20260622

# Biological / structural priors are added as fixed edge bias before softmax.
# Celltype is intentionally dominant because downstream MAGIC should diffuse
# mainly inside annotated cell types.
use_celltype_prior <- FALSE
prior_celltype_same_bonus <- 8
prior_celltype_cross_penalty <- 12
prior_wnn_similarity_weight <- 3
prior_spatial_distance_penalty <- 2
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

plot_df_heatmap <- function(df, x_col, y_col, fill_col, title, fill_label) {
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[fill_col]])) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient(low = "white", high = "#2c7fb8", na.value = "grey92") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    labs(title = title, x = NULL, y = NULL, fill = fill_label)
}

get_celltype_colors <- function(seu, group_col) {
  if (group_col != celltype_prior_col || !group_col %in% colnames(seu@meta.data)) return(NULL)
  group_levels <- unique(as.character(seu@meta.data[[group_col]]))
  missing_colors <- setdiff(group_levels, names(celltype_colors))
  color_values <- celltype_colors[intersect(names(celltype_colors), group_levels)]
  if (length(missing_colors) > 0) {
    extra_colors <- setNames(scales::hue_pal()(length(missing_colors)), missing_colors)
    color_values <- c(color_values, extra_colors)
  }
  color_values
}

make_spatial_meta <- function(seu, coord_cols = c("column", "row", "z")) {
  meta <- seu@meta.data
  meta$Cell <- rownames(meta)
  if (is.null(seu@misc$spatial_coords)) {
    stop("states@misc$spatial_coords is missing.")
  }
  coords <- as.data.frame(seu@misc$spatial_coords[rownames(meta), , drop = FALSE])
  coords$Cell <- rownames(coords)
  duplicated_coord_cols <- intersect(setdiff(colnames(coords), "Cell"), colnames(meta))
  if (length(duplicated_coord_cols) > 0) {
    meta <- meta[, setdiff(colnames(meta), duplicated_coord_cols), drop = FALSE]
  }
  spatial_meta <- left_join(meta, coords, by = "Cell")
  missing_coord_cols <- setdiff(c("column", "row"), colnames(spatial_meta))
  if (length(missing_coord_cols) > 0) {
    stop("Spatial metadata is missing required columns: ", paste(missing_coord_cols, collapse = ", "))
  }
  spatial_meta
}

save_gat_umap_plot <- function(seu, group_col, file, title, color_values = NULL, label = FALSE) {
  if (!group_col %in% colnames(seu@meta.data)) {
    warning("Skip GAT UMAP plot because column '", group_col, "' is missing.")
    return(invisible(NULL))
  }
  p <- DimPlot(
    seu,
    reduction = "gat.attention.umap",
    group.by = group_col,
    cols = color_values,
    label = label,
    repel = TRUE
  ) +
    ggtitle(title)
  ggsave(file, p, width = 8, height = 7, dpi = 300)
  invisible(p)
}

save_spatial_plot_by_sample <- function(meta_df, group_col, sample_col, file, title, point_size = 0.18, color_values = NULL) {
  if (!group_col %in% colnames(meta_df)) {
    warning("Skip spatial plot because column '", group_col, "' is missing.")
    return(invisible(NULL))
  }
  plot_df <- meta_df
  if (!is.null(sample_col) && sample_col %in% colnames(plot_df)) {
    plot_df$sample_value <- plot_df[[sample_col]]
  }
  p <- ggplot(plot_df, aes(x = column, y = row, color = .data[[group_col]])) +
    geom_point(size = point_size, alpha = 0.9) +
    coord_fixed() +
    scale_y_reverse() +
    theme_void() +
    labs(title = title, color = group_col)
  if (!is.null(color_values)) {
    p <- p + scale_color_manual(values = color_values, na.value = "grey80")
  }
  if (!is.null(sample_col) && sample_col %in% colnames(plot_df)) {
    p <- p + facet_wrap(~sample_value)
  }
  ggsave(file, p, width = 12, height = 8, dpi = 300)
  invisible(p)
}

weighted_median_num <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[which(cumsum(w) >= sum(w) / 2)[1]]
}

get_celltype_param <- function(celltype, param_map, default_value) {
  value <- unname(param_map[celltype])
  if (length(value) == 0 || is.na(value)) default_value else value
}

scale_to_range <- function(x, out_range, reverse = FALSE) {
  if (length(x) == 0) return(numeric())
  if (all(is.na(x))) return(rep(mean(out_range), length(x)))
  q <- stats::quantile(x, probs = c(0.05, 0.95), na.rm = TRUE, names = FALSE)
  if (isTRUE(all.equal(q[1], q[2]))) {
    scaled <- rep(0.5, length(x))
  } else {
    scaled <- (pmin(pmax(x, q[1]), q[2]) - q[1]) / (q[2] - q[1])
  }
  if (reverse) scaled <- 1 - scaled
  out_range[1] + scaled * diff(out_range)
}

compute_adaptive_celltype_spatial_params <- function(
    seu,
    spatial_coords,
    celltype_col,
    sample_col,
    aggregation_k = 10,
    min_cells_per_sample = 50,
    min_samples_per_celltype = 1,
    exclude_celltypes = c("Mix", "Unknown", ""),
    penalty_range = c(1.2, 3.8),
    candidate_k_range = c(20, 60),
    attention_k_range = c(10, 30),
    default_candidate_k = 30,
    default_attention_k = 20,
    default_penalty = 2,
    manual_candidate_k = c(),
    manual_attention_k = c(),
    manual_penalty = c()) {
  cell_names <- Cells(seu)
  spatial_coords <- spatial_coords[cell_names, , drop = FALSE]
  celltype <- as.character(seu@meta.data[[celltype_col]])
  sample <- if (!is.null(sample_col) && sample_col %in% colnames(seu@meta.data)) {
    as.character(seu@meta.data[[sample_col]])
  } else {
    rep("all_cells", length(cell_names))
  }
  sample[is.na(sample) | sample == ""] <- "UnknownSample"

  sample_chunks <- list()
  chunk_i <- 1
  for (sample_id in sort(unique(sample))) {
    sample_idx <- which(sample == sample_id)
    if (length(sample_idx) <= aggregation_k + 1) next

    all_k <- min(aggregation_k + 1, length(sample_idx))
    all_nn <- FindNeighbors(
      object = spatial_coords[sample_idx, , drop = FALSE],
      k.param = all_k,
      return.neighbor = TRUE,
      compute.SNN = FALSE,
      verbose = FALSE
    )
    all_ref_distance <- stats::median(all_nn@nn.dist[, all_k], na.rm = TRUE)
    if (!is.finite(all_ref_distance) || all_ref_distance <= 0) next

    for (ct in sort(unique(celltype[sample_idx]))) {
      if (is.na(ct) || ct %in% exclude_celltypes) next
      ct_idx <- sample_idx[celltype[sample_idx] == ct]
      if (length(ct_idx) < min_cells_per_sample || length(ct_idx) <= aggregation_k + 1) next

      same_k <- min(aggregation_k + 1, length(ct_idx))
      same_nn <- FindNeighbors(
        object = spatial_coords[ct_idx, , drop = FALSE],
        k.param = same_k,
        return.neighbor = TRUE,
        compute.SNN = FALSE,
        verbose = FALSE
      )
      same_ref_distance <- stats::median(same_nn@nn.dist[, same_k], na.rm = TRUE)
      if (!is.finite(same_ref_distance)) next

      sample_chunks[[chunk_i]] <- data.frame(
        Celltype = ct,
        SpatialSample = sample_id,
        NCells = length(ct_idx),
        SameCelltypeKnnDistance = same_ref_distance,
        AllCellKnnDistance = all_ref_distance,
        AggregationRatio = same_ref_distance / all_ref_distance,
        stringsAsFactors = FALSE
      )
      chunk_i <- chunk_i + 1
    }
  }

  sample_df <- bind_rows(sample_chunks)
  if (nrow(sample_df) == 0) {
    warning("No adaptive celltype spatial statistics were computed. Using global defaults.")
    return(list(
      summary = data.frame(),
      sample_stats = sample_df,
      candidate_spatial_k = c(),
      attention_k = c(),
      spatial_distance_penalty = c()
    ))
  }

  summary_df <- sample_df %>%
    mutate(WeightNCells = NCells) %>%
    group_by(Celltype) %>%
    summarise(
      NSamples = n(),
      TotalCells = sum(WeightNCells),
      MedianAggregationRatio = weighted_median_num(AggregationRatio, WeightNCells),
      MeanAggregationRatio = stats::weighted.mean(AggregationRatio, w = WeightNCells, na.rm = TRUE),
      MedianSameCelltypeKnnDistance = weighted_median_num(SameCelltypeKnnDistance, WeightNCells),
      MedianAllCellKnnDistance = weighted_median_num(AllCellKnnDistance, WeightNCells),
      .groups = "drop"
    ) %>%
    filter(NSamples >= min_samples_per_celltype) %>%
    mutate(NCells = TotalCells) %>%
    select(-TotalCells)

  if (nrow(summary_df) == 0) {
    warning("No celltypes passed adaptive spatial filters. Using global defaults.")
    return(list(
      summary = data.frame(),
      sample_stats = sample_df,
      candidate_spatial_k = c(),
      attention_k = c(),
      spatial_distance_penalty = c()
    ))
  }

  log_ratio <- log1p(summary_df$MedianAggregationRatio)
  scattered_score <- scale_to_range(log_ratio, c(0, 1), reverse = FALSE)
  summary_df$SpatialScatteredScore <- scattered_score
  summary_df$AdaptiveSpatialDistancePenalty <- scale_to_range(scattered_score, penalty_range, reverse = TRUE)
  summary_df$AdaptiveCandidateSpatialK <- as.integer(round(scale_to_range(scattered_score, candidate_k_range, reverse = FALSE)))
  summary_df$AdaptiveAttentionK <- as.integer(round(scale_to_range(scattered_score, attention_k_range, reverse = FALSE)))

  candidate_map <- setNames(summary_df$AdaptiveCandidateSpatialK, summary_df$Celltype)
  attention_map <- setNames(summary_df$AdaptiveAttentionK, summary_df$Celltype)
  penalty_map <- setNames(summary_df$AdaptiveSpatialDistancePenalty, summary_df$Celltype)

  candidate_map[names(manual_candidate_k)] <- manual_candidate_k
  attention_map[names(manual_attention_k)] <- manual_attention_k
  penalty_map[names(manual_penalty)] <- manual_penalty

  summary_df$FinalCandidateSpatialK <- as.integer(unname(candidate_map[summary_df$Celltype]))
  summary_df$FinalAttentionK <- as.integer(unname(attention_map[summary_df$Celltype]))
  summary_df$FinalSpatialDistancePenalty <- unname(penalty_map[summary_df$Celltype])

  list(
    summary = summary_df,
    sample_stats = sample_df,
    candidate_spatial_k = candidate_map,
    attention_k = attention_map,
    spatial_distance_penalty = penalty_map
  )
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
    spatial_sample_col = NULL,
    celltype_candidate_spatial_k = NULL,
    celltype_attention_k = NULL,
    celltype_spatial_distance_penalty = NULL,
    use_celltype_prior = TRUE,
    prior_spatial_distance_penalty = 2) {
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
  max_candidate_spatial_k <- max(candidate_spatial_k, unname(celltype_candidate_spatial_k), na.rm = TRUE)
  for (sample_id in sort(unique(spatial_sample))) {
    sample_idx <- which(spatial_sample == sample_id)
    if (length(sample_idx) <= 1) next

    sample_k <- min(max_candidate_spatial_k + 1, length(sample_idx))
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

      query_celltype <- celltype_label[i]
      candidate_spatial_k_i <- get_celltype_param(query_celltype, celltype_candidate_spatial_k, candidate_spatial_k)
      attention_k_i <- get_celltype_param(query_celltype, celltype_attention_k, attention_k)
      spatial_distance_penalty_i <- get_celltype_param(
        query_celltype,
        celltype_spatial_distance_penalty,
        prior_spatial_distance_penalty
      )

      candidate_keep <- seq_len(min(candidate_spatial_k_i, length(nbr_local_idx)))
      nbr_local_idx <- nbr_local_idx[candidate_keep]
      nbr_dist <- nbr_dist[candidate_keep]

      nbr_idx <- sample_idx[nbr_local_idx]
      same_raw_celltype <- safe_same(celltype_label, i, nbr_idx)
      real_query_celltype <- is_real_celltype(celltype_label[i])
      real_neighbor_celltype <- is_real_celltype(celltype_label[nbr_idx])
      same_celltype <- same_raw_celltype & real_query_celltype & real_neighbor_celltype
      same_available <- use_celltype_prior && real_query_celltype && sum(same_celltype) >= min_same_celltype_neighbors
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

      if (use_celltype_prior) {
        ord <- order(!same_celltype, nbr_dist)
      } else {
        ord <- order(nbr_dist)
      }
      ord <- ord[seq_len(min(attention_k_i, length(ord)))]
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
        SpatialCandidateK = candidate_spatial_k_i,
        AttentionK = attention_k_i,
        SpatialDistancePenalty = spatial_distance_penalty_i,
        SameCelltype = same_celltype,
        SameRawCelltype = same_raw_celltype,
        CelltypePriorEligible = use_celltype_prior & real_query_celltype & real_neighbor_celltype,
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

  cross_real_celltype <- use_celltype_prior & edges$CelltypePriorEligible & !edges$SameCelltype
  celltype_same_bonus_use <- if (use_celltype_prior) prior_celltype_same_bonus else 0
  celltype_cross_penalty_use <- if (use_celltype_prior) prior_celltype_cross_penalty else 0

  edges$PriorScore <- 0
  edges$PriorScore <- edges$PriorScore +
    celltype_same_bonus_use * as.numeric(edges$SameCelltype) -
    celltype_cross_penalty_use * as.numeric(cross_real_celltype) +
    prior_wnn_similarity_weight * edges$WNNSimilarityScaled -
    edges$SpatialDistancePenalty * edges$SpatialDistanceScaled

  if ("Same_z" %in% colnames(edges)) edges$PriorScore <- edges$PriorScore + prior_same_z_bonus * as.numeric(edges$Same_z)

  edges
}

gat_autoencoder <- nn_module(
  "GATAutoencoder",
  initialize = function(input_dim, hidden_dim, latent_dim, dropout = 0.05) {
    self$feature_proj <- nn_linear(input_dim, hidden_dim, bias = FALSE)
    self$attn_src <- nn_parameter(torch_randn(c(hidden_dim, 1)) * 0.01)
    self$attn_dst <- nn_parameter(torch_randn(c(hidden_dim, 1)) * 0.01)
    self$latent_proj <- nn_linear(hidden_dim, latent_dim)
    self$decoder <- nn_linear(latent_dim, input_dim)
    self$dropout <- dropout
  },
  forward = function(x, edge_src, edge_dst, edge_bias, n_nodes) {
    h <- self$feature_proj(x)
    h <- nnf_elu(h)
    src_score <- torch_matmul(h, self$attn_src)$squeeze(2)
    dst_score <- torch_matmul(h, self$attn_dst)$squeeze(2)
    edge_score <- src_score$index_select(1, edge_src) + dst_score$index_select(1, edge_dst)
    edge_score <- nnf_leaky_relu(edge_score, negative_slope = 0.2) + edge_bias
    edge_score <- torch_clamp(edge_score, min = -30, max = 30)

    exp_score <- torch_exp(edge_score)
    denom <- torch_zeros(c(n_nodes), device = x$device)
    denom$index_add_(1, edge_src, exp_score)
    alpha <- exp_score / (denom$index_select(1, edge_src) + 1e-12)
    alpha <- nnf_dropout(alpha, p = self$dropout, training = self$training)

    h_dst <- h$index_select(1, edge_dst)
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
if (!celltype_prior_col %in% colnames(states@meta.data)) {
  stop(sprintf("Celltype metadata column '%s' is missing from states@meta.data.", celltype_prior_col))
}
if (!spatial_sample_col %in% colnames(states@meta.data) &&
    !is.null(states@misc$harmony_batch_var) &&
    states@misc$harmony_batch_var %in% colnames(states@meta.data)) {
  spatial_sample_col <- states@misc$harmony_batch_var
}
if (!spatial_sample_col %in% colnames(states@meta.data)) {
  warning("Spatial sample column '", spatial_sample_col, "' is missing. Spatial KNN and QC plots will use all cells together.")
}

feature_mat <- get_reduction_features(
  states,
  reductions = c("rbrna.harmony", "ntrna.harmony"),
  dims = wnn_dims
)

if (adaptive_celltype_spatial && use_celltype_prior) {
  adaptive_spatial_params <- compute_adaptive_celltype_spatial_params(
    seu = states,
    spatial_coords = states@misc$spatial_coords,
    celltype_col = celltype_prior_col,
    sample_col = spatial_sample_col,
    aggregation_k = adaptive_aggregation_k,
    min_cells_per_sample = adaptive_min_cells_per_sample,
    min_samples_per_celltype = adaptive_min_samples_per_celltype,
    exclude_celltypes = adaptive_exclude_celltypes,
    penalty_range = adaptive_penalty_range,
    candidate_k_range = adaptive_candidate_spatial_k_range,
    attention_k_range = adaptive_attention_k_range,
    default_candidate_k = candidate_spatial_k,
    default_attention_k = attention_k,
    default_penalty = prior_spatial_distance_penalty,
    manual_candidate_k = manual_celltype_candidate_spatial_k,
    manual_attention_k = manual_celltype_attention_k,
    manual_penalty = manual_celltype_spatial_distance_penalty
  )
  celltype_candidate_spatial_k <- adaptive_spatial_params$candidate_spatial_k
  celltype_attention_k <- adaptive_spatial_params$attention_k
  celltype_spatial_distance_penalty <- adaptive_spatial_params$spatial_distance_penalty

  write.csv(
    adaptive_spatial_params$sample_stats,
    file.path(out_dir, "03_adaptive_celltype_spatial_stats_by_sample.csv"),
    row.names = FALSE
  )
  write.csv(
    adaptive_spatial_params$summary,
    file.path(out_dir, "03_adaptive_celltype_spatial_params.csv"),
    row.names = FALSE
  )
  if (nrow(adaptive_spatial_params$summary) > 0) {
    p_adaptive <- ggplot(
      adaptive_spatial_params$summary,
      aes(
        x = MedianAggregationRatio,
        y = FinalSpatialDistancePenalty,
        size = NCells,
        label = Celltype
      )
    ) +
      geom_point(color = "#2c7fb8", alpha = 0.75) +
      geom_text(size = 3, vjust = -0.6, check_overlap = FALSE) +
      scale_x_log10() +
      theme_minimal() +
      labs(
        title = "Adaptive spatial penalty by celltype aggregation",
        x = "Same-celltype / all-cell kth-neighbor distance ratio",
        y = "Final spatial distance penalty",
        size = "Cells"
      )
    ggsave(
      file.path(out_dir, "03_adaptive_celltype_spatial_params.png"),
      p_adaptive,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
} else {
  if (adaptive_celltype_spatial && !use_celltype_prior) {
    message("Skip adaptive celltype spatial parameters because use_celltype_prior = FALSE.")
  }
  adaptive_spatial_params <- NULL
  celltype_candidate_spatial_k <- manual_celltype_candidate_spatial_k
  celltype_attention_k <- manual_celltype_attention_k
  celltype_spatial_distance_penalty <- manual_celltype_spatial_distance_penalty
}

prior_edges <- build_prior_edges(
  seu = states,
  spatial_coords = states@misc$spatial_coords,
  celltype_prior_col = celltype_prior_col,
  candidate_spatial_k = candidate_spatial_k,
  attention_k = attention_k,
  min_same_celltype_neighbors = min_same_celltype_neighbors,
  wnn_graph_name = wnn_graph_name,
  spatial_sample_col = spatial_sample_col,
  celltype_candidate_spatial_k = celltype_candidate_spatial_k,
  celltype_attention_k = celltype_attention_k,
  celltype_spatial_distance_penalty = celltype_spatial_distance_penalty,
  use_celltype_prior = use_celltype_prior,
  prior_spatial_distance_penalty = prior_spatial_distance_penalty
)

message(sprintf(
  "Training GAT with %s cells, %s edges, mean %.2f edges/cell, hidden_dim = %s, device = %s",
  format(nrow(feature_mat), big.mark = ","),
  format(nrow(prior_edges), big.mark = ","),
  nrow(prior_edges) / nrow(feature_mat),
  gat_hidden_dim,
  gat_device
))

saveRDS(prior_edges, file.path(out_dir, "03_gat_attention_prior_edges.rds"))

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
if (use_celltype_prior) {
  attention_edges$CelltypeEdgeClass <- ifelse(
    !attention_edges$CelltypePriorEligible,
    "Mix/Unknown involved",
    ifelse(attention_edges$SameCelltype, "Same real celltype", "Cross real celltype")
  )
} else {
  attention_edges$CelltypeEdgeClass <- "Celltype prior disabled"
}

cell_meta <- data.frame(
  Cell = Cells(states),
  Celltype = as.character(states@meta.data[[celltype_prior_col]]),
  SpatialSample = if (spatial_sample_col %in% colnames(states@meta.data)) {
    as.character(states@meta.data[[spatial_sample_col]])
  } else {
    "all_cells"
  },
  stringsAsFactors = FALSE
)
cell_meta$Celltype[is.na(cell_meta$Celltype) | cell_meta$Celltype == ""] <- "Unknown"
cell_meta$SpatialSample[is.na(cell_meta$SpatialSample) | cell_meta$SpatialSample == ""] <- "UnknownSample"

effective_neighbors_by_cell <- attention_edges %>%
  group_by(Cell1) %>%
  summarise(
    AttentionEdges = n(),
    EffectiveNeighbors = 1 / sum(Weight^2),
    MaxWeight = max(Weight),
    .groups = "drop"
  )

cell_effective_meta <- cell_meta %>%
  left_join(effective_neighbors_by_cell, by = c("Cell" = "Cell1"))
cell_effective_meta$AttentionEdges[is.na(cell_effective_meta$AttentionEdges)] <- 0
cell_effective_meta$EffectiveNeighbors[is.na(cell_effective_meta$EffectiveNeighbors)] <- 0
cell_effective_meta$MaxWeight[is.na(cell_effective_meta$MaxWeight)] <- 0

celltype_effective_summary <- cell_effective_meta %>%
  group_by(Celltype) %>%
  summarise(
    NCells = n(),
    MeanAttentionEdges = mean(AttentionEdges),
    MeanEffectiveNeighbors = mean(EffectiveNeighbors),
    MedianEffectiveNeighbors = median(EffectiveNeighbors),
    SumEffectiveNeighbors = sum(EffectiveNeighbors),
    EffectiveNeighborFractionOfCelltype = MeanEffectiveNeighbors / NCells,
    MeanMaxWeight = mean(MaxWeight),
    .groups = "drop"
  ) %>%
  arrange(desc(NCells))
write.csv(
  celltype_effective_summary,
  file.path(out_dir, "03_qc_effective_neighbors_and_cell_counts_by_celltype.csv"),
  row.names = FALSE
)

celltype_sample_effective_summary <- cell_effective_meta %>%
  group_by(SpatialSample, Celltype) %>%
  summarise(
    NCells = n(),
    MeanAttentionEdges = mean(AttentionEdges),
    MeanEffectiveNeighbors = mean(EffectiveNeighbors),
    MedianEffectiveNeighbors = median(EffectiveNeighbors),
    SumEffectiveNeighbors = sum(EffectiveNeighbors),
    EffectiveNeighborFractionOfCelltype = MeanEffectiveNeighbors / NCells,
    MeanMaxWeight = mean(MaxWeight),
    .groups = "drop"
  ) %>%
  arrange(SpatialSample, desc(NCells))
write.csv(
  celltype_sample_effective_summary,
  file.path(out_dir, "03_qc_effective_neighbors_and_cell_counts_by_celltype_sample.csv"),
  row.names = FALSE
)

p_celltype_counts <- ggplot(
  celltype_effective_summary,
  aes(x = reorder(Celltype, NCells), y = NCells)
) +
  geom_col(fill = "#4c78a8") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Cell counts by celltype", x = NULL, y = "Cells")
ggsave(
  file.path(out_dir, "03_qc_cell_counts_by_celltype.png"),
  p_celltype_counts,
  width = 7,
  height = 5,
  dpi = 300
)

p_effective_ratio <- ggplot(
  celltype_effective_summary,
  aes(x = reorder(Celltype, EffectiveNeighborFractionOfCelltype), y = EffectiveNeighborFractionOfCelltype)
) +
  geom_col(fill = "#59a14f") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Effective neighbor fraction by celltype",
    x = NULL,
    y = "Mean effective neighbors / celltype cell count"
  )
ggsave(
  file.path(out_dir, "03_qc_effective_neighbor_fraction_by_celltype.png"),
  p_effective_ratio,
  width = 7,
  height = 5,
  dpi = 300
)

p_effective_ratio_sample <- ggplot(
  celltype_sample_effective_summary,
  aes(x = reorder(Celltype, EffectiveNeighborFractionOfCelltype), y = EffectiveNeighborFractionOfCelltype)
) +
  geom_col(fill = "#59a14f") +
  coord_flip() +
  facet_wrap(~SpatialSample, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Effective neighbor fraction by celltype and sample",
    x = NULL,
    y = "Mean effective neighbors / celltype cell count"
  )
ggsave(
  file.path(out_dir, "03_qc_effective_neighbor_fraction_by_celltype_sample.png"),
  p_effective_ratio_sample,
  width = 11,
  height = 7,
  dpi = 300
)

real_celltype_levels <- sort(unique(cell_meta$Celltype[is_real_celltype(cell_meta$Celltype)]))
celltype_mean_weight <- attention_edges %>%
  filter(CelltypePriorEligible) %>%
  group_by(Celltype1, Celltype2) %>%
  summarise(
    MeanWeight = mean(Weight),
    Edges = n(),
    .groups = "drop"
  )
celltype_heatmap_grid <- expand.grid(
  Celltype1 = real_celltype_levels,
  Celltype2 = real_celltype_levels,
  stringsAsFactors = FALSE
) %>%
  left_join(celltype_mean_weight, by = c("Celltype1", "Celltype2"))
celltype_heatmap_grid$MeanWeightFilled <- celltype_heatmap_grid$MeanWeight
celltype_heatmap_grid$MeanWeightFilled[is.na(celltype_heatmap_grid$MeanWeightFilled)] <- 0
celltype_heatmap_grid$Edges[is.na(celltype_heatmap_grid$Edges)] <- 0
write.csv(
  celltype_heatmap_grid,
  file.path(out_dir, "03_qc_celltype_annotation_mean_weight_heatmap_table.csv"),
  row.names = FALSE
)
p_celltype_heatmap <- plot_df_heatmap(
  celltype_heatmap_grid,
  x_col = "Celltype2",
  y_col = "Celltype1",
  fill_col = "MeanWeightFilled",
  title = "Mean attention weight by celltype pair",
  fill_label = "Mean weight"
)
ggsave(
  file.path(out_dir, "03_qc_celltype_annotation_mean_weight_heatmap_filled_zero.png"),
  p_celltype_heatmap,
  width = 8,
  height = 7,
  dpi = 300
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

states <- RunUMAP(
  states,
  reduction = "gat.attention",
  dims = seq_len(ncol(Embeddings(states, "gat.attention"))),
  reduction.name = "gat.attention.umap",
  reduction.key = "GATUMAP_",
  verbose = FALSE
)

#######attention embedding umap
states <- FindNeighbors(
  states,
  reduction = "gat.attention",
  dims = seq_len(ncol(Embeddings(states, "gat.attention"))),
  graph.name = gat_neighbor_graph_name,
  verbose = FALSE
)

states <- FindClusters(
  states,
  graph.name = gat_neighbor_graph_name,
  resolution = gat_cluster_resolution,
  verbose = FALSE
)

states[[gat_cluster_col]] <- as.character(Idents(states))

if (celltype_prior_col %in% colnames(states@meta.data)) {
  gat_cluster_celltype_table <- as.data.frame.matrix(
    table(states[[gat_cluster_col]][, 1], states[[celltype_prior_col]][, 1])
  )
  gat_cluster_celltype_table <- cbind(
    GATCluster = rownames(gat_cluster_celltype_table),
    gat_cluster_celltype_table
  )
  rownames(gat_cluster_celltype_table) <- NULL
  write.csv(
    gat_cluster_celltype_table,
    file.path(attention_plot_dir, "03_gat_cluster_by_celltype_annotation_table.csv"),
    row.names = FALSE
  )
}

plot_group_cols <- intersect(c(gat_cluster_col, celltype_prior_col), colnames(states@meta.data))
spatial_meta <- make_spatial_meta(states)

for (group_col in plot_group_cols) {
  group_label <- gsub("[^A-Za-z0-9_]+", "_", group_col)
  color_values <- get_celltype_colors(states, group_col)
  label_clusters <- group_col == gat_cluster_col
  
  save_gat_umap_plot(
    states,
    group_col = group_col,
    file = file.path(attention_plot_dir, paste0("03_gat_umap_by_", group_label, ".png")),
    title = paste("GAT attention UMAP by", group_col),
    color_values = color_values,
    label = label_clusters
  )
  
  save_spatial_plot_by_sample(
    spatial_meta,
    group_col = group_col,
    sample_col = spatial_sample_col,
    file = file.path(attention_plot_dir, paste0("03_spatial_by_", group_label, "_split_by_", spatial_sample_col, ".png")),
    title = paste("Spatial distribution by", group_col),
    color_values = color_values
  )
}

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
    use_celltype_prior = use_celltype_prior,
    spatial_sample_col = spatial_sample_col,
    adaptive_celltype_spatial = adaptive_celltype_spatial,
    adaptive_aggregation_k = adaptive_aggregation_k,
    adaptive_min_cells_per_sample = adaptive_min_cells_per_sample,
    adaptive_min_samples_per_celltype = adaptive_min_samples_per_celltype,
    adaptive_penalty_range = adaptive_penalty_range,
    adaptive_candidate_spatial_k_range = adaptive_candidate_spatial_k_range,
    adaptive_attention_k_range = adaptive_attention_k_range,
    adaptive_exclude_celltypes = adaptive_exclude_celltypes,
    celltype_candidate_spatial_k = celltype_candidate_spatial_k,
    celltype_attention_k = celltype_attention_k,
    celltype_spatial_distance_penalty = celltype_spatial_distance_penalty,
    manual_celltype_candidate_spatial_k = manual_celltype_candidate_spatial_k,
    manual_celltype_attention_k = manual_celltype_attention_k,
    manual_celltype_spatial_distance_penalty = manual_celltype_spatial_distance_penalty,
    prior_celltype_same_bonus = prior_celltype_same_bonus,
    prior_celltype_cross_penalty = prior_celltype_cross_penalty,
    prior_wnn_similarity_weight = prior_wnn_similarity_weight,
    prior_spatial_distance_penalty = prior_spatial_distance_penalty
  )
)
states@misc$magic_attention_graph <- magic_graph
states@misc$magic_transition <- magic_graph$transition

saveRDS(attention_edges, file.path(out_dir, "03_gat_attention_edges.rds"))
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

save(states, attention_edges, prior_edges, gat_result, magic_graph, adaptive_spatial_params, file = out_rdata)
saveRDS(states, file = out_rds)
message("Done. Saved GAT attention graph object to: ", out_rds)
