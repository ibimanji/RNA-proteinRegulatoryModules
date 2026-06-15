# MAGIC imputation with a STAGATE-like attention graph
getwd()
setwd("/media/zenglab/result/lingyuan/STEM/decon/SpaImpute")
rm(list = ls()); gc(); options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(matrixStats)
  library(anndata)
  library(ggplot2)
  library(ggrepel)
})

# =========================
# 0. Parameters
# =========================
h5ad_path <- "/media/zenglab/result/lingyuan/STEM/decon/celltype_annotation/preprocessh5ad/14mWT_20260527.h5ad"
spatial_coord_cols <- c("column", "row", "z")
wnn_dims <- 1:20
wnn_cluster_resolution <- 0.5
spatial_k <- 10
magic_k <- 10
expression_weight <- 1
spatial_weight <- 0.25
cluster_bonus <- 0.25
t_diff <- 4
t_max <- 10
output_te_imputed <- "/media/zenglab/result/lingyuan/STEM/decon/SpaImpute/14mWT_TE_imputed.csv"
output_te_rescaled <- "/media/zenglab/result/lingyuan/STEM/decon/SpaImpute/14mWT_TE_rescaled.csv"

# =========================
# 1. Data loading
# =========================
states <- read_h5ad(h5ad_path)
if (!all(c("rbRNA", "ntRNA") %in% names(states$layers))) {
  stop("The h5ad object must contain layers named 'rbRNA' and 'ntRNA'.")
}
if (!all(spatial_coord_cols %in% colnames(states$obs))) {
  stop(sprintf(
    "The h5ad obs must contain spatial coordinate columns: %s.",
    paste(spatial_coord_cols, collapse = ", ")
  ))
}

rb_raw_mat <- t(states$layers[["rbRNA"]])
nt_raw_mat <- t(states$layers[["ntRNA"]])
total_raw_mat <- t(states$layers[["ntRNA"]]) + t(states$layers[["rbRNA"]]) 
meta <- states$obs
spatial_mat <- as.matrix(meta[, spatial_coord_cols, drop = FALSE])
storage.mode(spatial_mat) <- "numeric"
rownames(spatial_mat) <- rownames(meta)
states_seu <- CreateSeuratObject(
  counts = total_raw_mat,
  meta.data = meta,
  project = "STEM_14mWT"
)
states_seu[["rbRNA"]] <- CreateAssayObject(
  counts = rb_raw_mat
)
states_seu[["ntRNA"]] <- CreateAssayObject(
  counts = nt_raw_mat
)
spatial_mat <- spatial_mat[Cells(states_seu), , drop = FALSE]
if (anyNA(spatial_mat)) {
  stop("Spatial coordinates could not be aligned to the filtered Seurat cells.")
}
if ("TE" %in% names(states$layers)) {
  states_seu[["TE_raw"]] <- CreateAssayObject(data = as(t(states$layers[["TE"]]), "dgCMatrix"))
}

#magic trial via R
#1.1 
#normalize should be carefully considered，here we shall use totalRNA layer to calculate spike-in-factor, to evaluate the cell sequencing-lib
#With this sequencing-lib, we can use them to normalize rbRNA layer (and ntRNA layer)
#!!!also seurat-innate normalize function still can be used.Normalization should be carefully considered since it is related to (1)PC calculation(2)D matrix in the data imputing process.

#totalRNA layer as spike in factor
rna_counts <- GetAssayData(states_seu, assay = "RNA", layer = "counts")
rna_lib_sizes <- colSums(rna_counts) #as spike-in factor for other layers

#totalRNA layer: which can be done on default seurat function pipeline
DefaultAssay(states_seu) <- "RNA"
states_seu <- NormalizeData(states_seu, normalization.method = "RC", scale.factor = median(rna_lib_sizes))
states_seu <- FindVariableFeatures(states_seu, selection.method = "vst", nfeatures = 1000)
states_seu <- ScaleData(states_seu)
states_seu <- RunPCA(states_seu, npcs = 50, reduction.name = "rna.pca", reduction.key = "rnaPC_")

#rbRNA layer
rb_counts <- GetAssayData(states_seu, assay = "rbRNA", layer = "counts")
scale_factor_ref_rb <- median(colSums(GetAssayData(states_seu, assay = "rbRNA", layer = "counts")))
rb_norm <- sweep(rb_counts, 2, rna_lib_sizes, "/") * scale_factor_ref_rb
states_seu <- SetAssayData(states_seu, assay = "rbRNA", layer = "data", new.data = as(rb_norm, "dgCMatrix"))
states_seu <- FindVariableFeatures(states_seu, selection.method = "vst", nfeatures = 1000, assay = "rbRNA")
states_seu <- ScaleData(states_seu, assay = "rbRNA")
states_seu <- RunPCA(states_seu, npcs = 50, assay = "rbRNA", reduction.name = "rbrna.pca", reduction.key = "rbrnaPC_")

#ntRNA layer
nt_counts <- GetAssayData(states_seu, assay = "ntRNA", layer = "counts")
scale_factor_ref_nt <- median(colSums(GetAssayData(states_seu, assay = "ntRNA", layer = "counts")))
nt_norm <- sweep(nt_counts, 2, rna_lib_sizes, "/") * scale_factor_ref_nt
states_seu <- SetAssayData(states_seu, assay = "ntRNA", layer = "data", new.data = as(nt_norm, "dgCMatrix"))
states_seu <- FindVariableFeatures(states_seu, selection.method = "vst", nfeatures = 1000, assay = "ntRNA")
states_seu <- ScaleData(states_seu, assay = "ntRNA")
states_seu <- RunPCA(states_seu, npcs = 50, assay = "ntRNA", reduction.name = "ntrna.pca", reduction.key = "ntrnaPC_")

#clustering based on expression: rbRNA + ntRNA WNN only for joint expression embedding/annotation
states_seu <- FindMultiModalNeighbors(
  states_seu,
  reduction.list = list("rbrna.pca", "ntrna.pca"), # here we shall mention that why ntRNA? totalRNA, is rb plus nt, if use total and rb as input, may lead to amplify the weight from rbRNA layer.
  dims.list = list(wnn_dims, wnn_dims),
  weighted.nn.name = "expression.weighted.nn",
  knn.graph.name = "expression_wknn",
  snn.graph.name = "expression_wsnn"
)
states_seu <- RunUMAP(
  states_seu,
  nn.name = "expression.weighted.nn",
  reduction.name = "expression.wnn.umap",
  reduction.key = "exprWNNUMAP_"
)
states_seu <- FindClusters(
  states_seu,
  graph.name = "expression_wsnn",
  resolution = wnn_cluster_resolution,
  algorithm = 3
)
######## 先给了一个占位式注释 ########
states_seu$expression_cluster <- as.character(Idents(states_seu))
states_seu$expression_celltype_annotation <- paste0("expression_cluster_", states_seu$expression_cluster)
if (!"celltype_annotation" %in% colnames(states_seu@meta.data)) {
  states_seu$celltype_annotation <- states_seu$expression_celltype_annotation
}

build_stagate_like_attention_edges <- function(
    seu,
    spatial_coords,
    reduction.list = c("rbrna.pca", "ntrna.pca"),
    dims.list = list(1:20, 1:20),
    spatial_k = 10,
    expression_weight = 1,
    spatial_weight = 0.25,
    cluster_bonus = 0.25) {
  cell_names <- Cells(seu)
  spatial_coords <- spatial_coords[cell_names, , drop = FALSE]
  if (anyNA(spatial_coords)) {
    stop("Spatial coordinates contain NA after aligning to Seurat cells.")
  }

  emb_list <- Map(function(reduction_name, dims_use) {
    emb <- Embeddings(seu, reduction = reduction_name)
    emb[cell_names, dims_use, drop = FALSE]
  }, reduction.list, dims.list)
  expression_emb <- scale(do.call(cbind, emb_list))
  expression_emb[is.na(expression_emb)] <- 0

  spatial_nn <- FindNeighbors(
    object = spatial_coords,
    k.param = spatial_k + 1,
    return.neighbor = TRUE,
    compute.SNN = FALSE
  )
  nn_idx_raw <- spatial_nn@nn.idx
  nn_dist_raw <- spatial_nn@nn.dist
  n_cells <- length(cell_names)
  cluster_label <- seu$expression_cluster
  names(cluster_label) <- cell_names

  edge_list <- vector("list", n_cells)
  nn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = spatial_k, dimnames = list(cell_names, NULL))
  nn_weight <- matrix(0, nrow = n_cells, ncol = spatial_k, dimnames = list(cell_names, NULL))
  nn_dist <- matrix(Inf, nrow = n_cells, ncol = spatial_k, dimnames = list(cell_names, NULL))

  for (i in seq_len(n_cells)) {
    keep <- nn_idx_raw[i, ] != i
    nbr_idx <- nn_idx_raw[i, keep]
    nbr_spatial_dist <- nn_dist_raw[i, keep]
    if (length(nbr_idx) > spatial_k) {
      nbr_idx <- nbr_idx[seq_len(spatial_k)]
      nbr_spatial_dist <- nbr_spatial_dist[seq_len(spatial_k)]
    }
    if (length(nbr_idx) == 0) {
      next
    }

    expr_delta <- sweep(expression_emb[nbr_idx, , drop = FALSE], 2, expression_emb[i, ], "-")
    expr_dist <- sqrt(rowSums(expr_delta^2))
    expr_sigma <- max(expr_dist)
    spatial_sigma <- max(nbr_spatial_dist)
    if (!is.finite(expr_sigma) || expr_sigma == 0) expr_sigma <- 1e-6
    if (!is.finite(spatial_sigma) || spatial_sigma == 0) spatial_sigma <- 1e-6

    same_cluster <- cluster_label[nbr_idx] == cluster_label[i]
    same_cluster[is.na(same_cluster)] <- FALSE
    attention_score <- -expression_weight * (expr_dist / expr_sigma)^2 -
      spatial_weight * (nbr_spatial_dist / spatial_sigma)^2 +
      cluster_bonus * as.numeric(same_cluster)
    attention_score <- attention_score - max(attention_score)
    attention_weight <- exp(attention_score)
    attention_weight <- attention_weight / sum(attention_weight)

    k_i <- length(nbr_idx)
    nn_idx[i, seq_len(k_i)] <- nbr_idx
    nn_weight[i, seq_len(k_i)] <- attention_weight
    nn_dist[i, seq_len(k_i)] <- nbr_spatial_dist
    edge_list[[i]] <- data.frame(
      Cell1 = cell_names[i],
      Cell2 = cell_names[nbr_idx],
      Weight = attention_weight,
      SpatialDistance = nbr_spatial_dist,
      ExpressionDistance = expr_dist
    )
  }

  list(
    edges = bind_rows(edge_list),
    nn = list(nn.idx = nn_idx, nn.weight = nn_weight, nn.dist = nn_dist)
  )
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
  nn_idx <- matrix(NA_integer_, nrow = length(cell_names), ncol = k_magic)
  nn_weight <- matrix(0, nrow = length(cell_names), ncol = k_magic)
  rownames(nn_idx) <- rownames(nn_weight) <- cell_names
  for (cell in names(split_edges)) {
    row_i <- cell_id[[cell]]
    cell_edges <- split_edges[[cell]]
    ord <- order(cell_edges$Weight, decreasing = TRUE)
    cell_edges <- cell_edges[ord[seq_len(min(k_magic, nrow(cell_edges)))], ]
    nn_idx[row_i, seq_len(nrow(cell_edges))] <- unname(cell_id[cell_edges$Cell2])
    nn_weight[row_i, seq_len(nrow(cell_edges))] <- cell_edges$Weight
  }

  i_idx <- unname(cell_id[attention_edges$Cell1])
  j_idx <- unname(cell_id[attention_edges$Cell2])
  A <- sparseMatrix(
    i = i_idx,
    j = j_idx,
    x = attention_edges$Weight,
    dims = c(length(cell_names), length(cell_names)),
    dimnames = list(cell_names, cell_names)
  )
  A_sym <- A + t(A)
  if (length(A_sym@x) == 0) {
    stop("The STAGATE-like attention graph did not return any positive edges for MAGIC.")
  }
  diag(A_sym) <- max(A_sym@x)
  row_sums <- rowSums(A_sym)
  row_sums[row_sums == 0] <- 1
  M <- Diagonal(x = 1 / row_sums) %*% A_sym
  list(
    nn = list(nn.idx = nn_idx, nn.weight = nn_weight, nn.dist = 1 / (nn_weight + 1e-8)),
    affinity = A_sym,
    transition = M
  )
}

stagate_like_graph <- build_stagate_like_attention_edges(
  seu = states_seu,
  spatial_coords = spatial_mat,
  reduction.list = c("rbrna.pca", "ntrna.pca"),
  dims.list = list(wnn_dims, wnn_dims),
  spatial_k = spatial_k,
  expression_weight = expression_weight,
  spatial_weight = spatial_weight,
  cluster_bonus = cluster_bonus
)
attention_magic_graph <- build_magic_transition_from_attention(
  attention_edges = stagate_like_graph$edges,
  cell_names = Cells(states_seu),
  k_magic = magic_k
)
nn_obj <- attention_magic_graph$nn
states_seu@misc$stagate_like.weighted.nn <- stagate_like_graph$nn
states_seu@misc$magic_attention.weighted.nn <- nn_obj
A_sym <- attention_magic_graph$affinity
M <- attention_magic_graph$transition #probability matrix, also transition matrix for single step
M_t <- M
if (t_diff > 1) {
  for (i in 2:t_diff) {
    M_t <- M_t %*% M
  }
}#K-C equation

##test as coefficient ：
delta_R2_list <- numeric(t_max)

D_prev <- t(GetAssayData(states_seu, assay = "RNA", layer = "data")) #cell x gene! MENTION HERE IS CEEL X GENE!,norm layer
dim(D_prev)
print("开始计算图扩散与 R^2 收敛率 ...")
for (t in 1:t_max) {
  D_curr <- M %*% D_prev
  
  # 2. 基因（列）归一化：为了防止个别高表达基因主导变化率
  # 使用 Matrix::Diagonal 进行稀疏矩阵友好的快速列除法
  col_sums_curr <- colSums(D_curr) #cell sum 
  col_sums_curr[col_sums_curr == 0] <- 1
  D_curr_norm <- D_curr %*% Diagonal(x = 1/col_sums_curr)
  
  col_sums_prev <- colSums(D_prev)
  col_sums_prev[col_sums_prev == 0] <- 1
  D_prev_norm <- D_prev %*% Diagonal(x = 1/col_sums_prev)
  
  # 3. 计算 SSE (误差平方和)
  SSE <- sum((D_curr_norm - D_prev_norm)^2)
  
  # 4. 计算 SST (总平方和)
  mean_curr <- mean(D_curr_norm)
  SST <- sum((D_curr_norm - mean_curr)^2)
  
  # 5. 计算变化率: 1 - R^2
  delta_r2 <- SSE / SST
  delta_R2_list[t] <- delta_r2
  
  cat(sprintf("t = %2d | 1 - R^2 = %.5f\n", t, delta_r2))
  
  # 更新 D_prev 以备下一次循环
  D_prev <- D_curr
}
df_plot <- data.frame(
  t = 1:t_max,
  Delta_R2 = delta_R2_list
)

# 按照 MAGIC 论文标准：寻找跌破 0.05 阈值后的“第二个” t 作为最佳 t,当然其实这个值在我们这里是4、5、6左右
below_threshold <- which(delta_R2_list < 0.05)
optimal_t <- if(length(below_threshold) >= 2) below_threshold[2] else NA
p <- ggplot(df_plot, aes(x = t, y = Delta_R2)) +
  geom_line(color = "#2c3e50", linewidth = 1) +
  geom_point(color = "#e74c3c", size = 3) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray50", linewidth = 1) +
  scale_x_continuous(breaks = 1:t_max) +
  labs(
    title = "MAGIC Data Diffusion Convergence",
    subtitle = "Finding the optimal t (Threshold = 0.05)",
    x = "Diffusion Time (t)",
    y = expression(1 - R^2 ~~ paste("(", Data[t], ", ", Data[t-1], ")"))
  )+
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )
if (!is.na(optimal_t)) {
  p <- p + geom_vline(xintercept = optimal_t, linetype = "dotted", color = "blue", linewidth = 1) +
    annotate("text", x = optimal_t + 0.2, y = max(delta_R2_list)/2, 
             label = paste("Optimal t =", optimal_t), color = "blue", hjust = 0, size = 5)
  cat(sprintf("\n=> recommanded t step (Optimal t): %d\n", optimal_t))
} else {
  cat("\n=> 提示：在当前 t_max 范围内未连续两次跌破 0.05 阈值，请考虑增大 t_max。\n")
}

print(p)
#the downstream should be sealed into a function
D <- t(GetAssayData(states_seu, assay = "RNA", layer = "data")) #totalRNA
D_imputed <- M_t %*% D #this function lead to the loss of rownames of D_imputed, we should get it back
rownames(D_imputed) = rownames(D)
p99_original <- apply(D, 2, function(x) quantile(x, probs = 0.99))#this is for the re-scaling process like that in the MAGIC 2018.cell
max_imputed <- apply(D_imputed, 2, max)
max_imputed[max_imputed == 0] <- 1#这里没啥用没一个东西赋到这个值
scale_factors <- p99_original / max_imputed
D_rescaled <- sweep(D_imputed, 2, scale_factors, FUN = "*")
rescaled_totalRNA = D_rescaled
imputed_totalRNA = D_imputed

D <- t(GetAssayData(states_seu, assay = "rbRNA", layer = "data")) #rbRNA
D_imputed <- M_t %*% D #this function lead to the loss of rownames of D_imputed...
rownames(D_imputed) = rownames(D)
p99_original <- apply(D, 2, function(x) quantile(x, probs = 0.99))
max_imputed <- apply(D_imputed, 2, max)
max_imputed[max_imputed == 0] <- 1#这里没啥用没一个东西赋到这个值
scale_factors <- p99_original / max_imputed
D_rescaled <- sweep(D_imputed, 2, scale_factors, FUN = "*")
rescaled_rbRNA = D_rescaled
imputed_rbRNA = D_imputed

all.equal(rownames(imputed_totalRNA), rownames(imputed_rbRNA))
all.equal(colnames(imputed_totalRNA), colnames(imputed_rbRNA))

TE_imputed = imputed_rbRNA / (imputed_totalRNA+1e-06)
TE_rescaled = rescaled_rbRNA / (rescaled_totalRNA+1e-06)

meta_imputed <- states_seu@meta.data[rownames(imputed_totalRNA), , drop = FALSE]
seurat_obj <- CreateSeuratObject(counts = t(imputed_totalRNA), meta.data = meta_imputed)
convert_to_seurat_matrix <- function(mat) {
  if (inherits(mat, "dgeMatrix")) {
    return(as(mat, "dgCMatrix"))
  }
  if (is.matrix(mat)) {
    return(as(mat, "dgCMatrix"))
  }
  return(mat)
}
seurat_obj[["rescaled_totalRNA"]] <- CreateAssayObject(data = convert_to_seurat_matrix(t(rescaled_totalRNA)))
seurat_obj[["imputed_rbRNA"]]      <- CreateAssayObject(data = convert_to_seurat_matrix(t(imputed_rbRNA)))
seurat_obj[["rescaled_rbRNA"]]     <- CreateAssayObject(data = convert_to_seurat_matrix(t(rescaled_rbRNA)))
seurat_obj[["TE_imputed"]]         <- CreateAssayObject(data = convert_to_seurat_matrix(t(TE_imputed)))
seurat_obj[["TE_rescaled"]]        <- CreateAssayObject(data = convert_to_seurat_matrix(t(TE_rescaled)))
print(seurat_obj)

#further scanning: 
max(TE_imputed)
max(TE_rescaled)
non_zero_count <- sum(TE_imputed > 0) 
total_elements <- length(TE_imputed)
non_zero_count / total_elements

non_zero_ratio <- (non_zero_count / total_elements) * 100
cat(sprintf("TE 矩阵中非 0 元素的个数: %d\n", non_zero_count))
cat(sprintf("TE 矩阵中非 0 元素的比例: %.2f%%\n", non_zero_ratio))


#Cnp - Cnp correlation cauculation validation
target_gene <- "Cnp"
if (target_gene %in% colnames(TE_imputed)) {
  ube2c_te <- TE_rescaled[, target_gene]
  ube2c_rna_imputed <- rescaled_totalRNA[, target_gene]
  cor_pearson_imp <- cor(ube2c_te, ube2c_rna_imputed, method = "pearson")
  cat(sprintf("TE 与 插值后 RNA 的 Pearson 相关性: %.4f\n", cor_pearson_imp))
  cat("=========================================================\n")
  plot_df <- data.frame(
    RNA_Expression = ube2c_rna_imputed,
    Translation_Efficiency = ube2c_te
  )
  ggplot(plot_df, aes(x = RNA_Expression, y = Translation_Efficiency)) +
    geom_point(alpha = 0.6, color = "#3498db") +
    geom_smooth(method = "lm", color = "#e74c3c", se = TRUE) +
    labs(
      title = paste("Correlation Analysis for", target_gene),
      x = "Imputed Total RNA Expression",
      y = "Imputed Translation Efficiency (TE)"
    ) +
    theme_minimal()
} else {
  cat(sprintf("Skip marker QC: gene %s was not found in TE_imputed.\n", target_gene))
}

te_vector <- as.vector(TE_rescaled)
cat("============== TE_imputed 全局数值统计 ==============\n")
cat(sprintf("最小值 (Min)   : %f\n", min(te_vector)))
cat(sprintf("最大值 (Max)   : %f\n", max(te_vector)))
cat(sprintf("均值   (Mean)  : %f\n", mean(te_vector)))
cat(sprintf("中位数 (Median): %f\n", median(te_vector)))

cat("\n--- 更详细的分位数分布 ---\n")
print(quantile(te_vector, probs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)))
cat("====================================================\n")

# =====================================================================
# 2. 检查极端值（如由于除以接近0的数导致的极大值）
# =====================================================================
# 统计大于 10 或者大于 100 的异常高值占比（根据你的生物学背景调整阈值）
high_val_10 <- sum(te_vector > 10)
high_val_100 <- sum(te_vector > 100)
total_elements <- length(te_vector)

cat(sprintf("TE > 10  的元素个数: %d (占比: %.4f%%)\n", high_val_10, (high_val_10/total_elements)*100))
cat(sprintf("TE > 100 的元素个数: %d (占比: %.4f%%)\n", high_val_100, (high_val_100/total_elements)*100))


# =====================================================================
# 3. 绘制整体数值分布图（直方图 + 密度曲线）
# =====================================================================
# 创建绘图数据框
plot_df <- data.frame(TE_Value = te_vector)
x_limit <- quantile(te_vector, probs = 0.99)

p <- ggplot(plot_df, aes(x = TE_Value)) +
  # 绘制直方图 (y = after_stat(density) 确保 y 轴为密度，方便与密度曲线叠加)
  geom_histogram(aes(y = after_stat(density)), bins = 100, fill = "#34495e", color = "white", alpha = 0.7) +
  # 绘制平滑密度曲线
  geom_density(color = "#e74c3c", linewidth = 1) +
  # 限制 X 轴范围，切掉 1% 的极端大值，让主体分布更清晰
  xlim(0, x_limit) + 
  labs(
    title = "Distribution of Imputed Translation Efficiency (TE)",
    subtitle = paste0("Showing 0 to 99th percentile (Cutoff: ", round(x_limit, 2), ")"),
    x = "TE Value (rbRNA / (RNA + 1e-06))",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )
print(p)

write.csv(TE_imputed, output_te_imputed)
write.csv(TE_rescaled, output_te_rescaled)
