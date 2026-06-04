#MAGIC for data impuating
getwd()
rm(list = ls()); gc(); options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(matrixStats)
  library(anndata)
  library(reticulate)
  library(ggplot2)
  library(ggrepel)
})
#0.data loading 
states = read_h5ad("/media/zenglab/data/jinpu/statesHeLa/cellline_filtered_data_20250614.h5ad")
rb_raw_mat <- t(states$layers[["rbRNA"]]) 
nt_raw_mat = t(states$layers[["ntRNA"]])
total_raw_mat <- t(states$layers[["ntRNA"]]) + t(states$layers[["rbRNA"]]) 
meta = states$obs
states_seu <- CreateSeuratObject(
  counts = total_raw_mat,
  meta.data = meta,
  project = "statesHeLa"
)
states_seu[["rbRNA"]] <- CreateAssayObject(
  counts = rb_raw_mat
)
states_seu[["ntRNA"]] <- CreateAssayObject(
  counts = nt_raw_mat
)
states_seu <- subset(
  states_seu,
  subset = sample == "C3control"
)
#magic trail via R
#1.1 we focus on the WT here, to see if it will return the linear biological events(cell cycle related ones， especially the UBE2C、TPX2 related ones)
#normalize should be carefully considered，here we shall use totalRNA layer to calculate spike-in-factor, to evaluate the cell sequencing-lib
#With this sequencing-lib, we can use them to normalize rbRNA layer (and ntRNA layer)

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

states_seu <- FindMultiModalNeighbors(
  states_seu, 
  reduction.list = list("rbrna.pca", "ntrna.pca"), # here we shall mention that why ntRNA? totalRNA, is rb plus nt, if use total and rb as input, may lead to amplify the weight from rbRNA layer.
  dims.list = list(1:20, 1:20), 
  knn.graph.name = "wknn", 
  snn.graph.name = "wsnn"
)
nn_obj <- states_seu@neighbors$weighted.nn
idx <- nn_obj@nn.idx   
dist <- nn_obj@nn.dist 
n_cells <- nrow(idx)
k_total <- ncol(idx)
ka <- 10 #这里可以自适应完善，对于选择的knn的k，选取第ka个邻居的距离，计算celli与它的距离作为sigma的adaptive估值
if(ka > k_total) stop("ka 不能大于 FindMultiModalNeighbors 中计算的邻居数(k.param)")
sigma <- dist[, ka] #sigma calculation in cell wise
sigma[sigma == 0] <- 1e-03#这里没啥用没一个东西赋到这个值
affinity_values <- exp(-(dist / sigma)^2)
i_idx <- as.vector(row(idx))       
j_idx <- as.vector(idx)           
x_val <- as.vector(affinity_values)

A <- sparseMatrix(i = i_idx, j = j_idx, x = x_val, dims = c(n_cells, n_cells)) #affinity matrix
A_sym <- A + t(A)
diag(A_sym) <- max(A_sym) #self loop handling
row_sums <- rowSums(A_sym)
M <- A_sym / row_sums #probability matrix, also transition matrix for single step
t_diff <- 4 #for t step,default we use 4
M_t <- M
if (t_diff > 1) {
  for (i in 2:t_diff) {
    M_t <- M_t %*% M
  }
}#M-C equation

##test as coefficient ：
t_max <- 10
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

seurat_obj <- CreateSeuratObject(counts = t(imputed_totalRNA), meta.data = meta)
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


#UBE2C-UBE2C correlation cauculation validation
target_gene <- "UBE2C"
if (!target_gene %in% colnames(TE_imputed)) {
  stop(paste("未在矩阵中找到基因:", target_gene))
}
ube2c_te <- TE_rescaled[, target_gene]
ube2c_rna_imputed  <- rescaled_totalRNA[, target_gene]
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

write.csv(TE_imputed, '/media/zenglab/data/jinpu/statesHeLa/hela_control_TE_imputed.csv')
write.csv(TE_rescaled, '/media/zenglab/data/jinpu/statesHeLa/hela_control_TE_rescaled.csv')
