# MAGIC imputation using the GAT attention transition matrix
rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(ggplot2)
})

# =========================
# 0. Parameters
# =========================
output_root <- "/media/zenglab/result/lingyuan/STEM/imputation/magic_stagate_like_pipeline"
input_dir <- file.path(output_root, "03_attention")
out_dir <- file.path(output_root, "04_magic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_rds <- file.path(input_dir, "03_attention_graph.rds")
out_rds <- file.path(out_dir, "04_magic_imputed_object.rds")
out_rdata <- file.path(out_dir, "04_magic_imputed_object.RData")

t_diff <- 4
t_max <- 10
target_gene <- "UBE2C"
output_te_imputed <- file.path(out_dir, "04_mousebrain_combine_TE_imputed.csv")
output_te_rescaled <- file.path(out_dir, "04_mousebrain_combine_TE_rescaled.csv")

convert_to_seurat_matrix <- function(mat) {
  if (inherits(mat, "dgeMatrix")) return(as(mat, "dgCMatrix"))
  if (is.matrix(mat)) return(as(mat, "dgCMatrix"))
  mat
}

magic_impute_and_rescale <- function(M_t, D) {
  D_imputed <- M_t %*% D
  rownames(D_imputed) <- rownames(D)
  colnames(D_imputed) <- colnames(D)

  p99_original <- apply(D, 2, function(x) quantile(x, probs = 0.99))
  max_imputed <- apply(D_imputed, 2, max)
  max_imputed[max_imputed == 0] <- 1
  scale_factors <- p99_original / max_imputed
  D_rescaled <- sweep(D_imputed, 2, scale_factors, FUN = "*")

  list(imputed = D_imputed, rescaled = D_rescaled)
}

# =========================
# 1. Load attention graph
# =========================
states <- readRDS(input_rds)
M <- states@misc$magic_transition
if (is.null(M)) {
  stop("states@misc$magic_transition is missing. Rerun 03_GAT_attention_graph.R.")
}

M_t <- M
if (t_diff > 1) {
  for (i in 2:t_diff) {
    M_t <- M_t %*% M
  }
}

# =========================
# 2. Check MAGIC diffusion convergence
# =========================
delta_R2_list <- numeric(t_max)
D_prev <- t(GetAssayData(states, assay = "RNA", layer = "data"))
cat("开始计算图扩散与 R^2 收敛率 ...\n")
for (t in seq_len(t_max)) {
  D_curr <- M %*% D_prev

  col_sums_curr <- Matrix::colSums(D_curr)
  col_sums_curr[col_sums_curr == 0] <- 1
  D_curr_norm <- D_curr %*% Diagonal(x = 1 / col_sums_curr)

  col_sums_prev <- Matrix::colSums(D_prev)
  col_sums_prev[col_sums_prev == 0] <- 1
  D_prev_norm <- D_prev %*% Diagonal(x = 1 / col_sums_prev)

  SSE <- sum((D_curr_norm - D_prev_norm)^2)
  mean_curr <- mean(D_curr_norm)
  SST <- sum((D_curr_norm - mean_curr)^2)
  delta_R2_list[t] <- SSE / SST
  cat(sprintf("t = %2d | 1 - R^2 = %.5f\n", t, delta_R2_list[t]))
  D_prev <- D_curr
}

df_plot <- data.frame(t = seq_len(t_max), Delta_R2 = delta_R2_list)
below_threshold <- which(delta_R2_list < 0.05)
optimal_t <- if (length(below_threshold) >= 2) below_threshold[2] else NA
p_convergence <- ggplot(df_plot, aes(x = t, y = Delta_R2)) +
  geom_line(color = "#2c3e50", linewidth = 1) +
  geom_point(color = "#e74c3c", size = 3) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray50", linewidth = 1) +
  scale_x_continuous(breaks = seq_len(t_max)) +
  labs(
    title = "MAGIC Data Diffusion Convergence",
    subtitle = "Finding the optimal t (Threshold = 0.05)",
    x = "Diffusion Time (t)",
    y = expression(1 - R^2 ~~ paste("(", Data[t], ", ", Data[t-1], ")"))
  ) +
  theme_minimal()
if (!is.na(optimal_t)) {
  p_convergence <- p_convergence +
    geom_vline(xintercept = optimal_t, linetype = "dotted", color = "blue", linewidth = 1) +
    annotate("text", x = optimal_t + 0.2, y = max(delta_R2_list) / 2,
             label = paste("Optimal t =", optimal_t), color = "blue", hjust = 0, size = 5)
}
ggsave(file.path(out_dir, "04_MAGIC_diffusion_convergence.png"), p_convergence, width = 7, height = 5, dpi = 300)

# =========================
# 3. MAGIC imputation
# =========================
total_result <- magic_impute_and_rescale(
  M_t = M_t,
  D = t(GetAssayData(states, assay = "RNA", layer = "data"))
)
rb_result <- magic_impute_and_rescale(
  M_t = M_t,
  D = t(GetAssayData(states, assay = "rbRNA", layer = "data"))
)

imputed_totalRNA <- total_result$imputed
rescaled_totalRNA <- total_result$rescaled
imputed_rbRNA <- rb_result$imputed
rescaled_rbRNA <- rb_result$rescaled

stopifnot(isTRUE(all.equal(rownames(imputed_totalRNA), rownames(imputed_rbRNA))))
stopifnot(isTRUE(all.equal(colnames(imputed_totalRNA), colnames(imputed_rbRNA))))

TE_imputed <- imputed_rbRNA / (imputed_totalRNA + 1e-6)
TE_rescaled <- rescaled_rbRNA / (rescaled_totalRNA + 1e-6)

meta_imputed <- states@meta.data[rownames(imputed_totalRNA), , drop = FALSE]
magic_seu <- CreateSeuratObject(counts = t(imputed_totalRNA), meta.data = meta_imputed, project = "mousebrain_combine_MAGIC")
magic_seu[["rescaled_totalRNA"]] <- CreateAssayObject(data = convert_to_seurat_matrix(t(rescaled_totalRNA)))
magic_seu[["imputed_rbRNA"]] <- CreateAssayObject(data = convert_to_seurat_matrix(t(imputed_rbRNA)))
magic_seu[["rescaled_rbRNA"]] <- CreateAssayObject(data = convert_to_seurat_matrix(t(rescaled_rbRNA)))
magic_seu[["TE_imputed"]] <- CreateAssayObject(data = convert_to_seurat_matrix(t(TE_imputed)))
magic_seu[["TE_rescaled"]] <- CreateAssayObject(data = convert_to_seurat_matrix(t(TE_rescaled)))

# =========================
# 4. QC plots and output
# =========================
te_vector <- as.vector(TE_rescaled)
cat("============== TE_rescaled 全局数值统计 ==============\n")
cat(sprintf("最小值 (Min)   : %f\n", min(te_vector)))
cat(sprintf("最大值 (Max)   : %f\n", max(te_vector)))
cat(sprintf("均值   (Mean)  : %f\n", mean(te_vector)))
cat(sprintf("中位数 (Median): %f\n", median(te_vector)))
print(quantile(te_vector, probs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)))

x_limit <- quantile(te_vector, probs = 0.99)
p_te <- ggplot(data.frame(TE_Value = te_vector), aes(x = TE_Value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 100, fill = "#34495e", color = "white", alpha = 0.7) +
  geom_density(color = "#e74c3c", linewidth = 1) +
  xlim(0, x_limit) +
  labs(
    title = "Distribution of Imputed Translation Efficiency (TE)",
    subtitle = paste0("Showing 0 to 99th percentile (Cutoff: ", round(x_limit, 2), ")"),
    x = "TE Value (rbRNA / (RNA + 1e-06))",
    y = "Density"
  ) +
  theme_minimal()
ggsave(file.path(out_dir, "04_TE_rescaled_distribution.png"), p_te, width = 7, height = 5, dpi = 300)

if (target_gene %in% colnames(TE_imputed)) {
  plot_df <- data.frame(
    RNA_Expression = rescaled_totalRNA[, target_gene],
    Translation_Efficiency = TE_rescaled[, target_gene]
  )
  cor_pearson_imp <- cor(plot_df$Translation_Efficiency, plot_df$RNA_Expression, method = "pearson")
  cat(sprintf("%s TE 与插值后 RNA 的 Pearson 相关性: %.4f\n", target_gene, cor_pearson_imp))
  p_gene <- ggplot(plot_df, aes(x = RNA_Expression, y = Translation_Efficiency)) +
    geom_point(alpha = 0.6, color = "#3498db") +
    geom_smooth(method = "lm", color = "#e74c3c", se = TRUE) +
    labs(
      title = paste("Correlation Analysis for", target_gene),
      x = "Imputed Total RNA Expression",
      y = "Imputed Translation Efficiency (TE)"
    ) +
    theme_minimal()
  ggsave(file.path(out_dir, paste0("04_", target_gene, "_TE_RNA_correlation.png")), p_gene, width = 6, height = 5, dpi = 300)
} else {
  cat(sprintf("Skip marker QC: gene %s was not found in TE_imputed.\n", target_gene))
}

write.csv(TE_imputed, output_te_imputed)
write.csv(TE_rescaled, output_te_rescaled)
save(states, magic_seu, imputed_totalRNA, rescaled_totalRNA, imputed_rbRNA, rescaled_rbRNA,
     TE_imputed, TE_rescaled, file = out_rdata)
saveRDS(magic_seu, file = out_rds)
message("Done. Saved MAGIC object to: ", out_rds)
