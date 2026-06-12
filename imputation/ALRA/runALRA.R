library(Seurat)
library(Matrix)
library(rsvd)

setwd("/media/zenglab/result/lingyuan/STEM/decon/ALRA")

source("alra.R")

obj <- readRDS("../states_with_plaque_info.rds")

target_type <- "14mWT"

cells_use <- rownames(obj@meta.data)[obj@meta.data$type == target_type]

obj_sub <- subset(obj, cells = cells_use)

counts <- GetAssayData(
  obj_sub,
  assay = "totalRNA",
  slot = "counts"
)

# counts 一般是 genes × cells
dim(counts)

# 转成 ALRA 要求的 cells × genes
A_counts <- t(as.matrix(counts))

# 做 ALRA 自带的 library-size normalize + log1p
A_norm <- normalize_data(A_counts)

# 自动选择 k
k_choice <- choose_k(A_norm)

k_choice$k

# 跑 ALRA
set.seed(1)
alra_out <- alra(A_norm, k = k_choice$k)

# 取最终补全结果
A_alra <- alra_out[[3]] # The actual adjusted, completed matrix, A_norm_rank_k_cor_sc=A_norm_rank_k_cor_sc

# 补回 cell names 和 gene names
rownames(A_alra) <- rownames(A_norm)
colnames(A_alra) <- colnames(A_norm)
#转回genes x cells
A_alra_gene_cell <- t(A_alra)
# 再确认一下
dim(A_alra_gene_cell)
head(rownames(A_alra_gene_cell))
head(colnames(A_alra_gene_cell))

alra_assay <- CreateAssayObject(data = A_alra_gene_cell)

obj_sub[["ALRA_totalRNA"]] <- alra_assay
DefaultAssay(obj_sub) <- "ALRA_totalRNA"

saveRDS(
  obj_sub,
  file = paste0("states_", target_type, "_ALRA.rds")
)


##################################### post check #########################################
cat("Before ALRA nonzero %:",
    100 * mean(A_norm > 0), "\n")

cat("After ALRA nonzero %:",
    100 * mean(A_alra > 0), "\n")
genes_per_cell_before <- rowSums(A_norm > 0)
genes_per_cell_after  <- rowSums(A_alra > 0)
summary(genes_per_cell_before)
summary(genes_per_cell_after)
boxplot(
  genes_per_cell_before,
  genes_per_cell_after,
  names = c("Before", "After"),
  ylab = "Detected genes per cell"
)

cells_per_gene_before <- colMeans(A_norm > 0)
cells_per_gene_after  <- colMeans(A_alra > 0)
summary(cells_per_gene_before)
summary(cells_per_gene_after)
plot(
  cells_per_gene_before,
  cells_per_gene_after,
  pch = 16,
  cex = 0.4,
  xlab = "Before ALRA: fraction of cells expressing gene",
  ylab = "After ALRA: fraction of cells expressing gene"
)
abline(0, 1, col = "red")

hist(
  A_norm[A_norm > 0],
  breaks = 100,
  main = "Before ALRA nonzero expression",
  xlab = "log-normalized expression"
)

hist(
  A_alra[A_alra > 0],
  breaks = 100,
  main = "After ALRA nonzero expression",
  xlab = "ALRA expression"
)
plot(
  density(A_norm[A_norm > 0]),
  main = "Nonzero expression distribution",
  xlab = "Expression"
)
lines(
  density(A_alra[A_alra > 0]),
  col = "blue"
)
legend(
  "topright",
  legend = c("Before", "After"),
  col = c("black", "blue"),
  lty = 1
)

