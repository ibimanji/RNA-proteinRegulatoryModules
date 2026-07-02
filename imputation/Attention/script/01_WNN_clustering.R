# WNN clustering before MixFind and manual annotation
rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# =========================
# 0. Parameters
# =========================
output_root <- "/media/zenglab/result/lingyuan/STEM/imputation/magic_stagate_like_pipeline"
input_dir <- file.path(output_root, "00_harmony")
out_dir <- file.path(output_root, "01_wnn_clustering")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_rds <- file.path(input_dir, "00_harmony_preprocessed.rds")
out_rds <- file.path(out_dir, "01_wnn_clustering.rds")
out_rdata <- file.path(out_dir, "01_wnn_clustering.RData")

wnn_dims <- 1:20
wnn_k_nn <- 50
wnn_prune_snn <- 0
wnn_cluster_resolution <- 2
wnn_algorithm <- 1
umap_min_dist <- 0.001
umap_spread <- 3

# =========================
# 1. Helpers
# =========================
save_dimplot_if_column_exists <- function(seu, group_col, file, title) {
  if (!group_col %in% colnames(seu@meta.data)) return(invisible(NULL))
  p <- DimPlot(seu, reduction = "expression.wnn.umap", group.by = group_col) +
    ggtitle(title)
  ggsave(file, p, width = 8, height = 6, dpi = 300)
  invisible(p)
}

# =========================
# 2. Load Harmony-preprocessed object
# =========================
states <- readRDS(input_rds)
if (!all(c("rbrna.harmony", "ntrna.harmony") %in% names(states@reductions))) {
  stop("Harmony reductions rbrna.harmony and ntrna.harmony are required. Rerun 00_Harmony_preprocess.R.")
}
if (is.null(states@misc$spatial_coords)) {
  stop("states@misc$spatial_coords is missing. Rerun 00_Harmony_preprocess.R.")
}
batch_var <- "protocol-replicate"
batch_var_seurat <- states@misc$harmony_batch_var
if (is.null(batch_var_seurat) || !batch_var_seurat %in% colnames(states@meta.data)) {
  batch_var_seurat <- make.names(batch_var)
}
if (!batch_var_seurat %in% colnames(states@meta.data)) {
  stop(sprintf("Batch variable '%s' is missing from states@meta.data. Rerun 00_Harmony_preprocess.R.", batch_var_seurat))
}

# =========================
# 3. WNN graph, UMAP and clustering
# =========================
states <- FindMultiModalNeighbors(
  states,
  reduction.list = list("rbrna.harmony", "ntrna.harmony"),
  dims.list = list(wnn_dims, wnn_dims),
  k.nn = wnn_k_nn,
  prune.SNN = wnn_prune_snn,
  weighted.nn.name = "expression.weighted.nn",
  knn.graph.name = "expression_wknn",
  snn.graph.name = "expression_wsnn",
  modality.weight.name = "rbRNA.weight",
  verbose = FALSE
)
states <- RunUMAP(
  states,
  nn.name = "expression.weighted.nn",
  reduction.name = "expression.wnn.umap",
  reduction.key = "exprWNNUMAP_",
  min.dist = umap_min_dist,
  spread = umap_spread,
  verbose = FALSE
)
states <- FindClusters(
  states,
  graph.name = "expression_wsnn",
  resolution = wnn_cluster_resolution,
  algorithm = wnn_algorithm,
  verbose = FALSE
)
states$expression_cluster <- as.character(Idents(states))

# =========================
# 4. Save WNN intermediate outputs
# =========================
p_cluster <- DimPlot(states, reduction = "expression.wnn.umap", group.by = "expression_cluster", label = TRUE) +
  ggtitle(paste0("Expression WNN clusters, resolution = ", wnn_cluster_resolution))
ggsave(file.path(out_dir, "01_expression_WNN_clusters.png"), p_cluster, width = 8, height = 6, dpi = 300)

save_dimplot_if_column_exists(
  states,
  group_col = batch_var_seurat,
  file = file.path(out_dir, paste0("01_expression_WNN_umap_by_", batch_var_seurat, ".png")),
  title = paste("Expression WNN UMAP by", batch_var)
)

weight_cols <- colnames(states@meta.data)[grepl("weight", colnames(states@meta.data))]
for (weight_col in weight_cols) {
  p_weight <- FeaturePlot(states, reduction = "expression.wnn.umap", features = weight_col) +
    ggtitle(paste("WNN modality weight:", weight_col))
  ggsave(file.path(out_dir, paste0("01_", weight_col, "_WNN_weight.png")), p_weight, width = 7, height = 6, dpi = 300)
}

umap_df <- as.data.frame(Embeddings(states, "expression.wnn.umap"))
umap_df$cell <- rownames(umap_df)
umap_df$expression_cluster <- states$expression_cluster
write.csv(umap_df, file.path(out_dir, "01_expression_WNN_umap_coordinates.csv"), row.names = FALSE)

cluster_table <- as.data.frame(table(states$expression_cluster))
colnames(cluster_table) <- c("expression_cluster", "n_cells")
write.csv(cluster_table, file.path(out_dir, "01_expression_WNN_cluster_sizes.csv"), row.names = FALSE)
write.csv(states@meta.data, file.path(out_dir, "01_expression_WNN_metadata.csv"))

save(states, file = out_rdata)
saveRDS(states, file = out_rds)
message("Done. Saved WNN clustering object to: ", out_rds)
