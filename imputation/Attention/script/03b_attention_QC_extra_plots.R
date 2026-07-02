rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(Seurat)
})

output_root <- "/media/zenglab/result/lingyuan/STEM/imputation/magic_stagate_like_pipeline"
attention_dir <- file.path(output_root, "03_attention")
out_dir <- file.path(attention_dir, "extra_qc_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

attention_rds <- file.path(attention_dir, "03_attention_graph.rds")
edges_csv <- file.path(attention_dir, "03_gat_attention_edges.csv")

states <- readRDS(attention_rds)
batch_var <- "protocol-replicate"
spatial_sample_col <- make.names(batch_var)
if (!is.null(states@misc$harmony_batch_var) &&
    states@misc$harmony_batch_var %in% colnames(states@meta.data)) {
  spatial_sample_col <- states@misc$harmony_batch_var
}
if (!spatial_sample_col %in% colnames(states@meta.data)) {
  warning("Spatial sample column '", spatial_sample_col, "' is missing. Sample-faceted QC plots will use all cells together.")
  states@meta.data[[spatial_sample_col]] <- "all_cells"
}

if (!is.null(states@misc$gat_attention$edges)) {
  edges <- states@misc$gat_attention$edges
} else {
  edges <- read.csv(edges_csv, stringsAsFactors = FALSE)
}

# =========================
# 1. Max weight / effective neighbor number
# =========================
cell_weight_qc <- edges %>%
  group_by(Cell1, SpatialSample, Celltype1) %>%
  summarise(
    n_edges = n(),
    max_weight = max(Weight),
    effective_neighbors = 1 / sum(Weight^2),
    same_real_weight_sum = sum(Weight[CelltypeEdgeClass == "Same real celltype"]),
    cross_real_weight_sum = sum(Weight[CelltypeEdgeClass == "Cross real celltype"]),
    mix_unknown_weight_sum = sum(Weight[CelltypeEdgeClass == "Mix/Unknown involved"]),
    .groups = "drop"
  )

write.csv(cell_weight_qc, file.path(out_dir, "03_qc_cell_attention_summary.csv"), row.names = FALSE)

p_max <- ggplot(cell_weight_qc, aes(x = max_weight)) +
  geom_histogram(bins = 80, fill = "#2c3e50", color = "white") +
  theme_minimal() +
  labs(title = "Max attention weight per cell", x = "Max attention weight", y = "Cells")
ggsave(file.path(out_dir, "03_qc_max_weight_per_cell.png"), p_max, width = 7, height = 5, dpi = 300)

p_eff <- ggplot(cell_weight_qc, aes(x = effective_neighbors)) +
  geom_histogram(bins = 80, fill = "#5b8fd9", color = "white") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, quantile(cell_weight_qc$effective_neighbors, 0.99, na.rm = TRUE))) +
  labs(title = "Effective attention neighbors per cell", x = "1 / sum(weight^2)", y = "Cells")
ggsave(file.path(out_dir, "03_qc_effective_neighbors_per_cell.png"), p_eff, width = 7, height = 5, dpi = 300)

p_eff_type <- ggplot(cell_weight_qc, aes(x = Celltype1, y = effective_neighbors)) +
  geom_boxplot(outlier.size = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, quantile(cell_weight_qc$effective_neighbors, 0.99, na.rm = TRUE))) +
  labs(title = "Effective attention neighbors by celltype", x = "Celltype", y = "Effective neighbors")
ggsave(file.path(out_dir, "03_qc_effective_neighbors_by_celltype.png"), p_eff_type, width = 11, height = 6, dpi = 300)

# =========================
# 2. Celltype-to-celltype attention heatmaps
# =========================
celltype_weight <- edges %>%
  group_by(Celltype1, Celltype2) %>%
  summarise(
    total_weight = sum(Weight),
    mean_weight = mean(Weight),
    n_edges = n(),
    .groups = "drop"
  )

write.csv(celltype_weight, file.path(out_dir, "03_qc_celltype_to_celltype_attention.csv"), row.names = FALSE)

p_ct_total <- ggplot(celltype_weight, aes(x = Celltype2, y = Celltype1, fill = log10(total_weight + 1e-6))) +
  geom_tile() +
  scale_fill_viridis_c(name = "log10 total weight") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Celltype-to-celltype total attention", x = "Neighbor celltype", y = "Source celltype")
ggsave(file.path(out_dir, "03_qc_celltype_attention_total_heatmap.png"), p_ct_total, width = 9, height = 8, dpi = 300)

p_ct_mean <- ggplot(celltype_weight, aes(x = Celltype2, y = Celltype1, fill = mean_weight)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean weight") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Celltype-to-celltype mean attention", x = "Neighbor celltype", y = "Source celltype")
ggsave(file.path(out_dir, "03_qc_celltype_attention_mean_heatmap.png"), p_ct_mean, width = 9, height = 8, dpi = 300)

# =========================
# 3. Edge class by sample
# =========================
sample_edge_class <- edges %>%
  group_by(SpatialSample, CelltypeEdgeClass) %>%
  summarise(
    n_edges = n(),
    total_weight = sum(Weight),
    mean_weight = mean(Weight),
    .groups = "drop"
  ) %>%
  group_by(SpatialSample) %>%
  mutate(
    edge_fraction = n_edges / sum(n_edges),
    weight_fraction = total_weight / sum(total_weight)
  ) %>%
  ungroup()

write.csv(sample_edge_class, file.path(out_dir, "03_qc_edge_class_by_sample.csv"), row.names = FALSE)

p_sample_edge <- ggplot(sample_edge_class, aes(x = SpatialSample, y = weight_fraction, fill = CelltypeEdgeClass)) +
  geom_col() +
  theme_bw() +
  labs(title = "Attention weight fraction by sample and edge class", x = batch_var, y = "Weight fraction")
ggsave(file.path(out_dir, "03_qc_edge_class_weight_fraction_by_sample.png"), p_sample_edge, width = 8, height = 5, dpi = 300)

# =========================
# 4. Weight vs spatial distance / WNN similarity
# =========================
set.seed(1)
plot_edges <- edges
if (nrow(plot_edges) > 300000) {
  plot_edges <- plot_edges[sample(seq_len(nrow(plot_edges)), 300000), ]
}

p_dist <- ggplot(plot_edges, aes(x = SpatialDistance, y = Weight, color = CelltypeEdgeClass)) +
  geom_point(alpha = 0.08, size = 0.2) +
  theme_minimal() +
  coord_cartesian(y = c(0, quantile(plot_edges$Weight, 0.995, na.rm = TRUE))) +
  labs(title = "Attention weight vs spatial distance", x = "Spatial distance", y = "Attention weight")
ggsave(file.path(out_dir, "03_qc_weight_vs_spatial_distance.png"), p_dist, width = 8, height = 5, dpi = 300)

p_wnn <- ggplot(plot_edges, aes(x = WNNSimilarity, y = Weight, color = CelltypeEdgeClass)) +
  geom_point(alpha = 0.08, size = 0.2) +
  theme_minimal() +
  coord_cartesian(y = c(0, quantile(plot_edges$Weight, 0.995, na.rm = TRUE))) +
  labs(title = "Attention weight vs WNN similarity", x = "WNN similarity", y = "Attention weight")
ggsave(file.path(out_dir, "03_qc_weight_vs_wnn_similarity.png"), p_wnn, width = 8, height = 5, dpi = 300)

# =========================
# 5. Spatial edge overlay
# =========================
make_spatial_meta <- function(seu) {
  meta <- seu@meta.data
  meta$Cell <- rownames(meta)
  
  coords <- as.data.frame(seu@misc$spatial_coords[rownames(meta), , drop = FALSE])
  coords$Cell <- rownames(coords)
  
  coord_cols <- setdiff(colnames(coords), "Cell")
  
  # 避免 column/row/z 重复：如果 meta 里已经有坐标列，就先删掉 meta 里的旧坐标列
  duplicated_coord_cols <- intersect(coord_cols, colnames(meta))
  if (length(duplicated_coord_cols) > 0) {
    meta <- meta[, setdiff(colnames(meta), duplicated_coord_cols), drop = FALSE]
  }
  
  spatial_meta <- dplyr::left_join(meta, coords, by = "Cell")
  
  if (!all(c("column", "row") %in% colnames(spatial_meta))) {
    stop("spatial_meta must contain columns: column and row.")
  }
  
  spatial_meta
}
spatial_meta <- make_spatial_meta(states)

plot_spatial_edges <- function(sample_id, max_edges = 3000, min_weight_quantile = 0.98) {
  sample_edges <- edges %>%
    filter(SpatialSample == sample_id)

  if (nrow(sample_edges) == 0) return(invisible(NULL))

  cutoff <- quantile(sample_edges$Weight, min_weight_quantile, na.rm = TRUE)
  sample_edges <- sample_edges %>%
    filter(Weight >= cutoff)

  if (nrow(sample_edges) > max_edges) {
    sample_edges <- sample_edges[sample(seq_len(nrow(sample_edges)), max_edges), ]
  }

  sample_meta <- spatial_meta %>%
    filter(.data[[spatial_sample_col]] == sample_id)

  edge_df <- sample_edges %>%
    left_join(
      sample_meta %>% select(Cell, column, row),
      by = c("Cell1" = "Cell")
    ) %>%
    rename(x1 = column, y1 = row) %>%
    left_join(
      sample_meta %>% select(Cell, column, row),
      by = c("Cell2" = "Cell")
    ) %>%
    rename(x2 = column, y2 = row)

  p <- ggplot() +
    geom_point(
      data = sample_meta,
      aes(x = column, y = row),
      color = "grey85",
      size = 0.15
    ) +
    geom_segment(
      data = edge_df,
      aes(x = x1, y = y1, xend = x2, yend = y2, color = CelltypeEdgeClass, alpha = Weight),
      linewidth = 0.15
    ) +
    coord_fixed() +
    scale_y_reverse() +
    scale_alpha(range = c(0.15, 0.9)) +
    theme_void() +
    labs(title = paste0("High-weight attention edges: ", sample_id), color = "Edge class")

  safe_sample <- gsub("[^A-Za-z0-9_]+", "_", sample_id)
  ggsave(file.path(out_dir, paste0("03_qc_spatial_edge_overlay_", safe_sample, ".png")),
         p, width = 8, height = 7, dpi = 300)
}

for (sample_id in sort(unique(edges$SpatialSample))) {
  plot_spatial_edges(sample_id)
}

# =========================
# 6. Optional: spatial map of max attention / effective neighbors
# =========================
cell_weight_qc2 <- cell_weight_qc %>%
  rename(Cell = Cell1)

plot_meta <- spatial_meta %>%
  left_join(cell_weight_qc2, by = "Cell")
plot_meta$sample_value <- plot_meta[[spatial_sample_col]]

p_max_spatial <- ggplot(plot_meta, aes(x = column, y = row, color = max_weight)) +
  geom_point(size = 0.2) +
  coord_fixed() +
  scale_y_reverse() +
  facet_wrap(~ sample_value) +
  scale_color_viridis_c() +
  theme_void() +
  labs(title = "Spatial map of max attention weight", color = "Max weight")
ggsave(file.path(out_dir, "03_qc_spatial_max_attention_weight_by_sample.png"),
       p_max_spatial, width = 14, height = 9, dpi = 300)

p_eff_spatial <- ggplot(plot_meta, aes(x = column, y = row, color = effective_neighbors)) +
  geom_point(size = 0.2) +
  coord_fixed() +
  scale_y_reverse() +
  facet_wrap(~ sample_value) +
  scale_color_viridis_c(option = "magma") +
  theme_void() +
  labs(title = "Spatial map of effective attention neighbors", color = "Effective neighbors")
ggsave(file.path(out_dir, "03_qc_spatial_effective_neighbors_by_sample.png"),
       p_eff_spatial, width = 14, height = 9, dpi = 300)

message("Done. Extra attention QC plots saved to: ", out_dir)
