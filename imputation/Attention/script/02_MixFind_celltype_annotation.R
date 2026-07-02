# MixFind and manual celltype annotation after WNN clustering
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
output_root <- "/media/zenglab/result/lingyuan/STEM/RNA-proteinRegulatoryModules/imputation/attention/magic_stagate_like_pipeline"
input_dir <- file.path(output_root, "01_wnn_clustering")
out_dir <- file.path(output_root, "02_mixfind_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
mixfind_qc_dir <- file.path(out_dir, "mixfind_qc")
annotation_evidence_dir <- file.path(out_dir, "annotation_evidence")
celltype_plot_dir <- file.path(out_dir, "celltype_annotation")
celltype_spatial_dir <- file.path(celltype_plot_dir, "spatial_by_celltype")
for (plot_dir in c(mixfind_qc_dir, annotation_evidence_dir, celltype_plot_dir, celltype_spatial_dir)) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
}

input_rds <- file.path(input_dir, "01_wnn_clustering.rds")
out_rds <- file.path(out_dir, "02_mixfind_celltype_annotation.rds")
out_rdata <- file.path(out_dir, "02_mixfind_celltype_annotation.RData")
out_markers <- file.path(annotation_evidence_dir, "02_wnn_mixfind_all_markers.csv")
out_top10 <- file.path(annotation_evidence_dir, "02_wnn_mixfind_top10_markers.csv")

mix_small_cluster_n <- 10
mix_auto_threshold_quantile <- 0.95
# Optional manual MixFind thresholds. Names must be expression_cluster labels.
# Example: manual_mix_threshold <- c("0" = 5, "1" = 8, "2" = 6)
manual_mix_threshold <- c("0" = 5, "2" = 5, "1" = 5, "5" = 5, "16" = 3,"9" = 4, "7" = 4,"4" = 5)

marker_panels <- list(
  `TEPN Layer` = c("Fezf2","Pde1a","Rgs4","Lamp5","Cux2"),
  `TEPN CA` = c("Chgb","Rgs14"),
  `TEPN DGGRC` = c("C1ql2"),
  `TEPN MSN` = c("Ppp1r1b","Adora2a"),
  VAS = c("Rgs5", "Ptgds", "Ly6e", "Ly6c1", "Flt1"),
  INH = c("Gad1", "Pvalb", "Sst", "Npy"),
  `DE/MEN` = c("Sparc", "Prkcd", "Ntng1"),
  `CHO/PEP` = c("Tac2", "Gng8", "Cadps2"),
  AC = c("Aldoc", "Gja1", "Atp1b2"),
  OLG = c("Plp1", "Mobp", "Fth1", "Cnp", "Trf"),
  OPC = c("Ptprz1", "Pdgfra", "Sox10", "Olig2", "Cacng4"),
  MLG = c("Hexb", "C1qa", "Csf1r", "C1qc", "Ctss"),
  CHOR = c("Ttr","Enpp2"),
  EPEN = c("Rarres2")
)
marker_genes <- unique(unlist(marker_panels, use.names = FALSE))

celltype_colors_l2 <- c(
  "TEPN CA"    = "#77ed8f",
  "TEPN Layer"   = "#40d102",
  "TEPN DGGRC"        = "#a6e8a6",
  "TEPN MSN"          = "#D96DA1",
  "VLMC"         = "#1f76b3",
  "VSMC"         = "#00aeef",
  "PeriVEC"     = "#d3a59c",
  "OPC"       = "#667872",
  "DE/MEN"    = "#b274e8",
  "CHO/PEP"   = "#97f4f7",
  "OLG"       = "#e6db17",
  "CHOR"         = "#7f52a9",
  "EPEN"         = "#c4b0d4",
  "MLG"       = "#8597c6",
  "INH"       = "#3182BD",
  "AC"        = "#FDC06F",
  "Mix"       = "#FAFAFA"
)

# =========================
# 1. Helpers
# =========================
save_spatial_plot <- function(meta_df, group_col, file, title, point_size = 0.25, color_values = NULL) {
  plot_df <- meta_df
  plot_df$group_value <- plot_df[[group_col]]
  p <- ggplot(plot_df, aes(x = column, y = row, color = group_value)) +
    geom_point(size = point_size) +
    coord_fixed() +
    scale_y_reverse() +
    theme_void() +
    labs(title = title, color = group_col)
  if (!is.null(color_values)) {
    p <- p + scale_color_manual(values = color_values, na.value = "grey80", drop = FALSE)
  }
  if ("z" %in% colnames(plot_df) && dplyr::n_distinct(plot_df$z) <= 24) {
    p <- p + facet_wrap(~ z)
  }
  ggsave(file, p, width = 8, height = 7, dpi = 300)
  invisible(p)
}

save_spatial_plot_by_batch <- function(meta_df, group_col, batch_col, batch_label, file, title, point_size = 0.18, color_values = NULL) {
  if (!batch_col %in% colnames(meta_df)) {
    warning("Skip batch-faceted spatial plot because column '", batch_col, "' is missing.")
    return(invisible(NULL))
  }
  plot_df <- meta_df
  plot_df$group_value <- plot_df[[group_col]]
  plot_df$batch_value <- plot_df[[batch_col]]
  p <- ggplot(plot_df, aes(x = column, y = row, color = group_value)) +
    geom_point(size = point_size) +
    coord_fixed() +
    scale_y_reverse() +
    facet_wrap(~ batch_value) +
    theme_void() +
    labs(title = title, color = group_col) +
    theme(strip.text = element_text(size = 11, face = "bold"), legend.position = "right")
  if (!is.null(color_values)) {
    p <- p + scale_color_manual(values = color_values, na.value = "grey80", drop = FALSE)
  }
  ggsave(file, p, width = 14, height = 9, dpi = 300)
  invisible(p)
}

save_spatial_highlight_plots <- function(meta_df, group_col, prefix, title_prefix, point_size = 0.18,
                                         file_ext = ".png", device = NULL, highlight_colors = NULL,
                                         batch_col = NULL, batch_label = NULL) {
  groups <- sort(unique(as.character(meta_df[[group_col]])))
  for (grp in groups) {
    safe_grp <- gsub("[^A-Za-z0-9_]+", "_", grp)
    plot_df <- meta_df
    plot_df$highlight <- ifelse(as.character(plot_df[[group_col]]) == grp, grp, "Other")
    if (!is.null(batch_col) && batch_col %in% colnames(plot_df)) {
      plot_df$batch_value <- plot_df[[batch_col]]
    }
    highlight_color <- if (!is.null(highlight_colors) && grp %in% names(highlight_colors)) highlight_colors[[grp]] else "#d73027"
    p <- ggplot(plot_df, aes(x = column, y = row)) +
      geom_point(color = "grey85", size = point_size) +
      geom_point(data = plot_df[plot_df$highlight == grp, , drop = FALSE], color = highlight_color, size = point_size) +
      coord_fixed() +
      scale_y_reverse() +
      theme_void() +
      labs(title = paste(title_prefix, grp))
    if (!is.null(batch_col) && batch_col %in% colnames(plot_df)) {
      p <- p +
        facet_wrap(~ batch_value) +
        theme(strip.text = element_text(size = 11, face = "bold"))
      ggsave(paste0(prefix, "_", safe_grp, file_ext), p, width = 14, height = 9, dpi = 300, device = device)
    } else {
      if ("z" %in% colnames(plot_df) && dplyr::n_distinct(plot_df$z) <= 24) {
        p <- p + facet_wrap(~ z)
      }
      ggsave(paste0(prefix, "_", safe_grp, file_ext), p, width = 8, height = 7, dpi = 300, device = device)
    }
  }
}

save_marker_dotplot <- function(seu, group_col, marker_panels, file, title) {
  marker_genes <- unique(unlist(marker_panels, use.names = FALSE))
  present_features <- intersect(marker_genes, rownames(seu))
  missing_features <- setdiff(marker_genes, rownames(seu))
  if (length(missing_features) > 0) {
    writeLines(missing_features, paste0(file, ".missing_genes.txt"))
  }
  if (length(present_features) == 0) {
    warning("No marker genes are present for ", file)
    return(invisible(NULL))
  }

  gene_to_panel <- stack(marker_panels)
  colnames(gene_to_panel) <- c("gene", "celltype_marker")
  gene_to_panel <- gene_to_panel[gene_to_panel$gene %in% present_features, , drop = FALSE]
  gene_to_panel <- gene_to_panel[!duplicated(gene_to_panel$gene), , drop = FALSE]
  gene_to_panel$gene_label <- paste0(gene_to_panel$gene, " (", gene_to_panel$celltype_marker, ")")

  group_ids <- levels(factor(seu@meta.data[[group_col]]))
  if (group_col == "celltype_annotation") {
    celltype_order <- c("TEPN Layer","TEPN CA","TEPN DGGRC","TEPN MSN", "PeriVEC","VLMC","VSMC","INH", "DE/MEN", "CHO/PEP", "AC", "OLG", "OPC", "MLG", "CHOR","EPEN")
    group_ids <- intersect(celltype_order, group_ids)
  }
  data_mat <- GetAssayData(seu, assay = DefaultAssay(seu), layer = "data")
  avg_expr <- sapply(group_ids, function(group_id) {
    cells_use <- rownames(seu@meta.data)[as.character(seu@meta.data[[group_col]]) == group_id]
    Matrix::rowMeans(data_mat[gene_to_panel$gene, cells_use, drop = FALSE])
  })
  avg_expr <- as.matrix(avg_expr)
  rownames(avg_expr) <- gene_to_panel$gene
  colnames(avg_expr) <- group_ids

  z_mat <- t(scale(t(as.matrix(avg_expr))))
  z_mat[is.na(z_mat)] <- 0

  pct_expr <- sapply(group_ids, function(group_id) {
    cells_use <- rownames(seu@meta.data)[as.character(seu@meta.data[[group_col]]) == group_id]
    Matrix::rowMeans(data_mat[gene_to_panel$gene, cells_use, drop = FALSE] > 0) * 100
  })
  pct_expr <- as.matrix(pct_expr)

  plot_df <- as.data.frame(as.table(z_mat))
  colnames(plot_df) <- c("gene", "group", "zscore")
  pct_df <- as.data.frame(as.table(pct_expr))
  colnames(pct_df) <- c("gene", "group", "pct_expr")
  plot_df <- merge(plot_df, pct_df, by = c("gene", "group"))
  plot_df <- merge(plot_df, gene_to_panel[, c("gene", "gene_label", "celltype_marker")], by = "gene")

  gene_label_order <- setNames(gene_to_panel$gene_label, gene_to_panel$gene)
  plot_df$gene <- factor(plot_df$gene, levels = rev(gene_to_panel$gene))
  plot_df$gene_label <- factor(plot_df$gene_label, levels = rev(gene_label_order))
  plot_df$group <- factor(plot_df$group, levels = group_ids)
  plot_df$zscore_plot <- pmax(pmin(plot_df$zscore, 2), -2)

  p <- ggplot(plot_df, aes(x = group, y = gene_label)) +
    geom_point(aes(size = pct_expr, color = zscore_plot)) +
    scale_color_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-2, 2),
      name = "Z-score"
    ) +
    scale_size(range = c(0.5, 6), name = "% expressed") +
    labs(
      title = title,
      x = group_col,
      y = "Marker gene (celltype marker)"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file,
    plot = p,
    width = max(9, length(group_ids) * 0.45),
    height = max(8, length(present_features) * 0.22),
    device = cairo_pdf
  )
  invisible(p)
}

make_spatial_meta <- function(seu, coord_cols = c("column", "row", "z")) {
  meta_df <- seu@meta.data
  missing_cols <- setdiff(coord_cols, colnames(meta_df))
  if (length(missing_cols) > 0) {
    coord_df <- as.data.frame(seu@misc$spatial_coords[Cells(seu), missing_cols, drop = FALSE])
    meta_df <- cbind(meta_df, coord_df)
  }
  meta_df
}

sort_cluster_labels <- function(x) {
  x <- unique(as.character(x))
  numeric_x <- suppressWarnings(as.numeric(x))
  x[order(is.na(numeric_x), numeric_x, x)]
}

# =========================
# 2. Load WNN clustering object
# =========================
states <- readRDS(input_rds)
if (!"expression_cluster" %in% colnames(states@meta.data)) {
  stop("expression_cluster is missing. Rerun 01_WNN_clustering.R.")
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
# 3. MixFind
# =========================
umap_coords <- Embeddings(states, "expression.wnn.umap")
umap_df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  expression_cluster = states$expression_cluster,
  row.names = Cells(states)
)
centroids <- umap_df %>%
  group_by(expression_cluster) %>%
  summarise(centroid_x = mean(UMAP1), centroid_y = mean(UMAP2), .groups = "drop")

states$distance2centroid_wnn <- NA_real_
for (cluster_id in centroids$expression_cluster) {
  idx <- which(states$expression_cluster == cluster_id)
  centroid <- centroids[centroids$expression_cluster == cluster_id, ]
  dm <- sqrt((umap_coords[idx, 1] - centroid$centroid_x)^2 + (umap_coords[idx, 2] - centroid$centroid_y)^2)
  states$distance2centroid_wnn[idx] <- dm
}

cluster_sizes <- table(states$expression_cluster)
small_clusters <- names(cluster_sizes[cluster_sizes < mix_small_cluster_n])
threshold_df <- states@meta.data %>%
  group_by(expression_cluster) %>%
  summarise(auto_threshold = quantile(distance2centroid_wnn, probs = mix_auto_threshold_quantile, na.rm = TRUE),
            n_cells = n(), .groups = "drop")
threshold_df$manual_threshold <- NA_real_
if (!is.null(manual_mix_threshold)) {
  threshold_df$manual_threshold <- unname(manual_mix_threshold[as.character(threshold_df$expression_cluster)])
}
threshold_df$threshold_used <- ifelse(is.na(threshold_df$manual_threshold), threshold_df$auto_threshold, threshold_df$manual_threshold)
threshold_df$is_small_cluster <- threshold_df$expression_cluster %in% small_clusters
write.csv(threshold_df, file.path(out_dir, "02_mixfind_thresholds.csv"), row.names = FALSE)

states$is_mix_wnn <- "False"
states$is_mix_wnn[states$expression_cluster %in% small_clusters] <- "True"
for (i in seq_len(nrow(threshold_df))) {
  cluster_id <- threshold_df$expression_cluster[i]
  idx <- which(states$expression_cluster == cluster_id)
  mix_idx <- idx[states$distance2centroid_wnn[idx] > threshold_df$threshold_used[i]]
  states$is_mix_wnn[mix_idx] <- "True"
}
states$is_mix_wnn <- factor(states$is_mix_wnn, levels = c("False", "True"))
states$wnn_cluster_mix <- as.character(states$expression_cluster)
states$wnn_cluster_mix[states$is_mix_wnn == "True"] <- "Mix"
states$wnn_cluster_mix <- factor(states$wnn_cluster_mix)

umap_df$is_mix_wnn <- states$is_mix_wnn
umap_df$wnn_cluster_mix <- states$wnn_cluster_mix
p_centroids <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = expression_cluster)) +
  geom_point(size = 0.4) +
  geom_point(data = centroids, aes(x = centroid_x, y = centroid_y), inherit.aes = FALSE, color = "red", size = 2) +
  theme_minimal() +
  ggtitle("WNN UMAP with cluster centroids")
ggsave(file.path(mixfind_qc_dir, "02_mixfind_centroids.png"), p_centroids, width = 7, height = 6, dpi = 300)

p_mix <- DimPlot(states, reduction = "expression.wnn.umap", group.by = "is_mix_wnn") +
  ggtitle("MixFind classification")
ggsave(file.path(mixfind_qc_dir, "02_mixfind_mix_classification.png"), p_mix, width = 7, height = 6, dpi = 300)

p_mix_cluster <- DimPlot(states, reduction = "expression.wnn.umap", group.by = "wnn_cluster_mix", label = TRUE) +
  ggtitle("WNN clusters after MixFind")
ggsave(file.path(mixfind_qc_dir, "02_mixfind_wnn_cluster_mix.png"), p_mix_cluster, width = 8, height = 7, dpi = 300)

for (cluster_id in sort(unique(states$expression_cluster))) {
  safe_cluster <- gsub("[^A-Za-z0-9_]+", "_", cluster_id)
  dist_df <- states@meta.data[states$expression_cluster == cluster_id, , drop = FALSE]
  threshold_used <- threshold_df$threshold_used[threshold_df$expression_cluster == cluster_id]
  p_dist <- ggplot(dist_df, aes(x = distance2centroid_wnn)) +
    geom_histogram(bins = 50, fill = "grey70", color = "white") +
    geom_vline(xintercept = threshold_used, color = "red") +
    theme_minimal() +
    ggtitle(paste("MixFind distance threshold, cluster", cluster_id))
  ggsave(file.path(mixfind_qc_dir, paste0("02_mixfind_cluster_", safe_cluster, "_distance.png")), p_dist, width = 6, height = 4, dpi = 300)
}

# =========================
# 4. Annotation evidence
# =========================
DefaultAssay(states) <- "RNA"
Idents(states) <- "wnn_cluster_mix"
markers <- FindAllMarkers(
  object = states,
  test.use = "wilcox",
  only.pos = TRUE,
  logfc.threshold = 0.1,
  min.pct = 0.25
)
all_markers <- markers %>%
  dplyr::select(gene, everything()) %>%
  dplyr::filter(p_val < 0.05)
top10 <- all_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE)
write.csv(all_markers, out_markers, row.names = FALSE)
write.csv(top10, out_top10, row.names = FALSE)

save_marker_dotplot(
  seu = states,
  group_col = "wnn_cluster_mix",
  marker_panels = marker_panels,
  file = file.path(annotation_evidence_dir, "02_marker_panel_dotplot_by_wnn_cluster_mix.pdf"),
  title = "Marker panels by WNN cluster after MixFind"
)

spatial_meta <- make_spatial_meta(states)
save_spatial_plot(spatial_meta, "wnn_cluster_mix", file.path(annotation_evidence_dir, "02_spatial_by_wnn_cluster_mix.png"),
                  "Spatial distribution by WNN cluster after MixFind")
save_spatial_plot_by_batch(spatial_meta, "wnn_cluster_mix", batch_var_seurat, batch_var,
                           file.path(annotation_evidence_dir, paste0("02_spatial_by_wnn_cluster_mix_split_by_", batch_var_seurat, ".png")),
                           paste("Spatial distribution by WNN cluster after MixFind, split by", batch_var))
save_spatial_highlight_plots(spatial_meta, "wnn_cluster_mix", file.path(annotation_evidence_dir, "02_spatial_cluster_highlight"),
                             "Spatial highlight cluster", batch_col = batch_var_seurat, batch_label = batch_var)

# =========================
# 5. Manual celltype annotation
# =========================
# Fill this mapping after checking:
# - annotation_evidence/02_wnn_mixfind_top10_markers.csv
# - annotation_evidence/02_marker_panel_dotplot_by_wnn_cluster_mix.pdf
# - annotation_evidence/02_spatial_by_wnn_cluster_mix.png
# - annotation_evidence/02_spatial_cluster_highlight_*.png, one cluster per PNG and split by protocol-replicate
cluster_levels <- sort_cluster_labels(states$wnn_cluster_mix)
cluster_to_celltype <- setNames(paste0("cluster_", cluster_levels), cluster_levels)
if ("Mix" %in% names(cluster_to_celltype)) cluster_to_celltype["Mix"] <- "Mix"

# Example:
cluster_to_celltype <- c(
  "0" = "TEPN Layer",
  "1" = "OLG",
  "2" = "AC",
  "3" = "TEPN Layer",
  "4" = "PeriVEC",
  "5" = "OLG",
  "6" = "MLG",
  "7" = "TEPN CA",
  "8" = "AC",
  "9" = "INH",
  "10" = "TEPN DGGRC",
  "11" = "OLG",
  "12" = "CHOR",
  "13" = "Mix",
  "14" = "TEPN MSN",
  "15" = "VLMC",
  "16" = "DE/MEN",
  "17" = "Mix",
  "18" = "EPEN",
  "19" = "CHO/PEP",
  "20" = "VSMC",
  "Mix" = "Mix"
)

states$celltype_annotation <- unname(cluster_to_celltype[as.character(states$wnn_cluster_mix)])
states$celltype_annotation[is.na(states$celltype_annotation)] <- "Unknown"
states$celltype_annotation <- factor(states$celltype_annotation)

# =========================
# 6. Post-annotation plots
# =========================
celltype_levels <- levels(states$celltype_annotation)
celltype_colors_use <- celltype_colors_l2[celltype_levels]
missing_celltype_colors <- setdiff(celltype_levels, names(celltype_colors_l2))
if (length(missing_celltype_colors) > 0) {
  warning("Missing colors for celltypes: ", paste(missing_celltype_colors, collapse = ", "))
  celltype_colors_use[missing_celltype_colors] <- grDevices::hcl.colors(length(missing_celltype_colors), palette = "Set 3")
}

p_celltype <- DimPlot(
  states,
  reduction = "expression.wnn.umap",
  group.by = "celltype_annotation",
  label = TRUE,
  cols = celltype_colors_use
) +
  ggtitle("Manual celltype annotation")
ggsave(file.path(celltype_plot_dir, "02_celltype_annotation_umap.png"), p_celltype, width = 8, height = 7, dpi = 300)

save_marker_dotplot(
  seu = states,
  group_col = "celltype_annotation",
  marker_panels = marker_panels,
  file = file.path(celltype_plot_dir, "02_marker_panel_dotplot_by_celltype.pdf"),
  title = "Marker panels by annotated celltype"
)

spatial_meta <- make_spatial_meta(states)
save_spatial_plot(
  spatial_meta,
  "celltype_annotation",
  file.path(celltype_plot_dir, "02_spatial_by_celltype.png"),
  "Spatial distribution by annotated celltype",
  color_values = celltype_colors_use
)
save_spatial_plot_by_batch(
  spatial_meta,
  "celltype_annotation",
  batch_var_seurat,
  batch_var,
  file.path(celltype_plot_dir, paste0("02_spatial_by_celltype_split_by_", batch_var_seurat, ".png")),
  paste("Spatial distribution by annotated celltype, split by", batch_var),
  color_values = celltype_colors_use
)
save_spatial_highlight_plots(
  spatial_meta,
  "celltype_annotation",
  file.path(celltype_spatial_dir, "02_spatial_celltype_highlight"),
  "Spatial highlight celltype",
  file_ext = ".pdf",
  device = cairo_pdf,
  highlight_colors = celltype_colors_use,
  batch_col = batch_var_seurat,
  batch_label = batch_var
)

save(states, file = out_rdata)
saveRDS(states, file = out_rds)
message("Done. Saved MixFind + annotation object to: ", out_rds)
