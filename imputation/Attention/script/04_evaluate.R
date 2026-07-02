rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

output_root <- "/media/zenglab/result/lingyuan/STEM/imputation/magic_stagate_like_pipeline"
input_csv <- file.path(output_root, "03_attention", "gat_umap_spatial_qc", "03_gat_cluster_by_celltype_annotation_table.csv")
out_dir <- file.path(output_root, "03_attention", "gat_umap_spatial_qc", "gat_cluster_celltype_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
  "VLMC" = "#1f76b3",
  "VSMC" = "#00aeef",
  "Peri/VEC" = "#d3a59c",
  "DE/MEN" = "#b274e8",
  "MSN" = "#D96DA1",
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
  "Mix" = "#F5F5F5",
  "Unknown" = "#BDBDBD"
)

save_plot <- function(p, name, width = 10, height = 7) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = height)
}

tab <- read.csv(input_csv, check.names = FALSE)
tab$GATCluster <- as.character(tab$GATCluster)

long_df_raw <- tab %>%
  pivot_longer(
    cols = -GATCluster,
    names_to = "Celltype",
    values_to = "N"
  )

celltype_order <- intersect(names(celltype_colors), unique(long_df_raw$Celltype))
extra_celltypes <- setdiff(unique(long_df_raw$Celltype), celltype_order)
celltype_order <- c(celltype_order, extra_celltypes)

dominant_df <- long_df_raw %>%
  group_by(GATCluster) %>%
  summarise(
    NCells = sum(N),
    DominantCelltype = Celltype[which.max(N)],
    DominantN = max(N),
    Purity = ifelse(NCells > 0, DominantN / NCells, 0),
    .groups = "drop"
  ) %>%
  mutate(
    DominantCelltype = factor(DominantCelltype, levels = celltype_order),
    cluster_num = suppressWarnings(as.numeric(GATCluster))
  ) %>%
  arrange(DominantCelltype, desc(Purity), is.na(cluster_num), cluster_num, GATCluster)

cluster_order <- dominant_df$GATCluster

long_df <- long_df_raw %>%
  mutate(
    GATCluster = factor(GATCluster, levels = cluster_order),
    Celltype = factor(Celltype, levels = rev(celltype_order))
  ) %>%
  group_by(GATCluster) %>%
  mutate(
    ClusterTotal = sum(N),
    ClusterFraction = ifelse(ClusterTotal > 0, N / ClusterTotal, 0)
  ) %>%
  ungroup() %>%
  group_by(Celltype) %>%
  mutate(
    CelltypeTotal = sum(N),
    CelltypeFraction = ifelse(CelltypeTotal > 0, N / CelltypeTotal, 0)
  ) %>%
  ungroup()

dominant_df_out <- dominant_df %>%
  select(GATCluster, DominantCelltype, NCells, DominantN, Purity)

write.csv(
  dominant_df_out,
  file.path(out_dir, "gat_cluster_dominant_celltype_purity.csv"),
  row.names = FALSE
)

p_count_heatmap <- ggplot(long_df, aes(x = GATCluster, y = Celltype, fill = log10(N + 1))) +
  geom_tile(color = "grey85", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = median(log10(long_df$N + 1), na.rm = TRUE),
    name = "log10(count + 1)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "GAT cluster by celltype annotation",
    x = "GAT attention cluster",
    y = "Celltype annotation"
  )

save_plot(p_count_heatmap, "01_gat_cluster_celltype_count_heatmap_red_white_blue", width = 13, height = 7)

p_fraction_heatmap <- ggplot(long_df, aes(x = GATCluster, y = Celltype, fill = ClusterFraction)) +
  geom_tile(color = "grey85", linewidth = 0.2) +
  scale_fill_gradient(
    low = "white",
    high = "#B2182B",
    labels = scales::percent,
    name = "Fraction in cluster"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Celltype composition within each GAT cluster",
    x = "GAT attention cluster",
    y = "Celltype annotation"
  )

save_plot(p_fraction_heatmap, "02_gat_cluster_celltype_fraction_heatmap", width = 13, height = 7)

p_stacked <- ggplot(long_df, aes(x = GATCluster, y = ClusterFraction, fill = Celltype)) +
  geom_col(width = 0.85, color = "grey90", linewidth = 0.1) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = celltype_colors, na.value = "grey80") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = "Celltype composition per GAT cluster",
    x = "GAT attention cluster",
    y = "Celltype fraction",
    fill = "Celltype"
  )

save_plot(p_stacked, "03_gat_cluster_celltype_stacked_bar", width = 13, height = 7)

p_purity <- ggplot(dominant_df_out, aes(x = factor(GATCluster, levels = cluster_order), y = Purity, fill = DominantCelltype)) +
  geom_col(width = 0.85, color = "grey30", linewidth = 0.15) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_fill_manual(values = celltype_colors, na.value = "grey80") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = "Dominant celltype purity of each GAT cluster",
    x = "GAT attention cluster",
    y = "Dominant celltype fraction",
    fill = "Dominant celltype"
  )

save_plot(p_purity, "04_gat_cluster_dominant_celltype_purity", width = 13, height = 7)

p_cluster_size <- ggplot(dominant_df_out, aes(x = factor(GATCluster, levels = cluster_order), y = NCells, fill = DominantCelltype)) +
  geom_col(width = 0.85, color = "grey30", linewidth = 0.15) +
  scale_fill_manual(values = celltype_colors, na.value = "grey80") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = "GAT cluster size colored by dominant celltype",
    x = "GAT attention cluster",
    y = "Cells",
    fill = "Dominant celltype"
  )

save_plot(p_cluster_size, "05_gat_cluster_size_by_dominant_celltype", width = 13, height = 7)

message("Done. Saved plots to: ", out_dir)
