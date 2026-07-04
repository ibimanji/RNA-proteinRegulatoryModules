source("R/scml_core.R")
source("R/spatial_plotting.R")
source("R/setup.R")
suppressPackageStartupMessages(library(patchwork))

heat_data <- readRDS(file.path(output_dir, "03_heat_states.rds"))
states_seu <- heat_data$states_seu
rb_heat <- heat_data$rb_heat
total_heat <- heat_data$total_heat

cat("Building raw and heat-kernel TE matrices...\n")

total_raw <- get_assay_matrix(states_seu, assay = "RNA", layer = "counts")
rb_raw <- get_assay_matrix(states_seu, assay = "rbRNA", layer = "counts")

te_raw <- rb_raw / (total_raw + 1e-6)
te_heat <- rb_heat / (total_heat + 1e-6)

tau <- apply(
  total_heat,
  1,
  function(x) quantile(x, 0.75)
)

weight <- total_heat
for (i in seq_len(nrow(total_heat))) {
  weight[i, ] <- total_heat[i, ] / (total_heat[i, ] + tau[i])
}
te_heat_weighted <- te_heat * weight

states_out <- states_seu
states_out[["total_raw"]] <- CreateAssayObject(counts = total_raw)
states_out[["te_raw"]] <- CreateAssayObject(counts = te_raw)
states_out[["te_heat"]] <- CreateAssayObject(counts = te_heat)
states_out[["te_heat_weighted"]] <- CreateAssayObject(counts = te_heat_weighted)

genes_use <- c(
  "Nrgn", "C1ql2", "Wfs1", "Chgb",
  "Mbp", "Plp1", "C1qa", "Hexb", "Pvalb",
  "Enpp2", "Ttr",
  "Ly6c1", "Rarres2", "Ptgds"
)

cat("Plotting Total RNA markers...\n")

plot_total <- list()
for (gene in genes_use) {
  vmax_raw <- GetVmax(total_raw, gene)
  vmax_heat <- GetVmax(total_heat, gene)

  plot_total[[length(plot_total) + 1]] <- PlotSpatialGene2(
    states_out,
    gene = gene,
    assay = "total_raw",
    vmax = vmax_raw,
    title = paste0(gene, "\nTotal Raw")
  )
  plot_total[[length(plot_total) + 1]] <- PlotSpatialGene2(
    states_out,
    gene = gene,
    assay = "total_heat",
    vmax = vmax_heat,
    title = paste0(gene, "\nTotal Heat")
  )
}

p_total <- wrap_plots(plot_total, ncol = 2) +
  plot_annotation(title = "Total RNA (ntRNA + rbRNA)")

cat("Plotting TE markers...\n")

plot_te <- list()
for (gene in genes_use) {
  vmax_raw <- GetVmax(te_raw, gene)
  vmax_heat <- GetVmax(te_heat_weighted, gene)

  plot_te[[length(plot_te) + 1]] <- PlotSpatialGene2(
    states_out,
    gene = gene,
    assay = "te_raw",
    vmax = vmax_raw,
    title = paste0(gene, "\nTE Raw")
  )
  plot_te[[length(plot_te) + 1]] <- PlotSpatialGene2(
    states_out,
    gene = gene,
    assay = "te_heat_weighted",
    vmax = vmax_heat,
    title = paste0(gene, "\nTE Heat Weighted")
  )
}

p_te <- wrap_plots(plot_te, ncol = 2) +
  plot_annotation(title = "Translation Efficiency")

ggsave(
  filename = file.path(output_dir, "TotalRNA_Heat2.pdf"),
  plot = p_total,
  width = 18,
  height = 49
)

ggsave(
  filename = file.path(output_dir, "TE_Heat2.pdf"),
  plot = p_te,
  width = 18,
  height = 49
)

states_out[["total_heat"]] <- NULL

saveRDS(
  list(
    states_out = states_out,
    te_raw = te_raw,
    te_heat = te_heat,
    te_heat_weighted = te_heat_weighted,
    p_total = p_total,
    p_te = p_te
  ),
  file.path(output_dir, "04_marker_plots.rds")
)

cat("Saved outputs/TotalRNA_Heat2.pdf\n")
cat("Saved outputs/TE_Heat2.pdf\n")
