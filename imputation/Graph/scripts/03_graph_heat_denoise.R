source("R/scml_core.R")
source("R/graph_filtering.R")
source("R/setup.R")

scml_data <- readRDS(file.path(output_dir, "02_scml_states.rds"))
states_seu <- scml_data$states_seu
SCML_embedding <- scml_data$SCML_embedding

cat("Computing SCML graph spectrum...\n")

W <- BuildSCMLGraph(SCML_embedding, k = 30)
spec <- ComputeGraphSpectrum(W, n_eigs = 400)
lambda <- spec$lambda
V <- spec$V

spectrum_df <- data.frame(
  idx = seq_along(lambda),
  lambda = sort(lambda)
)

pdf(file.path(output_dir, "SCML_spectrum.pdf"))
print(
  ggplot(spectrum_df, aes(idx, lambda)) +
    geom_line() +
    geom_point(size = 1) +
    theme_classic()
)
lambda_sort <- sort(lambda)
cum_energy <- cumsum(lambda_sort) / sum(lambda_sort)
plot(
  cum_energy,
  type = "l",
  lwd = 2,
  ylim = c(0, 1),
  ylab = "Cumulative Energy",
  xlab = "Eigenvalue Rank"
)
abline(h = c(0.8, 0.9, 0.95), lty = 2)
dev.off()

cat("Running graph heat-kernel denoising...\n")

g_heat <- GraphFilterKernel(lambda = lambda, type = "heat", beta = 3)

rb_raw <- get_assay_matrix(states_seu, assay = "rbRNA", layer = "counts")
nt_raw <- get_assay_matrix(states_seu, assay = "ntRNA", layer = "counts")
total_raw <- get_assay_matrix(states_seu, assay = "RNA", layer = "counts")

rb_heat <- GraphSpectralDenoise(expr = rb_raw, V = V, g = g_heat)
rownames(rb_heat) <- rownames(rb_raw)
colnames(rb_heat) <- colnames(rb_raw)

nt_heat <- GraphSpectralDenoise(expr = nt_raw, V = V, g = g_heat)
rownames(nt_heat) <- rownames(nt_raw)
colnames(nt_heat) <- colnames(nt_raw)

total_heat <- GraphSpectralDenoise(expr = total_raw, V = V, g = g_heat)
rownames(total_heat) <- rownames(total_raw)
colnames(total_heat) <- colnames(total_raw)

states_seu[["rb_heat"]] <- CreateAssayObject(counts = rb_heat)
states_seu[["nt_heat"]] <- CreateAssayObject(counts = nt_heat)
states_seu[["total_heat"]] <- CreateAssayObject(counts = total_heat)

heat_data <- scml_data
heat_data$states_seu <- states_seu
heat_data$graph_spectrum <- spec
heat_data$graph_filter <- g_heat
heat_data$rb_raw <- rb_raw
heat_data$nt_raw <- nt_raw
heat_data$total_raw <- total_raw
heat_data$rb_heat <- rb_heat
heat_data$nt_heat <- nt_heat
heat_data$total_heat <- total_heat

saveRDS(heat_data, file.path(output_dir, "03_heat_states.rds"))

cat("Saved outputs/03_heat_states.rds\n")
