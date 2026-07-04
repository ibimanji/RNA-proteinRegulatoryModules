source("R/scml_core.R")
source("R/setup.R")

preprocessed <- readRDS(file.path(output_dir, "01_preprocessed_states.rds"))
states_seu <- preprocessed$states_seu

cat("Extracting Harmony embeddings...\n")

Emb_rb <- Embeddings(states_seu, reduction = "rbrna.harmony")[, 1:30]
Emb_nt <- Embeddings(states_seu, reduction = "ntrna.harmony")[, 1:30]

cat("Building rbRNA graph...\n")

W_rb <- build_knn_graph(emb = Emb_rb, k = 30, symmetrize = "mean")

cat("Building ntRNA graph...\n")

W_nt <- build_knn_graph(emb = Emb_nt, k = 30, symmetrize = "mean")

cat("Building spatial graph...\n")

batch_ids <- unique(as.character(states_seu$type))
W_spatial_list <- vector("list", length(batch_ids))
names(W_spatial_list) <- batch_ids

for (sample_id in batch_ids) {
  cat("Spatial graph:", sample_id, "\n")

  cells_use <- colnames(states_seu)[states_seu$type == sample_id]
  coords <- states_seu@meta.data[cells_use, c("column", "row")]
  coords <- as.matrix(coords)
  rownames(coords) <- cells_use

  W_spatial_list[[sample_id]] <- build_knn_graph(emb = coords, k = 30)
}

W_spatial <- as(bdiag(W_spatial_list), "dgCMatrix")

stopifnot(nrow(W_spatial) == ncol(states_seu))
stopifnot(nrow(W_rb) == ncol(states_seu))
stopifnot(nrow(W_nt) == ncol(states_seu))

cat("Running SCML...\n")

scml_result <- run_scml_core(
  graph_list = list(
    ntRNA = W_nt,
    rbRNA = W_rb,
    spatial = W_spatial
  ),
  ndim = 30,
  default_alpha = 1,
  spatial_alpha = 0.0
)

SCML_embedding <- scml_result$U
rownames(SCML_embedding) <- colnames(states_seu)

states_seu[["scml"]] <- CreateDimReducObject(
  embeddings = SCML_embedding,
  key = "SCML_",
  assay = "RNA"
)

scml_data <- preprocessed
scml_data$states_seu <- states_seu
scml_data$scml_result <- scml_result
scml_data$SCML_embedding <- SCML_embedding
scml_data$graphs <- list(ntRNA = W_nt, rbRNA = W_rb, spatial = W_spatial)

saveRDS(scml_data, file.path(output_dir, "02_scml_states.rds"))

cat("Saved outputs/02_scml_states.rds\n")
