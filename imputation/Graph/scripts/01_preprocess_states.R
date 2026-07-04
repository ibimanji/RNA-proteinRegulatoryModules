setwd("/media/zenglab/result/lingyuan/STEM/imputation/graph")
source("R/scml_core.R")
source("R/setup.R")

cat("Loading STATES h5ad...\n")

states <- read_h5ad(states_h5ad_path)
rb_raw_mat <- t(states$layers[["rbRNA_raw"]])
total_raw_mat <- t(states$layers[["totalRNA_raw"]])
nt_raw_mat <- t(states$layers[["ntRNA_raw"]])
meta <- states$obs

states_seu <- CreateSeuratObject(
  counts = total_raw_mat,
  meta.data = meta,
  project = "statesmouse"
)
states_seu[["rbRNA"]] <- CreateAssayObject(counts = rb_raw_mat)
states_seu[["ntRNA"]] <- CreateAssayObject(counts = nt_raw_mat)
states_seu <- subset(states_seu, subset = type == "C1")

cat("Running RNA normalization...\n")

DefaultAssay(states_seu) <- "RNA"
rna_counts <- get_assay_matrix(states_seu, assay = "RNA", layer = "counts")
rna_lib_sizes <- colSums(rna_counts)

states_seu <- NormalizeData(
  states_seu,
  normalization.method = "RC",
  scale.factor = median(rna_lib_sizes)
)
states_seu <- FindVariableFeatures(states_seu, selection.method = "vst", nfeatures = 2000)
states_seu <- ScaleData(states_seu)
states_seu <- RunPCA(
  states_seu,
  npcs = 50,
  reduction.name = "rna.pca",
  reduction.key = "rnaPC_"
)

cat("Running rbRNA normalization...\n")

rb_counts <- get_assay_matrix(states_seu, assay = "rbRNA", layer = "counts")
scale_factor_ref_rb <- median(colSums(rb_counts))
rb_norm <- sweep(rb_counts, 2, rna_lib_sizes, "/") * scale_factor_ref_rb
states_seu <- set_assay_matrix(
  states_seu,
  assay = "rbRNA",
  layer = "data",
  value = as(rb_norm, "dgCMatrix")
)
states_seu <- FindVariableFeatures(
  states_seu,
  assay = "rbRNA",
  selection.method = "vst",
  nfeatures = 2000
)
states_seu <- ScaleData(states_seu, assay = "rbRNA")
states_seu <- RunPCA(
  states_seu,
  assay = "rbRNA",
  npcs = 50,
  reduction.name = "rbrna.pca",
  reduction.key = "rbrnaPC_"
)

cat("Running ntRNA normalization...\n")

nt_counts <- get_assay_matrix(states_seu, assay = "ntRNA", layer = "counts")
scale_factor_ref_nt <- median(colSums(nt_counts))
nt_norm <- sweep(nt_counts, 2, rna_lib_sizes, "/") * scale_factor_ref_nt
states_seu <- set_assay_matrix(
  states_seu,
  assay = "ntRNA",
  layer = "data",
  value = as(nt_norm, "dgCMatrix")
)
states_seu <- FindVariableFeatures(
  states_seu,
  assay = "ntRNA",
  selection.method = "vst",
  nfeatures = 2000
)
states_seu <- ScaleData(states_seu, assay = "ntRNA")
states_seu <- RunPCA(
  states_seu,
  assay = "ntRNA",
  npcs = 50,
  reduction.name = "ntrna.pca",
  reduction.key = "ntrnaPC_"
)

cat("Running Harmony on rbRNA...\n")

states_seu <- RunHarmony(
  object = states_seu,
  group.by.vars = "type",
  reduction.use = "rbrna.pca",
  reduction.save = "rbrna.harmony",
  verbose = TRUE
)

cat("Running Harmony on ntRNA...\n")

states_seu <- RunHarmony(
  object = states_seu,
  group.by.vars = "type",
  reduction.use = "ntrna.pca",
  reduction.save = "ntrna.harmony",
  verbose = TRUE
)

saveRDS(
  list(
    states_seu = states_seu,
    rb_norm = rb_norm,
    nt_norm = nt_norm,
    rna_lib_sizes = rna_lib_sizes
  ),
  file.path(output_dir, "01_preprocessed_states.rds")
)

cat("Saved outputs/01_preprocessed_states.rds\n")

