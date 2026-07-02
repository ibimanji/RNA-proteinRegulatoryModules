# Multi-sample normalization and Harmony preprocessing
rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(anndata)
  library(Seurat)
  library(harmony)
  library(ggplot2)
})

# =========================
# 0. Parameters
# =========================
h5ad_path <- "/media/zenglab/result/lingyuan/STEM/imputation/mousebrain_withannotation_WT.h5ad"
output_root <- "/media/zenglab/result/lingyuan/STEM/imputation/magic_stagate_like_pipeline"
out_dir <- file.path(output_root, "00_harmony")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

spatial_coord_cols <- c("column", "row", "z")
batch_var <- "protocol-replicate"
batch_var_seurat <- make.names(batch_var)
nfeatures <- 1000
npcs <- 50
harmony_dims <- 1:30
umap_dims <- 1:25

out_rds <- file.path(out_dir, "00_harmony_preprocessed.rds")
out_rdata <- file.path(out_dir, "00_harmony_preprocessed.RData")

# =========================
# 1. Read h5ad and build Seurat object
# =========================
states_h5ad <- read_h5ad(h5ad_path)
required_layers <- c("totalRNA_raw", "rbRNA_raw", "ntRNA_raw")
if (!all(required_layers %in% names(states_h5ad$layers))) {
  stop(sprintf(
    "The h5ad object must contain layers: %s.",
    paste(required_layers, collapse = ", ")
  ))
}
if (!all(spatial_coord_cols %in% colnames(states_h5ad$obs))) {
  stop(sprintf("The h5ad obs must contain spatial coordinate columns: %s.",
               paste(spatial_coord_cols, collapse = ", ")))
}

meta <- states_h5ad$obs
if (!batch_var %in% colnames(meta)) {
  stop(sprintf("batch_var '%s' is not present in h5ad obs.", batch_var))
}
message("Harmony batch variable: ", batch_var)
meta[[batch_var_seurat]] <- meta[[batch_var]]

total_raw_mat <- t(states_h5ad$layers[["totalRNA_raw"]])
rb_raw_mat <- t(states_h5ad$layers[["rbRNA_raw"]])
nt_raw_mat <- t(states_h5ad$layers[["ntRNA_raw"]])

states <- CreateSeuratObject(
  counts = total_raw_mat,
  meta.data = meta,
  assay = "RNA",
  project = "STEM_harmony"
)
states[["rbRNA"]] <- CreateAssayObject(counts = rb_raw_mat)
states[["ntRNA"]] <- CreateAssayObject(counts = nt_raw_mat)
if (!batch_var_seurat %in% colnames(states@meta.data)) {
  stop(sprintf("Harmony batch variable '%s' was not transferred to Seurat metadata.", batch_var_seurat))
}

spatial_mat <- as.matrix(meta[, spatial_coord_cols, drop = FALSE])
storage.mode(spatial_mat) <- "numeric"
rownames(spatial_mat) <- rownames(meta)
states@misc$spatial_coords <- spatial_mat[Cells(states), , drop = FALSE]
states@misc$harmony_batch_var <- batch_var_seurat

if ("TE" %in% names(states_h5ad$layers)) {
  states[["TE_raw"]] <- CreateAssayObject(data = as(t(states_h5ad$layers[["TE"]]), "dgCMatrix"))
}

# =========================
# 2. Normalize and PCA for each assay
# =========================
rna_counts <- GetAssayData(states, assay = "RNA", layer = "counts")
rna_lib_sizes <- Matrix::colSums(rna_counts)

DefaultAssay(states) <- "RNA"
states <- NormalizeData(states, normalization.method = "RC", scale.factor = median(rna_lib_sizes), verbose = FALSE)
states <- FindVariableFeatures(states, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
states <- ScaleData(states, verbose = FALSE)
states <- RunPCA(states, npcs = npcs, reduction.name = "rna.pca", reduction.key = "rnaPC_", verbose = FALSE)

rb_counts <- GetAssayData(states, assay = "rbRNA", layer = "counts")
scale_factor_ref_rb <- median(Matrix::colSums(rb_counts))
rb_norm <- sweep(rb_counts, 2, rna_lib_sizes, "/") * scale_factor_ref_rb
states <- SetAssayData(states, assay = "rbRNA", layer = "data", new.data = as(rb_norm, "dgCMatrix"))
states <- FindVariableFeatures(states, selection.method = "vst", nfeatures = nfeatures, assay = "rbRNA", verbose = FALSE)
states <- ScaleData(states, assay = "rbRNA", verbose = FALSE)
states <- RunPCA(states, npcs = npcs, assay = "rbRNA", reduction.name = "rbrna.pca", reduction.key = "rbrnaPC_", verbose = FALSE)

nt_counts <- GetAssayData(states, assay = "ntRNA", layer = "counts")
scale_factor_ref_nt <- median(Matrix::colSums(nt_counts))
nt_norm <- sweep(nt_counts, 2, rna_lib_sizes, "/") * scale_factor_ref_nt
states <- SetAssayData(states, assay = "ntRNA", layer = "data", new.data = as(nt_norm, "dgCMatrix"))
states <- FindVariableFeatures(states, selection.method = "vst", nfeatures = nfeatures, assay = "ntRNA", verbose = FALSE)
states <- ScaleData(states, assay = "ntRNA", verbose = FALSE)
states <- RunPCA(states, npcs = npcs, assay = "ntRNA", reduction.name = "ntrna.pca", reduction.key = "ntrnaPC_", verbose = FALSE)

# =========================
# 3. Harmony correction
# =========================
states <- RunHarmony(
  object = states,
  group.by.vars = batch_var_seurat,
  reduction.use = "rna.pca",
  dims.use = harmony_dims,
  reduction.save = "rna.harmony",
  verbose = TRUE
)
states <- RunHarmony(
  object = states,
  group.by.vars = batch_var_seurat,
  reduction.use = "rbrna.pca",
  dims.use = harmony_dims,
  reduction.save = "rbrna.harmony",
  verbose = TRUE
)
states <- RunHarmony(
  object = states,
  group.by.vars = batch_var_seurat,
  reduction.use = "ntrna.pca",
  dims.use = harmony_dims,
  reduction.save = "ntrna.harmony",
  verbose = TRUE
)

# QC UMAPs for checking batch mixing before WNN.
states <- RunUMAP(
  states,
  reduction = "rna.harmony",
  dims = umap_dims,
  reduction.name = "rna.harmony.umap",
  reduction.key = "rnaHarmonyUMAP_",
  n.neighbors = 50,
  min.dist = 0.001,
  spread = 3,
  seed.use = 42,
  verbose = FALSE
)
p_batch <- DimPlot(states, reduction = "rna.harmony.umap", group.by = batch_var_seurat) +
  ggtitle(paste("RNA Harmony UMAP by", batch_var))
ggsave(file.path(out_dir, "00_rna_harmony_umap_by_batch.png"), p_batch, width = 7, height = 6, dpi = 300)

save(states, file = out_rdata)
saveRDS(states, file = out_rds)
message("Done. Saved Harmony-preprocessed object to: ", out_rds)


harmony_umap_specs <- list(
  RNA = list(
    reduction = "rna.harmony",
    reduction.name = "rna.harmony.umap",
    reduction.key = "rnaHarmonyUMAP_"
  ),
  rbRNA = list(
    reduction = "rbrna.harmony",
    reduction.name = "rbrna.harmony.umap",
    reduction.key = "rbrnaHarmonyUMAP_"
  ),
  ntRNA = list(
    reduction = "ntrna.harmony",
    reduction.name = "ntrna.harmony.umap",
    reduction.key = "ntrnaHarmonyUMAP_"
  )
)

for (modality_name in names(harmony_umap_specs)) {
  spec <- harmony_umap_specs[[modality_name]]
  
  states <- RunUMAP(
    states,
    reduction = spec$reduction,
    dims = umap_dims,
    reduction.name = spec$reduction.name,
    reduction.key = spec$reduction.key,
    n.neighbors = 50,
    min.dist = 0.001,
    spread = 3,
    seed.use = 42,
    verbose = FALSE
  )
  
  for (group_col in c(batch_var_seurat)) {
    if (group_col %in% colnames(states@meta.data)) {
      p <- DimPlot(states, reduction = spec$reduction.name, group.by = group_col) +
        ggtitle(paste(modality_name, "Harmony UMAP by", batch_var))
      ggsave(
        file.path(
          out_dir,
          paste0("00_", tolower(modality_name), "_harmony_umap_by_", batch_var_seurat, ".png")
        ),
        p,
        width = 7,
        height = 6,
        dpi = 300
      )
    }
  }
}
