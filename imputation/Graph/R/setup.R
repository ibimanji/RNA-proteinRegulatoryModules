gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(Seurat)
  library(FNN)
  library(RSpectra)
  library(ggplot2)
  library(ggrepel)
  library(harmony)
})

if (!exists("read_h5ad") && requireNamespace("anndata", quietly = TRUE)) {
  read_h5ad <- anndata::read_h5ad
}

project_root <- normalizePath(getwd(), mustWork = FALSE)

output_dir <- file.path(project_root, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

states_h5ad_path <- Sys.getenv(
  "STATES_H5AD",
  "/media/zenglab/result/lingyuan/STEM/imputation/mousebrain_withannotation_WT.h5ad"
)
