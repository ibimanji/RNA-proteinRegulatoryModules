#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/run_knn_smooth_from_rds.R --type TYPE --k K [options]\n\n",
    "Required:\n",
    "  --type VALUE       Value in metadata column 'type' to keep\n",
    "  --k INT            Number of neighbors for knn_smooth.py\n\n",
    "Options:\n",
    "  --rds PATH         Input Seurat RDS [states_with_plaque_info.rds]\n",
    "  --meta-col NAME    Metadata column used for filtering [type]\n",
    "  --assay NAME       Assay to extract [totalRNA]\n",
    "  --slot NAME        Slot/layer to extract [counts]\n",
    "  --d INT            Number of PCs for knn_smooth.py [10]\n",
    "  --dither FLOAT     Dither for knn_smooth.py [0.03]\n",
    "  --seed INT         Random seed [0]\n",
    "  --python PATH      Python executable [python3]\n",
    "  --knn-script PATH  Path to knn_smooth.py [knn_smooth.py]\n",
    "  --outdir PATH      Output directory [knn_smooth_outputs]\n",
    "  --help             Show this help\n",
    sep = ""
  )
}

parse_args <- function(args) {
  cfg <- list(
    rds = "states_with_plaque_info.rds",
    type = NULL,
    meta_col = "type",
    assay = "totalRNA",
    slot = "counts",
    k = NULL,
    d = "10",
    dither = "0.03",
    seed = "0",
    python = "python3",
    knn_script = "knn_smooth.py",
    outdir = "knn_smooth_outputs"
  )

  key_map <- c(
    "--rds" = "rds",
    "--type" = "type",
    "--meta-col" = "meta_col",
    "--assay" = "assay",
    "--slot" = "slot",
    "--k" = "k",
    "--d" = "d",
    "--dither" = "dither",
    "--seed" = "seed",
    "--python" = "python",
    "--knn-script" = "knn_script",
    "--outdir" = "outdir"
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--help", "-h")) {
      usage()
      quit(status = 0)
    }
    if (!arg %in% names(key_map)) {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
    if (i == length(args)) {
      stop("Missing value for argument: ", arg, call. = FALSE)
    }
    cfg[[key_map[[arg]]]] <- args[[i + 1]]
    i <- i + 2
  }

  if (is.null(cfg$type)) stop("--type is required.", call. = FALSE)
  if (is.null(cfg$k)) stop("--k is required.", call. = FALSE)
  cfg
}

get_assay_matrix <- function(object, assay, slot_or_layer) {
  seurat_object_ok <- tryCatch({
    requireNamespace("SeuratObject", quietly = TRUE)
  }, error = function(e) {
    FALSE
  })

  if (isTRUE(seurat_object_ok)) {
    if (utils::packageVersion("SeuratObject") >= "5.0.0") {
      return(SeuratObject::GetAssayData(object, assay = assay, layer = slot_or_layer))
    }
    return(SeuratObject::GetAssayData(object, assay = assay, slot = slot_or_layer))
  }

  warning(
    "Could not load SeuratObject; trying direct slot extraction. ",
    "If this fails, fix the R SeuratObject/Matrix installation."
  )
  assay_obj <- object@assays[[assay]]
  if (is.null(assay_obj)) stop("Assay not found: ", assay, call. = FALSE)

  if (slot_or_layer %in% methods::slotNames(assay_obj)) {
    return(methods::slot(assay_obj, slot_or_layer))
  }
  if ("layers" %in% methods::slotNames(assay_obj)) {
    layers <- methods::slot(assay_obj, "layers")
    if (slot_or_layer %in% names(layers)) {
      return(layers[[slot_or_layer]])
    }
  }

  stop("Could not find slot/layer '", slot_or_layer, "' in assay '", assay, "'.", call. = FALSE)
}

sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

cfg <- parse_args(args)

if (!file.exists(cfg$rds)) stop("Input RDS does not exist: ", cfg$rds, call. = FALSE)
if (!file.exists(cfg$knn_script)) stop("knn_smooth.py does not exist: ", cfg$knn_script, call. = FALSE)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)

message("Reading Seurat RDS: ", cfg$rds)
object <- readRDS(cfg$rds)
if (!inherits(object, "Seurat")) stop("Input RDS is not a Seurat object.", call. = FALSE)

meta <- object@meta.data
if (!cfg$meta_col %in% colnames(meta)) {
  stop("Metadata column not found: ", cfg$meta_col, call. = FALSE)
}

keep_cells <- rownames(meta)[as.character(meta[[cfg$meta_col]]) == cfg$type]
if (length(keep_cells) == 0) {
  available <- paste(sort(unique(as.character(meta[[cfg$meta_col]]))), collapse = ", ")
  stop("No cells matched ", cfg$meta_col, " == '", cfg$type, "'. Available values: ", available, call. = FALSE)
}

message("Matched cells: ", length(keep_cells), " where ", cfg$meta_col, " == ", cfg$type)
message("Extracting assay=", cfg$assay, ", slot/layer=", cfg$slot)
counts <- get_assay_matrix(object, assay = cfg$assay, slot_or_layer = cfg$slot)

missing_cells <- setdiff(keep_cells, colnames(counts))
if (length(missing_cells) > 0) {
  stop("Some matched cells are not present in the assay matrix. First missing cell: ", missing_cells[[1]], call. = FALSE)
}

counts <- counts[, keep_cells, drop = FALSE]
message("Export matrix shape: ", nrow(counts), " genes x ", ncol(counts), " cells")

prefix <- paste0(
  sanitize_filename(tools::file_path_sans_ext(basename(cfg$rds))),
  "_",
  sanitize_filename(cfg$meta_col),
  "_",
  sanitize_filename(cfg$type),
  "_",
  sanitize_filename(cfg$assay),
  "_",
  sanitize_filename(cfg$slot)
)
input_tsv <- file.path(cfg$outdir, paste0(prefix, ".input.tsv"))
output_tsv <- file.path(cfg$outdir, paste0(prefix, ".knn_k", cfg$k, "_d", cfg$d, ".tsv"))

message("Writing dense TSV for knn_smooth.py: ", input_tsv)
counts_dense <- as.matrix(counts)
storage.mode(counts_dense) <- "double"
utils::write.table(
  counts_dense,
  file = input_tsv,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

cmd_args <- c(
  cfg$knn_script,
  "-f", input_tsv,
  "-o", output_tsv,
  "-k", cfg$k,
  "-d", cfg$d,
  "--dither", cfg$dither,
  "-s", cfg$seed
)

message("Running: ", paste(c(cfg$python, cmd_args), collapse = " "))
status <- system2(cfg$python, args = cmd_args)
if (!identical(status, 0L)) {
  stop("knn_smooth.py failed with exit status: ", status, call. = FALSE)
}

message("Done.")
message("Input TSV:  ", input_tsv)
message("Output TSV: ", output_tsv)
