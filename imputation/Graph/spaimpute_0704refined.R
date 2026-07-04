rm(list = ls()); gc(); options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(Seurat)
  library(FNN)
  library(RSpectra)
  library(ggplot2)
  library(ggrepel)
  library(harmony)
  library(anndata)
})

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

assert_named_list <- function(x, arg = "x") {
  if (!is.list(x) || is.null(names(x)) || any(names(x) == "")) {
    stop(arg, " must be a named list.", call. = FALSE)
  }
  invisible(TRUE)
}

row_normalize <- function(x, eps = 1e-10) {
  denom <- sqrt(rowSums(x^2))
  denom[denom < eps] <- 1 #maintaining the x while the weight is relatively small
  sweep(x, 1, denom, "/")
}

get_assay_matrix <- function(seurat_obj, assay, layer = "counts") {
  tryCatch(
    GetAssayData(seurat_obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(seurat_obj, assay = assay, slot = layer)
  )
}

set_assay_matrix <- function(seurat_obj, assay, layer = "data", value) {
  tryCatch(
    SetAssayData(seurat_obj, assay = assay, layer = layer, new.data = value),
    error = function(e) SetAssayData(seurat_obj, assay = assay, slot = layer, new.data = value)
  )
}

validate_embedding_list <- function(embeddings) {
  assert_named_list(embeddings, "embeddings")
  n_cells <- vapply(embeddings, nrow, integer(1))
  if (length(unique(n_cells)) != 1) {
    stop("All embedding layers must have the same number of rows/cells.", call. = FALSE)
  }
  
  row_ids <- lapply(embeddings, rownames)
  has_all_ids <- all(vapply(row_ids, function(z) !is.null(z), logical(1)))
  if (has_all_ids) {
    ref <- row_ids[[1]]
    ok <- vapply(row_ids, identical, logical(1), y = ref)
    if (!all(ok)) {
      stop("All embedding layers must use the same rownames in the same order.", call. = FALSE)
    }
  }
  
  invisible(TRUE)
}

build_knn_graph <- function(emb, k = 30, symmetrize = c("max", "mean"), sigma_floor = 1e-6) {
  symmetrize <- match.arg(symmetrize)
  emb <- as.matrix(emb)
  if (nrow(emb) <= k) {
    stop("k must be smaller than the number of cells/vertices.", call. = FALSE)
  }
  if (any(!is.finite(emb))) {
    stop("Embedding contains NA/NaN/Inf values.", call. = FALSE)
  }
  
  knn <- FNN::get.knn(emb, k = k)
  idx <- knn$nn.index
  dst <- knn$nn.dist
  n <- nrow(idx)
  
  sigma <- dst[, k]#here kth in the k-NN's K is borrowed as sigma, we can change this but this function is fast.
  sigma[sigma < sigma_floor] <- sigma_floor
  
  ii <- rep(seq_len(n), each = k)
  jj <- as.vector(t(idx))
  d <- as.vector(t(dst))
  
  sigma_i <- rep(sigma, each = k)
  sigma_j <- sigma[jj]
  w <- exp(-(d^2) / (sigma_i * sigma_j))
  
  W <- sparseMatrix(i = ii, j = jj, x = w, dims = c(n, n))
  W <- switch(
    symmetrize,
    max = pmax(W, t(W)),
    mean = (W + t(W)) / 2
  )
  diag(W) <- 0
  drop0(W)
}

compute_normalized_laplacian <- function(W, eps = 1e-10) {
  W <- as(W, "dgCMatrix")
  deg <- Matrix::rowSums(W)
  Dinv <- Diagonal(x = 1 / sqrt(deg + eps))
  Diagonal(nrow(W)) - Dinv %*% W %*% Dinv
}

get_laplacian_subspace <- function(L, ndim = 30) {
  if (ndim >= nrow(L)) {
    stop("ndim must be smaller than the number of cells/vertices.", call. = FALSE)
  }
  
  eig <- RSpectra::eigs_sym(L, k = ndim, which = "SM")
  U <- qr.Q(qr(eig$vectors))
  list(U = U, values = eig$values)
}

resolve_layer_values <- function(values, layer_names, default_value, spatial_value = NULL, arg = "values") {
  #this is for layer weight, especially for spatial layer
  if (is.null(values)) {
    out <- rep(default_value, length(layer_names))
    names(out) <- layer_names
  } else if (length(values) == 1 && is.null(names(values))) {
    out <- rep(as.numeric(values), length(layer_names))
    names(out) <- layer_names
  } else {
    out <- rep(default_value, length(layer_names))
    names(out) <- layer_names
    values <- unlist(values)
    if (is.null(names(values)) || any(names(values) == "")) {
      if (length(values) != length(layer_names)) {
        stop(arg, " must be length 1, length equal to layers, or named by layer.", call. = FALSE)
      }
      out[] <- as.numeric(values)
    } else {
      missing_names <- setdiff(names(values), layer_names)
      if (length(missing_names) > 0) {
        stop(arg, " contains unknown layer(s): ", paste(missing_names, collapse = ", "), call. = FALSE)
      }
      out[names(values)] <- as.numeric(values)
    }
  }
  
  if (!is.null(spatial_value)) {
    spatial_layers <- grepl("spatial|space|coord|location", layer_names, ignore.case = TRUE)
    implicit_values <- is.null(values) || (length(values) == 1 && is.null(names(values)))
    if (any(spatial_layers) && implicit_values) {
      out[spatial_layers] <- spatial_value
    }
  }
  
  if (any(!is.finite(out))) {
    stop(arg, " must contain finite numeric values.", call. = FALSE)
  }
  out
}
build_delaunay_graph <- function(
    coords,
    sigma_floor = 1e-6
){
  
  if(!requireNamespace("geometry", quietly = TRUE)){
    stop("Please install package 'geometry'.")
  }
  
  coords <- as.matrix(coords)
  
  if(ncol(coords) != 2){
    stop("coords must be an n × 2 matrix.")
  }
  
  if(any(!is.finite(coords))){
    stop("coords contains NA/NaN/Inf.")
  }
  
  n <- nrow(coords)
  
  ########################################################
  ## Delaunay triangulation
  ########################################################
  
  tri <- geometry::delaunayn(coords)
  
  edges <- rbind(
    tri[, c(1,2)],
    tri[, c(1,3)],
    tri[, c(2,3)]
  )
  
  edges <- t(apply(edges, 1, sort))
  edges <- unique(edges)
  
  ########################################################
  ## Euclidean distance
  ########################################################
  
  d <- sqrt(
    rowSums(
      (
        coords[edges[,1], , drop = FALSE] -
          coords[edges[,2], , drop = FALSE]
      )^2
    )
  )
  
  ########################################################
  ## Local sigma
  ########################################################
  
  sigma <- numeric(n)
  
  for(i in seq_len(n)){
    
    idx <- which(edges[,1] == i | edges[,2] == i)
    
    if(length(idx) == 0){
      
      sigma[i] <- sigma_floor
      
    }else{
      
      sigma[i] <- median(d[idx])
      
    }
    
  }
  
  sigma[sigma < sigma_floor] <- sigma_floor
  
  ########################################################
  ## Adaptive Gaussian weight
  ########################################################
  
  sigma_i <- sigma[edges[,1]]
  sigma_j <- sigma[edges[,2]]
  
  w <- exp(
    -(d^2) /
      (sigma_i * sigma_j)
  )
  
  ########################################################
  ## Sparse matrix
  ########################################################
  
  W <- Matrix::sparseMatrix(
    i = edges[,1],
    j = edges[,2],
    x = w,
    dims = c(n, n)
  )
  
  W <- pmax(W, Matrix::t(W))
  
  diag(W) <- 0
  
  Matrix::drop0(W)
}
build_radius_graph <- function(
    coords,
    radius,
    sigma_floor = 1e-6
){
  
  coords <- as.matrix(coords)
  
  if(ncol(coords) != 2){
    stop("coords must be an n × 2 matrix.")
  }
  
  if(any(!is.finite(coords))){
    stop("coords contains NA/NaN/Inf.")
  }
  
  n <- nrow(coords)
  
  ########################################################
  ## Pairwise distance
  ########################################################
  
  D <- as.matrix(dist(coords))
  
  ########################################################
  ## Radius neighbors
  ########################################################
  
  idx <- which(
    D <= radius &
      D > 0,
    arr.ind = TRUE
  )
  
  ## 保留上三角
  idx <- idx[idx[,1] < idx[,2], , drop = FALSE]
  
  if(nrow(idx) == 0){
    stop("No edges found. Increase radius.")
  }
  
  d <- D[idx]
  
  ########################################################
  ## Local sigma
  ########################################################
  
  sigma <- numeric(n)
  
  for(i in seq_len(n)){
    
    di <- D[i, ]
    
    di <- di[
      di > 0 &
        di <= radius
    ]
    
    if(length(di) == 0){
      
      sigma[i] <- sigma_floor
      
    }else{
      
      sigma[i] <- median(di)
      
    }
    
  }
  
  sigma[sigma < sigma_floor] <- sigma_floor
  
  ########################################################
  ## Adaptive Gaussian weight
  ########################################################
  
  sigma_i <- sigma[idx[,1]]
  sigma_j <- sigma[idx[,2]]
  
  w <- exp(
    -(d^2) /
      (sigma_i * sigma_j)
  )
  
  ########################################################
  ## Sparse adjacency
  ########################################################
  
  W <- Matrix::sparseMatrix(
    i = idx[,1],
    j = idx[,2],
    x = w,
    dims = c(n, n)
  )
  
  W <- pmax(W, Matrix::t(W))
  
  diag(W) <- 0
  
  Matrix::drop0(W)
}
run_scml_core <- function(
    graph_list,
    ndim = 30,
    alpha = NULL,
    layer_weight = NULL,
    default_alpha = 1.0,
    spatial_alpha = 0,
    return_laplacians = TRUE) {
  assert_named_list(graph_list, "graph_list")
  
  layer_names <- names(graph_list)
  n_vertices <- vapply(graph_list, nrow, integer(1))
  if (length(unique(n_vertices)) != 1 || any(vapply(graph_list, ncol, integer(1)) != n_vertices[1])) {
    stop("All graphs must be square matrices with the same dimensions.", call. = FALSE)
  }
  
  alpha_vec <- resolve_layer_values(
    alpha,
    layer_names,
    default_value = default_alpha,
    spatial_value = spatial_alpha,
    arg = "alpha"
  )
  weight_vec <- resolve_layer_values(
    layer_weight,
    layer_names,
    default_value = 1,
    spatial_value = NULL,
    arg = "layer_weight"
  )
  cat("\n========== Layer Parameters ==========\n")
  
  print(
    data.frame(
      layer = layer_names,
      weight = as.numeric(weight_vec[layer_names]),
      alpha = as.numeric(alpha_vec[layer_names]),
      row.names = NULL
    )
  )
  
  n <- n_vertices[1]
  L_sum <- Matrix(0, n, n, sparse = TRUE)
  P_sum <- Matrix(0, n, n, sparse = TRUE)
  
  U_layers <- list()
  eig_layers <- list()
  laplacians <- list()
  
  for (layer in layer_names) {
    L <- compute_normalized_laplacian(graph_list[[layer]])
    subspace <- get_laplacian_subspace(L, ndim = ndim)
    
    U_layers[[layer]] <- subspace$U
    eig_layers[[layer]] <- subspace$values
    if (return_laplacians) {
      laplacians[[layer]] <- L
    }
    
    L_sum <- L_sum + weight_vec[[layer]] * L
    P_sum <- P_sum + alpha_vec[[layer]] * tcrossprod(subspace$U)
  }
  
  L_modified <- L_sum - P_sum
  eig_consensus <- RSpectra::eigs_sym(L_modified, k = ndim, which = "SM")
  U_consensus <- row_normalize(eig_consensus$vectors)
  colnames(U_consensus) <- paste0("SCML_", seq_len(ncol(U_consensus)))
  
  list(
    U = U_consensus,
    eigenvalues = eig_consensus$values,
    U_layers = U_layers,
    eigenvalues_layers = eig_layers,
    L_modified = L_modified,
    laplacians = if (return_laplacians) laplacians else NULL,
    alpha = alpha_vec,
    layer_weight = weight_vec
  )
}

run_scml_from_embeddings <- function(
    embeddings,
    k = 30,
    ndim = 30,
    alpha = NULL,
    layer_weight = NULL,
    default_alpha = 1,
    spatial_alpha = 0,
    symmetrize = "max",
    return_graphs = TRUE) {
  validate_embedding_list(embeddings)
  
  graph_list <- lapply(
    embeddings,
    build_knn_graph,
    k = k,
    symmetrize = symmetrize
  )
  
  result <- run_scml_core(
    graph_list = graph_list,
    ndim = ndim,
    alpha = alpha,
    layer_weight = layer_weight,
    default_alpha = default_alpha,
    spatial_alpha = spatial_alpha
  )
  
  if (!is.null(rownames(embeddings[[1]]))) {
    rownames(result$U) <- rownames(embeddings[[1]])
  }
  result$graphs <- if (return_graphs) graph_list else NULL
  result$embeddings_input <- embeddings
  result
}
GrassmannDistance <- function(U1, U2){
  
  U1 <- qr.Q(qr(U1))
  U2 <- qr.Q(qr(U2))
  
  k <- ncol(U1)
  
  d <- sqrt(
    max(
      0,
      k - sum((crossprod(U1, U2))^2)
    )
  )
  
  return(d)
}
ComputeGrassmannDistance <- function(U_consensus, U_layers){
  
  emb <- c(list(consensus = U_consensus), U_layers)
  
  nm <- names(emb)
  n <- length(emb)
  
  D <- matrix(
    0,
    n,
    n,
    dimnames = list(nm, nm)
  )
  
  for(i in seq_len(n)){
    for(j in i:n){
      
      d <- GrassmannDistance(
        emb[[i]],
        emb[[j]]
      )
      
      D[i, j] <- d
      D[j, i] <- d
    }
  }
  
  D
}
################################
##1. STATES C1 Sample        m##
################################
states <- read_h5ad("/media/zenglab/data/keyao/20250224tissue_h5adtst/250228_output/0228-tissue-harmony.h5ad")#改了改了改了!
rb_raw_mat <- t(states$layers[["ribomap"]]) 
total_raw_mat <- t(states$layers[["raw"]])
nt_raw_mat <- total_raw_mat - rb_raw_mat
colnames(nt_raw_mat) = colnames(total_raw_mat)
rownames(nt_raw_mat) = rownames(total_raw_mat)
all.equal(rownames(nt_raw_mat), rownames(rb_raw_mat), rownames(total_raw_mat))
all.equal(colnames(nt_raw_mat), colnames(rb_raw_mat), colnames(total_raw_mat))
load("/media/zenglab/result/lingxuan/Github/01_mouse_brain_Downstream/states_celltypes_identification.RData")
meta <- states@meta.data
meta$row <- max(meta$row) + min(meta$row) - meta$row
meta$column <- max(meta$column) + min(meta$column) - meta$column

states_seu <- CreateSeuratObject(
  counts = total_raw_mat,
  meta.data = meta,
  project = "statesmouse"
)
states_seu[["rbRNA"]] <- CreateAssayObject(
  counts = rb_raw_mat
)
states_seu[["ntRNA"]] <- CreateAssayObject(
  counts = nt_raw_mat
)

states_seu <- subset(states_seu, subset = replicate == "C1")
states_seu

cat("Running RNA normalization...\n")

DefaultAssay(states_seu) <- "RNA"

rna_counts <- GetAssayData(
  states_seu,
  assay = "RNA",
  layer = "counts"
)

rna_lib_sizes <- colSums(rna_counts)

states_seu <- NormalizeData(
  states_seu,
  normalization.method = "RC",
  scale.factor = median(rna_lib_sizes)
)

states_seu <- FindVariableFeatures(
  states_seu,
  selection.method = "vst",
  nfeatures = 2000
)

states_seu <- ScaleData(states_seu)

states_seu <- RunPCA(
  states_seu,
  npcs = 50,
  reduction.name = "rna.pca",
  reduction.key = "rnaPC_"
)

############################################################
## 2. rbRNA normalization
############################################################

cat("Running rbRNA normalization...\n")

rb_counts <- GetAssayData(
  states_seu,
  assay = "rbRNA",
  layer = "counts"
)

scale_factor_ref_rb <- median(
  colSums(rb_counts)
)

rb_norm <- sweep(
  rb_counts,
  2,
  rna_lib_sizes,
  "/"
) * scale_factor_ref_rb

states_seu <- SetAssayData(
  states_seu,
  assay = "rbRNA",
  layer = "data",
  new.data = as(rb_norm, "dgCMatrix")
)

states_seu <- FindVariableFeatures(
  states_seu,
  assay = "rbRNA",
  selection.method = "vst",
  nfeatures = 2000
)

states_seu <- ScaleData(
  states_seu,
  assay = "rbRNA"
)

states_seu <- RunPCA(
  states_seu,
  assay = "rbRNA",
  npcs = 50,
  reduction.name = "rbrna.pca",
  reduction.key = "rbrnaPC_"
)

############################################################
## 3. ntRNA normalization
############################################################

cat("Running ntRNA normalization...\n")

nt_counts <- GetAssayData(
  states_seu,
  assay = "ntRNA",
  layer = "counts"
)

scale_factor_ref_nt <- median(
  colSums(nt_counts)
)

nt_norm <- sweep(
  nt_counts,
  2,
  rna_lib_sizes,
  "/"
) * scale_factor_ref_nt

states_seu <- SetAssayData(
  states_seu,
  assay = "ntRNA",
  layer = "data",
  new.data = as(nt_norm, "dgCMatrix")
)

states_seu <- FindVariableFeatures(
  states_seu,
  assay = "ntRNA",
  selection.method = "vst",
  nfeatures = 2000
)

states_seu <- ScaleData(
  states_seu,
  assay = "ntRNA"
)

states_seu <- RunPCA(
  states_seu,
  assay = "ntRNA",
  npcs = 50,
  reduction.name = "ntrna.pca",
  reduction.key = "ntrnaPC_"
)

############################################################
## 4. Harmony integration
############################################################

cat("Running Harmony on rbRNA...\n")

states_seu <- RunHarmony(
  object = states_seu,
  group.by.vars = "replicate",
  reduction.use = "rbrna.pca",
  reduction.save = "rbrna.harmony",
  verbose = TRUE
)

cat("Running Harmony on ntRNA...\n")

states_seu <- RunHarmony(
  object = states_seu,
  group.by.vars = "replicate",
  reduction.use = "ntrna.pca",
  reduction.save = "ntrna.harmony",
  verbose = TRUE
)

############################################################
## 5. Extract Harmony embeddings
############################################################

Emb_rb <- Embeddings(
  states_seu,
  reduction = "rbrna.harmony"
)[,1:30]

Emb_nt <- Embeddings(
  states_seu,
  reduction = "ntrna.harmony"
)[,1:30]

############################################################
## 6. Build expression graphs (global)
############################################################

cat("Building rbRNA graph...\n")

W_rb <- build_knn_graph(
  emb = Emb_rb,
  k = 30, symmetrize = 'max'
)

cat("Building ntRNA graph...\n")

W_nt <- build_knn_graph(
  emb = Emb_nt,
  k = 30, symmetrize = 'max'
)

############################################################
## 7. Build spatial graph
## sample-specific
############################################################

cat("Building spatial graph...\n")

batch_ids <- unique(
  as.character(states_seu$replicate)
)

W_spatial_list <- vector(
  "list",
  length(batch_ids)
)

names(W_spatial_list) <- batch_ids

for(sample_id in batch_ids){
  
  cat("Spatial graph:", sample_id, "\n")
  
  cells_use <- colnames(states_seu)[
    states_seu$replicate == sample_id
  ]
  
  coords <- states_seu@meta.data[
    cells_use,
    c("column","row")
  ]
  
  coords <- as.matrix(coords)
  
  rownames(coords) <- cells_use
  
  W_spatial_list[[sample_id]] <- build_knn_graph(
    emb = coords,
    k = 30,
    symmetrize = 'max'
  ) #!!!这里可以考虑更换构图方式或者改小k到5-10即可
}

############################################################
## block diagonal spatial graph
############################################################

W_spatial <- as(
  bdiag(W_spatial_list),
  "dgCMatrix"
)

############################################################
## sanity check
############################################################

stopifnot(
  nrow(W_spatial) == ncol(states_seu)
)

stopifnot(
  nrow(W_rb) == ncol(states_seu)
)

stopifnot(
  nrow(W_nt) == ncol(states_seu)
)

############################################################
## 8. Run SCML
############################################################

cat("Running SCML...\n")

scml_result <- run_scml_core(
  graph_list = list(
    ntRNA   = W_nt,
    rbRNA   = W_rb,
    spatial = W_spatial
  ),
  ndim = 30,
  default_alpha = 1,
  spatial_alpha = 0.0
)
gc()
############################################################
## 9. Store SCML embedding
############################################################

SCML_embedding <- scml_result$U
L_scml <- scml_result$L_modified 
rownames(SCML_embedding) <- colnames(states_seu)

states_seu[["scml"]] <- CreateDimReducObject(
  embeddings = SCML_embedding,
  key = "SCML_",
  assay = "RNA"
)

BuildSCMLGraph <- build_knn_graph
ComputeGraphSpectrum <- function(
    W,
    n_eigs = 300
){
  
  D <- Matrix::rowSums(W)
  
  D_inv <- Diagonal(
    x = 1/sqrt(D+1e-8)
  )
  
  L <- Diagonal(nrow(W)) -
    D_inv %*%
    W %*%
    D_inv
  
  eig <- eigs_sym(
    L,
    k=n_eigs,
    which="SM"
  )
  
  list(
    lambda=eig$values,
    V=eig$vectors
  )
}

W <- BuildSCMLGraph(SCML_embedding, k=30)

spec <- ComputeGraphSpectrum(
  W,
  n_eigs=200
)

lambda <- spec$lambda
V <- spec$V

GraphFilterKernel <- function(
    lambda,
    type = c(
      "heat",
      "tikhonov",
      "exponential",
      "wavelet",
      "adaptive_wavelet",
      "bandpass",
      "ideal",
      "butterworth",
      "hybrid"
    ),
    beta = 1,
    scale = 5,
    lambda_low = 0.1,
    lambda_high = 1,
    lambda_c = 0.4,
    order = 4,
    gamma = 0.25,
    peak_quantile = 0.15
){
  
  type <- match.arg(type)
  
  g <- switch(
    
    type,
    
    heat =
      exp(-beta * lambda),
    
    tikhonov =
      1 / (1 + beta * lambda),
    
    exponential =
      exp(-beta * sqrt(lambda)),
    
    wavelet =
      lambda * exp(-scale * lambda),
    
    adaptive_wavelet = {
      
      lambda2 <- sort(lambda)[2]
      
      lambda_peak <-
        quantile(
          lambda,
          peak_quantile
        )
      
      s <- 1 / lambda_peak
      
      pmax(
        lambda - lambda2,
        0
      ) *
        exp(
          -s * lambda
        )
      
    },
    
    bandpass =
      as.numeric(
        lambda > lambda_low &
          lambda < lambda_high
      ),
    
    ideal =
      as.numeric(
        lambda < beta
      ),
    
    butterworth =
      1 /
      (
        1 +
          (lambda / lambda_c)^(2 * order)
      ),
    
    hybrid =
      exp(-beta * lambda) +
      gamma *
      lambda *
      exp(-scale * lambda)
    
  )
  
  g <- g / max(g)
  
  g
}
lambda_sort <- sort(lambda)

g_heat <- GraphFilterKernel(
  lambda_sort,
  type = "heat",
  beta = 3
)

GraphWaveletDenoise <- function(
    expr,
    V,
    g
){
  
  g <- g
  
  X <- t(as.matrix(expr))
  
  X_hat <- t(V) %*% X
  
  X_hat_filt <- g * X_hat
  
  X_filt <- V %*% X_hat_filt
  
  X_filt[X_filt < 0] <- 0
  
  t(X_filt)
}

rb_raw_mat = GetAssayData(states_seu, layer = 'counts', assay = 'rbRNA')
total_raw_mat = GetAssayData(states_seu, layer  = 'counts', assay = 'RNA')
nt_raw_mat = GetAssayData(states_seu, layer = 'counts', assay = 'ntRNA')

total_norm <- GetAssayData(states_seu, layer = 'data',assay = 'RNA')
rb_norm <- GetAssayData(states_seu, layer = 'data',assay = 'rbRNA')
nt_norm <- GetAssayData(states_seu, layer = 'data',assay = 'ntRNA')

##rb filtered
rb_wavelet <- GraphWaveletDenoise(
  expr = rb_raw_mat,
  V = V,
  g = g_heat
)
colnames(rb_wavelet) <- colnames(rb_raw_mat)
all.equal(rownames(rb_wavelet),rownames(rb_raw_mat))

##nt filtered
nt_wavelet <- GraphWaveletDenoise(
  expr = nt_raw_mat,
  V = V,
  g = g_heat
)
colnames(nt_wavelet) <- colnames(nt_raw_mat)
all.equal(rownames(nt_wavelet),rownames(nt_raw_mat))

##totalRNA filtered
total_wavelet <- GraphWaveletDenoise(
  expr = total_raw_mat,
  V = V,
  g = g_heat
)
colnames(total_wavelet) <- colnames(total_raw_mat)

states_seu[["rb_wavelet"]] <-
  CreateAssayObject(
    counts =  rb_wavelet
  )
states_seu[["nt_wavelet"]] <-
  CreateAssayObject(
    counts = nt_wavelet
  )
states_seu[["total_wavelet"]] <-
  CreateAssayObject(
    counts = total_wavelet
  )

te_raw <- rb_raw_mat / (total_raw_mat + 1e-6)

te_wave <- rb_wavelet / (total_wavelet + 1e-6)

tau <- numeric(
  nrow(total_wavelet)
)
names(tau) <- rownames(total_wavelet)
for(i in seq_len(nrow(total_wavelet))){
  
  detected_cells <- total_raw_mat[i, ] > 0
  
  if(sum(detected_cells) == 0){
    
    tau[i] <- Inf
    
  } else {
    
    tau[i] <- quantile(
      total_wavelet[i, detected_cells],
      probs = 0.25,
      na.rm = TRUE
    )
    
  }
  
}
weight <- total_wavelet
for(i in seq_len(nrow(total_wavelet))){
  
  weight[i, ] <-
    total_wavelet[i, ] /
    (
      total_wavelet[i, ] +
        tau[i] +
        1e-8
    )
  
}
te_wave_weighted <- te_wave * weight

states_seu[["te_raw"]] <-
  CreateAssayObject(
    counts = te_raw
  )
states_seu[["te_wavelet"]] <-
  CreateAssayObject(
    counts = te_wave
  )

states_seu[["te_wavelet_weighted"]] <-
  CreateAssayObject(
    counts = te_wave_weighted
  )
PlotSpatialGene2 <- function(
    seu,
    gene,
    assay,
    layer = "data",
    vmax,
    title = NULL,
    point_size = 0.10,
    n_colors = 256
) {
  
  # 1. 获取表达数据并应用上限截断 (保留原本不共享色标的逻辑)
  expr <- GetAssayData(
    seu,
    assay = assay,
    layer = layer
  )[gene, ]
  
  expr_plot <- pmin(expr, vmax)
  
  # 2. 构建绘图数据框
  df <- data.frame(
    spatial_x = seu$column,
    spatial_y = -seu$row,
    Value = expr_plot
  )
  
  # 获取坐标极值用于等比例计算及坐标轴限制
  x_min <- min(df$spatial_x, na.rm = TRUE)
  x_max <- max(df$spatial_x, na.rm = TRUE)
  y_min <- min(df$spatial_y, na.rm = TRUE)
  y_max <- max(df$spatial_y, na.rm = TRUE)
  
  ar <- (y_max - y_min) / (x_max - x_min)
  
  # 3. 细腻的红色渐变色标 (同步目标函数)
  reds <- colorRampPalette(
    c(
      "white",
      "#fee8c8",
      "#fdbb84",
      "#fc8d59",
      "#d7301f",
      "#7f0000"
    )
  )(n_colors)
  
  # 默认标题
  if (is.null(title)) {
    title <- sprintf("Spatial expression of %s", gene)
  }
  
  # 4. ggplot 绘图
  p <- ggplot(df, aes(x = spatial_x, y = spatial_y, color = Value)) +
    geom_point(size = point_size) +
    
    # 坐标轴缩放与留白 (同步目标函数: expand = 0.02)
    scale_x_continuous(limits = c(x_min, x_max), expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(y_min, y_max), expand = expansion(mult = 0.02)) +
    
    # 独立渐变色标配置
    scale_color_gradientn(
      colours = reds,
      limits = c(0, vmax),
      breaks = c(0, vmax),
      labels = c(
        sprintf("%.2f", 0),
        sprintf("%.2f", vmax)
      ),
      guide = guide_colorbar(
        title = "Signal intensity",
        title.position = "top",
        barwidth = 0.6,   # 移除原本的 unit() 包裹，使用目标函数的简洁数字
        barheight = 6
      )
    ) +
    
    # 主题设定 (theme_bw + 细节微调)
    theme_bw(base_size = 12) +
    
    theme(
      aspect.ratio = ar,
      panel.grid = element_blank(),
      
      # 同步目标函数：轴标题增加恰当的 Margin，避免文字紧贴坐标轴
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 6)),
      
      # 保持原代码的清爽：去掉不必要的数字刻度，突出空间结构
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      
      plot.title = element_text(hjust = 0.5, face = "bold")
      # 注：已移除原代码中的 legend.title = element_blank()，以确保 Signal intensity 正常显示
    ) +
    
    # 坐标轴与图例名称 (同步目标函数)
    labs(
      title = title,
      x = "Column ", 
      y = "Row", 
      color = NULL
    )
  
  return(p)
}
GetVmax <- function(mat, gene, qmax = 0.99){
  
  quantile(
    mat[gene, ],
    probs = qmax,
    na.rm = TRUE
  )
}

genes_use <- c("Nrgn","C1ql2","Chgb", "Mbp","Plp1","C1qa","Ttr","Ly6c1","Rarres2","Ptgds","Nefl","Pcsk1n", "Kif5c","Olfm1")
genes_use <- c("Nrgn","C1ql2","Kif5c")
states_seu <- NormalizeData(
  states_seu,
  assay = "rb_wavelet",
  normalization.method = "RC",
  scale.factor = median(colSums(rb_wavelet)),  
  verbose = FALSE
)

states_seu <- NormalizeData(
  states_seu,
  assay = "total_wavelet",
  normalization.method = "RC",
  scale.factor = median(colSums(total_wavelet)),  
  verbose = FALSE
)

############################################################
## Total RNA
############################################################
total_wave_norm = GetAssayData(states_seu, assay = 'total_wavelet',layer = 'data')
plot_total <- list()

for(gene in genes_use){
  
  vmax_raw  <- GetVmax(total_norm, gene)
  vmax_wave <- GetVmax(GetAssayData(states_seu, assay = 'total_wavelet',layer = 'data'), gene)
  
  plot_total[[length(plot_total)+1]] <-
    PlotSpatialGene2(
      states_seu,
      gene = gene,
      assay = "RNA",layer = 'data',
      vmax = vmax_raw,
      title = paste0(gene, "\nTotal Raw")
    )
  
  plot_total[[length(plot_total)+1]] <-
    PlotSpatialGene2(
      states_seu,
      gene = gene,
      assay = "total_wavelet",layer = 'data',
      vmax = vmax_wave,
      title = paste0(gene, "\nTotal Filtered")
    )
}

p_total <- wrap_plots(plot_total, ncol = 2) +
  plot_annotation(title = "Total RNA (ntRNA + rbRNA)")
p_total

plot_rb <- list()

for(gene in genes_use){
  
  vmax_raw  <- GetVmax(rb_norm, gene)
  vmax_wave <- GetVmax(GetAssayData(states_seu, assay = 'rb_wavelet',layer = 'data'), gene)
  
  plot_rb[[length(plot_rb)+1]] <-
    PlotSpatialGene2(
      states_seu,
      gene = gene,
      assay = "rbRNA",layer = 'data',
      vmax = vmax_raw,
      title = paste0(gene, "\nRb Raw")
    )
  
  plot_rb[[length(plot_rb)+1]] <-
    PlotSpatialGene2(
      states_seu,
      gene = gene,
      assay = "rb_wavelet",layer = 'data',
      vmax = vmax_wave,
      title = paste0(gene, "\nRb Filtered")
    )
}

p_rb <- wrap_plots(plot_rb, ncol = 2) +
  plot_annotation(title = "rbRNA")
p_rb
############################################################
## TE
############################################################

plot_te <- list()

for(gene in genes_use){
  
  vmax_raw  <- GetVmax(te_raw, gene)
  vmax_wave <- GetVmax(te_wave_weighted, gene)
  
  plot_te[[length(plot_te)+1]] <-
    PlotSpatialGene2(
      states_seu,
      gene = gene,
      assay = "te_raw",layer = 'counts',
      vmax = vmax_raw,
      title = paste0(gene, "\nTE Raw")
    )
  
  plot_te[[length(plot_te)+1]] <-
    PlotSpatialGene2(
      states_seu,
      gene = gene,
      assay = "te_wavelet_weighted",layer = 'counts',
      vmax = vmax_wave,
      title = paste0(gene, "\nTE Filtered + Weight")
    )
}

p_te <- wrap_plots(plot_te, ncol = 2) +
  plot_annotation(title = "Translation Efficiency")

ggsave(
  filename = "0704_refined_total.pdf",
  plot = p_total,
  width = 12,
  height = 2 * length(genes_use),
  limitsize = FALSE,
  dpi = 1600
)

ggsave(
  filename = "0704_refined_rb.pdf",
  plot = p_rb,
  width = 12,
  height = 2 * length(genes_use),
  limitsize = FALSE,
  dpi = 1600
)

ggsave(
  filename = "0704_refined_TE.pdf",
  plot = p_te,
  width = 12,
  height = 2 * length(genes_use),
  limitsize = FALSE,
  dpi = 1600
)
D <- ComputeGrassmannDistance(
  U_consensus = scml_result$U,
  U_layers = scml_result$U_layers
)

D
