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
