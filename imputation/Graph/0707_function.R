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
  library(patchwork)
})

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

ComputeGraphSpectrum<- function(
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
assert_named_list <- function(x, arg = "x") {
  if (!is.list(x) || is.null(names(x)) || any(names(x) == "")) {
    stop(arg, " must be a named list.", call. = FALSE)
  }
  invisible(TRUE)
}
EstimateEnergyCutoff <- function(
    expr,
    V,
    energy = 0.95,
    lambda = NULL,
    plot = TRUE
){
  X <- t(as.matrix(expr))
  coef <- t(V) %*% X
  E <- rowMeans(coef^2)
  cumE <- cumsum(E)
  cumE <- cumE / sum(E)
  k <- which(cumE >= energy)[1]
  if(is.null(lambda)){
    lambda_cut <- NA
  }else{
    lambda_cut <- lambda[k]
  }
  if(plot){
    
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(
      bty = "l",
      las = 1,
      cex.lab = 1.2,
      cex.axis = 1.1,
      mar = c(5,5,3,2)
    )
    
    plot(
      cumE,
      type = "l",
      lwd = 3,
      col = "#2C7BB6",
      ylim = c(0,1),
      xlab = "Graph frequency index",
      ylab = "Cumulative spectral energy"
    )
    
    abline(
      h = energy,
      col = "#D7191C",
      lty = 2,
      lwd = 2
    )
    
    abline(
      v = k,
      col = "#FDAE61",
      lty = 2,
      lwd = 2
    )
    
    points(
      k,
      cumE[k],
      pch = 19,
      cex = 1.3,
      col = "#FDAE61"
    )
    
    text(
      x = k,
      y = min(cumE[k] + 0.06, 0.98),
      labels = paste0(
        "k = ",
        k,
        "\nEnergy = ",
        round(cumE[k],3)
      ),
      pos = 4,
      cex = 0.9
    )
    
  }
  
  list(
    energy = E,
    cumulative = cumE,
    k = k,
    lambda_cutoff = lambda_cut
  )
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

build_knn_graph <- function(
    emb,
    k = 30,
    sigma_nn = 20,
    symmetrize = c("max", "mean"),
    sigma_floor = 1e-6
){
  
  symmetrize <- match.arg(symmetrize)
  
  emb <- as.matrix(emb)
  
  if (nrow(emb) <= k)
    stop("k must be smaller than number of cells.")
  
  if(any(!is.finite(emb)))
    stop("Embedding contains NA/NaN/Inf.")
  
  ############################################################
  ## kNN
  ############################################################
  
  knn <- FNN::get.knn(emb, k = k)
  
  idx <- knn$nn.index
  dst <- knn$nn.dist
  
  n <- nrow(idx)
  
  ############################################################
  ## Neighbor list
  ############################################################
  
  nbr <- lapply(seq_len(n), function(i) idx[i, ])
  
  ############################################################
  ## Adaptive sigma
  ############################################################
  
  sigma <- numeric(n)
  
  for(i in seq_len(n)){
    
    ########################################################
    ## candidate cells:
    ## neighbors + neighbors-of-neighbors
    ########################################################
    
    cand <- unique(c(
      nbr[[i]],
      unlist(nbr[nbr[[i]]])
    ))
    
    cand <- setdiff(cand, i)
    
    if(length(cand)==0){
      
      sigma[i] <- dst[i,k]
      
      next
    }
    
    ########################################################
    ## Jaccard
    ########################################################
    
    Ni <- nbr[[i]]
    
    js <- numeric(length(cand))
    
    for(m in seq_along(cand)){
      
      j <- cand[m]
      
      Nj <- nbr[[j]]
      
      inter <- length(intersect(Ni, Nj))
      
      if(inter==0){
        
        js[m] <- 0
        
      }else{
        
        uni <- length(union(Ni, Nj))
        
        js[m] <- inter / uni
        
      }
    }
    
    ########################################################
    ## non-zero Jaccard
    ########################################################
    
    keep <- which(js > 0)
    
    if(length(keep)==0){
      
      sigma[i] <- dst[i,k]
      
      next
    }
    
    cand <- cand[keep]
    js <- js[keep]
    
    ########################################################
    ## Euclidean distance
    ########################################################
    
    d <- sqrt(
      rowSums(
        (emb[cand,,drop=FALSE] -
           matrix(
             emb[i,],
             nrow=length(cand),
             ncol=ncol(emb),
             byrow=TRUE
           ))^2
      )
    )
    
    ########################################################
    ## smallest Jaccard
    ########################################################
    
    ord <- order(js, -d)
    
    cand <- cand[ord]
    d <- d[ord]
    js <- js[ord]
    
    ########################################################
    ## choose sigma_nn cells
    ########################################################
    
    m <- min(sigma_nn, length(cand))
    
    sigma[i] <- mean(d[1:m])
    
  }
  
  sigma[sigma < sigma_floor] <- sigma_floor
  
  ############################################################
  ## Gaussian weight
  ############################################################
  
  ii <- rep(seq_len(n), each = k)
  jj <- as.vector(t(idx))
  d  <- as.vector(t(dst))
  
  sigma_i <- rep(sigma, each = k)
  sigma_j <- sigma[jj]
  
  w <- exp(-(d^2)/(sigma_i * sigma_j))
  
  ############################################################
  ## Sparse graph
  ############################################################
  
  W <- Matrix::sparseMatrix(
    i = ii,
    j = jj,
    x = w,
    dims = c(n,n)
  )
  
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

ReweightSpatialGraph <- function(
    W_spatial,
    Emb_rb,
    Emb_nt,
    center = 0.5,
    beta = 2,
    gate_thr = NULL,
    gate_quantile = 0.15,
    symmetrize = c("max", "mean")
){
  emb <- cbind(
    scale(Emb_rb),
    scale(Emb_nt)
  )
  
  ## L2 normalize
  emb <- emb / sqrt(rowSums(emb^2) + 1e-10)
  ## spatial Graph fetching
  ii <- W_spatial@i + 1
  
  jj <- rep(
    seq_len(ncol(W_spatial)),
    diff(W_spatial@p)
  )

  ## cosine similarity
  sim <- rowSums(
    emb[ii, , drop = FALSE] *
      emb[jj, , drop = FALSE]
  )
  
  ## sigmoid gate
  gate <- 1 /
    (
      1 +
        exp(
          -beta * (sim - center)
        )
    )
  if (is.null(gate_thr)) {
    
    gate_thr <- as.numeric(
      quantile(
        gate,
        probs = gate_quantile
      )
    )
    
  }
  
  keep <- gate > gate_thr
  
  W_new <- Matrix::sparseMatrix(
    i = ii[keep],
    j = jj[keep],
    x = W_spatial@x[keep] * gate[keep],
    dims = dim(W_spatial)
  )
  

  
  if(symmetrize == "max"){
    
    W_new <- pmax(
      W_new,
      Matrix::t(W_new)
    )
    
  } else if(symmetrize == "mean"){
    
    W_new <- (
      W_new +
        Matrix::t(W_new)
    ) / 2
    
  }
  
  diag(W_new) <- 0
  
  W_new <- Matrix::drop0(W_new)
  
  ########################################################
  ## QC summary
  ########################################################
  
  deg_old <- Matrix::rowSums(W_spatial > 0)
  deg_new <- Matrix::rowSums(W_new > 0)
  
  isolated <- which(deg_new == 0)
  
  removed_neighbor <- deg_old - deg_new
  
  cat("\n")
  cat("========================================\n")
  cat(" Spatial Graph Reweighting parameters Summary\n")
  cat("========================================\n")
  
  cat(sprintf("Center                : %.3f\n", center))
  cat(sprintf("Beta                  : %.3f\n", beta))
  cat(sprintf("Gate threshold        : %.3f\n", gate_thr))
  
  cat("\n")
  
  cat(sprintf("Original edges        : %d\n", length(gate)))
  cat(sprintf("Remaining edges       : %d\n", sum(keep)))
  cat(sprintf(
    "Removed edges         : %d (%.2f%%)\n",
    sum(!keep),
    100 * mean(!keep)
  ))
  
  cat("\n")
  
  cat(sprintf(
    "Mean degree           : %.2f -> %.2f\n",
    mean(deg_old),
    mean(deg_new)
  ))
  
  cat(sprintf(
    "Median degree         : %.2f -> %.2f\n",
    median(deg_old),
    median(deg_new)
  ))
  
  cat(sprintf(
    "Mean removed neighbor : %.2f\n",
    mean(removed_neighbor)
  ))
  
  cat(sprintf(
    "Max removed neighbor  : %d\n",
    max(removed_neighbor)
  ))
  
  cat(sprintf(
    "Isolated cells        : %d\n",
    length(isolated)
  ))
  
  if(length(isolated) > 0){
    cat(
      "First isolated IDs    : ",
      paste(head(isolated, 20), collapse = ", "),
      "\n",
      sep = ""
    )
  }
  
  cat("\n")
  
  cat(sprintf(
    "Mean cosine similarity : %.3f\n",
    mean(sim)
  ))
  
  cat(sprintf(
    "Mean gate (all edges)  : %.3f\n",
    mean(gate)
  ))
  
  cat(sprintf(
    "Mean gate (kept edges) : %.3f\n",
    mean(gate[keep])
  ))
  
  cat("\n")
  
  cat("Cosine similarity summary\n")
  print(summary(sim))
  
  cat("\n")
  
  cat("Gate summary\n")
  print(summary(gate))
  
  cat("\n")
  
  cat("Gate quantiles\n")
  print(
    quantile(
      gate,
      probs = c(
        0,
        0.01,
        0.05,
        0.10,
        0.25,
        0.50,
        0.75,
        0.90,
        0.95,
        0.99,
        1
      )
    )
  )
  cat(sprintf(
    "Gate threshold        : %.3f\n",
    gate_thr
  ))
  
  cat(sprintf(
    "Gate quantile         : %.2f%%\n",
    100 * gate_quantile
  ))
  cat("========================================\n\n")
  return(W_new)
}
GraphFilterKernel <- function(
    lambda = NULL,
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
    peak_quantile = 0.15,
    plot = FALSE
){
  
  ############################################################
  ## If only requesting kernel (no plotting)
  ############################################################
  
  if (is.null(lambda) && !plot) {
    return(function(lambda){
      GraphFilterKernel(
        lambda = lambda,
        type = type,
        beta = beta,
        scale = scale,
        lambda_low = lambda_low,
        lambda_high = lambda_high,
        lambda_c = lambda_c,
        order = order,
        gamma = gamma,
        peak_quantile = peak_quantile,
        plot = FALSE
      )
    })
  }
  
  ############################################################
  
  if (is.null(lambda))
    stop("`lambda` must be supplied when plot = TRUE.")
  
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
  
  ############################################################
  ## Plot
  ############################################################
  
  if (plot) {
    
    ord <- order(lambda)
    
    plot(
      lambda[ord],
      g[ord],
      type = "l",
      lwd = 3,
      xlab = expression(lambda),
      ylab = "Normalized filter response",
      main = paste("Graph Filter:", type)
    )
    
    grid()
    
  }
  
  return(g)
  
}

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
ComputeTEConfidence <- function(
    rb_mat,
    nt_mat,
    raw_total = NULL,
    tau_method = c("mad", "quantile"),
    tau_quantile = 0.75,
    evidence_method = c("linear", "log", "exp"),
    confidence_method = c("michaelis", "exp", "sigmoid"),
    lambda = 1,
    sigmoid_a = 1,
    sigmoid_b = 1,
    eps = 1e-6
){
  
  stopifnot(
    all(dim(rb_mat) == dim(nt_mat))
  )
  
  stopifnot(
    identical(
      rownames(rb_mat),
      rownames(nt_mat)
    )
  )
  
  stopifnot(
    identical(
      colnames(rb_mat),
      colnames(nt_mat)
    )
  )
  
  total_mat <- rb_mat + nt_mat
  te <- rb_mat / (total_mat + eps)
  
  tau_method <- match.arg(tau_method)
  
  tau <- numeric(nrow(total_mat))
  names(tau) <- rownames(total_mat)
  
  for(i in seq_len(nrow(total_mat))){
    
    if(is.null(raw_total)){
      
      x <- total_mat[i, ]
      
    }else{
      
      detected <- raw_total[i, ] > 0
      
      if(sum(detected) < 5){
        
        x <- total_mat[i, ]
        
      }else{
        
        x <- total_mat[i, detected]
        
      }
    }
    
    if(tau_method=="mad"){
      
      tau[i] <-
        median(x) +
        mad(x)
      
    }
    
    if(tau_method=="quantile"){
      
      tau[i] <-
        quantile(
          x,
          tau_quantile,
          names = FALSE
        )
      
    }
    
    tau[i] <- max(
      tau[i],
      eps
    )
  }
  evidence_method <- match.arg(evidence_method)
  ratio <-
    sweep(
      total_mat,
      1,
      tau,
      "/"
    )
  if(evidence_method=="linear"){
    evidence <- ratio
  }
  if(evidence_method=="log"){
    evidence <-
      log1p(ratio)
    
  }
  if(evidence_method=="exp"){
    evidence <-
      1 -
      exp(-ratio)
  }
  
  confidence_method <-
    match.arg(
      confidence_method
    )
  
  if(confidence_method=="michaelis"){
    confidence <-
      evidence /
      (
        evidence +
          lambda
      )
    
  }
  
  if(confidence_method=="exp"){
    
    confidence <-
      1 -
      exp(-evidence)
    
  }
  
  if(confidence_method=="sigmoid"){
    
    confidence <-
      1 /
      (
        1 +
          exp(
            -sigmoid_a *
              (
                evidence -
                  sigmoid_b
              )
          )
      )
    
  }
  te_weighted <-
    te *
    confidence
  
  return(
    list(
      te = te,
      tau = tau,
      evidence = evidence,
      confidence = confidence,
      te_weighted = te_weighted
    )
  )
}
BuildSCMLGraph <- build_knn_graph
