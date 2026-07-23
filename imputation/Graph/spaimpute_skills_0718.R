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


# 0. pre-part -------------------------------------------------------------

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
}#function related to vector normalization of U_consensus in run_scml_core eig_vectors


# 1.Graph_processing  ------------------------------------------------------

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

build_snn_graph <- function(
    emb,
    k = 30,
    prune = 1/15
){
  
  emb <- as.matrix(emb)
  
  if(nrow(emb) <= k)
    stop("k must be smaller than number of cells.")
  
  if(any(!is.finite(emb)))
    stop("Embedding contains NA/NaN/Inf.")
  
  
  ############################################################
  ## kNN
  ############################################################
  
  knn <- FNN::get.knn(
    emb,
    k = k
  )
  
  idx <- knn$nn.index
  
  n <- nrow(idx)
  
  
  ############################################################
  ## neighbor list
  ############################################################
  
  nbr <- lapply(
    seq_len(n),
    function(i) idx[i,]
  )
  
  
  ############################################################
  ## Compute SNN Jaccard weight
  ############################################################
  
  ii <- c()
  jj <- c()
  ww <- c()
  
  
  for(i in seq_len(n)){
    
    
    Ni <- nbr[[i]]
    
    
    ########################################################
    ## only calculate KNN neighbors
    ## same as Seurat sparse construction
    ########################################################
    
    for(j in Ni){
      
      
      Nj <- nbr[[j]]
      
      
      inter <- length(
        intersect(Ni,Nj)
      )
      
      
      if(inter==0){
        
        sij <- 0
        
      }else{
        
        uni <- length(
          union(Ni,Nj)
        )
        
        sij <- inter / uni
        
      }
      
      
      ####################################################
      ## prune.SNN
      ####################################################
      
      if(sij >= prune){
        
        ii <- c(ii,i)
        jj <- c(jj,j)
        ww <- c(ww,sij)
        
      }
      
    }
    
  }
  
  
  ############################################################
  ## sparse matrix
  ############################################################
  
  
  W <- Matrix::sparseMatrix(
    i = ii,
    j = jj,
    x = ww,
    dims=c(n,n)
  )
  
  
  
  diag(W)<-0
  
  
  drop0(W)
  return(W)
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

build_knn_graph_bandwidth_with_largest_jaccards <- function(
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
    ## LARGEST!!!!!!!!!!!! Jaccard!!!!!!!!!
    ########################################################
    
    ord <- order(-js, d) 
    
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

build_knn_graph_tradition <- function(emb, k = 30, symmetrize = c("max", "mean"), sigma_floor = 1e-6) {
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

NeighborExpressionEmbedding <- function(
    W,
    expr,
    normalize = TRUE,
    alpha = 0.5
){
  #现在主要的spatial利用来源
  deg <- Matrix::rowSums(W) #这个W是spatial original
  deg[deg==0] <- 1 #
  
  
  W_norm <-
    Diagonal(
      x=1/deg
    ) %*%
    W #D-1W, 类似于输出一个“celli neighbor中考量cell j表达多少的”权重矩阵
  
  
  
  ## neighborhood expression
  
  X_neighbor <-
    W_norm %*% expr 
  
  
  
  ## include self expression
  
  if(alpha>0){
    
    X_neighbor <-
      alpha * expr +
      (1-alpha)*X_neighbor #1-alpha为neighborhood的weight
    
  }
  
  
  return(
    X_neighbor
  )
  
}

BuildSCMLGraph <- build_knn_graph

BuildSCMLGraph_tradition <- build_knn_graph_tradition

BuildSCMLGraph_largest_Jaccard <- build_knn_graph_bandwidth_with_largest_jaccards

# 2.Laplacian Matrix/Spectrum related function -------------------------------------

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
    which="SA"
  )
  
  list(
    lambda=eig$values,
    V=eig$vectors,
    L = L
  )
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
  
  eig <- RSpectra::eigs_sym(L, k = ndim, which = "SA")
  U <- qr.Q(qr(eig$vectors))
  list(U = U, values = eig$values)
}

# 3.core function  --------------------------------------------------------

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
  eig_consensus <- RSpectra::eigs_sym(L_modified, k = ndim, which = "SA")
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
  ##OLD VERSION HERE
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

GraphSpectralDenoise <- function(
    expr,
    V,
    lambda,
    g,
    remove_dc = TRUE,
    sort_lambda = TRUE,
    add_dc = FALSE
){
  
  ################################################
  ## expr: gene x cell, seurat 
  ## V: cell x eigenvectors
  ## lambda: eigenvalues
  ## g:
  ##   vector(length=lambda) 
  ##   OR matrix(freq x gene)
  ## THIS IS OLD!!!!!!!
  ################################################
  
  
  X <- t(as.matrix(expr))   # cell x gene
  
  
  ################################################
  ## sort eigen system
  ################################################
  
  if(sort_lambda){
    
    ord <- order(lambda)
    
    lambda <- lambda[ord]
    V <- V[,ord,drop=FALSE]
    
  }
  
  
  ################################################
  ## remove DC component
  ################################################
  
  if(remove_dc){
    
    v0 <- V[,1,drop=FALSE]
    
    X_dc <- v0 %*% (t(v0)%*%X)
    
    X_spatial <- X - X_dc
    cat(
      "Spatial residual range:",
      range(X_spatial),
      "\n"
    )
    
    V_use <- V[,-1,drop=FALSE] #drop the smallest lambda and its egi-vector
    
    cat("drop lambda: ",lambda[1],"\n")
    
    lambda_use <- lambda[-1]
    
  }else{
    
    X_dc <- matrix(
      0,
      nrow=nrow(X),
      ncol=ncol(X)
    )
    
    X_spatial <- X
    
    V_use <- V
    lambda_use <- lambda
    
  }

  
  Xhat <- t(V_use)%*%X_spatial
  
  if(is.vector(g)){
    
    Xhat_filtered <- 
      g * Xhat
    
    
  }else{
    
    # gene-wise filter
    
    Xhat_filtered <-
      g * Xhat
    
  }
  
  X_spatial_filtered <-
    V_use %*% Xhat_filtered
  
  if(add_dc){
    
    X_filtered <-
      X_dc + X_spatial_filtered
    
  }else{
    
    X_filtered <-
      X_spatial_filtered
    
  }
  
  
  X_filtered[X_filtered<0] <- 0
  
  
  
  out <- t(X_filtered)
  
  rownames(out) <- rownames(expr)
  colnames(out) <- colnames(expr)
  
  
  return(
    list(
      expr_filtered = out,
      X_dc = t(X_dc),
      X_spatial = t(X_spatial),
      lambda = lambda_use,
      V = V_use,
      g = g
    )
  )
}

ComputeGeneWiseKernel <- function(
    X_spatial,
    L,
    lambda,
    beta_min = 1,
    beta_max = 10
){
  
  
  ################################################
  ## X_spatial:
  ## cell x gene
  ## when beta comes large, the high frequency gets low weight, like being killed
  ################################################

  
  roughness <- apply(
    X_spatial,
    2,
    function(x){
      
      as.numeric(
        t(x)%*%L%*%x
      ) /
        (sum(x^2)+1e-8)
      
    }
  )
  
  rough_scaled <-
    (roughness-min(roughness))/
    (max(roughness)-min(roughness)+1e-8)
  
  
  # high roughness -> low beta
  beta_gene <-
    beta_max -
    rough_scaled *
    (beta_max-beta_min)
  
  
  lambda_scaled <-
    lambda/max(lambda)
  
  G <-
    exp(
      -lambda_scaled %o% beta_gene
    )
  return(
    list(
      beta_gene = beta_gene,
      roughness = roughness,
      G = G
    )
  )
}


resolve_layer_values_pan <- function(
    values,
    layer_names,
    default_value = 1,
    arg = "values"
) {
  
  ## 默认值
  if (is.null(values)) {
    out <- rep(default_value, length(layer_names))
    names(out) <- layer_names
    
    ## 单个数字
  } else if (length(values) == 1 && is.null(names(values))) {
    
    out <- rep(as.numeric(values), length(layer_names))
    names(out) <- layer_names
    
    ## 否则开始解析
  } else {
    
    values <- unlist(values)
    
    ## 有重复名字
    if (!is.null(names(values))) {
      if (anyDuplicated(names(values)))
        stop(arg, " contains duplicated layer names.",
             call. = FALSE)
    }
    
    ## named vector
    if (!is.null(names(values)) &&
        all(names(values) != "")) {
      
      unknown <- setdiff(names(values), layer_names)
      if (length(unknown) > 0) {
        stop(
          arg,
          " contains unknown layer(s): ",
          paste(unknown, collapse = ", "),
          call. = FALSE
        )
      }
      
      missing <- setdiff(layer_names, names(values))
      if (length(missing) > 0) {
        stop(
          arg,
          " is missing layer(s): ",
          paste(missing, collapse = ", "),
          call. = FALSE
        )
      }
      
      out <- as.numeric(values[layer_names])
      names(out) <- layer_names
      
    } else {
      
      ## unnamed vector
      
      if (length(values) != length(layer_names)) {
        
        stop(
          arg,
          " must have length 1 or ",
          length(layer_names),
          ".",
          call. = FALSE
        )
        
      }
      
      out <- as.numeric(values)
      names(out) <- layer_names
      
    }
    
  }
  
  if (any(!is.finite(out)))
    stop(arg, " contains non-finite values.",
         call. = FALSE)
  
  out
}

run_scml_core_pan <- function(
    graph_list,
    ndim = 30,
    alpha = NULL,
    layer_weight = NULL,
    return_laplacians = TRUE
){
  
  assert_named_list(graph_list, "graph_list")
  
  layer_names <- names(graph_list)
  
  n_vertices <- vapply(
    graph_list,
    nrow,
    integer(1)
  )
  
  if (
    length(unique(n_vertices)) != 1 ||
    any(vapply(graph_list, ncol, integer(1)) != n_vertices[1])
  ){
    stop(
      "All graphs must be square matrices with identical dimensions.",
      call.=FALSE
    )
  }
  
  alpha_vec <- resolve_layer_values_pan(
    alpha,
    layer_names,
    default_value = 1,
    arg = "alpha"
  )
  
  weight_vec <- resolve_layer_values_pan(
    layer_weight,
    layer_names,
    default_value = 1,
    arg = "layer_weight"
  )
  
  cat("\n========== Layer Parameters ==========\n")
  
  print(
    data.frame(
      layer = layer_names,
      weight = as.numeric(weight_vec),
      alpha = as.numeric(alpha_vec),
      row.names = NULL
    )
  )
  
  n <- n_vertices[1]
  
  L_sum <- Matrix(0, n, n, sparse = TRUE)
  P_sum <- Matrix(0, n, n, sparse = TRUE)
  
  U_layers <- list()
  eig_layers <- list()
  laplacians <- list()
  
  for(layer in layer_names){
    
    L <- compute_normalized_laplacian(
      graph_list[[layer]]
    )
    
    subspace <- get_laplacian_subspace(
      L,
      ndim = ndim
    )
    
    U_layers[[layer]] <- subspace$U
    eig_layers[[layer]] <- subspace$values
    
    if(return_laplacians){
      laplacians[[layer]] <- L
    }
    
    L_sum <-
      L_sum +
      weight_vec[layer] * L
    
    P_sum <-
      P_sum +
      alpha_vec[layer] *
      tcrossprod(subspace$U)
    
  }
  
  L_modified <- L_sum - P_sum
  
  eig_consensus <-
    RSpectra::eigs_sym(
      L_modified,
      k = ndim,
      which = "SA"
    )
  
  U_consensus <-
    row_normalize(
      eig_consensus$vectors
    )
  
  colnames(U_consensus) <-
    paste0(
      "SCML_",
      seq_len(ncol(U_consensus))
    )
  
  list(
    U = U_consensus,
    eigenvalues = eig_consensus$values,
    U_layers = U_layers,
    eigenvalues_layers = eig_layers,
    L_modified = L_modified,
    laplacians =
      if(return_laplacians)
        laplacians
    else
      NULL,
    alpha = alpha_vec,
    layer_weight = weight_vec
  )
}#run_scml_core_pan(graph_list,alpha = c(A = i, ...), layer_weight = c(A = j, ...), ndims = ...)

RunGeneWiseGraphFilter <- function(
    expr,
    V,
    lambda,
    L,
    beta_min = 1,
    beta_max = 10
){
  
  X <- t(as.matrix(expr))
  v0 <- V[,1,drop=FALSE]
  X_dc <-
    v0 %*% (t(v0)%*%X) #对于一个图的谱分解，第一特征值理论为0，以及其对应的特征向量v0为一个元素完全一致的列向量；X_dc这里可以理解为全局的平均表达水平，一种library size，该矩阵没有区分cell的能力
  
  
  X_spatial <-
    X-X_dc  #可以立即为细胞间差异、空间结构、domain signal
  
  
  V_sp <- V[,-1,drop=FALSE]
  lambda_sp <- lambda[-1]
  
  kernel <- ComputeGeneWiseKernel( 
    X_spatial = X_spatial,
    L = L,
    lambda = lambda_sp,
    beta_min = beta_min,
    beta_max = beta_max
  )#滤波器计算，这里尝试了gene-wise的特异滤波，如果不考虑gene特异，直接设置betamin = betamax
  
  Xhat <-
    t(V_sp)%*%X_spatial
  
  
  
  Xhat_filtered <-
    kernel$G * Xhat
  
  
  X_spatial_filtered <-
    V_sp %*% Xhat_filtered
  
  
  
  X_filtered <-
    X_spatial_filtered
  
  
  X_filtered[X_filtered<0] <- 0
  
  
  
  out <- t(X_filtered)
  
  
  rownames(out)<-
    rownames(expr)
  
  colnames(out)<-
    colnames(expr)
  
  
  
  return(
    list(
      filtered = out,
      beta_gene = kernel$beta_gene,
      roughness = kernel$roughness,
      G = kernel$G,
      X_dc = t(X_dc),
      X_spatial = t(X_spatial)
    )
  )
}

QuantileRescaleGene <- function(
    expr_raw,#原来的矩阵，norm
    expr_denoise,#去噪，插补完成的矩阵
    q = 0.99
){
  
  ## both:
  ## gene x cell
  
  stopifnot(
    all(dim(expr_raw)==dim(expr_denoise))
  )
  
  
  q_raw <- apply(
    expr_raw,
    1,
    quantile,
    probs=q,
    na.rm=TRUE
  )
  
  
  q_denoise <- apply(
    expr_denoise,
    1,
    quantile,
    probs=q,
    na.rm=TRUE
  )
  
  
  scale_factor <-
    q_raw /
    (q_denoise + 1e-8)
  
  
  expr_rescale <-
    expr_denoise *
    scale_factor
  
  
  return(
    list(
      expr = expr_rescale,
      q_raw = q_raw,
      q_denoise = q_denoise,
      scale_factor = scale_factor
    )
  )
}

# 4. SpatialLayer_preprocessing -------------------------------------------

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
  print(quantile(
    sim,
    probs = seq(0, 1, 0.05),
    na.rm = TRUE
  ))
  
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

ReweightSpatialGraphAdaptive <- function(
    W_spatial,
    Emb_rb,
    Emb_nt,
    center = NULL,
    beta = NULL,
    gate_thr = NULL,
    gate_quantile = 0.15,
    symmetrize = c("max", "mean"),
    beta_scale = 4,
    eps = 1e-6
){
  
  symmetrize <- match.arg(symmetrize)
  
  ############################################################
  ## Joint embedding
  ############################################################
  
  emb <- cbind(
    scale(Emb_rb),
    scale(Emb_nt)
  )
  
  ############################################################
  ## L2 normalization
  ############################################################
  
  emb <- emb / sqrt(rowSums(emb^2) + 1e-10)
  
  ############################################################
  ## Fetch graph edges
  ############################################################
  
  ii <- W_spatial@i + 1
  
  jj <- rep(
    seq_len(ncol(W_spatial)),
    diff(W_spatial@p)
  )
  
  ############################################################
  ## Cosine similarity
  ############################################################
  
  sim <- rowSums(
    emb[ii, , drop = FALSE] *
      emb[jj, , drop = FALSE]
  )
  
  ############################################################
  ## Adaptive center
  ############################################################
  
  if(is.null(center)){
    
    center <- median(sim)
    
  }
  
  ############################################################
  ## Adaptive beta
  ############################################################
  
  if(is.null(beta)){
    
    sim_iqr <- IQR(sim)
    
    sim_iqr <- max(sim_iqr, eps)
    
    beta <- beta_scale / sim_iqr
    
  }
  
  ############################################################
  ## Sigmoid gate
  ############################################################
  
  gate <- 1 /
    (
      1 +
        exp(
          -beta * (sim - center)
        )
    )
  
  ############################################################
  ## Adaptive threshold
  ############################################################
  
  if(is.null(gate_thr)){
    
    gate_thr <- as.numeric(
      quantile(
        gate,
        probs = gate_quantile
      )
    )
    
  }
  
  ############################################################
  ## Edge pruning
  ############################################################
  
  keep <- gate > gate_thr
  
  W_new <- Matrix::sparseMatrix(
    
    i = ii[keep],
    
    j = jj[keep],
    
    x = W_spatial@x[keep] * gate[keep],
    
    dims = dim(W_spatial)
    
  )
  
  ############################################################
  ## Symmetrize
  ############################################################
  
  if(symmetrize == "max"){
    
    W_new <- pmax(
      W_new,
      Matrix::t(W_new)
    )
    
  }else{
    
    W_new <- (
      W_new +
        Matrix::t(W_new)
    ) / 2
    
  }
  
  diag(W_new) <- 0
  
  W_new <- Matrix::drop0(W_new)
  
  ############################################################
  ## Report
  ############################################################
  
  cat("Adaptive Spatial Graph Reweighting\n")
  cat("-----------------------------------\n")
  cat(sprintf("Center (median)      : %.3f\n", center))
  cat(sprintf("Beta (adaptive)      : %.3f\n", beta))
  cat(sprintf("Similarity median    : %.3f\n", median(sim)))
  cat(sprintf("Similarity IQR       : %.3f\n", IQR(sim)))
  cat(sprintf("Gate threshold       : %.3f\n", gate_thr))
  cat(sprintf("Original edges       : %d\n", length(W_spatial@x)))
  cat(sprintf("Remaining edges      : %d\n", length(W_new@x)))
  cat(sprintf(
    "Removed edges        : %d (%.2f%%)\n",
    length(W_spatial@x) - length(W_new@x),
    100 * (1 - length(W_new@x) / length(W_spatial@x))
  ))
  
  invisible(list(
    W = W_new,
    center = center,
    beta = beta,
    similarity = sim,
    gate = gate,
    gate_threshold = gate_thr
  ))
  
}
