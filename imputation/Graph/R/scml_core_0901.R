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

build_radius_capped_knn_graph <- function(
    coords,
    k = 15,
    radius = 500,
    sigma_floor = 1e-6,
    weight_mode = "adaptive_gaussian"
){
  
  coords <- as.matrix(coords)
  
  if(ncol(coords) != 2){
    stop("coords must be an n x 2 matrix.")
  }
  
  if(any(!is.finite(coords))){
    stop("coords contains NA/NaN/Inf.")
  }
  
  n <- nrow(coords)
  
  if(n <= 1){
    stop("At least two cells are required.")
  }
  
  ########################################################
  ## Pairwise distance
  ########################################################
  
  D <- as.matrix(dist(coords))
  
  diag(D) <- Inf
  
  
  ########################################################
  ## Construct radius-capped kNN
  ########################################################
  
  neighbor_list <- vector(
    "list",
    n
  )
  
  distance_list <- vector(
    "list",
    n
  )
  
  for(i in seq_len(n)){
    
    di <- D[i, ]
    
    ## only cells inside the maximum radius
    candidate <- which(
      is.finite(di) &
        di <= radius
    )
    
    if(length(candidate) == 0){
      
      neighbor_list[[i]] <- integer(0)
      distance_list[[i]] <- numeric(0)
      
      next
    }
    
    ## sort by distance
    candidate <- candidate[
      order(di[candidate])
    ]
    
    ## keep at most k nearest neighbors
    if(length(candidate) > k){
      candidate <- candidate[seq_len(k)]
    }
    
    neighbor_list[[i]] <-
      candidate
    
    distance_list[[i]] <-
      di[candidate]
  }
  
  
  ########################################################
  ## Local adaptive sigma
  ########################################################
  
  sigma <- numeric(n)
  
  for(i in seq_len(n)){
    
    di <- distance_list[[i]]
    
    if(length(di) == 0){
      
      sigma[i] <- sigma_floor
      
    } else {
      
      sigma[i] <- median(di)
      
    }
  }
  
  sigma[
    !is.finite(sigma) |
      sigma < sigma_floor
  ] <- sigma_floor
  
  
  ########################################################
  ## Build directed weights
  ########################################################
  
  ii <- integer(0)
  jj <- integer(0)
  ww <- numeric(0)
  
  for(i in seq_len(n)){
    
    nb <- neighbor_list[[i]]
    
    if(length(nb) == 0){
      next
    }
    
    d <- distance_list[[i]]
    
    if(weight_mode == "uniform"){
      
      w <- rep(
        1 / length(nb),
        length(nb)
      )
      
    } else if(weight_mode == "reciprocal"){
      
      w <- 1 / pmax(
        d,
        sigma_floor
      )
      
      w <- w / sum(w)
      
    } else if(weight_mode == "adaptive_gaussian"){
      
      sigma_i <- sigma[i]
      
      w <- exp(
        -d^2 /
          (2 * sigma_i^2)
      )
      
      w <- w / sum(w)
      
    } else {
      
      stop(
        "Unknown weight_mode. ",
        "Use 'uniform', 'reciprocal', ",
        "or 'adaptive_gaussian'."
      )
    }
    
    ii <- c(
      ii,
      rep(i, length(nb))
    )
    
    jj <- c(
      jj,
      nb
    )
    
    ww <- c(
      ww,
      w
    )
  }
  
  
  ########################################################
  ## Sparse row-normalized adjacency
  ########################################################
  
  W <- Matrix::sparseMatrix(
    i = ii,
    j = jj,
    x = ww,
    dims = c(n, n)
  )
  
  diag(W) <- 0
  
  W <- Matrix::drop0(W)
  
  
  ########################################################
  ## Neighbor number summary
  ########################################################
  
  neighbor_n <- lengths(neighbor_list)
  
  n_target <- sum(
    neighbor_n == k
  )
  
  n_less <- sum(
    neighbor_n < k
  )
  
  cat("\n")
  cat("========================================\n")
  cat("Radius-capped kNN neighbor summary\n")
  cat("========================================\n")
  
  cat(
    "Target neighbors per cell:",
    k,
    "\n"
  )
  
  cat(
    "Maximum radius:",
    radius,
    "\n"
  )
  
  cat(
    "Cells with target number of neighbors:",
    n_target,
    "/",
    n,
    sprintf(
      "(%.2f%%)",
      100 * n_target / n
    ),
    "\n"
  )
  
  cat(
    "Cells with fewer than target neighbors:",
    n_less,
    "/",
    n,
    sprintf(
      "(%.2f%%)",
      100 * n_less / n
    ),
    "\n"
  )
  
  if(n_less > 0){
    
    cat(
      "Neighbor number summary for cells with fewer than",
      k,
      "neighbors:\n"
    )
    
    print(
      summary(
        neighbor_n[
          neighbor_n < k
        ]
      )
    )
    
  } else {
    
    cat(
      "All cells reached the target number of neighbors.\n"
    )
    
  }
  
  cat(
    "Overall neighbor-number summary:\n"
  )
  
  print(
    summary(neighbor_n)
  )
  
  cat(
    "========================================\n\n"
  )
  
  
  return(W)
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

# 1. Neighbor_embedding ---------------------------------------------------

NeighborExpressionEmbedding <- function(
    W,
    expr,
    normalize = TRUE,
    alpha = 0.5
){
  #现在主要的spatial利用来源
  deg <- Matrix::rowSums(W) #this W should be the raw-spatial weight matrix, in which [i,j] stands for the weight in neighborhood space.
  deg[deg==0] <- 1 #孤立点
  W_norm <-
    Diagonal(
      x=1/deg
    ) %*%
    W #D-1W, 类似于输出一个“celli neighbor中考量cell j表达多少的”权重矩阵
  ## neighborhood expression
  X_neighbor <- W_norm %*% expr 

  if(alpha>0){ #控制neighborhood expression的权重
    X_neighbor <-
      alpha * expr +
      (1-alpha)*X_neighbor #1-alpha为neighborhood的weight
  }
  return(
    X_neighbor
  )
}


# 2. Laplace Matrix calculation and subspace spanning -------------------------------------------

#here we shall introduce 2 kinds of laplacian Matrix: random-walk(whose first eigen-vector equals to [1,1,1,1,1,1...], 
#and the normalized (whose first eigen-vector equals to [d1, d2, ...]))

compute_normalized_laplacian <- function(W, eps = 1e-10) {
  cat('symmetric normalized Laplacian: I-D-1/2WD-1/2 \n')
  W <- as(W, "dgCMatrix")
  deg <- Matrix::rowSums(W)
  Dinv <- Diagonal(x = 1 / sqrt(deg + eps))
  Diagonal(nrow(W)) - Dinv %*% W %*% Dinv
}

compute_random_walk_laplacian <- function(W, eps = 1e-10) {
  cat('random walk Laplacian:I-D-1W \n')
  W <- as(W, "dgCMatrix")
  deg <- Matrix::rowSums(W)
  Dinv <- Diagonal(
    x = 1 / (deg + eps)
  )
  L_rw <- Diagonal(nrow(W)) - Dinv %*% W 
  return(L_rw)
}

get_laplacian_subspace <- function(L, ndim = 30) {
  cat("subspace spanning for Laplacian Matrix")
  if (ndim >= nrow(L)) {
    stop("ndim must be smaller than the number of cells/vertices.", call. = FALSE)
  }
  eig <- RSpectra::eigs_sym(L, k = ndim, which = "SA")
  U <- qr.Q(qr(eig$vectors))
  list(U = U, values = eig$values)
}

GrassmannDistance <- function(U1, U2){
  cat("calculating Grassmann Distance for 2 spaces/n")
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


# 3. cell-specific alpha calculation --------------------------------------
#in this section we use embeddings and weight matrix for cell specific alpha calculation
#to evaluate whether Module A classify the cell better than Module B ultilizing dirichlet energy

NormalizeEmbedding <- function(Z){
  Z <- as.matrix(Z)
  Z / sqrt(rowSums(Z^2)+1e-10)
}

CheckCellConsistency <- function(
    graph_list,
    embedding_list
){
  layers <- names(graph_list)
  graph_cells <- lapply(
    graph_list,
    rownames
  )
  embedding_cells <- lapply(
    embedding_list,
    rownames
  )
  ref_cells <- graph_cells[[1]]
  for(nm in names(graph_cells)){
    
    if(!identical(
      graph_cells[[nm]],
      ref_cells
    )){
      
      stop(
        paste0(
          "Graph rownames mismatch: ",
          nm
        )
      )
    }
    
    if(!identical(
      colnames(graph_list[[nm]]),
      ref_cells
    )){
      
      stop(
        paste0(
          "Graph colnames mismatch: ",
          nm
        )
      )
    }
  }
  
  
  for(nm in names(embedding_cells)){
    
    if(!identical(
      embedding_cells[[nm]],
      ref_cells
    )){
      
      stop(
        paste0(
          "Embedding rownames mismatch: ",
          nm
        )
      )
    }
  }
  
  cat(
    "All graph and embedding cell identities are consistent."
  )
  
}

NormalizeEnergy <- function(E){
  E / (median(E) + 1e-10)
}

ComputeNodeDirichletEnergy <- function(
    W,
    embedding,
    normalize = TRUE
){
  stopifnot(
    nrow(W)==nrow(embedding)
  )
  n <- nrow(W)
  deg <- Matrix::rowSums(W)
  Z <- as.matrix(embedding)
  
  
  energy <- numeric(n)
  
  for(i in seq_len(n)){
    nei <- which(W[i,] > 0)
    if(length(nei)==0){
      energy[i] <- 0
      next
    }
    
    diff <- sweep(
      Z[nei,,drop=FALSE],
      2,
      Z[i,],
      "-"
    )
    
    dist2 <- rowSums(diff^2)
    
    energy[i] <-
      0.5 *
      sum(
        W[i,nei] * dist2
      ) #energy equals 1/2 Wij * dist2(i,j)
    
  }
  
  if(normalize){
    energy <-
      energy/(deg+1e-10)
    
  } #node wise normalization: based on degree
  return(energy)
}

ComputeCrossDirichletEnergy <- function(
    W_source,
    embedding_target,
    normalize = TRUE
){
  n <- nrow(W_source)
  deg <- Matrix::rowSums(W_source)
  Z <- as.matrix(
    embedding_target
  )
  
  
  energy <- numeric(n)
  
  
  for(i in seq_len(n)){
    
    nei <- which(
      W_source[i,] > 0
    )
    
    if(length(nei)==0){
      energy[i] <- 0
      next
    }
    
    
    diff <- sweep(
      Z[nei,,drop=FALSE],
      2,
      Z[i,],
      "-"
    )
    dist2 <- rowSums(diff^2)
    energy[i] <-
      0.5*
      sum(
        W_source[i,nei]*
          dist2
      )
    
  }
  if(normalize){
    
    energy <-
      energy/(deg+1e-10)
    
  }
  return(energy)
}


ComputeCellwiseAlpha <- function(
    graph_list,
    embedding_list,
    lambda_self = 1,
    beta = 1,
    normalize_alpha = TRUE
){
  
  cat("calculating cell-specific alpha\n")
  
  layers <- names(graph_list)
  n <- nrow(graph_list[[1]])
  M <- length(layers)
  
  
  ################################
  # self energy
  ################################
  
  energy_self <- lapply(
    layers,
    function(m){
      
      ComputeNodeDirichletEnergy(
        W = graph_list[[m]],
        embedding = embedding_list[[m]]
      )
      
    }
  )
  
  names(energy_self) <- layers
  
  
  energy_self_norm <- lapply(
    energy_self,
    NormalizeEnergy
  )
  
  
  ################################
  # cross energy
  ################################
  
  energy_cross <- vector(
    "list",
    M
  )
  
  energy_cross_norm <- vector(
    "list",
    M
  )
  
  names(energy_cross) <- layers
  names(energy_cross_norm) <- layers
  
  
  for(source in layers){
    
    energy_cross[[source]] <- list()
    energy_cross_norm[[source]] <- list()
    
    for(target in layers){
      
      if(source == target)
        next
      
      cross_energy <-
        ComputeCrossDirichletEnergy(
          W_source =
            graph_list[[source]],
          embedding_target =
            embedding_list[[target]]
        )
      
      energy_cross[[source]][[target]] <-
        cross_energy
      
      energy_cross_norm[[source]][[target]] <-
        NormalizeEnergy(
          cross_energy
        )
    }
  }
  
  
  ################################
  # pairwise score
  ################################
  
  pairwise_score <- vector(
    "list",
    M
  )
  
  names(pairwise_score) <- layers
  
  
  for(source in layers){
    
    pairwise_score[[source]] <- list()
    
    for(target in layers){
      
      if(source == target)
        next
      
      pairwise_score[[source]][[target]] <-
        energy_cross_norm[[source]][[target]] -
        lambda_self *
        energy_self_norm[[target]]
    }
  }
  
  
  ################################
  # pairwise score -> alpha score
  ################################
  
  alpha_score <- matrix(
    0,
    nrow = n,
    ncol = M,
    dimnames = list(
      NULL,
      layers
    )
  )
  
  
  for(source in layers){
    
    for(target in layers){
      
      if(source == target)
        next
      
      alpha_score[, target] <-
        alpha_score[, target] +
        pairwise_score[[source]][[target]]
    }
  }
  
  
  ################################
  # alpha
  ################################
  
  alpha <- exp(
    beta * alpha_score
  )
  
  alpha <-
    alpha /
    rowSums(alpha)
  
  
  ################################
  # minimum alpha
  ################################
  
  if(normalize_alpha){
    
    alpha[
      alpha < 1e-4
    ] <- 1e-4
    
    alpha <-
      alpha /
      rowSums(alpha)
  }
  
  alpha_layerwise <- colMeans(
    alpha,
    na.rm = TRUE
  )
  
  ################################
  # return
  ################################
  
  return(
    list(
      
      alpha =
        alpha,
      
      alpha_score =
        alpha_score,
      
      pairwise_score =
        pairwise_score,
      
      self_energy =
        energy_self,
      
      self_energy_norm =
        energy_self_norm,
      
      cross_energy =
        energy_cross,
      
      cross_energy_norm =
        energy_cross_norm,
      
      alpha_layerwise = 
        alpha_layerwise
    )
  )
}





# 4.SCML_core -------------------------------------------------------------


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


run_scml_core_pan_cell_specificalpha <- function(
    graph_list,
    ndim = 30,
    alpha = NULL,
    layer_weight = NULL,
    return_laplacians = TRUE
){
  
  assert_named_list(graph_list, "graph_list")
  
  layer_names <- names(graph_list)
  
  
  ####################################
  ## check graph dimensions
  ####################################
  
  n_vertices <- vapply(
    graph_list,
    nrow,
    integer(1)
  )
  
  if(
    length(unique(n_vertices)) != 1 ||
    any(vapply(graph_list,ncol,integer(1)) != n_vertices[1])
  ){
    stop(
      "All graphs must be square matrices with identical dimensions.",
      call.=FALSE
    )
  }
  
  
  n <- n_vertices[1]
  
  
  ####################################
  ## alpha parsing
  ####################################
  
  alpha_mode <- "global"
  
  
  if(is.null(alpha)){
    
    alpha_vec <- 
      resolve_layer_values_pan(
        NULL,
        layer_names,
        default_value = 1,
        arg="alpha"
      )
    
  }else{
    
    
    ################################
    ## cell-specific alpha
    ################################
    
    if(
      is.matrix(alpha) ||
      is.data.frame(alpha)
    ){
      
      alpha_cell <- as.matrix(alpha)
      
      if(
        nrow(alpha_cell)!=n
      ){
        stop(
          "alpha_cell must have rows equal to number of cells."
        )
      }
      
      if(
        !all(layer_names %in% colnames(alpha_cell))
      ){
        stop(
          "alpha_cell columns must match graph_list names."
        )
      }
      
      alpha_cell <-
        alpha_cell[,layer_names,drop=FALSE]
      
      
      ################################
      ## normalize per cell
      ################################
      
      alpha_cell <- alpha_cell 
      
      
      alpha_mode <- "cell_specific"
      
      alpha_vec <- NULL
      
    }else{
      
      alpha_vec <-
        resolve_layer_values_pan(
          alpha,
          layer_names,
          default_value=1,
          arg="alpha"
        )
    }
  }
  
  
  
  ####################################
  ## layer weight
  ####################################
  
  weight_vec <-
    resolve_layer_values_pan(
      layer_weight,
      layer_names,
      default_value=1,
      arg="layer_weight"
    )
  
  
  cat("\n========== Layer Parameters ==========\n")
  
  
  if(alpha_mode=="global"){
    
    print(
      data.frame(
        layer=layer_names,
        weight=as.numeric(weight_vec),
        alpha=as.numeric(alpha_vec)
      )
    )
    
  }else{
    
    cat(
      "Using cell-specific alpha matrix:",
      nrow(alpha_cell),
      "cells x",
      ncol(alpha_cell),
      "layers\n"
    )
    
    print(
      head(alpha_cell)
    )
  }
  
  
  
  ####################################
  ## initialize
  ####################################
  
  
  L_sum <-
    Matrix(
      0,
      n,
      n,
      sparse=TRUE
    )
  
  
  P_sum <-
    Matrix(
      0,
      n,
      n,
      sparse=TRUE
    )
  
  
  U_layers <- list()
  eig_layers <- list()
  laplacians <- list()
  
  
  
  ####################################
  ## layer loop
  ####################################
  
  
  for(layer in layer_names){
    
    
    cat(
      "Processing layer:",
      layer,
      "\n"
    )
    
    
    L <-
      compute_normalized_laplacian(
        graph_list[[layer]]
      )
    
    
    subspace <-
      get_laplacian_subspace(
        L,
        ndim=ndim
      )
    
    
    U_layers[[layer]] <-
      subspace$U
    
    eig_layers[[layer]] <-
      subspace$values
    
    
    if(return_laplacians){
      laplacians[[layer]] <- L
    }
    
    
    ################################
    ## Laplacian fusion
    ################################
    
    L_sum <-
      L_sum +
      weight_vec[layer] *
      L
    
    
    ################################
    ## Projection fusion
    ################################
    
    
    P_layer <-
      tcrossprod(
        subspace$U
      )
    
    
    if(alpha_mode=="global"){
      
      P_sum <-
        P_sum +
        alpha_vec[layer] *
        P_layer
      
      
    }else{
      
      ################################
      ## cell-specific alpha
      ################################
      
      a <-
        alpha_cell[,layer]
      
      D_half <-
        Matrix::Diagonal(
          x=sqrt(a)
        )
      
      
      P_sum <-
        P_sum +
        D_half %*%
        P_layer %*%
        D_half
      
    }
    
  }
  
  
  
  ####################################
  ## SCML eigenspace
  ####################################
  
  
  L_modified <-
    L_sum - P_sum
  
  
  
  eig_consensus <-
    RSpectra::eigs_sym(
      L_modified,
      k=ndim,
      which="SA"
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
    
    U=U_consensus,
    
    eigenvalues=
      eig_consensus$values,
    
    U_layers=
      U_layers,
    
    eigenvalues_layers=
      eig_layers,
    
    L_modified=
      L_modified,
    
    laplacians=
      if(return_laplacians)
        laplacians
    else
      NULL,
    
    alpha=
      if(alpha_mode=="global")
        alpha_vec
    else
      alpha_cell,
    
    layer_weight=
      weight_vec,
    
    alpha_mode=
      alpha_mode
  )
}
