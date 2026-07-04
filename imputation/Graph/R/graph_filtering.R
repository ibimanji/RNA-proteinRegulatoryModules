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
  
  eig <- RSpectra::eigs_sym(
    L,
    k=n_eigs,
    which="SM"
  )
  
  list(
    lambda=eig$values,
    V=eig$vectors
  )
}

GraphFilterKernel <- function(
    lambda,
    type = c(
      "heat",
      "tikhonov",
      "exponential",
      "bandpass",
      "ideal",
      "butterworth"
    ),
    beta = 1,
    lambda_low = 0.1,
    lambda_high = 1,
    lambda_c = 0.4,
    order = 4
){
  
  type <- match.arg(type)
  
  g <- switch(
    
    type,
    
    ########################################################
    ## Heat kernel
    ########################################################
    
    heat =
      exp(
        -beta * lambda
      ),
    
    ########################################################
    ## Tikhonov
    ########################################################
    
    tikhonov =
      1 /
      (
        1 +
          beta * lambda
      ),
    
    ########################################################
    ## Exponential
    ########################################################
    
    exponential =
      exp(
        -beta * sqrt(lambda)
      ),
    
    ########################################################
    ## Band-pass
    ########################################################
    
    bandpass =
      as.numeric(
        lambda > lambda_low &
          lambda < lambda_high
      ),
    
    ########################################################
    ## Ideal low-pass
    ########################################################
    
    ideal =
      as.numeric(
        lambda < beta
      ),
    
    ########################################################
    ## Butterworth low-pass
    ########################################################
    
    butterworth =
      1 /
      (
        1 +
          (lambda / lambda_c)^(2 * order)
      ),
    
  )
  
  ##########################################################
  ## Normalize
  ##########################################################
  
  g <- g / max(g)
  
  g
  
}
GraphSpectralDenoise <- function(
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
