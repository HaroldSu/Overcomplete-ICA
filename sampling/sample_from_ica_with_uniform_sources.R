if(!require("Matrix")){
  install.packages("Matrix")
  stopifnot(require("Matrix"))
}
if(!require("pracma")){
  install.packages("pracma")
  stopifnot(require("pracma"))
}

sample_from_ica_with_uniform_sources <- function(ds, N){
  ## ds: mxing matrix
  ## N: the number of observations, n = 1, 2, ..., N
  ##
  ## model: x[n] <- D * alpha[n]
  ##        alpha[n] ~ uniform(0.5, 1)
  
  
  p <- nrow(ds)
  k <- ncol(ds)
  
  Alpha <- rand(k,N) * abs(randn(k, N))
  X <- Matrix(ds %*% Alpha, sparse = T)
  output <- list()
  output$Alpha <- Alpha
  output$X <- X
  return(output)
}

