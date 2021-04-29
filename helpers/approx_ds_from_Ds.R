if(!require("Matrix")){
  install.packages("Matrix")
  stopifnot(require("Matrix"))
}
if(!require("pracma")){
  install.packages("pracma")
  stopifnot(require("pracma"))
}

approx_ds_from_Ds <- function(Ds){
  k <- ncol(Ds)
  p <- sqrt(nrow(Ds))
  ds <- zeros(p,k)
  eigmaxes <- zeros(k,1)
  
  for (i in 1:k) {
    D <- matrix(Ds[,i], nrow = p, ncol = p)
    extract_largest_eigenvector_result <- extract_largest_eigenvector(D)
    eigmaxes[i] <- extract_largest_eigenvector_result$e
    ds[,i] <- extract_largest_eigenvector_result$u
  }
  return(ds)
}