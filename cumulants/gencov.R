if(!require("Matrix")){
  install.packages("Matrix")
  stopifnot(require("Matrix"))
}


gencov <- function(X, omega){
  p <- nrow(X)
  n <- ncol(X)
  
  if(length(omega) == 1){
    omega <- omega * rep(1, p) / p
  }
  
  proj <- t(X) %*% omega
  eproj <- exp(proj)
  Eomega <- (X %*% eproj) / sum(eproj)
  
  C <- X %*% sparseMatrix(i = 1:n, j = 1:n, x = eproj) %*% t(X)
  C <- C / sum(eproj)
  C <- C - Eomega %*% t(Eomega)
  return(C)
}