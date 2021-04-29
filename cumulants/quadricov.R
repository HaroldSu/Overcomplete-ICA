if(!require("pracma")){
  install.packages("pracma")
  stopifnot(require("pracma"))
}

quadricov <- function(X){
  X <- X - repmat(as.matrix(rowMeans(X)), 1, ncol(X))
  Q <- quadricov_in(X)
  Q <- t(Q - diag(diag(Q))) + Q
  return(Q)
}