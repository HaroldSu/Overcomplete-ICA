if(!require("pracma")){
  install.packages("pracma")
  stopifnot(require("pracma"))
}

sample_orthogonal_matrix <- function(k){
  X <- matrix(runif(k^2), nrow = k, ncol = k)
  V <- orth(X)
  return(V)
}