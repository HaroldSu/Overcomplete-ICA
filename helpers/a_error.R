if(!require("Matrix")){
  install.packages("Matrix")
  stopifnot(require("Matrix"))
}
if(!require("pracma")){
  install.packages("pracma")
  stopifnot(require("pracma"))
}
if(!require("RcppHungarian")){
  install.packages("RcppHungarian")
  stopifnot(require("RcppHungarian"))
}

a_error <- function(ds_est, ds){
  k <- ncol(ds)
  F_matrix <- zeros(k)
  
  for (i in 1:k) {
    for (j in 1:k) {
      cosij <- abs(t(ds_est[,i]) %*% ds[,j]) / (norm(ds_est[,i], "2") * norm(ds[,j], "2"))
      loc <- acos(cosij)
      F_matrix[i,j] <- Re(loc)
    }
  }
  
  HungarianBipartiteMatching <- HungarianSolver(F_matrix)
  cost <- HungarianBipartiteMatching$cost
  perf <- 2 * cost / pi
  perf <- perf / k
  return(perf)
}
