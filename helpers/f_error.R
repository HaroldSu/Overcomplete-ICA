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

f_error <- function(ds_est, ds){
  k <- ncol(ds_est)
  signs <- c()
  for (i in 1:k) {
    signi = as.numeric(sign(t(ds_est[,i]) %*% ds[,i]))
    signs[i] <- signi
    ds_est[,i] <- ds_est[,i] * signi
  }
  perf <- (norm(ds_est - ds, "F")^2) / (norm(ds, "F")^2)
  return(perf)
}