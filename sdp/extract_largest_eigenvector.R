extract_largest_eigenvector <- function(D){
  D <- (D + t(D)) / 2
  eigen_result <- eigen(D)
  u <- eigen_result$vectors
  e <- eigen_result$values
  a <- max(abs(e))
  b <- which.max(abs(e))
  u <- Re(u[,b])
  e <- a
  output <- list()
  output$u <- u
  output$e <- e
  return(output)
}
