proj_simplex <- function(v){
  v <- (v > 0) * v
  u <- sort(v, decreasing = T)
  sv <- cumsum(u)
  
  rho <- which(u > ((sv - 1) / (1:length(u))))
  rho <- max(rho)
  theta <- max(0, (sv[rho] - 1) / rho)
  w <- v - theta
  for (k in 1:length(w))  w[k] <- max(w[k], 0)
  return(w)
}