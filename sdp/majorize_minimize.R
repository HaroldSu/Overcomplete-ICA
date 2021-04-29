if(!require("pracma")){
  install.packages("pracma")
  stopifnot(require("pracma"))
}

majorize_minimize <- function(G, Fs){
  p <- sqrt(nrow(Fs))
  Dinit <- eye(p) / p
  D <- Dinit
  mu <- 5
  maxiter <- 100
  tolerance <- 0.001
  nmmmax <- 100
  
  iter <- 1
  while (norm(D, "F") < 1) {
    u_result <- extract_largest_eigenvector(G)
    u <- as.matrix(u_result$u)
    Ginit <- u %*% t(u)
    D <- solve_relaxation_mezcal_approx_fista(Fs, Ginit, mu, Dinit, maxiter, tolerance)
    iter <- iter + 1
    if(iter > nmmmax){break}
  }
  
  return(D)
}