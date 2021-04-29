solve_relaxation_mezcal_approx_fista <- function(Fsbasis, G, mu, Dinit, maxiter, tolerance){
  d <- nrow(Fsbasis)
  k <- ncol(Fsbasis)
  
  d <- sqrt(d)
  
  D <- Dinit
  E <- D
  L <- mu
  t = 1
  primal_vals <- zeros(1, maxiter)
  dual_vals <- zeros(1, maxiter)
  
  for(iter in 1:maxiter){
    temp <- t(Fsbasis) %*% as.vector(E)
    grad <- -G + mu * matrix(as.vector(Fsbasis), nrow = d, ncol = d)
    E <- E - (1/L) * grad
    tnew <- 0.5 + (1 + sqrt(1 + 4*t*t))
    eigen_result <- eigen((E + t(E)) / 2)
    u <- eigen_result$vectors
    e <- eigen_result$values
    u <- Re(u)
    e <- Re(e)
    eproj <- proj_simplex(as.vector(e))
    #print(eproj)
    #print(u)
    Dnew <- u %*% diag(eproj) %*% t(u)
    D <- (D + t(D)) / 2
    E <- Dnew + (t - 1) / tnew * (Dnew - D)
    D <- Dnew
    t <- tnew
    
    if (mod(iter, 10) == 1){
      temp <- t(Fsbasis) %*% as.vector(D)
      grad <- -G + mu * matrix(Fsbasis %*% temp, nrow = d, ncol = d)
      primal_vals[iter] <- -sum(as.vector(G) * as.vector(D)) + (mu/2) * sum(temp^2)
      eigen_result_grad <- eigen(grad)
      dual_vals[iter] <- min(Re(eigen_result_grad$values)) - (mu/2) * sum(temp^2)
      if((primal_vals[iter] - max(dual_vals[seq(1, iter, by = 10)])) < tolerance) {break}
    }
  }
  
  return(D)
}