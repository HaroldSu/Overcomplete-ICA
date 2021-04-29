if(!require("pracma")){
  install.packages("pracma")
  stopifnot(require("pracma"))
}

sdp_cluster <- function(Hs, k, ctype = "h"){
  p <- sqrt(nrow(Hs))
  
  nclust <- 3*k
  
  extract_basis_result <- extract_basis(Hs, k)
  Fs <- extract_basis_result$Fsbasis
  Dss <- zeros(p^2, nclust)
  
  for(irep in 1:nclust){
    u <- randn(p, 1)
    u <- u / norm(u, "F")
    G <- u %*% t(u)
    D <- majorize_minimize(G, Fs)
    Dss[,irep] <- as.vector(D)
  }
  
  Ds_est <- cluster_Dss(Dss, k, ctype)
  ds_est <- approx_ds_from_Ds(Ds_est)
  
  output <- list()
  output$Ds_est <- Ds_est
  output$ds_est <- ds_est
  return(output)
}