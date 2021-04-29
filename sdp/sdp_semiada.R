sdp_semiada <- function(Hs, k, ctype){
  sdp_cluster_result <- sdp_cluster(Hs, k, ctype = "h")
  Ds_clust <- sdp_cluster_result$Ds_est
  
  G <- abs(t(Ds_clust) %*% Ds_clust) - eye(k)
  G_col_max <- c()
  for(j in 1:ncol(G)){G_col_max[j] <- max(G[,j])}
  mind <- which.min(G_col_max)
  Ds_est <- Ds_clust[,mind]
  Ds_est <- matrix(Ds_est, ncol = 1)
  
  ind <- 2
  for (i in 1:k) {
    if(i == mind){next}
    
    if(max(G[,ind]) < 0.8){
      Ds_est <- cbind(Ds_est, Ds_clust[,i])
      ind <- ind + 1
    }
  }
  
  extract_basis_result <- extract_basis(Hs, k)
  Fs <- extract_basis_result$Fsbasis
  Ds_est <- adaptive_deflation(Fs, k, Ds_est)
  ds_est <- approx_ds_from_Ds(Ds_est)
  
  #output <- list()
  #output$ds_est <- ds_est
  #output$Ds_est <- Ds_est
  return(ds_est)
}