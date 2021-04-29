sdp_adaptive <- function(Hs, k){
  extract_basis_result <- extract_basis(Hs, k)
  Fs <- extract_basis_result$Fsbasis
  Ds_est <- adaptive_deflation(Fs, k)
  ds_est <- approx_ds_from_Ds(Ds_est)
  
  return(Ds_est)
}