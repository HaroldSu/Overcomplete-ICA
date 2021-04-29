sample_mixing_matrix <- function(p, k){
  ds <- matrix(rnorm(p*k), nrow = p, ncol = k)
  
  for(i in 1:k){
    ds[,i] <- ds[,i] / norm(ds[,i], "2")
  }
  return(ds)
}