if(!require("pracma")){
  install.packages("pracma")
  stopifnot(require("pracma"))
}


adaptive_deflation <- function(Fs, k, Ds_est){
  p <- sqrt(nrow(Fs))
  if (ncol(Ds_est) < k){
    kloc <- k - ncol(Ds_est)
    Fsloc <- cbind(Fs, Ds_est)
    
    for(i in 1:kloc){
      Matlab$startServer(port = 9999)
      matlab <- Matlab(port = 9999)
      open(matlab)
      setVariable(matlab, Fsloc = Fsloc)
      evaluate(matlab, "[uuu,~,~] = svd(Fsloc);")
      matlab_result <- getVariable(matlab, "uuu")
      uuu <- matlab_result$uuu
      close(matlab)
      # svd_result <- svd(Fsloc)
      # uuu <- svd_result$u
      Fsloc <- uuu[,1:(ncol(uuu) - (kloc - i + 1))]
      Esloc <- uuu[,(ncol(uuu) - (kloc - i)):ncol(uuu)]
      if(i == kloc){
        g <- Esloc
      }else{
        g <- Esloc[,1]
      }
      G <- matrix(g, nrow = p, ncol = p)
      D <- majorize_minimize(G, Fsloc)
      Ds_est <- cbind(Ds_est, as.vector(D))
      Fsloc <- cbind(Fsloc, as.vector(D))
    }
  }
  return(Ds_est)
}