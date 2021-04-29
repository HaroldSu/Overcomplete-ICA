if(!require("R.matlab")){
  install.packages("R.matlab")
  stopifnot(require("R.matlab"))
}

extract_basis <- function(Es, k){
  Matlab$startServer(port = 9999)
  matlab <- Matlab(port = 9999)
  open(matlab)
  setVariable(matlab, Es = Es)
  evaluate(matlab, "[Q,~] = qr(Es);")
  matlab_result <- getVariable(matlab, "Q")
  Q <- matlab_result$Q
  close(matlab)
  # qr_result <- qr(Es)
  # Q <- qr.Q(qr_result)
  Esbasis <- Q[,1:k]
  Fsbasis <- Q[,(k+1):ncol(Q)]
  output <- list()
  output$Esbasis <- Esbasis
  output$Fsbasis <- Fsbasis
  return(output)
}