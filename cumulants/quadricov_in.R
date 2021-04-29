if(!require("Matrix")){
  install.packages("Matrix")
  stopifnot(require("Matrix"))
}
if(!require("pracma")){
  install.packages("pracma")
  stopifnot(require("pracma"))
}

quadricov_in <- function(data){
  p <- nrow(data)
  n <- ncol(data)
  C <- zeros(p)
  Q <- zeros(p*p)
  temp <- c()
  
  for (i in 1:n){
    for (a in 1:p) {
      for (b in 1:p) {
        tt <- data[a,i] * data[b,i]
        C[a,b] <- C[a,b] + tt
        temp[a+(b-1)*p] <- tt
      }
    }
    
    for (a in 1:(p*p)) {
      for (b in 1:(p*p)) {
        if (b >= a) {
          # debug: print(c(a,b))
          Q[a,b] <- Q[a,b] + temp[a] * temp[b]
        }
        
      }
    }
  }
  
  C <- C / n
  Q <- Q / n
  
  for (a in 1:p) {
    for (b in 1:p) {
      for (c in 1:p) {
        for (d in 1:p) {
          if ((c + (d-1)*p) >= (a + (b-1)*p)) {
            Q[(a + (b-1)*p),(c + (d-1)*p)] <- Q[(a + (b-1)*p),(c + (d-1)*p)] - C[a,b]*C[c,d] - C[a,c]*C[b,d] - C[a,d]*C[c,d]
          }
        }
      }
    }
  }
  return(Q)
}