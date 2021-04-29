if(!require("pracma")){
  install.packages("pracma")
  stopifnot(require("pracma"))
}

genquadricov <- function(X, u){
######################################################################
# X: p x N, p dimension, N number of samples
# u: p, processing point
# Q: p*p x p*p, matricized 4-th order generalized cumulant
#
# The matricization rule:
# Q( (i1-1)*p + i2, (i3-1)*p + i4 ) = CUM(i1,i2,i3,i4);
# (i.e. row-wise for the the first two and the last two dimensions)
#
# THE OUTPUT IS A COMPLEX MATRIX
######################################################################
  p <- nrow(X)
  N <- ncol(X)
  
  Xu <- t(X) %*% u
  expXu <- exp(1i * Xu)
  M0 <- sum(expXu) / N
  
  M1 <- (X %*% expXu) / (N*M0)
  
  XC <- X - repmat(M1, 1, N)
  
  M2 <- matrix(rep(0, p*p), nrow = p, ncol = p)
  M4 <- matrix(rep(0, p^4), nrow = p*p, ncol = p*p)
  temp <- matrix(rep(0, p*p), nrow = p, ncol = p)
  
  for (n in 1:N) {
    xn <- XC[,n]
    temp <- xn %*% t(xn)
    M2 <- M2 + expXu[n] * temp
    M4 <- M4 + (expXu[n] * as.vector(temp)) %*% t(expXu[n] * as.vector(temp))
  }
  
  M2 <- M2 / (N*M0)
  M4 <- M4 / (N*M0)
  Q <- matrix(rep(0, p^4), nrow = p*p, ncol = p*p)
  
  for (i3 in 1:p){
    for (i4 in 1:p) {
      icol <- (i3 - 1) * p + i4
      temp <- matrix(M4[,icol], ncol = p, nrow = p, byrow = F)
      temp <- temp - M2[i3, i4] * t(M2) - M2[,i3] %*% t(M2[i4,]) - M2[,i4] %*% t(M2[i3,])
      Q[,icol] <- as.vector(temp)
    }
  }
  
  return(Q)
}

