if(!require("pracma")){
  install.packages("pracma")
  stopifnot(require("pracma"))
}
if(!require("R.matlab")){
  install.packages("R.matlab")
  stopifnot(require("R.matlab"))
}

fourier_pca <- function(X, k){
  p <- nrow(X)
  
  u <- randn(p, 1)
  u <- u / norm(u, "2")
  
  Q <- genquadricov(X, u)
  Q1 <- Re(Q)
  Q2 <- Im(Q)
  
  Matlab$startServer(port = 9999)
  matlab <- Matlab(port = 9999)
  open(matlab)
  setVariable(matlab, Q1 = Q1)
  evaluate(matlab, "[W,S,~] = svd(Q1);")
  matlab_result <- getVariable(matlab, c("W", "S"))
  W <- matlab_result$W
  S <- matlab_result$S
  close(matlab)
  
  inds <- order(diag(S), decreasing = T)
  W <- W[,inds]
  W <- W[,1:k]
  
  q1 <- t(W) %*% Q1 %*% W
  q2 <- t(W) %*% Q2 %*% W
  
  Matlab$startServer(port = 9999)
  matlab <- Matlab(port = 9999)
  open(matlab)
  setVariable(matlab, q1 = q1, q2 = q2)
  evaluate(matlab, "M = q1/q2;", "[V,~] = eig(M);")
  matlab_result <- getVariable(matlab, "V")
  V <- matlab_result$V
  close(matlab)
  
  C <- W %*% V
  
  for (j in 1:k) {
    c <- C[,j]
    
    a <- Re(c)
    b <- Im(c)
    
    theta <- atan(-(2 * sum(a * b)) / sum(a^2 - b^2)) / 2
    
    while(theta < 0) theta <- theta + pi
    while (theta > 2 *pi) theta <- theta - pi
    
    temp <- Re(exp(1i * theta) * c)
    C[,j] <- temp / norm(temp, "2")
  }
  
  ds <- zeros(p,k)
  
  for (j in 1:k) {
    cc <- C[,j]
    Matlab$startServer(port = 9999)
    matlab <- Matlab(port = 9999)
    open(matlab)
    setVariable(matlab, cc = c)
    evaluate(matlab, "p = sqrt(length(cc));","[v,s,~] = svd( reshape(cc,p,p) );")
    matlab_result <- getVariable(matlab, c("v", "s"))
    v <- matlab_result$v
    s <- matlab_result$s
    close(matlab)
    inds <- order(diag(s), decreasing = T)
    v <- v[,inds]
    ds[,j] <- v[,1]
  }
  
  return(ds)
}