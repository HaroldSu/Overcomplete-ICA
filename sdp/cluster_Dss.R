if(!require("pracma")){
  install.packages("pracma")
  stopifnot(require("pracma"))
}

cluster_Dss <- function(Dss, k, ctype){
  D1 <- Dss[,1]
  
  for (i in 2:ncol(Dss)) {
    Di <- Dss[,i]
    Dss[,i] <- as.numeric(sign(t(D1) %*% Di)) * Di
  }
  
  Ds_est <- extract_clusters(DD = Dss, nclust = k, ctype = ctype)
  return(Ds_est)
}

extract_clusters <- function(DD, nclust, ctype){
  p <- sqrt(nrow(DD))
  
  if (ctype == "h"){
    hc <- hclust(dist(t(DD)))
    cc <- cutree(hc, k = nclust)
  }
  
  if(ctype == "km"){
    km <- kmeans(t(DD), nclust)
    cc <- km$cluster
  }
  
  Ds_temp <- zeros(round(p^2), nclust)
  
  for (i in 1:nclust) {
    DDi <- as.matrix(DD[,(cc == i)])
    Ds_temp[,i] <- DDi[,1]
  }
  return(Ds_temp)
}