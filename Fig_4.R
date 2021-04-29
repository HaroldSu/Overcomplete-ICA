data <- readMat("data_batch_1.mat")
data_matrix <- data$data
p <- 7
sr <- 32
sc <- 32
N <- 10000

d <- floor(p/2) + mod(p, 2) - 1
step <- (sr - 2*d - 1) * (sc - 2*d - 1)
X <- zeros(p*p, N*step)
ind <- 1

img_array <- array(dim = c(32,32,1,3))

for(n in 1:N){
  data_row <- data_matrix[n,]
  img_array[,,1,1] <- matrix(data_row[1:1024], nrow = 32, ncol = 32)
  img_array[,,1,2] <- matrix(data_row[1025:2048], nrow = 32, ncol = 32)
  img_array[,,1,3] <- matrix(data_row[2049:3072], nrow = 32, ncol = 32)
  # plot(cimg(img_array))
  img_gray_array <- grayscale(cimg(img_array))
  img_gray_matrix <- as.matrix(img_gray_array)
  
  #################
  ## Matlab part ##
  #################
  
  Matlab$startServer(port = 9999)
  matlab <- Matlab(port = 9999)
  open(matlab)
  setVariable(matlab, gg = img_gray_matrix)
  evaluate(matlab, "pp = nlfilter(gg, [7 7], @(block) {block});")
  matlab_result <- getVariable(matlab, "pp")
  pp <- matlab_result$pp
  close(matlab)
  
  for(i in 1:32){
    for (j in 1:32) {
      patch <- pp[[(i + (j - 1)*32)]]
      X[,ind] <- as.vector(patch)
      ind <- ind + 1
    }
  }
}

XC <- X - repmat(as.matrix(rowMeans(X)), 1, ncol(X))

k <- 150

C <- quadricov(X)
svd_result <- svd(C)
CU <- svd_result$u
Hs <- CU[,1:k]
ds_est <- sdp_semiada(Hs, k)

## Here we can obtain the estimated mixing matrix. 
## However, due to the limit of R, we can not plot the atoms out. 

