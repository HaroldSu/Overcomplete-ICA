############
## Figure 3
############
source("./sampling/sample_mixing_matrix.R")
source("./sampling/sample_from_ica_with_uniform_sources.R")
source("./cumulants/quadricov_in.R")
source("./cumulants/quadricov.R")
source("./cumulants/genquadricov.R")
source("./sdp/extract_basis.R")
source("./sdp/extract_largest_eigenvector.R")
source("./sdp/cluster_Dss.R")
source("./sdp/majorize_minimize.R")
source("./sdp/sdp_cluster.R")
source("./sdp/sdp_semiada.R")
source("./sdp/adaptive_deflation.R")
source("./sdp/solve_relaxation_mezcal_approx_fista.R")
source("./sdp/proj_simplex.R")
source("./helpers/approx_ds_from_Ds.R")
source("./helpers/a_error.R")
source("./helpers/f_error.R")
source("./comparisons/fourier_pca.R")
## Fix the latent dimension k
k <- 30
## Fix the observed dimension p
p <- 15
## Fix the number of replicates nrep
nrep <- 10


## Vary the latent dimension from 5 to 60
output_n <- sapply(c(seq(1000, 10000, by = 1000), seq(10000, 210000, by = 10000)), function(ns){
  output <- sapply(1:nrep, function(replicate){
    ds <- sample_mixing_matrix(p, k)
    sample_from_ica_with_uniform_sources_result <- sample_from_ica_with_uniform_sources(ds, ns)
    X <- sample_from_ica_with_uniform_sources_result$X
    X <- as.matrix(X)
    C <- quadricov(X)
    svd_result <- svd(C)
    CU <- svd_result$u
    Hs <- CU[,1:k]
    ## The ds estimates
    ds_est_oica_semiada <- sdp_semiada(Hs, k)
    ds_est_fpca <- fourier_pca(X, k)
    ds_est_rand <- sample_mixing_matrix(p,k)
    ## A-Error
    a_error_oica_semiada <- a_error(ds_est_oica_semiada, ds)
    a_error_fpca <- a_error(ds_est_fpca, ds)
    a_error_rand <- a_error(ds_est_rand, ds)
    ## F-Error
    f_error_oica_semiada <- f_error(ds_est_oica_semiada, ds)
    f_error_fpca <- f_error(ds_est_fpca, ds)
    f_error_rand <- f_error(ds_est_rand, ds)
    ## Perfect Recovery
    perf_recovery_oica_semiada <- evaluation_recovery(ds_est_oica_semiada, ds, th = 8)
    perf_recovery_rand <- evaluation_recovery(ds_est_rand, ds, th = 8)
    return(c(a_error_oica_semiada, a_error_rand, a_error_fpca, 
             f_error_oica_semiada, f_error_rand, f_error_fpca))
  })
  output <- t(as.matrix(output))
  output <- colMeans(output)
  return(output)
})

output_n <- t(as.matrix(output_n))
colnames(output_n) <- c("a_error_oica_semiada", "a_error_rand", 'a_error_fpca',
                        "f_error_oica_semiada", "f_error_rand", "f_error_fpca")
output_n <- as.data.frame(output_n)

## A-Error-1
plot(y = output_n$a_error_oica_semiada[1:10], x = seq(1000, 10000, by = 1000),
     col = "red", type = "b",
     ylab = "sample size", ylab = 'A-Error', ylim = c(0,1))
points(y = output_n$a_error_rand[1:10], x = seq(1000, 10000, by = 1000),
       col = "orange", type = "b")
points(y = output_n$a_error_fpca[1:10], x = seq(1000, 10000, by = 1000),
       col = "blue", type = "b")

## A-Error-2
plot(y = output_n$a_error_oica_semiada[11:31], x = seq(10000, 210000, by = 10000),
     col = "red", type = "b",
     ylab = "sample size", ylab = 'A-Error', ylim = c(0,1))
points(y = output_n$a_error_rand[11:31], x = seq(10000, 210000, by = 10000),
       col = "orange", type = "b")
points(y = output_n$a_error_fpca[11:31], x = seq(10000, 210000, by = 10000),
       col = "blue", type = "b")

## F-Error-1
plot(y = output_n$f_error_oica_semiada[1:10], x = seq(1000, 10000, by = 1000),
     col = "red", type = "b",
     ylab = "sample size", ylab = "F-Error", ylim = c(0,1))
points(y = output_n$f_error_rand[1:10], x = seq(1000, 10000, by = 1000),
       col = "orange", type = "b")
points(y = output_n$f_error_fpca[1:10], x = seq(1000, 10000, by = 1000),
       col = "blue", type = "b")

## F-Error-2
plot(y = output_n$f_error_oica_semiada[11:31], x = seq(10000, 210000, by = 10000),
     col = "red", type = "b",
     ylab = "sample size", ylab = 'F-Error', ylim = c(0,1))
points(y = output_n$f_error_rand[11:31], x = seq(10000, 210000, by = 10000),
       col = "orange", type = "b")
points(y = output_n$f_error_rand[11:31], x = seq(10000, 210000, by = 10000),
       col = "blue", type = "b")
