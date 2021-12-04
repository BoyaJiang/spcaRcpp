#devtools::use_package("testthat")

#context("Sparse PCA")

#Load library
library(spcaRcpp)     

#Set seed
set.seed(322)

#--------------------------------------------------------------------
# Generate Some Data in R
#--------------------------------------------------------------------
m = 10000
V1 = rnorm(m, 0, 290)
V2 = rnorm(m, 0, 300)
V3 = -0.1*V1 + 0.1*V2 + rnorm(m,0,100)

X = cbind(V1,V1,V1,V1, V2,V2,V2,V2, V3,V3)
X = X + matrix(rnorm(length(X),0,1), ncol = ncol(X), nrow = nrow(X))


#*************************************************************************************
# Test: spcaRcpp - center = TRUE
#*************************************************************************************

pca_out <- prcomp(X, rank = 3,center = TRUE,) #Sparse PCA
spca_out <- spcaRcpp(X, k=3, alpha=0, beta=0, center = TRUE) #Sparse PCA with Rcpp

#Test1: SPCA recovers PCA for alpha = 0 and beta = 0.
testthat::test_that("Test cov; alpha = 0 and beta = 0", {
  testthat::expect_equal(pca_out$sdev[1:3], spca_out$sdev[1:3])
  testthat::expect_equal(sum(diag(1, 3, 3) - t(spca_out$loadings)%*%spca_out$loadings), 0 )
})


#*************************************************************************************
# Test: spcaRcpp - center = FALSE, k = NULL
#*************************************************************************************

pca_out <- prcomp(X, center = FALSE,) #Sparse PCA
spca_out <- spcaRcpp(X, alpha=0, beta=0, center = FALSE) #Sparse PCA with Rcpp

#Test2: SPCA recovers PCA for center = FALSE, k = NULL.
testthat::test_that("Test cov; alpha = 0 and beta = 0", {
  testthat::expect_equal(pca_out$sdev, spca_out$sdev)
  testthat::expect_equal(sum(diag(1, 10, 10) - t(spca_out$loadings)%*%spca_out$loadings), 0 )
})

#*************************************************************************************
# Test: spcaRcpp - k > min(n, p)
#*************************************************************************************

pca_out <- prcomp(X,center = FALSE,) #Sparse PCA
spca_out <- spcaRcpp(X, k = 11, alpha=0, beta=0, center = FALSE) #Sparse PCA with Rcpp

#Test3: SPCA recovers PCA for k > min(n, p).
testthat::test_that("Test cov; alpha = 0 and beta = 0", {
  testthat::expect_equal(pca_out$sdev, spca_out$sdev)
  testthat::expect_equal(sum(diag(1,10,10) - t(spca_out$loadings)%*%spca_out$loadings), 0 )
})
