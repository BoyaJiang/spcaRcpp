#devtools::use_package("testthat")

context("Sparse PCA")

#Load rsvd library
library(spcaRcpp)     

#Set seed
set.seed(1234)

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
# Test: SPARSE PCA - center = TRUE
#*************************************************************************************

pca_out <- prcomp(X, rank = 3,center = TRUE,) #Sparse PCA
spca_out <- spcaRcpp(X, k=3, alpha=0, beta=0, center = TRUE) #Sparse PCA with Rcpp

#Test1: SPCA recovers PCA for alpha = 0 and beta = 0.
testthat::test_that("Test cov; alpha = 0 and beta = 0", {
  testthat::expect_equal(pca_out$sdev[1:3], spca_out$sdev[1:3])
  testthat::expect_equal(sum(diag(1,3,3) - t(spca_out$loadings)%*%spca_out$loadings), 0 )
})


