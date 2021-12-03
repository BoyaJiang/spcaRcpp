## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(spcaRcpp)

## -----------------------------------------------------------------------------
m <- 10000
V1 <- rnorm(m, -100,200)
V2 <- rnorm(m, -100, 300)
V3 <- -0.1*V1 + 0.1*V2 + rnorm(m, 0, 100)

X <- cbind(V1,V1,V1,V1, V2,V2,V2,V2, V3,V3)
X <- X + matrix(rnorm(length(X),0,1), ncol = ncol(X), nrow = nrow(X))

## -----------------------------------------------------------------------------
rcpp_out = spcaRcpp(X, k = 3, alpha = 0, beta = 0, max_iter = 10000,center = TRUE)
print (rcpp_out)

## -----------------------------------------------------------------------------
library(sparsepca)
spca_out = sparsepca::spca(X, k = 3, alpha = 0, beta = 0, center = TRUE, scale = FALSE, verbose = FALSE)

## -----------------------------------------------------------------------------
all.equal(spca_out$sdev, rcpp_out$sdev)
all.equal(spca_out$var, rcpp_out$var)
all.equal(spca_out$loadings, rcpp_out$loadings)
all.equal(spca_out$eigenvalues, rcpp_out$eigenvalues)
all.equal(spca_out$center, rcpp_out$center)

## -----------------------------------------------------------------------------
library(microbenchmark)
library(ggplot2)
runtime = microbenchmark::microbenchmark(spca_out, rcpp_out, times = 200)
autoplot(runtime)

## -----------------------------------------------------------------------------
prcomp_out = prcomp(X, center = TRUE, rank. = 3)

all.equal(prcomp_out$sdev[1:3], rcpp_out$sdev)
all.equal(prcomp_out$rotation, rcpp_out$loadings, check.attributes = FALSE)
all.equal(spca_out$center, rcpp_out$center)

## -----------------------------------------------------------------------------
library(tidyverse, quietly = TRUE)
data("USArrests")
head(USArrests)

## -----------------------------------------------------------------------------
rcpp_out2 = spcaRcpp(USArrests,k = 4, center = TRUE, alpha = 0, beta = 0)
spca_out2 = sparsepca::spca(USArrests,k = 4, center = TRUE, alpha = 0, beta = 0)
rcpp_out2
spca_out2

## -----------------------------------------------------------------------------
all.equal(spca_out2$sdev, rcpp_out2$sdev)
all.equal(spca_out2$var, rcpp_out2$var)
all.equal(spca_out2$loadings, rcpp_out2$loadings)
all.equal(spca_out2$eigenvalues, rcpp_out2$eigenvalues)
all.equal(spca_out2$center, rcpp_out2$center)

## -----------------------------------------------------------------------------
runtime = microbenchmark::microbenchmark(spca_out2, rcpp_out2, times = 200)
autoplot(runtime)

