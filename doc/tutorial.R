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
rcpp_out = spcaRcpp(X, k = 3, alpha = 1e-5, beta = 1e-5, max_iter = 10000,center = TRUE)
print (rcpp_out)

## -----------------------------------------------------------------------------
rcpp_out2 = spcaRcpp(X, k = 3, alpha = 1e-3, beta = 1e-3, max_iter = 10000,center = TRUE)
print (rcpp_out2)

## -----------------------------------------------------------------------------
rcpp_out3 = spcaRcpp(X, k = 3, alpha = 0, beta = 0, max_iter = 10000,center = TRUE)
prcomp_out = prcomp(X, center = TRUE, rank. = 3)

all.equal(prcomp_out$sdev[1:3], rcpp_out3$sdev)
all.equal(prcomp_out$rotation, rcpp_out3$loadings, check.attributes = FALSE)
all.equal(prcomp_out$center, rcpp_out3$center)

## -----------------------------------------------------------------------------
library(tidyverse, quietly = TRUE)
data("USArrests")
head(USArrests)

## -----------------------------------------------------------------------------
rcpp_out_rd = spcaRcpp(USArrests,k = 4, center = TRUE, alpha = 0, beta = 0)
prcomp_out_rd = prcomp(USArrests,center = TRUE, rank. = 4)
rcpp_out_rd
prcomp_out_rd

## -----------------------------------------------------------------------------
all.equal(prcomp_out_rd$sdev, rcpp_out_rd$sdev)
all.equal(prcomp_out_rd$rotation, rcpp_out_rd$loadings,
          check.attributes = FALSE)
all.equal(prcomp_out_rd$center, rcpp_out_rd$center)

## -----------------------------------------------------------------------------
runtime = microbenchmark::microbenchmark(prcomp_out_rd, rcpp_out_rd, times = 200)
autoplot(runtime)

