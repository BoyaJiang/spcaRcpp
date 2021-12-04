<!-- badges: start -->
[![R-CMD-check](https://github.com/BoyaJiang/spcaRcpp/workflows/R-CMD-check/badge.svg)](https://github.com/BoyaJiang/spcaRcpp/actions)
<!-- badges: end -->
# Rcpp Integration of SPCA via Variable Projection

Sparse principal componenet analysis is a specialized variant of PCA. Specifically, SPCA promotes sparsity
in the modes, i.e., the sparse modes have only a few active (nonzero) coefficients, while the majority of coefficients
are constrained to be zero. This approach leads to a improved localization and interpretability of the model
compared to the global PCA modes obtained from traditional PCA. In addition, SPCA avoids overfitting in a 
high-dimensional data setting where the number of variables ``p`` is greater than the number of observations ``n``.

This package provides SPCA routines in R/Rcpp:
 
* Sparse PCA: ``spca()``.

## Problem Formulation

Given a data matrix ``X`` with shape ``(m,p)``, SPCA attemps to minimize the following
optimization problem:

```
minimize f(A,B) = 1/2⋅‖X - X⋅B⋅Aᵀ‖² + α⋅‖B‖₁ + 1/2⋅β‖B‖², subject to AᵀA = I.
```

The matrix ``B`` is the sparse weight (loadings) matrix and ``A`` is an orthonormal matrix.

Then, the principal components ``Z`` are then formed as

```
Z = X %*% B.
```

Specifically, the interface of the SPCA function is:

```R
spca(X, k, alpha=1e-4, beta=1e-4, center=TRUE,max_iter=1000, tol=1e-4)
```
The description of the arguments is listed in the following:

* ``X`` is a real ``n`` by ``p`` data matrix, where ``n`` denotes the number of observations and ``p`` the number of variables.

* ``k`` specifies the target rank, i.e., number of components to be computed.

* ``alpha`` is a sparsity controlling parameter. Higher values lead to sparser components.

* `` beta`` is the amount of ridge shrinkage to apply in order to improve conditioning.

* ``center`` logical value which indicates whether the variables should be zero centered (TRUE by default).

* ``max_iter`` maximum number of iterations to perform before exiting (default is 1000).

* ``tol`` stopping criteria for convergence (default is ``1e-5``).

A list with the following components is returned:

* ``loadings`` sparse loadings (weight) vector.
* ``standard deviations`` the approximated standard deviations; ``k`` dimensional array.
* ``eigenvalues`` the approximated eigenvalues; 

## Installation

Install the developer version of sparsepca package via github
```R
#install.packages("devtools")
library(devtools)
devtools::install_github("BoyaJiang/spcaRcpp")
library(spcaRcpp)
```

## References

* [N. Benjamin Erichson, et al. "Sparse Principal Component Analysis via Variable Projection." (2018)](https://arxiv.org/abs/1804.00341)

* [N. Benjamin Erichson, et al. "Randomized Matrix Decompositions using R." (2016)](http://arxiv.org/abs/1608.02148)

* [N. B. Erichson, P. Zheng, S. Aravkin, sparsepca, (2018), GitHub repository](https://github.com/erichson/spca) 
