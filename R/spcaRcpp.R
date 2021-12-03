#' @title  Rcpp Integration for Sparse Principal Component Analysis (spca).
#
#' @description Implementation of SPCA, using variable projection as an optimization strategy.
#
#' @details
#' Sparse principal componenet analysis is a specialized variant of PCA. Specifically, SPCA promotes sparsity
#' in the modes, i.e., the sparse modes have only a few active (nonzero) coefficients, while the majority of coefficients
#' are constrained to be zero. This approach leads to a improved localization and interpretability of the model
#' compared to the global PCA modes obtained from traditional PCA. In addition, SPCA avoids overfitting in a 
#' high-dimensional data setting where the number of variables \eqn{p} is greater than the number of
#' observations \eqn{n}.  
#' 
#' Given an \eqn{(n,p)} data matrix \eqn{X}, SPCA attemps to minimize the following
#' objective function:
#'
#' minimize f(A,B) = 1/2⋅‖X - X⋅B⋅Aᵀ‖² + α⋅‖B‖₁ + 1/2⋅β‖B‖², subject to AᵀA = I.
#'
#' where \eqn{B} is the sparse weight matrix and \eqn{A} is an orthonormal matrix.
#' The principal components \eqn{Z} are formed as
#'
#' \deqn{ Z = X B }{Z = X * B}
#'
#' and the data can be approximately rotated back as
#'
#' \deqn{ \tilde{X} = Z A^\top }{ X = Z t(A)}
#'
#' The print method can be used to present the results in a nice format.
#'
#' @param X       array_like; \cr
#'                a real \eqn{(n, p)} input matrix.
#'
#' @param k       integer; \cr
#'                specifies the target rank, i.e., the number of components to be computed.
#'
#' @param alpha   float; \cr
#'                Sparsity controlling parameter. Higher values means sparser components.
#'
#' @param beta    float; \cr
#'                Amount of ridge shrinkage to apply in order to improve conditioning.
#'                
#' @param center  bool; \cr
#'                logical value indicating if the variables should be zero centered
#'                (TRUE by default).
#'                
#' @param max_iter integer; \cr
#'                 maximum number of iterations to perform.
#'
#' @param tol     float; \cr
#'                stopping criteria for the convergence.
#'
#'@return \code{spca} returns a list containing the following four components:
#'\item{loadings}{  array_like; \cr
#'           sparse weight vector;  \eqn{(p, k)} dimensional array.
#'}
#'
#'\item{standard deviations}{  array_like; \cr
#'          the approximated standard deviations; \eqn{(k)} dimensional array.
#'}
#'
#'\item{eigenvalues}{  array_like; \cr
#'          the approximated eigenvalues; \eqn{(k)} dimensional array.
#'}
#'
#'\item{center}{  array_like; \cr
#'          the centering used.
#'}
#'
#'\item{var}{  float; \cr
#'          the variance.
#'}
#'
#'\code{\link{iteration}}.
#'
#' @references
#' \itemize{
#'
#'  \item [1] N. B. Erichson, P. Zheng, K. Manohar, S. Brunton, J. N. Kutz, A. Y. Aravkin.
#'  "Sparse Principal Component Analysis via Variable Projection."
#'  SIAM Journal on Applied Mathematics 2020 80:2, 977-1002
#'  (available at `arXiv \url{https://arxiv.org/abs/1804.00341}).
#' 
#' 
#' \item [2] N. B. Erichson, P. Zheng, S. Aravkin, sparsepca, (2018), GitHub repository,
#' \url{https://github.com/erichson/spca}. 
#'}
#'
#' @author Boya Jiang
#' 
#' @examples
#'
#' # Create artifical data
#' m <- 10000
#' V1 <- rnorm(m, -100,200)
#' V2 <- rnorm(m, -100, 300)
#' V3 <- -0.1*V1 + 0.1*V2 + rnorm(m, 0, 100)
#'
#' X <- cbind(V1,V1,V1,V1, V2,V2,V2,V2, V3,V3)
#' X <- X + matrix(rnorm(length(X),0,1), ncol = ncol(X), nrow = nrow(X))
#'
#' # Compute SPCA
#' out <- spcaRcpp(X, k=3, alpha=1e-3, beta=1e-3, center = TRUE)
#' print(out)
#'
#'
#' @import Rcpp
#' @import RcppArmadillo
#' @export
spcaRcpp <- function(X, k=NULL, alpha=1e-4, beta=1e-4, center=TRUE, max_iter=1000, tol=1e-5) UseMethod("spcaRcpp")

#' @export
spcaRcpp.default <- function(X, k=NULL, alpha=1e-4, beta=1e-4, center=TRUE, max_iter=1000, tol=1e-5) {
  options(warn=-1)
  library(Rcpp)
  library(RcppArmadillo)

  X = as.matrix(X)
  #initiate spca object
  spcaObj = list(loadings = NULL,
                 eigenvalues = NULL,
                 center = center)

  n <- nrow(X)
  p <- ncol(X)

  #set target rank, k must be less or equal to the minimum of (n,p)
  if(is.null(k)) {
    k <- min(n,p)
  }else if(k > min(n,p)) {
    k <- min(n,p)
  }

  #Center Data
  if(center == TRUE) {
    spcaObj$center <- colMeans(X)
    X <- sweep(X, MARGIN = 2, STATS = spcaObj$center, FUN = "-", check.margin = TRUE)
  } else { 
    spcaObj$center <- FALSE 
  }

  #initialization step, calculate svd
  svd_X <- svd(X)
  Dmax  <- svd_X$d[1] # l2 norm

  B <- svd_X$v[,1:k]
  V <- svd_X$v
  
  VD  = sweep(V, MARGIN = 2, STATS = svd_X$d, FUN = "*")
  VD2 = sweep(V, MARGIN = 2, STATS = svd_X$d**2, FUN = "*")

  #set tuning parameters
  alpha <- alpha * Dmax^2
  beta  <- beta  * Dmax^2
  nu    <- 1.0 / (Dmax^2 + beta)
  kappa <- nu * alpha

  #apply variable projection solver using RcppArmadillo iterations
  ret = iter(max_iter, tol, VD, VD2, V, B, alpha, beta, nu, kappa)

  #update spca objects
  spcaObj$loadings    <- ret$B
  spcaObj$eigenvalues <- as.vector(t(ret$d / (n - 1)))

  #explained variance
  spcaObj$sdev <- sqrt( spcaObj$eigenvalues )
  spcaObj$var  <- sum( apply( Re(X) , 2, stats::var ) )
  
  #return the real part for complex values
  if (is.complex(X)){
    spcaObj$var <- Re(spcaObj$var + sum( apply( Im(X) , 2, stats::var ) ))
  }
  
  class (spcaObj) <- "spca"
  return (spcaObj)
}

#' @export
print.spcaRcpp <- function(x , ...) {
  #print spca
  cat("Standard deviations:\n")
  print (round(x$sdev,4))
  cat("\nEigenvalues:\n")
  print (round(x$eigenvalues,4))
  cat("\nSparse loadings:\n")
  print (round(x$loadings,4))
}

