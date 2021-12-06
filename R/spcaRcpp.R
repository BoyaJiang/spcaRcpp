#' @title  Rcpp Integration for Sparse Principal Component Analysis (spca).
#
#' @description Implementation of SPCA, using variable projection as an optimization strategy.
#
#' @details
#' Sparse principal component analysis is a specialized variant of PCA. Specifically, SPCA promotes sparsity
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
#' @param X       a numeric matrix or data.frame which provides the data for 
#'                the sparse principal components analysis.
#'
#' @param k       optional, a number specifying the maximal rank.
#' 
#' @param alpha   Sparsity controlling parameter. Higher values means sparser components.
#'
#' @param beta    Amount of ridge shrinkage to apply in order to improve conditioning.
#'                
#' @param center  a logical value indicating whether the variables should be 
#'                shifted to be zero centered.  
#'                
#' @param max_iter maximum number of iterations to perform.
#'
#' @param tol     stopping criteria for the convergence.
#'
#'@return \code{spca} returns a list containing the following six components:
#'\item{loadings}{  the matrix of variable loadings. 
#'}
#'
#'\item{standard deviations}{  the approximated standard deviations.
#'}
#'
#'\item{eigenvalues}{  the approximated eigenvalues.
#'}
#'
#'\item{center}{  the centering used.
#'}
#'
#'\item{var}{  the variance.
#'}
#'
#'\item{scores}{  the principal component scores.
#'}
#'
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
#' V1 <- rnorm(m, -100, 200)
#' V2 <- rnorm(m, -100, 300)
#' V3 <- -0.1 * V1 + 0.1 * V2 + rnorm(m, 0, 100)
#'
#' X <- cbind(V1, V1, V1, V1, V2, V2, V2, V2, V3, V3)
#' X <- X + matrix( rnorm( length(X), 0, 1 ), ncol = ncol(X), nrow = nrow(X) )
#'
#' # Compute SPCA
#' out <- spcaRcpp(X, k=3, alpha=1e-4, beta=1e-4, center = TRUE)
#' print(out)
#'
#'
#' @import Rcpp
#' @import RcppArmadillo
#' @rawNamespace useDynLib(spcaRcpp)
#' @export
spcaRcpp <- function(X, k=NULL, alpha=1e-4, beta=1e-4, center=TRUE, max_iter=1000, tol=1e-5) UseMethod("spcaRcpp")

#' @export
spcaRcpp.default <- function(X, k=NULL, alpha=1e-4, beta=1e-4, center=TRUE, max_iter=1000, tol=1e-5) {
  options(warn = -1)

  #make sure X is a matrix
  X = as.matrix(X)
  
  #initiate spca object
  spcaObj = list(loadings = NULL,
                 eigenvalues = NULL,
                 center = center,
                 scores = NULL)

  n <- nrow(X)
  p <- ncol(X)

  #set target rank, k must be less or equal to the minimum of (n,p)
  if( is.null(k) ) {
    k <- min(n, p)
  }else if( k > min(n, p) ) {
    k <- min(n, p)
  } 

  #Center Data
  if(center == TRUE) {
    spcaObj$center <- colMeans(X)
    X <- sweep(X, MARGIN = 2, STATS = spcaObj$center, FUN = "-", check.margin = TRUE)
  }else{ 
    spcaObj$center <- FALSE 
  }

  #initialization step, calculate svd
  svd_X <- svd(X)
  Dmax  <- svd_X$d[1] # l2 norm

  B <- svd_X$v[,1:k]
  V <- svd_X$v
  
  VD  = sweep(V, MARGIN = 2, STATS = svd_X$d, FUN = "*", check.margin = TRUE)
  VD2 = sweep(V, MARGIN = 2, STATS = svd_X$d ^ 2, FUN = "*", check.margin = TRUE)

  #set tuning parameters
  alpha <- alpha * Dmax^2
  beta  <- beta  * Dmax^2
  nu    <- 1.0 / (Dmax^2 + beta)
  kappa <- nu * alpha

  #apply variable projection solver using RcppArmadillo iterations
  ret = iter(max_iter, tol, VD, VD2, V, B, alpha, beta, nu, kappa)

  #update spca objects
  spcaObj$loadings    <- ret$B
  spcaObj$eigenvalues <- as.vector( t( ret$d / (n - 1) ) )
  spcaObj$scores      <- X %*% B

  #explained variance
  spcaObj$sdev <- sqrt(spcaObj$eigenvalues)
  spcaObj$var  <- sum( apply( Re(X) , 2, stats::var ) )
  
  class (spcaObj) <- "spca"
  return (spcaObj)
}


