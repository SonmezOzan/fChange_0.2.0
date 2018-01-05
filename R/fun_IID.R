#' Simulate Independent Functional Data
#'
#' This function generates independent functional observations of sample size n with a desired eigenvalue decay structure of the covariance operator.
#'
#'
#' @param n Sample size of generated functional data. A strictly positive integer
#' @param nbasis Number of basis functions used to represent functional observations
#' @param Sigma  Eigen value decay of the covariance operator of the functional data. The eigenvalues of
#' the covariance operator of the generated functional sample are given by \code{Sigma}.
#' The length of \code{Sigma} must match number of basis. By default it is set as \code{(1:nbasis)^-1}
#' @param basis A functional basis object defining the basis. It can be the class of
#' "\code{basisfd, fd, fdPar}". As a default it is set to be a Fourier basis
#' @param rangeval A vector of length 2 containing the initial and final values of the
#' interval over which the functional data object can be evaluated. As a default it is
#' set to be [0,1].
#' @param ... Further arguments to pass
#' @import stats
#' @import fda
#' @export
#'
#' @return An independent functional data sample (class \code{fd}) containing:
#'\item{coefs}{
#' The coefficient array
#'}
#'\item{basis}{
#' A basis object
#'}
#'\item{fdnames}{
#' A list containing names for the arguments, function values
#' and variables
#'}
#'
#'@details Independent functional sample is generated based on a linear combination of basis functions where the \eqn{i}-th linear
#'combination coefficient is normally distributed with mean zero and standard deviation \eqn{\sigma[i]}.
#'
#'
#'@references Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
#' Data Analysis, 2nd ed.}, Springer, New York.
#' @references Aue A., Rice G., Sonmez O. (2017+), \emph{Detecting and dating structural breaks in
#' functional data without dimension reduction} (https://arxiv.org/pdf/1511.04020.pdf)
#' @seealso \code{\link{Data2fd}}
#' @examples # Functional data with 21 fourier basis with a geometric eigenvalue decay
#' fun_IID(n=100, nbasis=21)
#'
#' @examples # Define eigenvalue decay
#' Sigma1=2^-(1:21)
#' # Then generate functional data
#' fun_IID(n=100, nbasis=21, Sigma=Sigma1)
#'
#' @examples # Define eigenvalue decay, and basis function
#' library(fda)
#' basis1 = create.bspline.basis(rangeval = c(0,1), nbasis=21)
#' Sigma1=2^-(1:21)
#' # Then generate functional data
#' fun_IID(n=100, nbasis=21, Sigma=Sigma1, basis=basis1)
#'



fun_IID <- function(n, nbasis, Sigma = NULL, basis=NULL, rangeval=c(0,1), ...){
   if (is.null(Sigma)){
      Sigma = (1:nbasis)^-1}
   if (length(Sigma)!=nbasis){
      stop("Length of Sigma must be equal to nbasis")
   } else {
      if (is.null(basis)){
         basis = create.fourier.basis(rangeval = rangeval, nbasis=nbasis)
      }
      mdata = matrix(0, nbasis, n)
      for (j in 1:n){
         mdata[,j] = rnorm(nbasis, 0, Sigma)
      }
      fd(mdata, basis)
   }
}
