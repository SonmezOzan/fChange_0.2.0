#' Simulate Functional Auto-Regressive Process
#'
#'
#'This generates a functional Auto-Regressive, FAR, process of sample size n  with a specific
#' number of basis functions where the eigenvalue decay of the covariance operator is given
#' by the defined vector Sigma. The norm of the FAR
#' operators are defined by the vector \code{kappa}. The generic function uses Fourier basis in [0,1],
#' however one can define a different basis and different range values. If the order or kappa is not defined
#' then the function generates iid functional data by default.
#' @param n Sample size of generated functional data. A strictly positive integer
#' @param nbasis Number of basis functions used to represent functional observations
#' @param order Order of the FAR process
#' @param kappa Vector of norm of the FAR operators. The length of this vector must be same as the
#' FAR \code{order}
#' @param Sigma  Eigen value decay of the covariance operator of the functional data. The eigenvalues of
#' the covariance operator of the generated functional sample are given by \code{Sigma}.
#' The length of \code{Sigma} must match number of basis. By default it is set as \code{(1:nbasis)^-1}
#' @param basis A functional basis object defining the basis. It can be the class of
#' \code{basisfd}, \code{fd}, \code{fdPar}. As a default it is set to be a Fourier basis
#' @param rangeval A vector of length 2 containing the initial and final values of the
#' interval over which the functional data object can be evaluated. As a default it is
#' set to be [0,1].
#' @param ... Further arguments to pass
#' @import stats
#' @import fda
#' @export
#'
#' @return Functional Auto-Regressive data sample (class \code{fd}) containing:
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
#'@details This function should be used for a simple FAR data generation
#'for a desired eigenvalue decay of covariance operator. The j-th FAR operator \eqn{\Psi[j]} is generated
#' by \eqn{\Psi[j] = \kappa[j]\Psi}, where \eqn{\Psi} has a unit norm with \eqn{\Psi[i,j] = N(0, \sigma[i]\sigma[j])}.
#' For more details see Aue A., Rice G., Sonmez O. (2017+).
#'
#'@references Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
#' Data Analysis, 2nd ed.}, Springer, New York.
#' @references Aue A., Rice G., Sonmez O. (2017+), \emph{Detecting and dating structural breaks in
#' functional data without dimension reduction} (https://arxiv.org/pdf/1511.04020.pdf)
#' @seealso \code{\link{Data2fd}, \link{fun_IID}, \link{fun_MA}}
#
#' @examples # FAR(1) data with 21 fourier basis with a geometric eigenvalue decay
#' fun_AR(n=100, nbasis=21, order=1, kappa=0.8)
#'
#'
#' @examples # Define eigenvalue decay
#' Sigma1 = 2^-(1:21)
#' # Then generate FAR(2) data
#' fun_AR(n=100, nbasis=21, order=2, kappa= c(0.5, 0.3), Sigma=Sigma1)
#'
#' @examples # Define eigenvalue decay, and basis function
#' library(fda)
#' basis1 = create.bspline.basis(rangeval = c(0,1), nbasis=21)
#' Sigma1 = 2^-(1:21)
#' # Then generate FAR(1)
#' fun_AR(n=100, nbasis=21, order=1, kappa= 0.3,Sigma=Sigma1, basis=basis1)
#'
#' @examples # Not defining order will result in generating IID functions
#' fun_AR(n=100, nbasis=21) # same as fun_IID(n=100, nbasis=21)



fun_AR = function(n , nbasis, order=NULL, kappa=NULL, Sigma = NULL,  basis = NULL, rangeval=c(0,1), ...){
   if (is.null(Sigma)){
      Sigma = (1:nbasis)^-1}
   if (length(Sigma)!=nbasis){
      stop("Length of Sigma must be equal to nbasis")
   }
   if (is.null(order) || order==0){
      fun_IID(n = n , nbasis = nbasis, Sigma = Sigma)
   } else {
      if(order!=0 & is.null(kappa)==FALSE & length(kappa)!=order){
         stop("Order of FAR does not match with kappa")
      }
      if(order!=0 & is.null(kappa)){
         stop("Please Define kappa")
      }else{
         if (is.null(basis)){
            basis = create.fourier.basis(rangeval = rangeval, nbasis=nbasis)
         }
         #burnin
         l_burnin = n/2
         edata = matrix(0, nbasis, (n+l_burnin))
         for (i in 1:(n+l_burnin)){
            edata[,i] = rnorm(nbasis, 0, Sigma)
         }
         # Psi operator
         Psi = matrix(0, nbasis, nbasis)
         for (i in 1:nbasis){
            for (j in 1:nbasis){
               Psi[i, j] = rnorm(1, 0, t(Sigma[i])%*%Sigma[j])
            }
         }
         #adjust the norm
         Psi = Psi/norm(Psi, type="F")
         Psi.i = list()
         for (j in 1:order){
            Psi.i[[j]] = kappa[j] * Psi
         }
         coef = matrix(0, nbasis, (n+l_burnin))
         coef[, 1:order] = edata[, 1:order]
         for (i in (order+1):(n+l_burnin)){
            coef[,i] = rowSums(sapply(1:order, function(j) Psi.i[[j]]%*% coef[, i-j])) + edata[,i]
         }
         dat = coef[, (l_burnin+1):(n+l_burnin)]
         fd(dat, basis)
      }
   }
}
