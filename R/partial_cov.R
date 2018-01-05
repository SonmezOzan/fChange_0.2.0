
#' Partial Sample Estimates of the Covariance Function in Functional Data Analysis
#'
#' This function computes the partial sum estimate of the covariance function in functional data analysis.
#' It also generate the eigenvalues and eigenfunctions of the parial sum estimate, along with the coefficient matrix
#' of the estimated covariance function.
#'
#'
#' @param fdobj A functional data object of class '\code{fd}'
#' @param x Fraction of the sample size that the partial sum is computed. This input must be in (0,1]. The default is \code{x=1}
#' which corresponds to the regular covariance function estimate in functional data analysis.
#'
#' @import fda
#' @export
#'
#'
#'
#' @return
#'\item{\code{eigen_val}}{
#' Eigenvalues of the partial sum estimate of the covariance function
#'}
#'\item{\code{eigen_fun}}{
#' Eigenfunctions of the partial sum estimate of the covariance function
#'}
#'\item{\code{coef_matrix}}{
#' Coefficient matrix of the partial sum estimate of the covariance function
#'}
#'
#'@details This function simply estimates the covariance function based on the partial sum of the centered functional observations.
#'The length of the sum is determined by \code{x}, and when \code{x=1} this estimate corresponds to the regular covariance function
#'estimates using the whole sample.
#'
#' @seealso \code{\link{pca.fd}}
#'
#' @examples
#' library(fda)
#' # generate functional data
#' fdata = fun_IID(n=100, nbasis=21)
#' # Estimated eigenvalues
#' e1 = partial_cov(fdata)$eigen_val
#' e2 = pca.fd(fdata, nharm = 21, centerfns = TRUE)$values
#' # e1 and e2 will both estimate the eigenvalues of the covariance
#' # operator based on the whole sample
#' # estimates using only 90% of the data
#' Cov = partial_cov(fdata, 0.9)
#'
partial_cov = function(fdobj, x = NULL){
  if (class(fdobj)!="fd"){stop("Must insert functional data")}
  n = ncol(fdobj$coefs)
  D = nrow(fdobj$coefs)
  if (is.null(x)){x=1}
  if (x>1 | x<=0){stop("x should be in (0,1]")}
  k = floor(n*x)
  cdata = center.fd(fdobj)
  dat = cdata$coefs
  C = matrix(0, D, D)
  for (i in 1:k){
    C = C + dat[,i]%*%t(dat[,i])
  }
  E = eigen(C/n)
  e_val = E$value
  e_fun = fd(E$vector, fdobj$basis)
  list(eigen_val = e_val, eigen_fun = e_fun, coef_matrix = C/n)
}
