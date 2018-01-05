
#' Confidence Intervals For Change Point Estimator
#'
#' A change point in the mean function of functional data is estimated, and for a user selected level
#' aconfidence interval is computed using the fully functional procedure introduced in Aue, Rice and Sonmez (2017+).
#'
#' @param fdobj Functional data object, class of \code{"fd"}
#' @param h Bandwidth parameter to estimate the long run covariance operator.
#' @param alpha Nominal level used to construct the confidence intervals
#' @param kern_type Kernel to be used for the estimation of the long run covariance function
#' @param ... Further arguments to pass
#'
#' @import stats
#' @import fda
#' @export
#'
#'@return
#' \item{\code{Lower CI}}{
#' Lower Confidence Band
#'}
#'\item{\code{Estimate}}{
#' Estimated change point location
#'}
#'\item{\code{Upper CI}}{
#' Upper Confidence Band
#'}
#'
#'@references Aue, A., Rice, G. & O. Sonmez   (2017+), \emph{ Detecting and dating structural breaks ¨
#' in functional data without dimension reduction}, (https://arxiv.org/pdf/1511.04020.pdf)
#'
#' @seealso \code{\link{LongRun}} \code{\link{opt_bandwidth}} \code{\link{change_FF}}
#'
#'@examples
#'fdata1 = fun_AR(n=100, nbasis=21, order=1, kappa=0.5)
#'Conf_int(fdata1, h=2)


Conf_int = function(fdobj, h=0, kern_type = "BT", alpha=0.05, ...){
  samp = fdobj$coefs
  N = ncol(fdobj$coefs)
  D =  nrow(fdobj$coefs)
  Sn=(1:N)
  Sn[1]=0
  for(j in (2:N)){
    Sn[j]= sum((rowSums(samp[,1:j]) - (j/N)*rowSums(samp[,1:N]))^2)/N
  }
  k.star = min(which(Sn==max(Sn)))
  theta.hat = k.star/N
  delta.h = mean(fdobj[(k.star+1):N]) - mean(fdobj[1:k.star])
  norm.d = inprod(delta.h, delta.h)
  LongRunC = LongRun(fdobj=fdobj, h=h, kern_type = kern_type)
  lambda.hat = LongRunC$e_val
  phi.hat = LongRunC$e_fun
  sigma.h = (norm.d)^-1*sum(sapply(1:D, function(j) lambda.hat[j]*(inprod(delta.h,phi.hat[j]))^2))
  sigma.hat = sqrt(sigma.h)
  q1 = quant((1-theta.hat), theta.hat, sigma.hat, sigma.hat, alpha/2)
  q2 = quant((1-theta.hat), theta.hat, sigma.hat, sigma.hat, 1-alpha/2)
  # upper and lower 1-alpha CI
  upper = k.star + q2/norm.d
  lower = k.star + q1/norm.d
  if(lower<0){lower=0}
  if(upper>N){upper=N}
  out = c(lower, k.star, upper)
  names(out) = c("Lower CI", "Estimate", "Upper CI")
  return(out)
}

