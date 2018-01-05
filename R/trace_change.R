

#' Testing changes in the trace of the covariance operator in functional data
#'
#' This function tests and detects changes in the trace of the covariance operator.
#'
#'
#' @param fdobj A functional data object of class '\code{fd}'
#' @param M Number of monte carlo simulations used to get the critical values. The default value is \code{M=1000}
#' value is \code{h=2}.
#' @param mean_change If \code{TRUE} then the data is centered considering the change in the mean function.
#' @param delta Trimming parameter to estimate the covariance function using partial sum estimates.
#'
#' @import sde
#' @import fda
#' @import sandwich
#' @export
#'
#'
#'
#' @return
#'\item{\code{pvalue}}{
#' Approximate p value for testing whether there is a significant change in the desired eigenvalue of the covariance operator
#'}
#'\item{\code{change}}{
#' Estimated change location
#'}
#'\item{\code{trace_before}}{
#' Estimated trace before the change
#'}
#'\item{\code{trace_after}}{
#' Estimated trace after the change
#'}
#'@details This function dates and detects changes in trace of the covariance function. This can be interpreted as the changes in
#'the total variation of the the functional data. Trace is defined as the infinite sum of the eigenvalues of the covariance operator
#' and for the sake of implementation purpose, the sum is truncated up to the total number of basis functions that defines the functional
#' data at hand. The critical values are approximated via \code{M} Monte Carlo simulations.
#'
#'
#'
#' @examples
#' # generate functional data
#' fdata = fun_IID(n=100, nbasis=21)
#' trace_change(fdata)


trace_change = function(fdobj, mean_change = FALSE,  delta = 0.1, M = 1000){
  N = ncol(fdobj$coefs)
  D = nrow(fdobj$coefs)
  cdata = center_data(fdobj, change = mean_change)
  Cov_op = partial_cov(cdata, 1)
  lambda= Cov_op$eigen_val
  T_1 = sum(lambda)
  Xi = sapply(1:N, function(i) inprod(cdata[i], cdata[i]))
  sigma_sq = lrvar(Xi, prewhite=F)
  sigma = sqrt(sigma_sq)
  Values = sapply(1:M, function(k) max(abs(BBridge(0, 0, 0, 1, N)[(floor(delta*N)+1):N])))
  s = floor(delta*N)
  Tn = c(rep(0,s))
  for (k in (s+1):N){
    T_x = sum(Xi[1:k])/N
    Tn[k] = (1/sigma)*abs(T_x - (k)/N*T_1)
  }
  Sn = max(Tn)
  k_star = min(which.max(Tn))
  z = Sn <= Values
  p = length(z[z==TRUE])/length(z)
  tr_before = sum(pca.fd(fdobj[1:k_star], nharm = D, centerfns = T)$values)
  tr_after = sum(pca.fd(fdobj[(1+k_star):N], nharm = D, centerfns = T)$values)
  list(change = k_star/N,
       pvalue = p,
       trace_before = tr_before,
       trace_after = tr_after)
}




