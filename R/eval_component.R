

#' Detecting Changes in the Eigenvalues of the Covariance Operator of the Functional Data
#'
#' This function tests and detects changes in the specific eigenvalue of the covariance operator.
#'
#'
#' @param fdobj A functional data object of class '\code{fd}'
#' @param component The eigenvalue that the componentwise test is applied to.
#' @param M Number of monte carlo simulations used to get the critical values. The default value is \code{M=1000}
#' @param  h The window parameter for the estimation of the long run covariance matrix. The default
#' value is \code{h=2}.
#' @param mean_change If \code{TRUE} then the data is centered considering the change in the mean function.
#' @param delta Trimming parameter to estimate the covariance function using partial sum estimates.
#'
#' @import sde
#' @import fda
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
#'@details This function dates and detects changes in the defined eigenvalue of the covariance function. The critical values are
#'approximated via \code{M} Monte Carlo simulations.
#'
#'
#'
#' @examples
#' # generate functional data
#' fdata = fun_IID(n=100, nbasis=21)
#' eval_component(fdata, 2)
#'
#'
#'


eval_component = function(fdobj, component, h =2, mean_change = FALSE, delta = 0.1, M = 1000){
  if (class(fdobj)!="fd"){stop("Must insert functional data")}
  N = ncol(fdobj$coefs)
  D = nrow(fdobj$coefs)
  cdata = center_data(fdobj, change = mean_change)
  Cov_op = partial_cov(cdata,1)
  PSI = Cov_op$coef_matrix
  Phi = Cov_op$eigen_fun
  lambda = Cov_op$eigen_val

  Projections = inprod(Phi, cdata)
  Proj_sq = Projections^2
  Psi_diag = diag(PSI)
  THETA = Proj_sq - Psi_diag
  theta = matrix(THETA[component,], ncol = N, nrow = 1)
  Sigma_d = LongRunCovMatrix(theta, h=2)
  Values = sapply(1:M, function(k) max(BBridge(0, 0, 0, 1, N)[(floor(delta*N)+1):N]^2))

  s = floor(delta*N)
  Tn = c(rep(0,s))
  for (k in (s+1):N){
    lam_i = partial_cov(cdata, k/N)$eigen_val[component]
    Tn[k] = (N*(lam_i - (k/N)*lambda[component])^2)/Sigma_d
  }
  Sn = max(Tn)
  k_star = min(which.max(Tn))
  z = Sn <= Values
  p = length(z[z==TRUE])/length(z)
  list(change = k_star/N, pvalue = p)
}
