
#' Detecting Changes Jointly in the Eigenvalues of the Covariance Operator of the Functional Data
#'
#' This function tests and detects changes jointly in the eigenvalue of the covariance operator.
#'
#'
#' @param fdobj A functional data object of class '\code{fd}'
#' @param d Number of eigenvalues to include in testing.
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
#'\item{\code{eval_before}}{
#' Estimated eigenvalues before the change
#'}
#'\item{\code{eval_after}}{
#' Estimated eigenvalues after the change
#'}
#'@details This function dates and detects changes in the joint eigenvalues that is defined by \code{d} of the covariance function.
#'The critical values are approximated via \code{M} Monte Carlo simulations.
#'
#'
#'
#' @examples
#' # generate functional data
#' fdata = fun_IID(n=100, nbasis=21)
#' eval_joint(fdata, 2)



eval_joint = function(fdobj, d, h =2, mean_change = FALSE, delta = 0.1, M = 1000){
  N = ncol(fdobj$coefs)
  D = nrow(fdobj$coefs)
  cdata = center_data(fdobj, change = mean_change)
  Cov_op = partial_cov(cdata, 1)
  PSI = Cov_op$coef_matrix
  Phi = Cov_op$eigen_fun
  lambda = Cov_op$eigen_val

  Projections = inprod(Phi, cdata)
  Proj_sq = Projections^2
  Psi_diag = diag(PSI)
  THETA = Proj_sq - Psi_diag
  theta = matrix(THETA[1:d,], ncol = N, nrow = d)
  Sigma_d = LongRunCovMatrix(theta, h=2)

  asymp = function(){
    s = floor(N * delta)
    WW = matrix(0, d, N-s)
    for (i in 1:d){
      WW[i, ] = BBridge(0, 0, 0, 1, N)[(s+1):N]
    }
    WW_2 = WW^2
    val = colSums(WW_2)
    max(val)
  }
  Values_j = sapply(1:M, function(i) asymp())

  s = floor(delta*N)
  Tn_j =  c(rep(0,s))
  for (k in (s+1):N){
    lam_i = partial_cov(cdata, k/N)$eigen_val[1:d]
    kapa = lam_i - (k/N)*lambda[1:d]
    Tn_j[k] = N*t(kapa)%*%solve(Sigma_d)%*%kapa
  }
  Sn_j = max(Tn_j)
  k_star = min(which.max(Tn_j))
  z_j = Sn_j <= Values_j
  p_j = length(z_j[z_j==TRUE])/length(z_j)
  l1 = pca.fd(fdobj[1:k_star], nharm = D, centerfns = T)$values
  l2 = pca.fd(fdobj[(1+k_star):N], nharm = D, centerfns = T)$values
  list(change = k_star/N,
       pvalue = p_j,
       eval_before = l1,
       eval_after = l2)
}

