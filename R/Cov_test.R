
#' Testing the Equality of Covariance Operators in Functional Samples
#'
#' This function tests  for the equality of the covariance structures in two functional samples.
#'
#'
#' @param fdobj1 A functional data object of class '\code{fd}'
#' @param fdobj2 A functional data object of class '\code{fd}'
#' @param d Level ofdimension reduction needed to represent the data. One can use \code{pick_dim}
#'
#' @import sde
#' @import fda
#' @export
#'
#'
#'
#' @return
#'\item{\code{pvalue}}{
#' Approximate p value for testing equality of the covariance structures in two functional samples
#'}
#'
#'@details  test for the equality of the covariance structures
#'in two functional samples. The test statistic has a chi-square asymptotic distribution
#'with a known number of degrees of freedom, which depends on the level of
#'dimension reduction needed, \code{d}, to represent the data.
#'
#' @seealso \code{\link{change_fPCA}}
#'
#' @references Fremdt S., Horvath L., Kokoszka P., Steinebach J. (2017+), \emph{Testing the Equality of Covariance Operators in
#' Functional Samples} Scandinavian Journal of Statistics, 2013.
#'
#' @examples
#' # generate functional data
#' fdata1 = fun_IID(n=100, nbasis=21)
#' fdata2 = fun_IID(n=150, nbasis=21)
#' Cov_test(fdata1, fdata2, d=5)
#'
#'
#'

Cov_test = function(fdobj1, fdobj2, d){
  if (class(fdobj1)!="fd"){stop("Must insert functional data")}
  if (class(fdobj2)!="fd"){stop("Must insert functional data")}
  D1 = nrow(fdobj1$coefs)
  D2 = nrow(fdobj2$coefs)
  if(D1!=D2){stop("Fundtional data have different number of basis functions")}
  N = ncol(fdobj1$coefs)
  M = ncol(fdobj2$coefs)
  cdata1 = center.fd(fdobj1)
  cdata2 = center.fd(fdobj2)
  D = D1
  th = N/(N+M)
  m_pooled = matrix(0, D, N+M)
  for (j in 1:N){
    m_pooled[, j] = cdata1[j]$coefs
  }
  for (j in (1+N):(N+M)){
    m_pooled[, j] = cdata2[j-N]$coefs
  }
  f_pooled = fd(m_pooled, fdobj1$basis)
  phi = pca.fd(f_pooled, nharm = D, centerfns = T)$harmonics
  Lam1 = inprod(cdata1, phi)^2
  l1 = colMeans(Lam1)
  Lam2 = inprod(cdata2, phi)^2
  l2 = colMeans(Lam2)
  Inner1 = inprod(cdata1, phi)
  Inner2 = inprod(cdata2, phi)
  res = matrix(0, d, d)
  for (i in 1:d){
    for (j in 1:d){
      num = mean(sapply(1:N, function(l) Inner1[l,i]*Inner1[l,j])) - mean(sapply(1:M, function(l) Inner2[l,i]*Inner2[l,j]))
      denom = (th * l1[i] + (1-th) * l2[i]) * (th * l1[j] + (1-th) * l2[j])
      res[i, j] = num^2/denom
    }
  }
  T_hat = (N+M)/2 * th * (1-th) * sum(res)
  1 - pchisq(T_hat, df = (d*(d+1))/2)
}
