
#' Long Run Covariance Operator Estimation for Functional Time Series
#'
#' This function estimates the long run covariance operator of a given functional data sample and its estimated eigenelements.
#'
#' @param fdobj A functional data object
#' @param h The bandwidth parameter. It is strictly non-zero. Choosing the bandwidth parameter to be zero is identical
#' to estimating covariance operator assuming iid data.
#' @param kern_type Kernel function to be used for the estimation of the long run covariance
#' function. The choices are \code{c("BT", "PR", "SP", "FT")} which are respectively, bartlett, parzen, simple and flat-top kernels.
#' By default the function uses a \code{"barlett"} kernel.
#' @param is_change If \code{TRUE} then the data is centered considering the change in the mean function.
#' @param ... Further arguments to pass
#' @import stats
#' @import fda
#' @import lattice
#' @export
#'
#' @return
#'\item{\code{e_fun}}{
#' Eigenfunctions of the estimated long run covariance function
#'}
#'\item{\code{e_val}}{
#' Eigenvalues of the estimated long run covariance function
#'}
#'\item{\code{covm}}{
#' Coefficient matrix of the estimated long run covariance operator.
#'}
#'\item{\code{contour_plot}}{
#' The estimated covariance function \eqn{C(t,s)} surface plot if \code{plot=TRUE}
#'}
#'
#'
#' @seealso \code{\link{opt_bandwidth}}
#'
#' @references Aue A., Rice G., Sonmez O. (2017+), \emph{Detecting and dating structural breaks in
#' functional data without dimension reduction} (https://arxiv.org/pdf/1511.04020.pdf)
#' @references Rice G. and Shang H. L. (2017), \emph{A plug-in bandwidth selection procedure for
#' long run covariance estimation with stationary functional time series}, Journal of Time Series Analysis, 38(4), 591-609
#'
#' @examples
#' # Generate FAR(1) process
#' fdata = fun_AR(n=100, nbasis=31, order=1, kappa=0.9)
#' # Estimate the Long run covrariance
#' C_hat = LongRun(fdata, h=2)
#' C_hat$e_fun # eigenfunctions of Long Run Cov
#' C_hat$e_val # eigenvalues of Long Run Cov
#' C_hat$covm # Estimated covariance matrix
#'


LongRun <- function(fdobj, h, kern_type = "BT", is_change = TRUE, ...){
  if (class(fdobj)!="fd"){stop("Must input Functional data")}
  if (h<0) {stop("h must be non-negative")}
  kerneltype = switch(kern_type, BT = "Bartlett", PR = "Parzen",
                      FT = "flat_top", SP = "Simple")
  N = ncol(fdobj$coefs)
  D = nrow(fdobj$coefs)
  basis = fdobj$basis
  Kernel <- function(i, h){
    x = i/h
    if (kerneltype == "flat"){
      return(1)
    }
    if (kerneltype == "Simple"){
      return(0)
    }
    if (kerneltype == "Bartlett"){
      return(1 - x)
    }

    if (kerneltype == "flat_top"){
      if (x < 0.1){
        return(1)
      } else {
        if (x >= 0.1 & x < 1.1){
          return(1.1 - x)
        }else{
          return(0)
        }
      }
    }
    if (kerneltype == "Parzen"){
      if (x < 1/2){
        return(1 - 6 * x^2 + 6 * abs(x)^3)
      }else{
        return(2 * (1 - abs(x))^3)
      }
    }
  }
  D_mat = matrix(0, D, D)
  fdobj_centered = center_data(fdobj, change = is_change)
  # Long Run Cov Est
  for (k in 1:D) {
    for (r in k:D) {
      s = fdobj_centered$coefs[k, 1:N] %*% fdobj_centered$coefs[r, 1:N]
      if (h > 0) {
        for (i in 1:h) {
          a = fdobj_centered$coefs[k, 1:(N - i)] %*% fdobj_centered$coefs[r, (i + 1):N]
          a = a + fdobj_centered$coefs[r, 1:(N - i)] %*% fdobj_centered$coefs[k, (i + 1):N]
          s = s + Kernel(i, h) * a
        }
      }
      D_mat[k, r] = s
      D_mat[r, k] = D_mat[k, r]
    }
  }
  D_mat = D_mat/N
  eigen_struct = eigen(D_mat, symmetric = TRUE)
  eigenfunc = fd(eigen_struct$vectors, basisobj = basis)
    a = seq(0, 1, length.out = D)
    Psi = eval.basis(a, fdobj$basis)
    C_ts = matrix(0, D, D)
    for (t in 1:D){
      for (s in 1:D){
        C_ts[t,s] = sum(sapply(1:D, function(i) sapply(1:D, function(j) D_mat[i,j]*Psi[t,i]*Psi[s,j])))
      }
    }
    colnames(C_ts) = rownames(C_ts) = a
    CC = data.frame(melt(C_ts))
    colnames(CC) = c("t", "s", "C(t,s)")
    p = wireframe(CC[,3] ~ t * s, CC, shade = TRUE, aspect = c(1, 1),
              light.source = c(1,5,10), main = "C(t,s)", zlab="C(t,s)",
              scales = list(z.ticks=5,arrows=FALSE, col="black", font=10, tck=0.5),
              screen = list(z = 40, x = -65, y = 0))
  list(e_fun = eigenfunc, e_val = abs(eigen_struct$values), covm = D_mat, contour_plot = p)

}
