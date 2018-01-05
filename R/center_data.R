#' Center Functional Data With Change
#'
#' This function centers the functional data by subtracting the pointwise mean from each of the functions in a functional data
#' by taking into account a potential change in the mean function. If there is a change in the mean function, the location of the change
#' is estimated using a fully functional estimator implemented in \code{change_FF}, and the mean before and after the change is computed and subtracted from the respective
#' part of the functional data.
#'
#' @param fdobj A functional data object
#' @param change If \code{TRUE} centering is done by considering the mean change, if \code{FALSE}, the global mean function is subtracted
#' from the functional data. The default is \code{change=TRUE}.
#' @import stats
#' @import fda
#' @export
#'
#' @return Centered functional data sample (class \code{fd}) containing:
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
#' @seealso \code{\link{center.fd}}
#'
#' @examples # Generate FAR(1) process with change in the mean
#' f_AR = fun_AR(n=100, nbasis=21, kappa=0.9)
#' f_AR_change = insert_change(f_AR, k=20, change_location = 0.5, SNR=5)
#' fdata = f_AR_change$fundata
#' c_fdata = center_data(fdata)
#' par(mfrow=c(1,2))
#' plot(fdata, main="Functional Data")
#' plot(c_fdata, main="Centered Functional Data")




center_data = function(fdobj, change=TRUE){
  if(change==TRUE){
    fdobj = center.fd(fdobj)
    basis = fdobj$basis
    samp = fdobj$coefs
    N = ncol(samp)
    D = nrow(samp)
    Sn=(1:N)
    Sn[1]=0
    for(j in (2:N)){
      Sn[j]= sum((rowSums(samp[,1:j]) - (j/N)*rowSums(samp[,1:N]))^2)/N
    }
    k.star = min(which(Sn==max(Sn)))

    dat_b = center.fd(fdobj[1:k.star])
    dat_a = center.fd(fdobj[(1+k.star):N])

    c_dat = cbind(dat_b$coefs, dat_a$coefs)
    c_fdata = fd(c_dat, basis)
  }
  if (change==FALSE){
    c_fdata = center.fd(fdobj)
  }
  c_fdata
}
