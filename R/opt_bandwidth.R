
#' Optimal Bandwidth Selection for the Long Run Covariance Estimation
#'
#' This function estimates an optimal window parameter for long run covariance operator estimation  in functional time series using the method of Rice G. and Shang H. L. (2017)
#'
#' @param fdobj Functional data object, class of \code{"fd"}
#' @param kern_type Kernel that is used for the long run covariance estimation. The available options are
#' \code{c("BT", "PR", "TH", "QS")} where \code{"BT"} is Bartlett, \code{"PR"} is Parzen, \code{"TH"} is Tukey-Hanning, and \code{"QS"} is Quadratic Spectral kernel.
#' @param kern_type_ini Initial Kernel function to start the optimal bandwidth search
#' @param is_change If \code{TRUE} then the data is centered considering the change in the mean function
#' @param ... Further arguments to pass
#'
#' @import sandwich
#' @import fda
#' @export
#'
#'@return
#' \item{\code{hat_h_opt}}{
#' Estimated optimal bandwidth
#'}
#'\item{\code{C_0_est}}{
#' Estimated Long run covariance kernel using the optimal bandwidth  \code{hat_h_opt}
#'}
#'@references Rice G. and Shang H. L. (2017), \emph{A plug-in bandwidth selection procedure for
#' long run covariance estimation with stationary functional time series}, Journal of Time Series Analysis, 38(4), 591-609
#'
#' @seealso \code{\link{LongRun}}
#'
#'@examples
#'fdata1 = fun_AR(n=100, nbasis=21, order=1, kappa=0.8)
#'opt_bandwidth(fdata1, "PR", "BT")


opt_bandwidth = function(fdobj, kern_type, kern_type_ini, is_change = TRUE, ...){
  if(kern_type%!in% c("BT", "PR", "TH", "QS")){stop("Provide a valid kernel!")}
  if (class(fdobj)!="fd"){stop("Must input Functional data")}
  center_dat = center_data(fdobj, change = is_change)$coefs
  N = ncol(center_dat)
  band_ini = N^(1/3)
  kern_name_kweights_ini = switch(kern_type_ini, BT = "Bartlett", PR = "Parzen",
                              TH = "Tukey-Hanning", QS = "Quadratic Spectral")

  gamma_l=function(lag){
    gamma_lag_sum = 0
    if(lag >= 0){
      for(ij in 1:(N-lag)){
        gamma_lag_sum = gamma_lag_sum + t(as.matrix(center_dat[,ij]) %*% t(as.matrix(center_dat[,(ij+lag)])))
      }
    } else {
      for(ij in 1:(N+lag)){
        gamma_lag_sum = gamma_lag_sum + t(as.matrix(center_dat[,ij-lag]) %*% t(as.matrix(center_dat[,ij])))
      }
    }
    return(gamma_lag_sum/N)
  }

  cov_l = function(porder, band_ini, kern_type){
    cov_sum = gamma_l(0)
    for(ik in 1:(N-1)){
      cov_sum = cov_sum + kweights(ik/band_ini, kernel = kern_type) * abs(ik)^(porder) * (gamma_l(ik) + t(gamma_l(ik)))
    }
    return(cov_sum)
  }

  kern_name_kweights = switch(kern_type, BT = "Bartlett", PR = "Parzen",
                              TH = "Tukey-Hanning", QS = "Quadratic Spectral")

  w_weights = switch(kern_type, BT = 1, PR = 6, TH = pi^2/4, QS = 18*pi*pi/125)
  q = switch(kern_type, BT = 1, PR = 2, TH = 2, QS = 2)

  C_0 = cov_l(porder = 0, band_ini=band_ini, kern_type = kern_name_kweights_ini)
  C_2 = w_weights * cov_l(porder = q, band_ini=band_ini, kern_type = kern_name_kweights_ini)
  first_part = (2 * q * sum(C_2^2))^(1/(1 + 2*q))
  kernel_square_int = switch(kern_type, BT = 2/3, PR = 0.539285, TH = 3/4, QS = 1, FT = 4/3) #(QS kernel p 822 Andrews(1991))
  second_part = ((sum(C_0^2) + sum(diag(C_0))^2) * kernel_square_int)^(-1/(1+2*q))
  c_0 = first_part * second_part
  hat_h_opt =  c_0 * (N^(1/(1+2*q)))
  C_0_est = cov_l(porder = 0, band_ini = hat_h_opt, kern_type = kern_name_kweights)
  return(list(hat_h_opt = hat_h_opt, C_0_est = C_0_est))
}



