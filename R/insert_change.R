#' Insert A Change In the Mean Function Of Functional Data
#'
#' This function inserts a change in the mean function to a given functional data sample. The change function
#' can either be directly defined by the user or it can be generated based on the sum of first
#' \code{k} basis functions defined by the \code{fdobj}. Once the change function is defined the change is inserted at the
#' defined change location with a signal magnitude defined by signal to noise ratio, \code{SNR}.
#' For more details on how these quantities are defined. See Aue, Rice and Sonmez (2017+).
#'
#' @param fdobj Functional data object of class \code{'fd'}
#' @param change_fun Self defined change function. It has to be a functional data object having the same
#' number of basis functions.
#' @param k Number of basis functions to be summed to construct the change function. It should
#' be used when \code{change_fun} is not defined. It has to be less than number of basis functions.
#' @param change_location Location of the change to be inserted. It is scaled to be in [0,1].
#' @param SNR Signal to Noise Ratio to determine the magnitude of the change function that is being
#' inserted.
#' @param plot Plots the functional data before (blue) and after (red) the change.
#' @param ... Further information to pass
#'
#' @import stats
#' @import fda
#' @export
#'
#'@return \code{fundata:} functional data with an inserted change in the mean function
#'@return \code{change_fun:} inserted change function
#'@return \code{plot:} of the functional data with inserted change
#'
#'@details This function should only be used to artificially insert a change function to the mean of
#'the functional data set either by defining a specific change function or generating the change function
#'based from the basis functions.
#'
#'@references Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
#' Data Analysis, 2nd ed.}, Springer, New York.
#' @references Aue A., Rice G., Sonmez O. (2017+), \emph{Detecting and dating structural breaks in
#' functional data without dimension reduction} (https://arxiv.org/pdf/1511.04020.pdf)
#' @seealso \code{\link{Data2fd}, \link{fun_IID}, \link{fun_MA}, \link{fun_AR}}
#'
#' @examples ####################################
#' #first generate FAR(1) process
#' fdata = fun_AR(n=100, nbasis=25, Sigma=2^-(1:25))
#' # insert the change which is the sum of first 3 basis functions
#' # in the middle of the data with SNR=2
#' insert_change(fdata, k=3, change_location=0.5, SNR=2)
#'
#' @examples ####################################
#' #first generate FAR(1) process
#' fdata = fun_AR(n=100, nbasis=25, Sigma=2^-(1:25))
#' # insert the change which is the 20th onservation
#' # in the middle of the data with SNR=2
#' insert_change(fdata, change_fun = fdata[20], change_location=0.5, SNR=2)

insert_change = function(fdobj, change_fun=NULL, k=NULL, change_location, SNR, plot=TRUE, ...){
   if (class(fdobj)!="fd"){stop("Data is not a 'fd' class ")}
   n = ncol(fdobj$coefs)
   nbasis = nrow(fdobj$coefs)
   basis = fdobj$basis
   fun_mu = center.fd(fdobj)
   cdata = center.fd(fdobj)
   dat = cdata$coefs
   tr = sum(diag(cov(t(dat))))
   theta = change_location
   c = SNR*tr/(theta*(1-theta)*sqrt(nbasis))
   if (is.null(k) & is.null(change_fun)){
      stop("Define the change function")
   }
   if (is.null(change_fun) & is.null(k)==FALSE){
      change_coef = c(rep(1, k), rep(0, nbasis-k))
      v_hat = fd(change_coef, fdobj$basis)
      Change = v_hat$coefs*sqrt(c/k)
   }else{
      v_hat = change_fun
      kk = inprod(v_hat, v_hat)[1,1]
      Change = v_hat$coefs*sqrt(c/kk)
   }
   newdata = dat
   if (theta>1 || theta<0){stop("Change location must be in [0,1]")}
   x = theta*n
   for (i in (x+1):n){
      newdata[ ,i] = dat[ ,i] + Change
   }
   fdata = fd(newdata, basis) + fun_mu
   if (plot==FALSE){
      list(fundata = fdata, change_fun = v_hat)
   }else{
      plot(fdata, col="grey")
      lines(fdata[(x+1):n], col="red")
      lines(fdata[1:x], col="blue")
      list(fundata = fdata, change_fun = v_hat)
   }
}
