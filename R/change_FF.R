
#' @importFrom graphics par legend lines
NULL

#' Change Point Analysis Of Functional Data Without Dimension Reduction (Fully Functional)
#'
#' This function tests whether there is a significant change in the mean function of the functional data, and it will
#' give an estimate for the location of the change. The procedure is based on the standard L-2 norm and hence does not depend on any dimension
#' reduction technique such as fPCA.
#'
#'
#' @param fdobj A functional data object of class '\code{fd}'
#' @param M Number of monte carlo simulations used to get the critical values. The default value is \code{M=1000}
#' @param  h The window parameter for the estimation of the long run covariance kernel. The default
#' value is \code{h=0}, i.e., it assumes iid data
#' @param plot If \code{TRUE} plot of the functional data before and after the estimated change and plot of the
#' estimated change function is given
#' @param ... Further arguments to pass
#'
#' @import sde
#' @import fda
#' @export
#'
#'
#'
#' @return
#'\item{\code{pvalue}}{
#' Approximate p value for testing whether there is a significant change in the mean function
#'}
#'\item{\code{change}}{
#' Estimated change location
#'}
#'\item{\code{DataBefore}}{
#' Data before the estimated change
#'}
#'\item{\code{DataAfter}}{
#' Data after the estimated change
#'}
#'\item{\code{MeanBefore}}{
#' Mean function before the estimated change
#'}
#'\item{\code{MeanAfter}}{
#' Mean function after the estimated change
#'}
#'\item{\code{change_fun}}{
#' Estimated change function
#'}
#'
#'@details This function dates and detects changes in the mean function of functional data using a fully functional technique that does not
#'dependent on dimension reduction. For more details, see Aue, Rice, Sonmez (2017+)
#'
#' @seealso \code{\link{change_fPCA}}
#'
#' @references Aue A., Rice G., Sonmez O. (2017+), \emph{Detecting and dating structural breaks in
#' functional data without dimension reduction} (https://arxiv.org/pdf/1511.04020.pdf)
#'
#' @examples
#' # generate functional data
#' fdata = fun_IID(n=100, nbasis=21)
#' # insert an artifiical change
#' data_c = insert_change(fdata, k=21, change_location = 0.5, SNR=1)$fundata
#' change_FF(data_c)$change
#'
#'
#'


change_FF = function(fdobj, M=1000, h=0, plot=FALSE, ...){
   if (class(fdobj)!="fd"){stop("Must insert functional data")}
   fdata = center.fd(fdobj)
   mean_fun = mean(fdobj)
   basis = fdata$basis
   samp = fdata$coefs
   N = ncol(samp)
   D = nrow(samp)
   Sn=(1:N)
   Sn[1]=0
   for(j in (2:N)){
      Sn[j]= sum((rowSums(samp[,1:j]) - (j/N)*rowSums(samp[,1:N]))^2)/N
   }
   k.star = min(which(Sn==max(Sn)))
   Tn = max(Sn)
   LongRunC = LongRun(fdobj=fdata, h=h)
   lambda = LongRunC$e_val
   asymp <- function(N){
      BridgeLam= matrix(0,D,N)
      for(j in (1:D)){
         BridgeLam[j,]=lambda[j]*(BBridge(0,0,0,1,N-1)^2)
      }
      max(colSums(BridgeLam))
   }
   Values = sapply(1:M, function(k) asymp(N))
   z = Tn<=Values
   p = length(z[z==TRUE])/length(z)
   dat.b = fdobj[1:k.star]
   dat.a = fdobj[(k.star+1):N]
   mean.b = mean(dat.b)
   mean.a = mean(dat.a)
   delta = mean.a - mean.b
   if (plot == TRUE){
      par(mfrow=c(1,2))
      plot(fdobj, col="grey", main="Functional Data", ylab="values")
      lines(dat.a, col="pink")
      lines(dat.b, col="lightblue")
      lines(mean.b, col="blue")
      lines(mean.a, col="red")
      legend("topleft", c("before", "after"), col=c("blue", "red"), lty=c(1,1), cex=0.5)
      plot(delta, main="Estimated Change Function", ylab="values")
      list(pvalue = p , change = k.star,
           DataBefore = dat.b, DataAfter = dat.a,
           chanfe_fun = delta)
   }else{
      list(pvalue = p , change = k.star,
           DataBefore = dat.b, DataAfter = dat.a,
           MeanBefore = mean.b, MeanAfter = mean.a,
           chanfe_fun = delta)
   }
}
