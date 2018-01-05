#' @importFrom graphics par legend lines
NULL

#' Change Point Analysis Of Functional Data Via Dimension Reduction
#'
#' This function tests whether there is a significant change in the mean function of functional data, and it gives an estimate
#' of the location of the change. The procedure will reduce the dimension of the functional data using functional
#' principal component analysis and will use \code{d} leading principal curves to carry out the change point
#' analysis. The projection dimension \code{d} can be chosen via total variation explained (TVE) using the function
#' \code{\link{pick_dim}}.
#'
#' @param fdobj A functional data object of class '\code{fd}'
#' @param d Number of principal components
#' @param M Number of monte carlo simulations to get the critical values. The default value is \code{M=1000}
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
#' @return
#'\item{\code{pvalue}}{
#' An approximate p value for testing whether there is a significant change in the mean function
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
#' @seealso \code{\link{change_FF}}
#'
#'@details This functions performs structural break analysis for the functional data using an fPCA based initial dimension reduction. It
#'is recommended that the dimension of the subspace, \code{d}, that the functional observations are projected onto should be selected based on
#'TVE using \code{\link{pick_dim}}.
#'
#'@references Berkes, I., Gabrys, R.,Hovarth, L. & P. Kokoszka (2009)., \emph{Detecting changes in the mean of functional observations}
#' Journal of the Royal Statistical Society, Series B 71, 927–946
#' @references Aue, A., Gabrys, R.,Hovarth, L. & P. Kokoszka (2009)., \emph{Estimation of a change-point in the mean function
#' of functional data}Journal of Multivariate Analysis 100, 2254–2269.
#'
#' @examples
#' # generate functional data
#' fdata = fun_IID(n=100, nbasis=21)
#' # insert an artifiical change
#' data_c = insert_change(fdata, k=21, change_location = 0.5, SNR=1)$fundata
#' d.hat = pick_dim(data_c, 0.9)$d
#' change_fPCA(data_c, d=d.hat)$change



change_fPCA = function(fdobj, d, M=1000, h=0, plot=FALSE, ...){
   if (class(fdobj)!="fd"){stop("Must insert functional data")}
   D = nrow(fdobj$coefs)
   n = ncol(fdobj$coefs)
   fdata = center.fd(fdobj)
   mean_fun = mean(fdobj)
   FpCa <- pca.fd(fdata, nharm=D, centerfns=TRUE)
   eta.hat = matrix(t(FpCa$scores)[1:d, ], nrow=d)
   S_n = function(k){
      if(d==1){
         eta.bar = sum(eta.hat)/n
         out = sum(eta.hat[,1:k]) - k*eta.bar
      } else {
         eta.bar = as.matrix(rowSums(eta.hat)/n)
         out = rowSums(as.matrix(eta.hat[, 1:k])) - k*eta.bar
      }
      out
   }
   Sigma.hat = LongRun(fdobj=fdata, h=h)$covm[1:d, 1:d]
   TT = sapply(1:n, function(k) 1/n *(t(S_n(k))%*%solve(Sigma.hat)%*%S_n(k)))
   Tn = max(TT)
   k.star = min(which(TT==max(TT)))
   asymp <- function(N){
      B.Bridges= matrix(0,d,(N))
      for(j in (1:d)){
         B.Bridges[j,]=BBridge(0,0,0,1,N-1)^2
      }
      max(colSums(B.Bridges))
   }
   Values = sapply(1:M, function(k) asymp(n))
   z = Tn<=Values
   p = length(z[z==TRUE])/length(z)
   dat.b = fdobj[1:k.star]
   dat.a = fdobj[(k.star+1):n]
   mean.b = mean(dat.b)
   mean.a = mean(dat.a)
   delta = mean.a - mean.b
   if (plot == TRUE){
      par(mfrow=c(1,2))
      plot(fdobj, col="grey", main="Functional Data")
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
