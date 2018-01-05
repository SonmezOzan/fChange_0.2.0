
#' Number Of Principal Component Selection Based On Variation
#'
#' This function selects number of principal components based on the total variation explained (TVE). The functional data is
#' projected onto a smaller number of principal curves via functional Principal Component Analysis (fPCA). This function
#' picks the dimension of the projection space based on the desired total variation explained.
#'
#' @param fdobj Functional data object, class of \code{"fd"}
#' @param TVE Total Variation Explained. It must be in [0,1].
#'
#' @import stats
#' @import fda
#' @export
#'
#' @return
#'\item{\code{d}}{
#' Minimum number of Principle components needed in order to reach desired \code{TVE}
#'}
#'\item{\code{TVEs}}{
#' Vector of TVEs. This has the same length of number of basis that the functional data is represented
#'}
#'
#'@details This function is used to determine the dimension of the space that the functional data is projected onto
#'based on the variation. One of the common treatments of the functional data is to transform it to multivariate
#'objects using so called score vectors and utilize the multivariate techniques. This function will enable users
#'to pick the dimension of the score vectors.
#'
#'@references Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
#' Data Analysis, 2nd ed.}, Springer, New York
#'
#' @seealso \code{\link{pca.fd}}
#'
#'@examples
#'fdata1 = fun_IID(n=100, nbasis=21)
#'pick_dim(fdata1, 0.95)
#'@examples
#'fdata2 = fun_IID(n=100, nbasis=21, Sigma=3^-(1:21))
#'pick_dim(fdata2, 0.95)
#'
#'
pick_dim = function(fdobj, TVE){
   if (class(fdobj)!="fd"){stop("Must input Functional data")}
   if (TVE > 1 || TVE<0){stop("TVE must be in [0,1]")}
   D = nrow(fdobj$coefs)
   Fpca <- pca.fd(fdobj, nharm=D, centerfns=TRUE)
   Lambda = Fpca$values
   tves = sapply(1:D, function(j) sum(Lambda[1:j])/sum(Lambda))
   d = length(which(tves<TVE))+1
   list(d=d, TVEs=tves)
}
