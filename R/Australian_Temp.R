#'  Australian Climate Data
#'
#'
#' Australian daily minimum temperature climate data for 8 different stations are provided. The data is taken from Australian Government Bureau of Meteorology.
#'
#' @docType data
#' @name Australian_Temp
#'
#' @import fda
#' @import reshape2
#'
#' @return Australian daily minimum temperature climate data for 8 different stations:
#'\item{Sydney}{
#' Sydney (Observatory Hill), taken from 1859 to 2012
#'}
#'\item{Melbourne}{
#' Melbourne (Regional Office), taken from 1855 to 2012
#'}
#'\item{Boulia}{
#' Boulia Airport, taken from 1888 to 2012
#'}
#'\item{Cape_Otway}{
#' Cape Otway Lighthouse, taken from 1864 to 2012
#'}
#'\item{Gayndah}{
#' Gayndah Post Office, taken from 1893 to 2009
#'}
#'\item{Gunnedah}{
#' Gunnedah Pool, taken from 1876 to 2011
#'}
#'\item{Hobart}{
#' Hobart (Ellerslie Road), taken from 1882 to 2012
#'}
#'\item{Robe}{
#' Robe Comparison, taken from 1884 to 2012
#'}
#'
#'
#' @format An object of class \code{data.table} or \code{data.frame}.
#'
#'@keywords datasets
#'
#' @references Australian Government Bureau of Meteorology
#'
#' @source The daily observations are available from \href{http://www.bom.gov.au/climate/data}{http://www.bom.gov.au/climate/data}.
#' Copyright Commonwealth of Australia 2010, Bureau of Meteorology. Definitions adapted from
#' \href{http://www.bom.gov.au/climate/dwo/IDCJDW0000.shtml}{http://www.bom.gov.au/climate/dwo/IDCJDW0000.shtml}
#'
#' @examples
#'library(fda)
#'library(reshape2)
#'fun_data_S = Australian_Temp$Sydney
#'D = 21
#'basis = create.fourier.basis(rangeval = c(0, 1), nbasis = D)
#'nas = which(is.na(fun_data_S$Days.of.accumulation.of.minimum.temperature))
#'fun_data_S = fun_data_S[-nas, ]
#'yy = unique(fun_data_S$Year)
#'mat.S = matrix(0, D, length(yy))
#'for (i in 1:length(yy)){
#'  aa = subset(fun_data_S, Year==yy[i])
#'  cc = aa$Minimum.temperature..Degree.C.
#'  a = which(is.na(cc))
#'  if (length(a)>0){
#'  cc = cc[-which(is.na(cc))]
#'  }else{
#'    cc = cc
#'  }
#'  f_Obs = Data2fd(argvals=seq(0, 1, length = length(cc)) , cc, basisobj = basis)
#'  mat.S[, i] = f_Obs$coefs
#'}
#'fdata = fd(mat.S, basis)
#'# note that the last year, has data only up to 6 months
#'# therefore we remove it
#'fdata = fdata[-length(yy)]
#'plot(fdata)


"Australian_Temp"
