#' Exclude function
#'
#' this function is used in \code{opt_bandwidth} function
#' @name exclude
#' @keywords internal

'%!in%' <- function(x,y){!('%in%'(x,y))}
