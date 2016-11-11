#' Check is class
#'
#' Performs a check to see inheritance
#' @param x A \code{diag_resid} object
#' @return A \code{boolean} indicating \code{TRUE} or \code{FALSE}
is.diag_resid = function(x){
  inherits(x, "diag_resid")
}

#' @importFrom ggplot2 autoplot
NULL
