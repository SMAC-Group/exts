#' EDA on Time Series
#'
#' Computes P/ACF and packages it alongside the TS
#' @param xt     A data set
#' @param robust A \code{bool} indicating whether to robustly estimate only the
#' ACF (not PACF).
#' @export
eda_ts = function(xt, robust = FALSE, ...){

  if(!is.gts(xt)){ xt = gts(xt) }

  structure(list(data = xt, acf = emp_acf(xt, robust = robust), pacf = emp_pacf(xt)),
            class = "eda_ts")
}


#' Visualize EDA on Time Series
#'
#' Generates a time series EDA plot
#' @param x,object A \code{eda_ts} object
#' @rdname plot.eda_ts
#' @export
plot.eda_ts = function(x, ...){
  autoplot.eda_ts(object = x, ...)
}

#' @rdname plot.eda_ts
#' @export
autoplot.eda_ts = function(object, ...){
  grid.arrange(plot(object$data), plot(object$acf), plot(object$pacf), nrow = 1)
}


#' Calculate both Classical and Robust ACF
#'
#' Computes Classical & Robust ACF and packages it alongside the TS
#' @param xt A data set
#' @export
compare_acf = function(xt, ...){

  if(!is.gts(xt)){ xt = gts(xt) }

  structure(list(data = xt,
                 acf = emp_acf(xt, robust = FALSE),
                 racf = emp_acf(xt, robust = TRUE)),
            class = "eda_ts")
}


#' Visualize EDA on Time Series
#'
#' Generates a time series EDA plot
#' @param x,object A \code{eda_ts} object
#' @rdname plot.eda_ts
#' @export
plot.eda_ts = function(x, ...){
  autoplot.eda_ts(object = x, ...)
}

#' @rdname plot.eda_ts
#' @export
autoplot.eda_ts = function(object, ...){
  grid.arrange(plot(object[[1]]), plot(object[[2]]), plot(object[[3]]), nrow = 1)
}
