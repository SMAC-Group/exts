#' EDA on Time Series
#'
#' Computes P/ACF and packages it alongside the TS
#' @param xt     A data set
#' @param robust A \code{bool} indicating whether to robustly estimate only the
#' ACF (not PACF).
#' @export
eda_ts = function(xt, lag.max = 30, robust = FALSE, ...){

  if(!is.gts(xt)){ xt = gts(xt) }

  structure(list(data = xt,
                 acf = emp_acf(xt, lag.max = lag.max, robust = robust),
                 pacf = emp_pacf(xt, lag.max = lag.max)),
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

#' Viable Transforms
#'
#' Provides the ability to graph multiple time series under different
#' conditions
#' @param x     A data object
#' @param t1,t2 A \code{function} to transform the data
#' @param both  A \code{boolean} to apply both transforms to the data
#' e.g \code{t1(t2(x))}  (\code{TRUE}) or not e.g. \code{t2(x)} (\code{FALSE}).
#' @return A \code{eda_change} object with
#' \itemize{
#' \item{original}{Untouched data}
#' \item{transform1}{First Transform}
#' \item{transform2}{Second Transform}
#' }
#' @export
eda_change = function(x, t1 = diff, t2 = log, both = FALSE){
  o = gts(x, start = 1947, freq = 4)

  val1 = gts(t1(x), start = 1947, freq = 4)

  if(both){
    val2 = gts(t1(t2(o)))
  } else{
    val2 = gts(t2(o))
  }

  # Extend to add in type of transforms...
  structure(list(original = o, transform1 = val1, transform2 = val2),
            both = both, class = c("eda_change","list"))
}

#' Plot Different Time Series Transforms
#'
#' Graphs time series transforms together with original data
#' @param x,object A \code{\link{eda_change}} object
#' @export
#' @rdname plot.eda_change
plot.eda_change = function(x, ...){
  autoplot.eda_change(object = x, ...)
}


#' @export
#' @rdname plot.eda_change
autoplot.eda_change = function(object, ...) {

  p1 = autoplot(object$original)

  p2 = autoplot(object$transform1)

  p3 = autoplot(object$transform2)

  grid.arrange(p1, p2, p3, nrow = 3)
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
