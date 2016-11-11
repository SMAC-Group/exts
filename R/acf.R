# hidden casting function
cast_acf = function(object, n, name_ = "Empirical", type = "Autocorrelation",
                    class = "theo_arma"){

  # Force to array
  if(is.null(nrow(object))){ dim(object) = c(length(object), 1, 1)}

  # Make pretty names
  ids = seq_len(nrow(object))
  if(type == "Autocorrelation" || type == "Autocovariance"){
     ids = ids - 1
  }
  dimnames(object)  = list(ids, name_, name_)

  structure(object, type = type, n = n, class = c(class,"ACF","array"))
}


#' Plot ACF-like graphs
#'
#' Generic ACF plot function for empirical and theoretical graphs
#' @param x,object  A object from \code{\link{emp_acf}}, \code{\link{emp_pacf}},
#' \code{\link{theo_acf}}, or \code{\link{theo_pacf}}.
#' @param show.ci   A \code{bool} indicating whether to show confidence region
#' @param ci        A \code{double} containing the 1-alpha level. Default is 0.95
#' @param ...       Additional parameters
#' @rdname plot.acf_plot
#' @export
#' @examples
#' # Calculate the Autocorrelation
#' m = emp_acf(datasets::AirPassengers)
#'
#' # Plot with 95% CI
#' plot(m)
#'
#' # Plot with 90% CI
#' plot(m, ci = 0.90)
#'
#' # Plot without 95% CI
#' plot(m, show.ci = FALSE)
plot.acf_plot = function(x, show.ci = TRUE, ci = 0.95, ...){
  autoplot.acf_plot(object = x, show.ci = show.ci, ci = ci, ...)
}

#' @rdname plot.acf_plot
#' @export
autoplot.acf_plot = function(object, show.ci = TRUE, ci = 0.95, ...){

  # Quiet the warnings...
  Lag = xmin = xmax = ymin = ymax = NULL

  # Wide to long array transform
  x2 = as.data.frame.table(object, responseName = "ACF", stringsAsFactors = FALSE)

  colnames(x2) = c("Lag", "Signal X", "Signal Y", "ACF")

  # Remove character cast
  x2$Lag = as.numeric(x2$Lag)

  # Create key
  x2$key = ifelse(
    x2$`Signal X` == x2$`Signal Y`,
    paste0(x2$`Signal X`),
    paste0(x2$`Signal X`, " & ", x2$`Signal Y`)
  )

  # Plot with facetting
  g = ggplot(data = x2, mapping = aes(x = Lag, y = ACF)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = Lag, yend = 0)) +
    facet_wrap(~ key) + theme_bw()


  if(show.ci){

    clim0 = qnorm( (1 + ci)/2 ) / sqrt(attr(object,'n'))

    ci.region = data.frame(xmin = -Inf, xmax = Inf, ymin = -clim0, ymax = clim0)

    g = g + geom_rect(data = ci.region,
                      aes(xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax),
                      fill = "blue", alpha = 0.10,
                      inherit.aes = FALSE)
  }

  g = g + ylab(attr(object,'type'))

  g
}



############################################
# Theoretical P/ACF
############################################

#' Computes the Theoretical Autocorrelation (ACF) of an ARMA process
#'
#' Acts as a wrapper around `ARMAacf` to calculate the theoretical ACF process
#' of an ARMA process.
#' @param ar A \code{vector} containing the AR coefficients
#' @param ma A \code{vector} containing the MA coefficients
#' @param lag.max A \code{integer} indicating the max length.
#' @examples
#' # Computes the theoretical ACF for an ARMA(1,0) or better known as an AR(1)
#' theo_acf(ARMA(ar = -0.25, ma = NULL))
#' # Computes the theoretical ACF for an ARMA(2, 1)
#' theo_acf(ARMA(ar = c(.50, -0.25), ma = 0.20), lag.max = 10)
#' @importFrom gmwm is.ts.model ARMA
#' @export
theo_acf = function(model = ARMA(ar = c(.50, -0.25), ma = .20), lag.max = 20){
  theo_arma_(model = model, lag.max = lag.max, pacf = FALSE)
}

#' Computes the Theoretical Partial Autocorrelation (PACF) of an ARMA process
#'
#' Acts as a wrapper around `ARMAacf` to calculate the theoretical PACF process
#' of an ARMA process.
#' @inheritParams  theo_acf
#' @export
#' @examples
#' # Computes the theoretical ACF for an ARMA(1,0) or better known as an AR(1)
#' theo_pacf(ARMA(ar = -0.25, ma = NULL), lag.max = 7)
#' # Computes the theoretical ACF for an ARMA(2, 1)
#' theo_pacf(ARMA(ar = ARMA(ar = c(.50, -0.25), ma = .20), lag.max = 10)
theo_pacf = function(model = ARMA(ar = c(.50, -0.25), ma = .20), lag.max = 20){
  theo_arma_(model = model, lag.max = lag.max, pacf = TRUE)
}

# Work horse of the above two functions
theo_arma_ = function(model, lag.max = 20, pacf = FALSE){

  if(!is.ts.model(model)){ stop("`model` must be a `ts.model` object.")}
  if(model$starting){ stop("`model` must have specific parameter values.")}
  if(length(model$desc) != 1 && model$desc != "SARIMA"){ stop("`model` must contain only 1 process.")}


  objdesc = model$objdesc[[1]]

  ar = model$theta[model$process.desc == "AR"]
  ma = model$theta[model$process.desc == "MA"]

  if(is.null(lag.max)){ lag.max = max(length(ar), length(ma) + 1) + 2}

  cast_acf(ARMAacf(ar = ar, ma = ma, lag.max = lag.max, pacf = pacf),
           n = lag.max,
           name_ = paste("Theoretical", if(pacf){ "PACF"} else{"ACF"}),
           type = "Autocorrelation",
           class = "theo_arma")
}


#' Plot Theoretical Autocorrelation (ACF) for ARMA Models
#'
#' Displays the theoretical autocorrelation for ARMA Models
#' @param x,object  An \code{"theo_arma"} object from \code{\link{theo_acf}}
#' or  \code{\link{theo_pacf}}.
#' @param ...       Additional parameters
#' @return An \code{array} of dimensions \eqn{N \times S \times S}{N x S x S}.
#' @rdname plot.theo_arma
#' @export
#' @examples
#' # Compute Theoretical ACF
#' m = theo_acf(ARMA(ar = -0.25, ma = NULL), lag.max = 7)
#'
#' # Compute Theoretical PACF
#' m2 = theo_pacf(ARMA(ar = -0.25, ma = NULL), lag.max = 7)
#'
#' # Plot either the theoretical ACF or PACF
#' plot(m); plot(m2)
plot.theo_arma = function(x, ...){
  autoplot.acf_plot(x, show.ci = FALSE, ...)
}

#' @rdname plot.theo_arma
#' @export
autoplot.theo_arma = function(object, ...){
  autoplot.acf_plot(object, show.ci = FALSE)
}


#' Compute Theoretical ACF and PACF
#'
#' Computes both the ACF and PACF for a given process.
#' @inheritParams  theo_acf
#' @export
#' @examples
#' # Computes the theoretical ACF for an ARMA(2, 1)
#' theo_corr(ARMA(ar = c(.50, -0.25), ma = 0.20), lag.max = 10)
theo_corr = function(model = ARMA(ar = c(.50, -0.25), ma = .20), lag.max = 20){

  structure(list(acf = theo_acf(model, lag.max = lag.max),
                 pacf = theo_pacf(model, lag.max = lag.max)), class = "theo_corr")

}

#' Plot Theoretical Autocorrelation (ACF) for ARMA Models
#'
#' Displays the theoretical autocorrelation for ARMA Models
#' @param x,object  An \code{"theo_arma"} object from \code{\link{theo_acf_arma}}
#' or  \code{\link{theo_pacf_arma}}.
#' @param ...       Additional parameters
#' @return An \code{array} of dimensions \eqn{N \times S \times S}{N x S x S}.
#' @rdname plot.theo_corr
#' @export
plot.theo_corr = function(x, ...){
  autoplot.theo_corr(object = x, ...)
}

#' @rdname plot.theo_corr
#' @export
#' @importFrom gridExtra grid.arrange
autoplot.theo_corr = function(object, ...){
  grid.arrange(autoplot(object[[1]]),autoplot(object[[2]]), nrow = 1)
}

############################################
# Empirical P/ACF
############################################

#' Empirical ACF Plot
#'
#' @param x       A data set, \code{\link{arima}} or \code{\link{lm}}
#' @param lag.max A \code{integer} indicating the length of the data.
#' @param type    A \code{string} with either \code{"correlation"} or \code{"covariance"}
#' @param robust  A \code{boolean} indicating whether to use a robust estimation.
#' @rdname emp_acf
#' @export
emp_acf = function(x, lag.max = 30L, type = "correlation", robust = FALSE, ...){
  UseMethod("emp_acf")
}

#' @rdname emp_acf
#' @export
emp_acf.Arima = function(x, lag.max = 30L, type = "correlation", robust = FALSE, ...){
  emp_acf.default(x = x$residuals, lag.max = lag.max, type = type, robust = robust)
}

#' @rdname emp_acf
#' @export
emp_acf.lm = function(x, lag.max = 30L, type = "correlation", robust = FALSE, ...){
  emp_acf.default(x = x$residuals, lag.max = lag.max, type = type, robust = robust)
}

#' @rdname emp_acf
#' @export
emp_acf.default = function(x, lag.max = 30L, type = "correlation", robust = FALSE, ...){

  # Force to matrix form
  if(is.ts(x) || is.atomic(x)){
    x = data.matrix(x)
  }

  if(robust){
    o = robcor::robacf(x, lag.max = lag.max, type = type, plot=FALSE)$acf
  }else{
    o = acf(x, lag.max = lag.max, type = type, plot=FALSE)$acf
  }

  cast_acf(o, nrow(x),
           name_ = paste0("Empirical ",
                          if(robust){ "Robust"} else { "Classical"} ,
                          " ACF"),
           type = if(type == "correlation") {"Autocorrelation"} else {"Autocovariance"},
           class = "emp_acf")

}

#' @title Auto-Covariance and Correlation Functions
#' @description The acf function computes the estimated
#' autocovariance or autocorrelation for both univariate and multivariate cases.
#' @param x,object  An \code{"ACF"} object from \code{\link{ACF}}.
#' @param show.ci   A \code{bool} indicating whether to show confidence region
#' @param ci        A \code{double} containing the 1-alpha level. Default is 0.95
#' @param ...       Additional parameters
#' @return An \code{array} of dimensions \eqn{N \times S \times S}{N x S x S}.
#' @rdname plot.emp_acf
#' @export
#' @examples
#' # Calculate the Autocorrelation
#' m = emp_acf(datasets::AirPassengers)
#'
#' # Plot with 95% CI
#' plot(m)
#'
#' # Plot with 90% CI
#' plot(m, ci = 0.90)
#'
#' # Plot without 95% CI
#' plot(m, show.ci = FALSE)
plot.emp_acf = function(x, show.ci = TRUE, ci = 0.95, ...){
  autoplot.emp_acf(object = x, show.ci = show.ci, ci = ci, ...)
}

#' @rdname plot.emp_acf
#' @export
autoplot.emp_acf = function(object, show.ci = TRUE, ci = 0.95, ...){
  autoplot.acf_plot(object = object, show.ci = show.ci, ci = ci, ...)
}

#' Empirical PACF
#'
#' Computes the empirical partial autocorrelation (PACF).
#' @param x       A data set, \code{\link{arima}} or \code{\link{lm}}
#' @param lag.max A \code{integer} indicating the maximum lag to compute.
#' @export
emp_pacf = function(x, lag.max = 30, ...){
  UseMethod("emp_pacf")
}

#' @rdname emp_pacf
#' @export
emp_pacf.Arima = function(x, lag.max = 30, ...){
  emp_pacf.default(x = x$residuals, lag.max = lag.max, ...)
}

#' @rdname emp_pacf
#' @export
emp_pacf.lm = function(x, lag.max = 30, ...){
  emp_pacf.default(x = x$residuals, lag.max = lag.max, ...)
}

#' @rdname emp_pacf
#' @export
emp_pacf.default = function(x, lag.max = 30, ...){

  # Force to matrix form
  if(is.ts(x) || is.atomic(x)){
    x = data.matrix(x)
  }

  cast_acf(pacf(x, lag.max = lag.max, plot=FALSE)$acf, nrow(x),
           name_ = "Empirical PACF",
           type = "Partial Autocorrelation",
           class = "emp_acf")

}

#' Compute Empirical ACF and PACF
#'
#' Computes both the ACF and PACF for a given dataset.
#' @inheritParams  emp_acf
#' @export
#' @examples
#' # Computes the empirical correlation
#' emp_corr(gnp)
emp_corr = function(x, lag.max = 30, robust = FALSE){
  structure(list(acf = emp_acf(x, lag.max = lag.max, robust = robust),
                 pacf = emp_pacf(x, lag.max = lag.max)), class = "emp_corr")

}

#' Plot Empirical Autocorrelation (ACF) for ARMA Models
#'
#' Displays the empirical autocorrelation for ARMA Models
#' @param x,object  An \code{\link{emp_corr}} using \code{\link{emp_acf}}
#' and  \code{\link{emp_pacf}}.
#' @param ...       Additional parameters
#' @return An \code{array} of dimensions \eqn{N \times S \times S}{N x S x S}.
#' @rdname plot.emp_corr
#' @export
plot.emp_corr = function(x, ...){
  autoplot.emp_corr(object = x, ...)
}

#' @rdname plot.emp_corr
#' @export
#' @importFrom gridExtra grid.arrange
autoplot.emp_corr = function(object, ...){
  grid.arrange(autoplot(object[[1]]),autoplot(object[[2]]), nrow = 1)
}
