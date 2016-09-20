
###############################
# Helper Functions
###############################

#' @importFrom scales percent
resid_hist = function(object){
  ggplot(object, aes(x = residuals)) +
    geom_histogram(aes(y=..density..), binwidth=.5) +
    stat_density(geom="line", position="identity", aes(color = "Kernel")) +
    stat_function(fun = dnorm, aes(color = "Normal")) +
    #stat_function(fun = dt, args = list(df = df), aes(color = "Student's t")) +
    theme_bw() + theme(legend.position =c(0.08, 0.85)) +
    labs(
      y = "Percent",
      x = "Residuals",
      colour = "Density"
    ) +
    scale_colour_manual(values=c("grey", "blue", "orange")) +
    scale_y_continuous(labels = scales::percent)

}


# Hack to use the gmwm gts plot object
resid_plot = function(object, ...) {

  res = gts(object)

  type = if(attr(object, "std")){ "Standardized Residuals" } else { "Residuals"}

  plot(res) + ylab(type) + xlab("Observation Number")
}


#' Residuals of the Process
#'
#' Extracts the residuals of the process
#' @param x   A \code{\link{arima}}, \code{\link{lm}}, or a data object.
#' @param std A \code{boolean} indicating whether the residuals should be
#' standardized.
#' @export
#' @rdname diag_resid
diag_resid = function(x, std = FALSE) {
  UseMethod("diag_resid")
}

#' @export
#' @rdname diag_resid
diag_resid.Arima = function(x, std = FALSE) {
  diag_resid.default(resid(x), std = std)
}

#' @export
#' @rdname diag_resid
diag_resid.lm = function(x, std = FALSE) {
  diag_resid.default(resid(x), std = std)
}

#' @export
#' @rdname diag_resid
diag_resid.default = function(x, std = FALSE) {

  if(is.ts(x)){
    x = as.numeric(x)
  }

  if(std){
    x = x/sd(x)
  }

  structure(data.frame(residuals = x, stringsAsFactors = FALSE),
            std = std,
            class = c("diag_resid","data.frame"))
}

#' Graph the Distribution of Residuals
#'
#' Makes a histogram of the distribution of residuals
#' @param x,object A \code{\link{diag_resid}} object.
#' @param type     A \code{string} indicating either:
#' \code{"hist"} (residual histogram), \code{"resid"} (residual plot),
#' or \code{"both"}
#' @seealso \code{\link[gmwm]{gmwm}} and \code{\link[gmwm]{plot.gmwm}}
#' @return A \code{ggplot2} object.
#' @export
#' @rdname plot.diag_resid
plot.diag_resid = function(x, type = "hist", ...){
  autoplot.diag_resid(object = x, type = type, ...)
}

#' @export
#' @rdname plot.diag_resid
autoplot.diag_resid = function(object, type = "hist", ...){
  type = tolower(type)

  if(type == "hist"){
    resid_hist(object, ...)
  } else if(type == "resid") {
    resid_plot(object, ...)
  } else{
    grid.arrange(
      resid_plot(object, ...),
      resid_hist(object, ...),
      nrow = 1
      )
  }


}


#' Wavelet Variance (WV) fit for White Noise
#'
#' Computes whether the residuals match with a
#' White Noise by using the Generalized Method of Wavelet Moments (GMWM)
#' @param x A \code{\link{arima}}, \code{\link{lm}}, or data set object.
#' @return A \code{diag_wv} object that inherits \code{gmwm} object information.
#' @export
#' @rdname diag_wv
diag_wv = function(x){
  UseMethod("diag_wv")
}

#' @export
#' @rdname diag_wv
diag_wv.lm = function(x){
  diag_wv.default(resid(x))
}

#' @export
#' @rdname diag_wv
diag_wv.Arima = function(x){
  diag_wv.default(resid(x))
}

#' @export
#' @rdname diag_wv
diag_wv.default = function(x){
  structure(gmwm(WN(),x), class = c("diag_wv","gmwm"))
}

#' Graph the Wavelet Variance (WV) Estimation of the White Noise
#'
#' Creates a plot containing the empirical WV and the implied WV for an
#' assumed white noise process.
#' @param x,object A \code{\link{diag_wv}} object.
#' @seealso \code{\link[gmwm]{gmwm}} and \code{\link[gmwm]{plot.gmwm}}
#' @return A \code{ggplot2} object.
#' @export
#' @rdname plot.diag_wv
plot.diag_wv = function(x, ...){
  autoplot(x, ...)
}

#' @export
#' @rdname plot.diag_wv
autoplot.diag_wv = function(object, ...){

  class(object) = "gmwm"

  autoplot(object, title = "Haar WV") + theme(legend.position = "none")
}

#' QQ Normal plot functionality
#'
#' Provides the calculations behind a qq normal plot.
#' @inheritParams diag_resid
#' @return A \code{data.frame} with structure
#' \itemize{
#' \item{residuals}{Residuals}
#' \item{quantiles}{Quantiles}
#' \item{mirror}{Line y-values}
#' }
#' @export
#' @rdname diag_qq
diag_qq = function(x, std = FALSE){
  UseMethod("diag_qq")
}

#' @export
#' @rdname diag_qq
#' @importFrom stats resid
diag_qq.Arima = function(x, std = FALSE){
  diag_qq.default(resid(x), std = std)
}

#' @export
#' @rdname diag_qq
diag_qq.lm = function(x, std = FALSE){
  diag_qq.default(resid(x), std = std)
}

#' @export
#' @rdname diag_qq
diag_qq.ts = function(x, std = FALSE){
  diag_qq.default(data.matrix(x), std = std)
}

#' @export
#' @rdname diag_qq
diag_qq.default = function(x, std = FALSE){

  if(is.ts(x)){
    x = data.matrix(x)
  }

  if(std){
    x = x/sd(x)
  }

  df_x = data.frame(residuals = sort(x))
  n = nrow(df_x)

  df_x$quantiles = qnorm((1:n)/(n+1))

  m = (quantile(x, 0.75)-quantile(x,0.25))/
    (qnorm(0.75)-qnorm(0.25))
  b = quantile(x,0.25) - m*qnorm(0.25)

  structure(df_x, m = m, b = b, std = std,
            class = c("diag_qq","data.frame"))
}

#' QQ Normal Plot
#'
#' Generates a QQ Normal plot
#' @param x,object A \code{\link{diag_qq}} object
#' @param ... Additional parameters (not used)
#' @export
#' @rdname plot.diag_qq
plot.diag_qq = function(x, ...){
  autoplot.diag_qq(object = x, ...)
}

#' @export
#' @rdname plot.diag_qq
autoplot.diag_qq = function(object, ...){

  b = attr(object, "b")
  m = attr(object, "m")
  std = attr(object, "std")

  class(object) = "data.frame"

  ggplot(object, aes(x = quantiles)) +
    geom_point(aes(y = residuals)) +
    geom_abline(intercept = b,
                slope = m,
                color = "blue") +
    labs(
      x = "Theoretical Quantiles",
      y = paste0(if(std){"Standardized "} else { "" }, "Sample Quantiles"),
      title = "Normal Q-Q") +
    theme_bw()
}


#' Bind Fitted and Residual Values
#'
#' Creates a \code{data.frame} containing fitted and residual values from a
#' model.
#' @param fitted A \code{numeric vector} containing fitted values from the model.
#' @param resid  A \code{numeric vector} containing residuals.
#' @return A \code{data.frame} with structure
#' \itemize{
#' \item{fitted}{Fitted Values}
#' \item{residuals}{Residuals}
#' }
#' @details
#' This function will be reworked at a later time to provide model support.
#' @export
#' @rdname diag_fitted
diag_fitted = function(fitted, resid) {
  structure(data.frame(fitted = as.numeric(fitted),
                       residuals = as.numeric(resid), stringsAsFactors = FALSE),
            class = c("diag_fitted","data.frame"))
}

#' Plot Fitted vs. Residual
#'
#' Generates a fitted vs. residual graph
#' @param x,object A \code{\link{diag_fitted}} object
#' @export
#' @rdname plot.diag_fitted
plot.diag_fitted = function(x, ...) {
  autoplot.diag_fitted(object = x, ...)
}


#' @export
#' @rdname plot.diag_fitted
autoplot.diag_fitted = function(object, ...){
  ggplot(object, aes(x = fitted, y = residuals)) +
    geom_point() +
    geom_hline(yintercept = 0, color = "red", linetype="dashed") +
    stat_smooth(method = "loess") +
    labs(
      x = "Fitted values",
      y = "Residuals",
      title = "Residuals vs Fitted Plot"
    ) + theme_bw()
}


#' Time Series Diagnostics
#'
#' Provides a list structure containing time series diagnostic information
#' @param model  An \code{\link{arima}} object.
#' @param xt     The data used to construct said model.
#' @param std    A \code{boolean} indicating whether residuals should be standardized.
#' @param test   A \code{string} indicating the type of portmanteau test to run.
#' @param robust A \code{boolean} indicating whether the ACF should be robust (\code{TRUE}) or
#' not (\code{FALSE}).
#' @return A \code{diag_ts} object
#' @export
diag_ts = function(model, xt, test = "ljung-box", robust = FALSE){

  res = resid(model)
  fitted = as.numeric(xt - res) # fitted() no go for arima...

  structure(list(
    prt = if(tolower(test) == "ljung-box"){
      diag_ljungbox(model, stop_lag = 30)
    } else{
      diag_boxpierce(model, stop_lag = 30)
    },
    res = diag_resid(model),
    fit = diag_fitted(fitted,res),
    acf = emp_acf(model),
    pacf = emp_pacf(model),
    wv = diag_wv(model),
    qq = diag_qq(model)), class = c("diag_ts","list"))

}

#' Graph Time Series Diagnostic
#'
#' Creates a time series diagnostic plot
#' @param x,object A \code{\link{diag_ts}} object.
#' @return A \code{grid.arrange} object.
#' @export
#' @rdname plot.diag_ts
plot.diag_ts = function(x, ...){
  autoplot.diag_ts(object = x, ...)
}

#' @export
#' @rdname plot.diag_ts
autoplot.diag_ts = function(object, ...){

  with(object,
       grid.arrange(autoplot(res, type = "resid"), autoplot(fit), autoplot(res, type = "hist"), autoplot(qq),
                    autoplot(acf), autoplot(pacf),autoplot(wv), autoplot(prt), nrow = 2)
  )
}
