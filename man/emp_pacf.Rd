% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/acf.R
\name{emp_pacf}
\alias{emp_pacf}
\alias{emp_pacf.Arima}
\alias{emp_pacf.default}
\alias{emp_pacf.lm}
\title{Empirical PACF}
\usage{
emp_pacf(x, lag.max = 30, ...)

\method{emp_pacf}{Arima}(x, lag.max = 30, ...)

\method{emp_pacf}{lm}(x, lag.max = 30, ...)

\method{emp_pacf}{default}(x, lag.max = 30, ...)
}
\arguments{
\item{x}{A data set, \code{\link{arima}} or \code{\link{lm}}}

\item{lag.max}{A \code{integer} indicating the maximum lag to compute.}
}
\description{
Computes the empirical partial autocorrelation (PACF).
}

