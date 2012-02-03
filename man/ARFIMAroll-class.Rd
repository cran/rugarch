\name{ARFIMAroll-class}
\Rdversion{1.1}
\docType{class}
\alias{ARFIMAroll-class}
\alias{as.ARFIMAforecast,ARFIMAroll-method}
\alias{as.ARFIMAforecast,ANY-method}
\alias{as.ARFIMAforecast}
\alias{as.data.frame,ARFIMAroll-method}
\alias{report,ARFIMAroll-method}
\alias{fpm,ARFIMAroll-method}
\title{class: ARFIMA Rolling Forecast Class}
\description{
Class for the ARFIMA rolling forecast.}
\section{Slots}{
  \describe{
    \item{\code{roll}:}{Object of class \code{"vector"} }
    \item{\code{forecast}:}{Object of class \code{"vector"} }
    \item{\code{model}:}{Object of class \code{"vector"} }
  }
}
\section{Extends}{
Class \code{"\linkS4class{ARFIMA}"}, directly.
Class \code{"\linkS4class{rGARCH}"}, by class "ARFIMA", distance 2.
}
\section{Methods}{
  \describe{
    \item{as.ARFIMAforecast}{\code{signature(object = "ARFIMAroll")}: 
 	 extracts  and converts the forecast object contained in the roll object to
    one of \code{\linkS4class{ARFIMAforecast}} given the refit number supplied 
    by additional argument \sQuote{refit} (defaults to 1).}
    \item{as.data.frame}{\code{signature(x = "ARFIMAroll")}: extracts various 
    values from object (see note). }
   \item{fpm}{\code{signature(object = "ARFIMAroll")}: 
   Forecast performance measures.}
    \item{report}{\code{signature(object = "ARFIMAroll")}: roll backtest reports 
    (see note).}
	 }
}
\note{
The \code{as.data.frame} extractor method allows the extraction of a variety of 
values from the object. Additional arguments are:\cr
\emph{which} indicates the type of value to return. Valid values are \dQuote{coefs} 
returning the parameter coefficients for all refits, \dQuote{density} for the 
parametric density, \dQuote{coefmat} for the parameter coefficients with their 
respective standard errors and t- and p- values, \dQuote{LLH} for the likelihood
across the refits, and \dQuote{VaR} for the Value At Risk measure if it was 
requested in the roll function call.\cr
\emph{n.ahead} for the n.ahead forecast horizon to return if \code{which} was 
used with arguments \dQuote{density} or \sQuote{VaR}.\cr
\emph{refit} indicates which refit window to return the \dQuote{coefmat} 
if that is chosen. If \dQuote{series} is chosen under via the \emph{which} 
argument, then the forecast series is returned for a particular refit, else 
when \dQuote{all} is used it returns the complete forecasted series across all 
refits.\cr
The \code{report} method takes the following additional arguments:\cr
\emph{type} for the report type. Valid values are \dQuote{VaR} for the Value at 
Risk report based on the unconditional and conditional coverage tests for VaR 
exceedances (discussed below) and \dQuote{fpm} for forecast performance measures.\cr
\emph{n.ahead} for the rolling n.ahead forecasts (defaults to 1).\cr
\emph{VaR.alpha} for the Value at Risk backtest report, this is the tail 
probability and defaults to 0.01.\cr
\emph{conf.level} the confidence level upon which the conditional coverage 
hypothesis test will be based on (defaults to 0.95).\cr
Kupiec's unconditional coverage test looks at whether the amount of expected 
versus actual exceedances given the tail probability of VaR actually occur as 
predicted, while the conditional coverage test of Christoffersen is a joint test 
of the unconditional coverage and the independence of the exceedances. Both the 
joint and the separate unconditional test are reported since it is always 
possible that the joint test passes while failing either the independence or 
unconditional coverage test.\cr
The \dQuote{fpm} does not take any additional arguments, but instead returns the
forecast performance measures for all \dQuote{n.ahead} values.
}
\author{Alexios Ghalanos}
\keyword{classes}