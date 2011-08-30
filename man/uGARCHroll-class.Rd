\name{uGARCHroll-class}
\docType{class}
\alias{uGARCHroll-class}
\alias{as.data.frame,uGARCHroll-method}
\alias{plot,uGARCHroll,missing-method}
\alias{report,uGARCHroll-method}
\alias{report,ANY-method}
\alias{report}
\alias{as.uGARCHforecast,uGARCHroll-method}
\alias{as.uGARCHforecast,ANY-method}
\alias{as.uGARCHforecast}
\alias{fpm,uGARCHroll-method}
\title{class: Univariate GARCH Rolling Forecast Class}
\description{Class for the univariate GARCH rolling forecast.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Extends}{
Class \code{"\linkS4class{GARCHroll}"}, directly.
Class \code{"\linkS4class{rGARCH}"}, by class "GARCHroll", distance 2.
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "uGARCHroll")}: 
    Extracts various values from object (see note).}
    \item{plot}{\code{signature(x = "uGARCHroll", y = "missing")}:
    Roll result backtest plots (see note).}
    \item{report}{\code{signature(object = "uGARCHroll")}:
	Roll backtest reports (see note).}
	\item{fpm}{\code{signature(object = "uGARCHroll")}: 
   Forecast performance measures.}
    \item{as.uGARCHforecast}{\code{signature(object = "uGARCHroll")}:
	Extracts  and converts the forecast object contained in the roll object to
    one of \code{\linkS4class{uGARCHforecast}} given the refit number supplied by
    additional argument \sQuote{refit} (defaults to 1).}
	 }
}
\note{
The \code{as.data.frame} extractor method allows the extraction of a variety of 
values from the object. Additional arguments are:
\code{which} indicates the type of value to return. Valid values are \dQuote{coefs} 
returning the parameter coefficients for all refits, \dQuote{density} for the 
parametric density, \dQuote{coefmat} for the parameter coefficients with their 
respective standard errors and t- and p- values, \dQuote{LLH} for the likelihood
across the refits, and \dQuote{VaR} for the Value At Risk measure if it was 
requested in the roll function call.\cr
\code{n.ahead} for the n.ahead forecast horizon to return if \code{which} was 
used with arguments \dQuote{density} or \dQuote{VaR}.\cr
\code{refit} indicates which refit window to return the \dQuote{coefmat} if that 
was chosen.\cr
The plot method takes the following additional arguments:\cr
\emph{which} allows for either a numeric value of 1:4, else will default to 
\dQuote{ask} for interactive printing of the options in the command windows.
Additionally, the value of \dQuote{all} wil create a 2x2 chart with all plots.\cr
\emph{n.ahead} for the rolling n.ahead forecasts (defaults to 1).\cr
\emph{VaR.alpha} for the Value at Risk backtest plot, this is the tail probability
and defaults to 0.01.\cr
\emph{density.support} the support for the time varying density plot density, 
defaults to c(-0.15, 0.15) but you should change this to something more 
appropriate for your data and period under consideration.\cr
The report method takes the following additional arguments:\cr
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
unconditional coverage test.
}
\author{Alexios Ghalanos}
\keyword{classes}