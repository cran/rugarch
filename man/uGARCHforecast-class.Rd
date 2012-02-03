\name{uGARCHforecast-class}
\docType{class}
\alias{uGARCHforecast-class}
\alias{as.data.frame,uGARCHforecast-method}
\alias{as.array,uGARCHforecast-method}
\alias{as.list,uGARCHforecast-method}
\alias{plot,uGARCHforecast,missing-method}
\alias{show,uGARCHforecast-method}
\alias{fpm}
\alias{fpm,ANY-method}
\alias{fpm,uGARCHforecast-method}
\title{class: Univariate GARCH Forecast Class}
\description{Class for the univariate GARCH forecast.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Extends}{
Class \code{\linkS4class{GARCHforecast}}, directly.
Class \code{\linkS4class{rGARCH}}, by class \code{\linkS4class{GARCHforecast}}, distance 2.
}
\section{Methods}{
  \describe{
   \item{as.array}{\code{signature(x = "uGARCHforecast")}: 
   Extracts the forecast array with matrix column dimensions equal to the 
   n.ahead value and row dimension 2 (sigma and series forecast), and array 
   dimension equal to the number of rolling forecasts chosen.}
    \item{as.data.frame}{\code{signature(x = "uGARCHforecast")}: 
    Extracts the forecasts. Takes many additional arguments (see note below).}
     \item{as.list}{\code{signature(x = "uGARCHforecast")}: 
     Extracts the forecast list with all rollframes.}
    \item{plot}{\code{signature(x = "uGARCHforecast", y = "missing")}:
    Forecast plots with \code{n.roll} optional argument indicating the rolling 
    sequence to plot. }
    \item{fpm}{\code{signature(object = "uGARCHforecast")}: 
   Forecast performance measures.}
    \item{show}{\code{signature(object = "uGARCHforecast")}: 
    Forecast summary returning the 0-roll frame only.}
	 }
}
\note{
There are 3 main extractor functions for the uGARCHforecast object which is 
admittedly the most complex in the package as a result of allowing for rolling 
forecasts. The \code{as.array} extracts an array object where each page of the 
array represents a roll. The \code{as.list} method works similarly returns 
instead a list object. There are no additional arguments to these extractor 
functions and they will return all the forecasts. The \code{as.data.frame} 
method on the other hand provides for 5 additional arguments. The argument 
\code{which} indicates the type of forecast value to return(with valid valued 
being \dQuote{sigma} and \dQuote{series}). The \code{rollframe} option is for 
the rolling frame to return (with 0 being the default no-roll) and allows either 
a valid numeric value or alternatively the character value \dQuote{all} for which
additional options then come into play. When \dQuote{all} is chosen in the 
\code{rollframe} argument, the data.frame returned may be time aligned (logical 
option \code{aligned}) in which case the logical option \code{prepad} indicates 
whether to pad the values prior to the forecast start time with actual values or 
NA (value \code{FALSE}). Finally, the \code{type} option controls whether to 
return all forecasts (value 0, default), return only those forecasts which have 
in sample equivalent data (value 1) or return only those values which are truly 
forecasts without in sample data (value 2). Depending on the intended usage of 
the forecasts, some or all these options may be useful to the user when 
extracting data from the forecast object.\cr
The plot method takes additional arguments \code{which} and \code{n.roll} 
indicating which roll frame to plot.
}
\seealso{
Classes \code{\linkS4class{uGARCHfit}}, \code{\linkS4class{uGARCHsim}} and 
\code{\linkS4class{uGARCHspec}}.
}
\author{Alexios Ghalanos}
\examples{
\dontrun{
# Basic GARCH(1,1) Spec
data(dmbp)
spec = ugarchspec()
fit = ugarchfit(data = dmbp[,1], spec = spec)
forc = ugarchforecast(fit, n.ahead=20)
forc
head(as.data.frame(forc))
#plot(forc, which = "all")
}}
\keyword{classes}