\name{uGARCHmultiforecast-class}
\docType{class}
\alias{uGARCHmultiforecast-class}
\alias{as.list,uGARCHmultiforecast-method}
\alias{as.array,uGARCHmultiforecast-method}
\alias{show,uGARCHmultiforecast-method}
\title{class: Univariate GARCH Multiple Forecast Class}
\description{Class for the univariate GARCH Multiple forecast.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Extends}{
Class \code{\linkS4class{GARCHforecast}}, directly.
Class \code{\linkS4class{rGARCH}}, by class \code{\linkS4class{GARCHforecast}}, distance 3.
}
\section{Methods}{
\describe{
  \item{as.array}{\code{signature(x = "uGARCHmultiforecast")}: 
   		extracts the forecast array with matrix column dimensions equal to
   		the number of assets, row dimension the n.ahead and
   		array dimension equal to the number of rolling forecasts chosen. The optional argument
   		\dQuote{which} allows to choose from \dQuote{sigma} and \dQuote{series} to return the
   		forecasts for.}
     \item{as.list}{\code{signature(x = "uGARCHforecast")}: extracts the forecast list
     of length equal to the number of assets, sublists equal to n.roll, row dimension
     of each sublist equal to n.ahead and column dimension equal to 2 (sigma and series forecasts).}
    \item{show}{\code{signature(object = "uGARCHforecast")}: forecast summary.}
  }
}
\author{Alexios Ghalanos}
\seealso{
Classes \code{\linkS4class{uGARCHmultifilter}}, \code{\linkS4class{uGARCHmultifit}} and 
\code{\linkS4class{uGARCHmultispec}}.
}
\keyword{classes}