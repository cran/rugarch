\name{ARFIMAmultiforecast-class}
\Rdversion{1.1}
\docType{class}
\alias{ARFIMAmultiforecast-class}
\alias{as.array,ARFIMAmultiforecast-method}
\alias{as.list,ARFIMAmultiforecast-method}
\alias{show,ARFIMAmultiforecast-method}

\title{class: ARFIMA Multiple Forecast Class}
\description{
Class for the ARFIMA Multiple forecast.}
\section{Slots}{
  \describe{
    \item{\code{forecast}:}{Object of class \code{"vector"} }
    \item{\code{desc}:}{Object of class \code{"vector"} }
  }
}
\section{Extends}{
Class \code{"\linkS4class{ARFIMA}"}, directly.
Class \code{"\linkS4class{rGARCH}"}, by class "ARFIMA", distance 2.
}
\section{Methods}{
  \describe{
    \item{as.array}{\code{signature(x = "ARFIMAmultiforecast")}: 
    extracts the forecast array with matrix column dimensions equal to the number of assets, 
    row dimension the n.ahead and array dimension equal to the number of rolling forecasts 
    chosen.}
    \item{as.list}{\code{signature(x = "ARFIMAmultiforecast")}: 
    extracts the forecast list of length equal to the number of assets, sublists equal to n.roll, 
    row dimension of each sublist equal to n.ahead and column dimension equal to 1 (series forecasts).}
    \item{show}{\code{signature(object = "ARFIMAmultiforecast")}: forecast summary. }
	 }
}
\author{Alexios Ghalanos}
\keyword{classes}