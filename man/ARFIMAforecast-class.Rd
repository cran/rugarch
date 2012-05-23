\name{ARFIMAforecast-class}
\Rdversion{1.1}
\docType{class}
\alias{ARFIMAforecast-class}
\alias{as.array,ARFIMAforecast-method}
\alias{as.data.frame,ARFIMAforecast-method}
\alias{as.list,ARFIMAforecast-method}
\alias{show,ARFIMAforecast-method}
\alias{fpm,ARFIMAforecast-method}
\title{class: ARFIMA Forecast Class}
\description{
Class for the ARFIMA forecast.
}
\section{Slots}{
  \describe{
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
   \item{as.array}{\code{signature(x = "ARFIMAforecast")}: 
   Extracts the forecast array with matrix column dimensions equal to the 
   n.ahead value and row dimension 1 (series forecast), and array dimension 
   equal to the number of rolling forecasts chosen.}
   \item{as.data.frame}{\code{signature(x = "ARFIMAforecast")}: 
   Extracts the forecasts. Takes many additional arguments (see note below).}
   \item{as.list}{\code{signature(x = "ARFIMAforecast")}: 
   Extracts the forecast list with all rollframes.}
   \item{fpm}{\code{signature(object = "ARFIMAforecast")}: 
   Forecast performance measures.}
   \item{show}{\code{signature(object = "ARFIMAforecast")}: 
   Forecast summary returning the 0-roll frame only.}
	}
}
\note{
There are 3 main extractor functions for the ARFIMA object which is admittedly 
the most complex in the package as a result of allowing for rolling forecasts. 
The \code{as.array} extracts an array object where each page of the array 
represents a roll. The \code{as.list} method works similarly returns instead a 
list object. There are no additional arguments to these extractor functions and 
they will return all the forecasts. The \code{as.data.frame} method on the other 
hand provides for 4 additional arguments. The \code{rollframe} option is for the 
rolling frame to return (with 0 being the default no-roll) and allows either a 
valid numeric value or alternatively the character value \dQuote{all} for which 
additional options then come into play. When \dQuote{all} is chosen in the 
\code{rollframe} argument, the data.frame returned may be time aligned 
(logical option \code{aligned}) in which case the logical option \code{prepad} 
indicates whether to pad the values prior to the forecast start time with actual 
values or NA (value \code{FALSE}). Finally, the \code{type} option controls 
whether to return all forecasts (value 0, default), return only those forecasts 
which have in sample equivalent data (value 1) or return only those values which 
are truly forecasts without in sample data (value 2). Depending on the intended 
usage of the forecasts, some or all these options may be useful to the user when 
extracting data from the forecast object.\cr
}
\author{Alexios Ghalanos}

\keyword{classes}
