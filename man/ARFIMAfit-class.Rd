\name{ARFIMAfit-class}
\Rdversion{1.1}
\docType{class}
\alias{ARFIMAfit-class}
\alias{as.data.frame,ARFIMAfit-method}
\alias{coef,ARFIMAfit-method}
\alias{fitted,ARFIMAfit-method}
\alias{getspec,ARFIMAfit-method}
\alias{infocriteria,ARFIMAfit-method}
\alias{likelihood,ARFIMAfit-method}
\alias{residuals,ARFIMAfit-method}
\alias{show,ARFIMAfit-method}
\alias{uncmean,ARFIMAfit-method}
\title{class: ARFIMA Fit Class}
\description{
Class for the ARFIMA fit.
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Object of class \code{"vector"} }
    \item{\code{model}:}{Object of class \code{"vector"} }
  }
}
\section{Extends}{
Class \code{"\linkS4class{ARFIMA}"}, directly.
Class \code{"\linkS4class{rGARCH}"}, by class "ARFIMA", distance 2.
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "ARFIMAfit")}: 
    Extracts the position (dates), data, fitted values and residuals.}
    \item{coef}{\code{signature(object = "ARFIMAfit")}: 
    Extracts the coefficients.}
    \item{fitted}{\code{signature(object = "ARFIMAfit")}: 
    Extracts the fitted values.}
    \item{infocriteria}{\code{signature(object = "ARFIMAfit")}: 
    Calculates and returns various information criteria.}
    \item{likelihood}{\code{signature(object = "ARFIMAfit")}: 
    Extracts the likelihood.}
    \item{residuals}{\code{signature(object = "ARFIMAfit")}: 
    Extracts the residuals.}
    \item{show}{\code{signature(object = "ARFIMAfit")}: 
    Fit summary.}
    \item{uncmean}{\code{signature(object = "ARFIMAfit")}: 
    Calculates and returns the unconditional mean. Takes additional arguments 
    \sQuote{method} with option for \dQuote{analytical} or \dQuote{simulation}, 
    \sQuote{n.sim} for the number of simulations (if that method was chosen, 
    and defaults to 100000) and \sQuote{rseed} for the simulation random 
    generator initialization seed.}
	}
}
\author{Alexios Ghalanos}
\examples{
showClass("ARFIMAfit")
}
\keyword{classes}
