\name{ARFIMApath-class}
\Rdversion{1.1}
\docType{class}
\alias{ARFIMApath-class}
\alias{as.data.frame,ARFIMApath-method}
\alias{show,ARFIMApath-method}

\title{class: ARFIMA Path Simulation Class}
\description{
Class for the ARFIMA Path simulation.}
\section{Slots}{
  \describe{
    \item{\code{path}:}{Object of class \code{"vector"} }
    \item{\code{model}:}{Object of class \code{"vector"} }
    \item{\code{seed}:}{Object of class \code{"integer"} }
  }
}
\section{Extends}{
Class \code{"\linkS4class{ARFIMA}"}, directly.
Class \code{"\linkS4class{rGARCH}"}, by class "ARFIMA", distance 2.
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "ARFIMApath")}: 
    Extracts the simulated path values (see note).}
    \item{show}{\code{signature(object = "ARFIMApath")}: path simulation summary. }
	 }
}
\note{
The \code{as.data.frame} function takes optionally 1 additional arguments, namely \code{which},
indicating the type of simulation path series to extract. Valid values \dQuote{series} for the 
simulated series and \dQuote{residuals} for the simulated residuals. The dimension of the 
\code{data.frame} will be \code{n.sim} by \code{m.sim}.
}
\author{Alexios Ghalanos}
\keyword{classes}
