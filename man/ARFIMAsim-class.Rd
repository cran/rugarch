\name{ARFIMAsim-class}
\Rdversion{1.1}
\docType{class}
\alias{ARFIMAsim-class}
\alias{as.data.frame,ARFIMAsim-method}
\alias{show,ARFIMAsim-method}
\title{class: ARFIMA Simulation Class}
\description{
Class for the ARFIMA simulation.
}
\section{Slots}{
  \describe{
    \item{\code{simulation}:}{Object of class \code{"vector"} }
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
    \item{as.data.frame}{\code{signature(x = "ARFIMAsim")}: extracts the simulated values (see note).}
    \item{show}{\code{signature(object = "ARFIMAsim")}: simulation summary. }
	 }
}
\note{
The \code{as.data.frame} function takes optionally 1 additional arguments, namely \code{which},
indicating the type of simulation series to extract. Valid values are \dQuote{series} for the 
simulated series and \dQuote{residuals} for the simulated residuals. The dimension of the 
\code{data.frame} will be \code{n.sim} by \code{m.sim}.
}
\author{Alexios Ghalanos}
\keyword{classes}
