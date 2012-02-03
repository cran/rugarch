\name{uGARCHpath-class}
\docType{class}
\alias{uGARCHpath-class}
\alias{as.data.frame,uGARCHpath-method}
\alias{plot,uGARCHpath,missing-method}
\alias{show,uGARCHpath-method}
\title{class: Univariate GARCH Path Simulation Class}
\description{Class for the univariate GARCH Path simulation.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Extends}{
Class \code{"\linkS4class{uGARCHpath}"}, directly.
Class \code{"\linkS4class{rGARCH}"}, by class "GARCHpath", distance 2.
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "uGARCHpath")}: extracts the simulated path values (see note).}
    \item{plot}{\code{signature(x = "uGARCHpath", y = "missing")}: path simulation plots. }
    \item{show}{\code{signature(object = "uGARCHpath")}: path simulation summary. }
	 }
}
\note{
The \code{as.data.frame} function takes optionally 1 additional arguments, namely \code{which},
indicating the type of simulation path series to extract. Valid values are \dQuote{sigma} for the simulated 
sigma, \dQuote{series} for the simulated series and \dQuote{residuals} for the simulated residuals.
The dimension of the \code{data.frame} will be \code{n.sim} by \code{m.sim}.
}
\author{Alexios Ghalanos}
\seealso{
Classes \code{\linkS4class{uGARCHsim}}, \code{\linkS4class{uGARCHfit}} and 
\code{\linkS4class{uGARCHspec}}.
}
\keyword{classes}