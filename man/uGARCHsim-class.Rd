\name{uGARCHsim-class}
\docType{class}
\alias{uGARCHsim-class}
\alias{as.data.frame,uGARCHsim-method}
\alias{plot,uGARCHsim,missing-method}
\alias{show,uGARCHsim-method}
\title{class: Univariate GARCH Simulation Class}
\description{Class for the univariate GARCH simulation.}
\section{Extends}{
Class \code{"\linkS4class{GARCHsim}"}, directly.
Class \code{"\linkS4class{rGARCH}"}, by class "GARCHsim", distance 2.
}
\section{Slots}{
  \describe{
    \item{\code{simulation}:}{Object of class \code{"vector"}
    Holds data on the simulation.}
    \item{\code{model}:}{Object of class \code{"vector"}
    The model specification common to all objects.}
    \item{\code{seed}:}{Object of class \code{"integer"}
    The random seed used.}
  }
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "uGARCHsim")}: 
    Extracts the simulated values (see note).}
    \item{plot}{\code{signature(x = "uGARCHsim", y = "missing")}: 
    Simulation plots. }
    \item{show}{\code{signature(object = "uGARCHsim")}: 
    Simulation summary. }
	 }
}
\note{
The \code{as.data.frame} function takes optionally 1 additional arguments, 
namely \code{which}, indicating the type of simulation series to extract. 
Valid values are \dQuote{sigma} for the simulated sigma, \dQuote{series} for the 
simulated series and \dQuote{residuals} for the simulated residuals. The 
dimension of the \code{data.frame} will be \code{n.sim} by \code{m.sim}.
}
\author{Alexios Ghalanos}
\seealso{
Classes \code{\linkS4class{uGARCHforecast}}, \code{\linkS4class{uGARCHfit}} and 
\code{\linkS4class{uGARCHspec}}.
}
\examples{
\dontrun{
# Basic GARCH(1,1) Spec
data(dmbp)
spec = ugarchspec()
fit = ugarchfit(data = dmbp[,1], spec = spec)
sim = ugarchsim(fit,n.sim=1000, n.start=1, m.sim=1, startMethod="sample")
sim
# plot(sim, which="all")
# as.data.frame takes an extra argument which
# indicating one of "sigma", "series" and "residuals"
head(as.data.frame(sim, which = "sigma"))
}}
\keyword{classes}
