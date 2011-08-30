\name{ugarchfit-methods}
\docType{methods}
\alias{ugarchfit}
\alias{ugarchfit-methods}
\alias{ugarchfit,ANY-method}
\alias{ugarchfit,uGARCHspec-method}
\title{function: Univariate GARCH Fitting}
\description{
Method for fitting a variety of univariate GARCH models.
}
\usage{
ugarchfit(spec, data, out.sample = 0, solver = "solnp", solver.control = list(), 
fit.control = list(stationarity = 1, fixed.se = 0, scale = 0), ...)
}
\arguments{
\item{data}{
A univariate data object. Can be a numeric vector, matrix, data.frame, zoo, xts, 
timeSeries, ts or irts object.}
\item{spec}{
A univariate GARCH spec object of class \code{\linkS4class{uGARCHspec}}.}
\item{out.sample}{
A positive integer indicating the number of periods before the last to keep for 
out of sample forecasting (see details).}
\item{solver}{
One of either \dQuote{nlminb}, \dQuote{solnp} or \dQuote{gosolnp}. }
\item{solver.control}{
Control arguments list passed to optimizer.}
\item{fit.control}{
Control arguments passed to the fitting routine. Stationarity explicitly imposes 
the variance stationarity constraint during optimization. The fixed.se argument 
controls whether standard errors should be calculated for those parameters which 
were fixed (through the fixed.pars argument of the \code{\link{ugarchspec} 
function). The scale parameter controls whether the data should be scaled before 
being submitted to the optimizer.}}
\item{...}{ . }
}
\value{
A \code{\linkS4class{uGARCHfit}} object containing details of the GARCH fit.
}
\details{
The GARCH optimization routine first calculates a set of feasible starting 
points which are used to initiate the GARCH recursion. The main part of the 
likelihood calculation is performed in C-code for speed.\cr
The out.sample option is provided in order to carry out forecast performance 
testing against actual data. A minimum of 5 data points are required for these 
tests. If the out.sample option is positive, then the routine will fit only 
N - out.sample (where N is the total data length) data points, leaving 
out.sample points for forecasting and testing using the forecast performance 
measures. In the \code{\link{ugarchforecast}} routine the n.ahead may also be 
greater than the out.sample number resulting in a combination of out of sample 
data points matched against actual data and some without, which the forecast 
performance tests will ignore.\cr
The \dQuote{gosolnp} solver allows for the initialization of multiple restarts 
of the solnp solver with randomly generated parameters (see documentation in the 
Rsolnp-package for details of the strategy used). The solver.control list then 
accepts the following additional (to the solnp) arguments: \dQuote{n.restarts} 
is the number of solver restarts required (defaults to 1), \dQuote{parallel} and 
\dQuote{parallel.control} for use of the parallel functionality, \dQuote{rseed} 
is the seed to initialize the random number generator, and \dQuote{n.sim} is the 
number of simulated parameter vectors to generate per n.restarts.
}
\author{Alexios Ghalanos}
\seealso{
For specification \code{\link{ugarchspec}},filtering \code{\link{ugarchfilter}}, 
forecasting \code{\link{ugarchforecast}}, simulation \code{\link{ugarchsim}}, 
rolling forecast and estimation \code{\link{ugarchroll}}, parameter distribution 
and uncertainty \code{\link{ugarchdistribution}}, bootstrap forecast 
\code{\link{ugarchboot}}.
}
\examples{
# Basic GARCH(1,1) Spec
data(dmbp)
spec = ugarchspec()
fit = ugarchfit(data = dmbp[,1], spec = spec)
fit
coef(fit)
head(as.data.frame(fit))
#plot(fit,which="all")
# in order to use fpm (forecast performance measure function)
# you need to select a subsample of the data:
spec = ugarchspec()
fit = ugarchfit(data = dmbp[,1], spec = spec, out.sample=100)
forc = ugarchforecast(fit, n.ahead=100)
# this means that 100 data points are left from the end with which to
# make inference on the forecasts
fpm(forc)
}
\keyword{methods}