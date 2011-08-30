\name{ugarchroll-methods}
\docType{methods}
\alias{ugarchroll}
\alias{ugarchroll-methods}
\alias{ugarchroll,ANY-method}
\alias{ugarchroll,uGARCHspec-method}
\title{function: Univariate GARCH Rolling Density Forecast and Backtesting}
\description{
Method for creating rolling density forecast from ARMA-GARCH models with option
for refitting every n periods and some multicore parallel functionality.
}
\usage{
ugarchroll(spec, data, n.ahead = 1, forecast.length = 500, refit.every = 25, 
refit.window = c("recursive", "moving"), parallel = FALSE, 
parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2), solver = "solnp", 
fit.control = list(), solver.control = list(), calculate.VaR = TRUE, 
VaR.alpha = c(0.01, 0.05), ...)
}
\arguments{
\item{spec}{
A univariate GARCH spec object specifiying the desired model for testing.}
\item{data}{
A univariate dataset.}
\item{n.ahead}{
The number of periods to forecast.}
\item{forecast.length}{
The length of the total forecast for which out of sample data 
from the dataset will be excluded for testing.}
\item{refit.every}{
Determines every how many periods the model is re-estimated.}
\item{refit.window}{
Whether the refit is done on an expanding window including all the previous data 
or a moving window, the length of the window determined by the argument above 
(\dQuote{refit.every}).}
\item{parallel}{
Whether to make use of parallel processing on multicore systems.}
\item{parallel.control}{
The parallel control options including the type of package for performing the 
parallel calculations (\sQuote{multicore} for non-windows O/S and 
\sQuote{snowfall} for all O/S), and the number of cores to make use of.}
\item{solver}{
The solver to use. }
\item{fit.control}{
Control parameters parameters passed to the fitting function.}
\item{solver.control}{
Control parameters passed to the solver.}
\item{calculate.VaR}{
Whether to calculate forecast Value at Risk during the estimation.}
\item{VaR.alpha}{
The Value at Risk tail level to calculate.}
\item{...}{ . }
}
\details{
GARCH models generate a partially time varying density based on the variation in 
the conditional sigma and mean values (skewness and shape are usually not time 
varying in GARCH models unless the underlying distribution has an interaction 
with the conditional sigma). The function first generates rolling forecasts of 
the ARMA-GARCH model and then rescales the density from a standardized 
(0, 1, skew, shape) to the one representing the underlying return process 
(mu, sigma, skew, shape). Given this information it is then a simple matter to 
generate any measure of risk through the analytical evaluation of some type of 
function of the time varying density. The function calculates one such measure 
(VaR), but since the full time varying density parameters are returned, the user
can calculate many others (see for example partial moments based measures or the
Pedersen-Satchell family of measures).\cr
The argument \code{refit.every} determines every how many periods the fit is
recalculated and the total forecast length actually calculated. For example, for
a forecast length of 500 and \code{refit.every} of 25, this is 20 windows of 25 
periods each for a total actual forecast length of 500. However, for a 
\code{refit.every} of 30, we take the floor of the division of 500 by 30 which 
is 16 windows of 30 periods each for a total actual forecast length of 480 
(16 x 30). The important thing to remember about the refit.every is that it acts 
like the n.roll argument in the \code{\link{ugarchforecast}} function as it 
determines the number of rolls to perform. For example for n.ahead of 1 and 
refit.every of 25, the forecast is rolled every day using the filtered (actual) 
data of the previous period while for n.ahead of 1 and refit.every of 1 we will 
get 1 n.ahead forecasts for every day after which the model is refitted and 
reforecast for a total of 500 refits (when length.forecast is 500)!\cr
The function has 2 main methods for viewing the data, a standard plot method and
a new report methods (see class \code{\linkS4class{uGARCHroll}} for details on 
how to use these methods).
}
\value{
An object of class \code{\linkS4class{uGARCHroll}}.
}
\author{Alexios Ghalanos}
\seealso{
For specification \code{\link{ugarchspec}}, fitting \code{\link{ugarchfit}}, 
filtering \code{\link{ugarchfilter}}, forecasting \code{\link{ugarchforecast}}, 
simulation \code{\link{ugarchsim}}, parameter distribution and uncertainty 
\code{\link{ugarchdistribution}}, bootstrap forecast \code{\link{ugarchboot}}.
}
\examples{
\dontrun{
data(sp500ret)
ctrl = list(rho = 1, delta = 1e-9, outer.iter = 100, tol = 1e-7)
spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
		mean.model = list(armaOrder = c(0,0), include.mean = TRUE), 
		distribution.model = "std")

sp500.bktest = ugarchroll(spec, data = sp500ret, n.ahead = 1, 
forecast.length = 100, refit.every = 25, refit.window = "recursive", 
solver = "solnp", fit.control = list(), solver.control = ctrl, 
calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.025, 0.05))
report(sp500.bktest, type="VaR", n.ahead = 1, VaR.alpha = 0.01, 
conf.level = 0.95) 
report(sp500.bktest, type="fpm")
}}
\keyword{methods}