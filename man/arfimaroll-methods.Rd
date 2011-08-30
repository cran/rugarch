\name{arfimaroll-methods}
\docType{methods}
\alias{arfimaroll-methods}
\alias{arfimaroll,ANY-method}
\alias{arfimaroll,ARFIMAspec-method}
\alias{arfimaroll}
\title{function: ARFIMA Rolling Density Forecast and Backtesting}
\description{
Method for creating rolling density forecast from ARFIMA models with option
for refitting every n periods and some multicore parallel functionality.}
\usage{
arfimaroll(spec, data, n.ahead = 1, forecast.length = 500, refit.every = 25, 
refit.window = c("recursive", "moving"), parallel = FALSE, 
parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2), 
solver = "solnp", fit.control = list(), solver.control = list(), 
calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), ...)
}
\arguments{
  \item{spec}{An ARFIMA spec object specifiying the desired model for testing.}
  \item{data}{ A univariate dataset.}
  \item{n.ahead}{ The number of periods to forecast.}
  \item{forecast.length}{ The length of the total forecast for which out of 
  sample data from the dataset will be excluded for testing.}
  \item{refit.every}{ Determines every how many periods the model is 
  re-estimated.}
  \item{refit.window}{ Whether the refit is done on an expanding window 
  including all the previous data or a moving window, the length of the window 
  determined by the argument above (\emph{refit.every}).}
  \item{parallel}{ Whether to make use of parallel processing on multicore 
  systems.}
  \item{parallel.control}{ The parallel control options including the type of 
  package for performing the parallel calculations (\sQuote{multicore} for 
  non-windows O/S and \sQuote{snowfall} for all O/S), and the number of cores to 
  make use of.}
  \item{solver}{ The solver to use. }
  \item{fit.control}{ Control parameters parameters passed to the fitting 
  function.}
  \item{solver.control}{ Control parameters passed to the solver.}
  \item{calculate.VaR}{ Whether to calculate forecast Value at Risk during the 
  estimation.}
  \item{VaR.alpha}{ The Value at Risk tail level to calculate. }
  \item{...}{ . }
}
\details{
ARFIMA models generate a partially time varying density based on the variation 
in the conditional mean values (sigma, skewness and shape are not time varying). 
The function first generates rolling forecasts of the ARFIMA model and then 
rescales the density from a standardized (0, 1, skew, shape) to the one 
representing the underlying return process (mu, sigma, skew, shape). Given this 
information it is then a simple matter to generate any measure of risk through 
the analytical evaluation of some type of function of the density. The function 
calculates one such measure (VaR), but since the full density parameters are 
returned, the user can calculate many others.\cr
The argument \emph{refit.every} determines every how many periods the fit is 
recalculated and the total forecast length actually calculated. For example, for 
a forecast length of 500 and \emph{refit.every} of 25, this is 20 windows of 25 
periods each for a total actual forecast length of 500. However, for a 
\emph{refit.every} of 30, we take the floor of the division of 500 by 30 which 
is 16 windows of 30 periods each for a total actual forecast length of 480 
(16 x 30). The important thing to remember about the refit.every is that it acts 
like the n.roll argument in the  \code{\link{arfimaforecast}} function as it 
determines the number of rolls to perform. For example for n.ahead of 1 and 
refit.every of 25, the forecast is rolled every day using the filtered (actual) 
data of the previous period while for n.ahead of 1 and refit.every of 1 we will 
get 1 n.ahead forecasts for every day after which the model is refitted and 
reforecast for a total of 500 refits (when length.forecast is 500)!\cr 
}
\value{
An object of class \code{\linkS4class{ARFIMAroll}}.
}
\author{Alexios Ghalanos}
\keyword{methods}
