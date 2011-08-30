\name{arfimadistribution-methods}
\docType{methods}
\alias{arfimadistribution-methods}
\alias{arfimadistribution,ANY-method}
\alias{arfimadistribution,ARFIMAfit-method}
\alias{arfimadistribution,ARFIMAspec-method}
\alias{arfimadistribution}
\title{function: ARFIMA Parameter Distribution via Simulation}
\description{
Method for simulating and estimating the parameter distribution from an ARFIMA 
models as well as the simulation based consistency of the estimators given the 
data size.}
\usage{
arfimadistribution(fitORspec, n.sim = 2000, n.start = 1, m.sim = 100, 
recursive = FALSE, recursive.length = 6000, recursive.window = 1000, 
prereturns = NA, preresiduals = NA, rseed = NA, 
custom.dist = list(name = NA,  distfit = NA, type = "z"), mexsimdata = NULL, 
fit.control = list(), solver = "solnp", solver.control = list(), 
parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), 
cores = 2), ...)
}
\arguments{
\item{fitORspec}{ Either an ARFIMA fit object of class 
\code{\linkS4class{ARFIMAfit}} or alternatively an ARFIMA specification object 
of class \code{\linkS4class{ARFIMAspec}} with valid parameters supplied via the 
\code{fixed.pars} argument in the specification.}
\item{n.sim}{ The simulation horizon.}
\item{n.start}{ The burn-in sample.}
\item{m.sim}{ The number of simulations.}
\item{recursive}{ Whether to perform a recursive simulation on an expanding 
window.}
\item{recursive.length}{ If \code{recursive} is TRUE, this indicates the final 
length of the simulation horizon, with starting length \code{n.sim}.}
\item{recursive.window}{ If \code{recursive} is TRUE, this indicates the 
increment to the expanding window. Together with \code{recursive.length}, it 
determines the total number of separate and increasing length windows which will 
be simulated and fitted.}
\item{prereturns}{ Allows the starting return data to be provided by the user.}
\item{preresiduals}{ Allows the starting residuals to be provided by the user.}
\item{rseed}{ Optional seeding value(s) for the random number generator.}
\item{custom.dist}{ Optional density with fitted object from which to simulate.}
\item{mexsimdata}{ Matrix of simulated external regressor-in-mean data. If the 
fit object contains external regressors in the mean equation, this must be 
provided.}
\item{solver}{One of either \dQuote{nlminb} or \dQuote{solnp}. }
\item{solver.control}{Control arguments list passed to optimizer.}
\item{fit.control}{Control arguments passed to the fitting routine (as in the 
\code{arfimafit} method).}
\item{parallel}{
Whether to make use of parallel processing on multicore systems.
}
\item{parallel.control}{
The parallel control options including the type of package for performing the 
parallel calculations (\sQuote{multicore} for non-windows O/S and 
\sQuote{snowfall} for all O/S), and the number of cores to make use of.
}
\item{...}{ . }
}
\value{
A \code{\linkS4class{ARFIMAdistribution}} object containing details of the 
ARFIMA simulated parameters distribution.
}
\details{
This method facilitates the simulation and evaluation of the uncertainty of 
ARFIMA model parameters. The recursive option also allows the evaluation of the 
simulation based consistency (in terms of sqrt(N) ) of the parameters as the 
length (n.sim) of the data increases, in the sense of the root mean square error 
(rmse) of the difference  between the simulated and true (hypothesized) 
parameters.\cr
This is an expensive function, particularly if using the \code{recursive} 
option, both on memory and CPU resources, performing many re-fits of the 
simulated data in order to generate the parameter distribution.
}
\examples{
\dontrun{
spec = arfimaspec( mean.model = list(armaOrder = c(2,2), include.mean = TRUE, 
arfima = FALSE), distribution.model = "norm", fixed.pars = list(ar1=0.6, 
ar2=0.21, ma1=-0.7, ma2=0.3, mu = 0.02, sigma = 0.02))
dist = arfimadistribution(spec, n.sim = 2000, n.start = 100, m.sim = 100, 
recursive = TRUE, recursive.length = 10000, recursive.window = 1000)
# slots:
slotNames(dist)
# methods:
# summary
show(dist)
# as.data.frame(...., window, which=c("rmse", "stats", "coef", "coefse"))
# default
as.data.frame(dist)

as.data.frame(dist, window = 1, which = "rmse")
as.data.frame(dist, window = 1, which = "stats")
as.data.frame(dist, window = 1, which = "coef")
as.data.frame(dist, window = 1, which = "coefse")


as.data.frame(dist, window = 8, which = "rmse")
as.data.frame(dist, window = 8, which = "stats")
as.data.frame(dist, window = 8, which = "coef")
as.data.frame(dist, window = 8, which = "coefse")


# create some plots
# 
nwindows = dist@dist$details$nwindows
# 2000/3000/4000/5000/6000/7000/8000/9000/10000

# expected reduction factor in RMSE for sqrt(N) consistency
expexcted.rmsegr = sqrt(2000/seq(3000,10000,by=1000))

# actual RMSE reduction
actual.rmsegr = matrix(NA, ncol = 8, nrow = 6)
rownames(actual.rmsegr) = c("mu", "ar1", "ar2", "ma2", "ma2", "sigma")
# start at 2000 (window 1)
rmse.start = as.data.frame(dist, window = 1, which = "rmse")
for(i in 2:nwindows) actual.rmsegr[,i-1] = as.numeric(as.data.frame(dist, 
window = i, which = "rmse")/rmse.start)
par(mfrow = c(2,3))
for(i in 1:6){
	plot(seq(3000,10000,by=1000),actual.rmsegr[i,], type = "l", lty = 2, 
	ylab = "RMSE Reduction", xlab = "N (sim)",main = rownames(actual.rmsegr)[i])
	lines(seq(3000,10000,by=1000), expexcted.rmsegr, col = 2)
	legend("topright", legend = c("Actual", "Expected"), col = 1:2, bty = "m", 
	lty = c(2,1))
	}
}
}
\author{Alexios Ghalanos}
\keyword{methods}