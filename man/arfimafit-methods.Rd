\name{arfimafit-methods}
\docType{methods}
\alias{arfimafit-methods}
\alias{arfimafit,ANY-method}
\alias{arfimafit,ARFIMAspec-method}
\alias{arfimafit}
\title{function: ARFIMA Fit}
\description{
Method for fitting an ARFIMA models.
}
\usage{
arfimafit(spec, data, out.sample = 0, solver = "solnp", solver.control = list(), 
fit.control = list(fixed.se = 0, scale = 0), ...)
}
\arguments{
\item{data}{ A univariate data object. Can be a numeric vector, matrix, 
data.frame, zoo, xts, timeSeries, ts or irts object.}
\item{spec}{ An ARFIMA spec object of class \code{\linkS4class{ARFIMAspec}}.}
\item{out.sample}{ A positive integer indicating the number of periods before 
the last to keep for out of sample forecasting (see details).}
\item{solver}{One of either \dQuote{nlminb}, \dQuote{solnp}, \dQuote{gosolnp}
or \dQuote{nloptr}.}
\item{solver.control}{Control arguments list passed to optimizer.}
\item{fit.control}{Control arguments passed to the fitting routine. The 
fixed.se argument controls whether standard errors should be calculated for 
those parameters which were fixed (through the fixed.pars argument of the 
\code{\link{arfimaspec} function). The scale parameter controls whether the
data should be scaled before being submitted to the optimizer.}}
\item{...}{ . }
}
\value{
A \code{\linkS4class{ARFIMAfit}} object containing details of the ARFIMA fit.
}
\details{
The ARFIMA optimization routine first calculates a set of feasible starting 
points which are used to initiate the ARFIMA Maximum Likelihood recursion. The 
main part of the likelihood calculation is performed in C-code for speed.\cr
The out.sample option is provided in order to carry out forecast performance 
testing against actual data. A minimum of 5 data points are required for these 
tests. If the out.sample option is positive, then the routine will fit only 
N - out.sample (where N is the total data length) data points, leaving 
out.sample points for forecasting and testing using the forecast performance 
measures. In the \code{\link{arfimaforecast}} routine the 
n.ahead may also  be greater than the out.sample number resulting in a 
combination of out of sample data points matched against actual data and some 
without, which the forecast performance tests will ignore.\cr
The \dQuote{gosolnp} solver allows for the initialization of multiple restarts 
of the solnp solver with randomly generated parameters (see documentation in 
the Rsolnp-package for details of the strategy used). The solver.control list 
then accepts the following additional (to the solnp) arguments: \dQuote{n.restarts} 
is the number of solver restarts required (defaults to 1), \dQuote{parallel} and 
\dQuote{parallel.control} for use of the parallel functionality, \dQuote{rseed} 
is the seed to initialize the random number generator, and \dQuote{n.sim} is the 
number of simulated parameter vectors to generate per n.restarts.
}
\examples{
\dontrun{
data(sp500ret)
fit = vector(mode = "list", length = 9)
dist = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "jsu")
for(i in 1:9){
	spec = arfimaspec(mean.model = list(armaOrder = c(1,1), include.mean = TRUE, 
	arfima = FALSE), distribution.model = dist[i])
	fit[[i]] = arfimafit(spec = spec, data = sp500ret, solver = "solnp", 
	fit.control = list(scale = 1))
}
cfmatrix = matrix(NA, nrow = 9, ncol = 7)
colnames(cfmatrix) = c("mu", "ar1", "ma1", "sigma", "skew", "shape", "ghlambda")
rownames(cfmatrix) = dist

for(i in 1:9){
	cf = coef(fit[[i]])
	cfmatrix[i, match(names(cf), colnames(cfmatrix))] =  cf
}
sk = ku = rep(0, 9)
for(i in 1:9){
	cf = coef(fit[[i]])
	if(fit[[i]]@model$modelinc[16]>0) sk[i] = dskewness(distribution = dist[i], 
				skew = cf["skew"], shape = cf["shape"], lambda = cf["ghlambda"])		
	if(fit[[i]]@model$modelinc[17]>0) ku[i] = dkurtosis(distribution = dist[i], 
				skew = cf["skew"], shape = cf["shape"], lambda = cf["ghlambda"])
}
hq = sapply(fit, FUN = function(x) infocriteria(x)[4])
cfmatrix = cbind(cfmatrix, sk, ku, hq)
colnames(cfmatrix)=c(colnames(cfmatrix[,1:7]), "skewness", "ex.kurtosis","HQIC")

print(cfmatrix, digits = 4)
# notice that for the student distribution kurtosis is NA since shape (dof) < 4.
cat("\nas.data.frame method:\n")
head(as.data.frame(fit[[1]]))
cat("\ncoef method:\n")
coef(fit[[1]])
cat("\nfitted method:\n")
head(fitted(fit[[1]]))
cat("\ninfocriteria method:\n")
infocriteria(fit[[1]])
cat("\nlikelihood method:\n")
likelihood(fit[[1]])
cat("\nresiduals method:\n")
# Note that we the package will always return the full length residuals and 
# fitted object irrespective of the lags (i.e. since this is an ARMA(1,1) 
# i.e. max lag = 1, the first row is zero and should be discarded).
head(residuals(fit[[1]]))
cat("\nuncmean method:\n")
uncmean(fit[[1]])
cat("\nuncmean method (by simulation):\n")
uncmean(fit[[1]], method = "simulation", n.sim = 100000, rseed = 100)
cat("\nsummary method:\n")
show(fit[[1]])
}
}
\author{Alexios Ghalanos}
\keyword{methods}