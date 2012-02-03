\name{arfimaforecast-methods}
\docType{methods}
\alias{arfimaforecast-methods}
\alias{arfimaforecast,ANY-method}
\alias{arfimaforecast,ARFIMAfit-method}
\alias{arfimaforecast,ARFIMAspec-method}
\alias{arfimaforecast}
\title{function: ARFIMA Forecasting}
\description{
Method for forecasting from an ARFIMA model.
}
\usage{
arfimaforecast(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0, 
external.forecasts = list(mregfor = NULL), ...)
}
\arguments{
\item{fitORspec}{ Either an ARFIMA fit object of class 
\code{\linkS4class{ARFIMAfit}} or alternatively an ARFIMA specification object 
of class \code{\linkS4class{ARFIMAspec}} with valid parameters supplied via the 
\code{fixed.pars} argument in the specification.}
\item{data}{ Required if a specification rather than a fit object is supplied.}
\item{n.ahead}{ The forecast horizon.}
\item{n.roll}{ The no. of rolling forecasts to create beyond the first one (see 
details).}
\item{out.sample}{ Optional. If a specification object is supplied, indicates 
how many data points to keep for out of sample testing.}
\item{external.forecasts}{ A list with a matrix of forecasts for the external 
regressors in the mean.}
\item{...}{ . }
}
\value{
A \code{\linkS4class{ARFIMAforecast}} object containing details of the ARFIMA 
forecast. See the class for details on the returned object and methods for 
accessing it and performing some tests.
}
\details{
The forecast function has two dispatch methods allowing the user to call it with 
either a fitted object (in which case the data argument is ignored), or a 
specification object (in which case the data is required) with the parameters 
entered via the \code{set.fixed<-} methods on an \code{\linkS4class{ARFIMAspec}} 
object.\cr
One step ahead forecasts are based on the value of the previous data, while 
n-step ahead (n>1) are based on the unconditional mean of the model.\cr
The ability to roll the forecast 1 step at a time is implemented with the 
\code{n.roll} argument which controls how many times to roll the n.ahead 
forecast. The default argument of n.roll = 0 denotes no rolling beyond the first 
forecast and returns the standard n.ahead forecast. Critically, since n.roll 
depends on data being available from which to base the rolling forecast, the 
\code{\link{arfimafit}} function needs to be called with the argument 
\code{out.sample} being at least as large as the n.roll argument, or in the case 
of a specification being used instead of a fit object, the \code{out.sample}
argument directly in the forecast function.
}
\examples{
\dontrun{
# Long Horizon Forecast
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

umean = rep(0, 9)
for(i in 1:9){
	umean[i] = uncmean(fit[[i]])
}

forc = vector(mode = "list", length = 9)
for(i in 1:9){
	forc[[i]] = arfimaforecast(fit[[i]], n.ahead = 100)
}

lmean40 = sapply(forc, FUN = function(x) as.numeric(as.data.frame(x)[40,1]))
cfmatrix1 = cbind(cfmatrix, umean, lmean40)
colnames(cfmatrix1) = c(colnames(cfmatrix1[,1:7]), "uncmean", "forecast40")

# forecast with spec to check results
forc2 = vector(mode = "list", length = 9)
for(i in 1:9){
	spec = arfimaspec(mean.model = list(armaOrder = c(1,1), include.mean = TRUE, 
	arfima = FALSE), distribution.model = dist[i])
	setfixed(spec) = as.list(coef(fit[[i]]))
	forc2[[i]] = arfimaforecast(spec, data = sp500ret, n.ahead = 100)
}
lmean240 = sapply(forc2, FUN = function(x) as.numeric(as.data.frame(x)[40,1]))
cfmatrix2 = cbind(cfmatrix, umean, lmean240)
colnames(cfmatrix2) = c(colnames(cfmatrix2[,1:7]), "uncmean", "forecast40")

cat("\nARFIMAforecast from ARFIMAfit and ARFIMAspec check:")
cat("\nFit\n")	
print(cfmatrix1, digits = 4)
cat("\nSpec\n")	
print(cfmatrix2, digits = 4)
# methods and slots
slotNames(forc[[1]])
showMethods(classes="ARFIMAforecast")
# summary
show(forc[[1]])
# Extractor Functions
# as array (array dimension [3] is 1 since n.roll = 0 i.e. no rolling beyond 
# the first)
as.array(forc[[1]])
# as.data.frame
as.data.frame(forc[[1]])
# as.list
as.list(forc[[1]])

# Rolling Forecast
data(sp500ret)
fit = vector(mode = "list", length = 9)
dist = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "jsu")
for(i in 1:9){
	spec = arfimaspec(mean.model = list(armaOrder = c(1,1), include.mean = TRUE, 
			arfima = FALSE), distribution.model = dist[i])
	fit[[i]] = arfimafit(spec = spec, data = sp500ret, solver = "solnp", 
			out.sample = 1000, fit.control = list(scale = 1))
}
cfmatrix = matrix(NA, nrow = 9, ncol = 7)
colnames(cfmatrix) = c("mu", "ar1", "ma1", "sigma", "skew", "shape", "ghlambda")
rownames(cfmatrix) = dist

for(i in 1:9){
	cf = coef(fit[[i]])
	cfmatrix[i, match(names(cf), colnames(cfmatrix))] =  cf
}


forc = vector(mode = "list", length = 9)
for(i in 1:9){
	forc[[i]] = arfimaforecast(fit[[i]], n.ahead = 1, n.roll = 999)
}
rollforc = sapply(forc, FUN = function(x) t(unlist(as.data.frame(x, 
	rollframe = "all", aligned = FALSE))))

# forecast performance measures:
fpmlist = vector(mode = "list", length = 9)
for(i in 1:9){
	fpmlist[[i]] = fpm(forc[[i]], summary = FALSE)
}

par(mfrow = c(1,2))
dd = rownames(tail(sp500ret, 1250))
clrs = rainbow(9, alpha = 1, start = 0.4, end = 0.95)
plot(as.Date(dd), tail(sp500ret[,1], 1250), type = "l", 
	ylim = c(-0.02, 0.02), col = "lightgrey", ylab = "", xlab = "", 
	main = "Rolling 1-ahead Forecasts\nvs Actual")
for(i in 1:9){
	tmp = tail(sp500ret[,1], 1250)
	tmp[251:1250] = rollforc[1:1000,i]
	lines(as.Date(dd), c(rep(NA, 250), tmp[-(1:250)]), col = clrs[i])
}
legend("topleft", legend = dist, col = clrs, fill = clrs, bty = "n")

# plot deviation measures and range
tmp = vector(mode = "list", length = 9)
for(i in 1:9){
	tmp[[i]] = fpmlist[[i]][,"AE"]
	names(tmp[[i]]) = dist[i]
}
boxplot(tmp, col = clrs, names = dist, range  = 6, notch = TRUE, 
		main = "Rolling 1-ahead Forecasts\nAbsolute Deviation Loss")

# fpm comparison
compm = matrix(NA, nrow = 3, ncol = 9)
compm = sapply(fpmlist, FUN = function(x) c(mean(x[,"SE"]), mean(x[,"AE"]), 
	mean(x[,"DAC"])))
colnames(compm) = dist
rownames(compm) = c("MSE", "MAD", "DAC")

cat("\nRolling Forecast FPM\n")
print(compm, digits = 4)
cat("\nMethods Check\n")
as.data.frame(forc[[1]], rollframe = 0)
as.data.frame(forc[[1]], rollframe = 999)
t(as.data.frame(forc[[1]], rollframe = "all", aligned = FALSE))
fpm(forc[[1]], summary = TRUE)
show(forc[[1]])
}}
\author{Alexios Ghalanos}
\keyword{methods}