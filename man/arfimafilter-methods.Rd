\name{arfimafilter-methods}
\docType{methods}
\alias{arfimafilter-methods}
\alias{arfimafilter,ANY-method}
\alias{arfimafilter,ARFIMAspec-method}
\alias{arfimafilter}
\title{function: ARFIMA Filtering}
\description{
Method for filtering an ARFIMA model.
}
\usage{
arfimafilter(spec, data, out.sample = 0, n.old=NULL, ...)
}
\arguments{
\item{data}{ A univariate data object. Can be a numeric vector, matrix, 
data.frame, zoo, xts, timeSeries, ts or irts object.}
\item{spec}{ An ARFIMA spec object of class \code{\linkS4class{ARFIMAspec}} 
with the fixed.pars argument having the model parameters on which the filtering 
is to take place.}
\item{out.sample}{ A positive integer indicating the number of periods before 
the last to keep for out of sample forecasting (as in \code{\link{arfimafit}} 
function).}
\item{n.old}{ For comparison with ARFIMA models using the out.sample argument,
this is the length of the original dataset (see details).}
\item{...}{ . }
}
\value{
A \code{\linkS4class{ARFIMAfilter}} object containing details of the ARFIMA 
filter.
}
\details{
The n.old argument is optional and indicates the length of the original data
(in cases when this represents a dataseries augmented by newer data). The reason
for using this is so that the old and new datasets agree since the original 
recursion uses the sum of the residuals to start the recursion and therefore is 
influenced by new data. For a small augmentation the values converge after x 
periods, but it is sometimes preferable to have this option so that there is no 
forward looking information contaminating the study.
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


# filter the data to check results
filt = vector(mode = "list", length = 9)
for(i in 1:9){
	spec = arfimaspec(mean.model = list(armaOrder = c(1,1), include.mean = TRUE, 
	arfima = FALSE), distribution.model = dist[i])
	setfixed(spec) = as.list(coef(fit[[i]]))
	filt[[i]] = arfimafilter(spec = spec, data = sp500ret)
}
print(cfmatrix, digits = 4)
cat("\nARFIMAfit and ARFIMAfilter residuals check:\n")
print(head(sapply(filt, FUN = function(x) residuals(x))) == head(sapply(fit, 
FUN = function(x) residuals(x))))
cat("\nas.data.frame method:\n")
print(cbind(head(as.data.frame(filt[[1]])), head(as.data.frame(fit[[1]]))))
cat("\ncoef method:\n")
print(cbind(coef(filt[[1]]), coef(fit[[1]])))
cat("\nfitted method:\n")
print(cbind(head(fitted(filt[[1]])), head(fitted(fit[[1]]))))
cat("\ninfocriteria method:\n")
# For filter, it assumes estimation of parameters else does not make sense!
print(cbind(infocriteria(filt[[1]]), infocriteria(fit[[1]])))
cat("\nlikelihood method:\n")
print(cbind(likelihood(filt[[1]]), likelihood(fit[[1]])))
cat("\nresiduals method:\n")
# Note that we the package will always return the full length residuals and 
# fitted object irrespective of the lags (i.e. since this is an ARMA(1,1) 
# i.e. max lag = 1, the first row is zero and should be discarded.
print(cbind(head(residuals(filt[[1]])), head(residuals(fit[[1]]))))
cat("\nuncmean method:\n")
print(cbind(uncmean(filt[[1]]), uncmean(fit[[1]])))
cat("\nuncmean method (by simulation):\n")
# For spec and fit
spec = arfimaspec( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, 
arfima = FALSE), distribution.model = dist[1])
setfixed(spec) = as.list(coef(fit[[1]]))
print(cbind(uncmean(spec, method = "simulation", n.sim = 100000, rseed = 100), 
uncmean(fit[[1]], method = "simulation", n.sim = 100000, rseed = 100)))
cat("\nsummary method:\n")
show(filt[[1]])
show(fit[[1]])
}
}
\author{Alexios Ghalanos}
\keyword{methods}
