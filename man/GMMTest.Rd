\name{GMMTest}
\alias{GMMTest}
\title{
The GMM Orthogonality Test of Hansen
}
\description{
Implements the GMM Orthogonality Test of Hansen.

}
\usage{
GMMTest(z, lags = 1, skew=0, kurt=3, conf.level = 0.95)
}
\arguments{
\item{z}{
A numeric vector the standardized residuals.}
\item{lags}{
The number of lags to test for.}
\item{skew}{
The skewness of the standardized residuals (derived from the estimated model). 
This can be either a scalar or numeric vector the same size as z.}
\item{kurt}{
The kurtosis (not excess) of the standardized residuals (derived from the 
estimated model). This can be either a scalar or numeric vector the same size 
as z.}
\item{conf.level}{
The confidence level at which the Null Hypothesis is evaluated.}
}
\details{
This is a mispecification test based on Hansen's GMM procedure. Under a correctly 
specified model, certain population moment conditions should be satisfied and 
hold in the sample using the standardized residuals. The moment conditions can 
be tested both individually using a t-test or jointly using a Wald test 
(the vignette gives more details). The test returns a matrix containing the 
first 4 moments statistics, their standard error and t-values, \sQuote{M1} to 
\sQuote{M4}, the t-value of the test of the the joint squared residuals under 
the specified number of lags, \sQuote{Q2} , and the t-value of the test of joint 
nullness, \sQuote{J}. The joint tests (\sQuote{Q2} and \sQuote{J}) are 
distributed chi-squared with n.lag and 4 + n.lag d.o.f. respectively, while the 
individual moment conditions may be tested for  significance using a standard 
t-test.
}
\value{
A list with the following items:\cr
\item{statmat}{ The matrix of the statistics (see details).}
\item{H0}{ The Null Hypothesis.}
\item{critical.value}{ The critical value for each statistic.}
\item{Decision}{Whether to reject or not the Null given the conf.level.}
}
\references{
Hansen, L. (1982), Large Sample Properties of Generalized Method of Moments 
Estimators, \emph{Econometrica}, \bold{50(4)}, 1029--1054.
}
\author{Alexios Ghalanos}
\examples{
\dontrun{
data(dji30ret)
spec = ugarchspec(mean.model = list(armaOrder = c(1,1), include.mean = TRUE),
variance.model = list(model = "gjrGARCH"), distribution.model = "sstd")
fit = ugarchfit(spec, data = dji30ret[, 1, drop = FALSE])
z = residuals(fit)\/sigma(fit)
skew = dskewness("sstd",skew = coef(fit)["skew"], shape= coef(fit)["shape"])
# add back 3 since dkurtosis returns the excess kurtosis
kurt = 3+dkurtosis("sstd",skew = coef(fit)["skew"], shape= coef(fit)["shape"])
print(GMMTest(z, lags = 1, skew=skew, kurt=kurt))
}}