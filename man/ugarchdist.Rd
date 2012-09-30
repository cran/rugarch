\name{rgarchdist}
\alias{rgarchdist}
\alias{ddist}
\alias{pdist}
\alias{qdist}
\alias{rdist}
\alias{fitdist}
\alias{dskewness}
\alias{dkurtosis}
\title{Distribution: rugarch distribution functions}
\description{
Density, distribution function, quantile function, random generation and fitting 
from the univariate distributions implemented in the rugarch package, with 
functions for skewness and excess kurtosis given density skew and shape 
parameters.    
    \tabular{ll}{
    \code{rgarchdist} \tab rugarch univariate distributions, \cr
	\code{fitdist} \tab MLE parameter fit for the rugarch univariate distributions, \cr
    }
}
\usage{
ddist(distribution = "norm", y, mu = 0, sigma = 1, lambda = -0.5, skew = 1, 
shape = 5)
pdist(distribution = "norm", q, mu = 0, sigma = 1, lambda = -0.5, skew = 1, 
shape = 5)
qdist(distribution = "norm", p, mu = 0, sigma = 1, lambda = -0.5, skew = 1, 
shape = 5)
rdist(distribution = "norm", n, mu = 0, sigma = 1, lambda = -0.5, skew = 1, 
shape = 5)
fitdist(distribution = "norm", x, control=list())
dskewness(distribution = "norm", skew = 1, shape = 5, lambda = -0.5)
dkurtosis(distribution = "norm", skew = 1, shape = 5, lambda = -0.5)
}
\arguments{
\item{distribution}{
The distribution name. Valid choices are \dQuote{norm}, \dQuote{snorm}, 
\dQuote{std}, \dQuote{sstd}, \dQuote{ged}, \dQuote{sged}, \dQuote{nig}, 
\dQuote{jsu}.}
\item{mu, sigma, skew, shape}{
location, scale and skewness and shape parameters (see details).}
\item{lambda}{
The additional shape parameter  for the Generalized Hyperbolic and NIG 
distributions.}
\item{n}{
The number of observations.} 
\item{p}{
A numeric vector of probabilities.}
\item{y, q}{
A numeric vector of quantiles.}
\item{x}{
A univariate dataset (for fitting routine).}
\item{control}{
Control parameters passed to the \code{solnp} solver.}   
}
\value{
    \code{d*} returns the density,
    \code{p*} returns the distribution function,
    \code{q*} returns the quantile function, and
    \code{r*} generates random deviates, \cr
    all values are numeric vectors.
    
    \code{fitdist} returns a list with the following components: 
    \cr
    \item{par}{
        The best set of parameters found. 
        }
    \item{value}{
        The likelihood values of the optimization (vector whose length represents 
        the number of major iterations).
        }
    \item{convergence}{
        An integer code. 0 indicates successful convergence.
        }
    \item{lagrange}{
    	The lagrange multiplier value at convergence.
    }
    \item{h}{
    	The hessian at the solution.
    }
	\item{xineq0}{
    	The value of the inequality constraint multiplier (NULL for the 
    	distribution fit problems).
    }
    \code{dskewness} returns the skewness of the distribution. 
    \code{dkurtosis} returns the excess kurtosis of the distribution. 
}
\details{
For the dQuote{nig} and \dQuote{ghyp} distributions, the shape, skew and lambda
are transformed from the \sQuote{zeta-rho} to the \sQuote{alpha-beta} 
parametrization and then scaled by the mean and standard deviation. The fitting 
routines use the \code{solnp} solver and minimize the negative of the 
log-likelihood.
The \dQuote{dskewness} and \dQuote{dkurtosis} functions take as inputs the 
distribution name, skew and shape parameters and return the skewneness and 
excess kurtosis of the distribution. The functions are not at present 
vectorized.
}
\references{
Johnson, N. L. 1954, Systems of frequency curves derived from the first law of 
Laplace, \emph{Trabajos de Estadistica}, \bold{5}, 283--291.\cr
Barndorff-Nielsen, O. E. 1995, Normal inverse Gaussian processes and the 
modeling of stock returns, \emph{mimeo: Univ.of Aarhus Denmark}.\cr
Fernandez C. and Steel, M.F.J. 1998, On Bayesian Modelling of Fat Tails and 
Skewness, \emph{Journal of the American Statistical Association}, 359--371.\cr
}
\author{
Diethelm Wuertz for the Rmetrics \R-port of the \dQuote{norm}, \dQuote{snorm}, 
\dQuote{std}, \dQuote{sstd}, \dQuote{ged}, \dQuote{sged} and \dQuote{nig} 
distrbutions.\cr
Rigby, R. A. and Stasinopoulos D. M for the JSU distribution in the 
\code{gamlss} package.\cr
Alexios Ghalanos for rugarch implementation and higher moment distribution 
functions.
}