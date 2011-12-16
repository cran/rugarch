\name{rugarch-package}
\alias{rugarch-package}
\alias{rugarch}
\title{The rugarch package}
\description{
The rugarch package aims to provide a flexible and rich univariate GARCH 
modelling and testing environment. Modelling is a simple process of defining a 
specification and fitting the data. Inference can be made from summary, various 
tests and plot methods, while the forecasting, filtering and simulation methods 
complete the modelling environment. Finally, specialized methods are implemented 
for simulating parameter distributions and evaluating parameter consistency, 
and a bootstrap forecast method which takes into account both parameter and 
predictive distribution uncertainty.\cr 
The testing environment is based on a rolling backtest function which considers 
the more general context in which GARCH models are based, namely the conditional 
time varying estimation of density parameters and the implication for their use 
in analytical risk management measures.\cr
The mean equation allows for AR(FI)MA, arch-in-mean and external regressors, 
while the variance equation implements a wide variety of univariate GARCH models 
as well as the possibility of including external regressors. Finally, a set of 
rich distributions from the \dQuote{fBasics} package and Johnson's reparametrized 
SU from the \dQuote{gamlss} package are used for modelling innovations.\cr
This package is part of what used to be the rgarch package, which was split 
into univariate (rugarch) and multivariate (rmgarch) models for easier 
maintenance and use. The rmgarch package is still under re-write so the old 
rgarch package should be used in the meantime for multivariate models (and
hosted on r-forge).
}
\details{
\tabular{ll}{
Package: \tab rugarch\cr
Type: \tab Package\cr
Version: \tab 1.00-8\cr
Date: \tab 2011-12-16\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
Depends: \tab R (>= 2.10.0), Rcpp (>= 0.8.5), RcppArmadillo (>= 0.2.5), numDeriv, chron, Rsolnp\cr }
While the package has implemented some safeguards, both during pre-estimation 
as well as the estimation phase, there is no guarantee of convergence in the 
fitting procedure. As a result, the fit method allows the user to input starting 
parameters as well as keep any parameters from the spec as fixed (including 
the case of all parameters fixed).\cr
The functionality of the packages is contained in the main methods for defining 
a specification \code{\link{ugarchspec}}, fitting \code{\link{ugarchfit}}, 
forecasting \code{\link{ugarchforecast}}, simulation from fit object 
\code{\link{ugarchsim}}, path simulation from specification object 
\code{\link{ugarchpath}}, parameter distribution by simulation 
\code{\link{ugarchdistribution}}, bootstrap forecast \code{\link{ugarchboot}} 
and rolling estimation and forecast \code{\link{ugarchroll}}. There are also 
some functions which enable multiple fitting of assets in an easy to use wrapper 
with the option of multicore functionality, namely \code{\link{multispec}}, 
\code{\link{multifit}}, \code{\link{multifilter}} and \code{\link{multiforecast}}. 
Explanations on the available methods for the returned classes can be found in 
the documentation for those classes.\cr
A separate subset of methods and classes has been included to calculate pure 
ARFIMA models with constant variance. This subset includes similar functionality 
as with the GARCH methods, with the exception that no plots are yet implemented, 
and neither is a forecast based on the bootstrap. These may be added in the 
future. While there are limited examples in the documentation on the ARFIMA 
methods, the interested user can search the rugarch.tests folder of the source 
installation for some tests using ARFIMA models as well as equivalence to the 
base R arima methods (particularly replication of simulation). Finally, no 
representation is made about the adequacy of ARFIMA models, particularly the 
statistical properties of parameters when using distributions which go beyond 
the Gaussian.\cr
The conditional distributions used in the package are also exposed for the 
benefit of the user through the \code{\link{rgarchdist}} functions which contain 
methods for density, distribution, quantile, sampling and fitting. Additionally, 
\code{\link{ghyptransform}} function provides the necessary parameter 
transformation and scaling methods for moving from the location scale invariant 
\sQuote{rho-zeta} parametrization with mean and standard deviation, to the 
standard \sQuote{alpha-beta-delta-mu} parametrization of the Generalized 
Hyperbolic Distribution family.\cr
The type of data handled by the package is quite varied, accepting \dQuote{timeSeries}, 
\dQuote{xts}, \dQuote{zoo}, \dQuote{zooreg}, \dQuote{data.frame} with dates as 
rownames, \dQuote{matrix} and \dQuote{numeric} vector with dates as names. 
For the \dQuote{numeric} vector and \dQuote{data.frame} with characterdates in 
names or rownames, the package tries a variety of methods to try to recognize 
the type and format of the date else will index the data numerically. The 
package holds dates internally as class \code{\link{Date}}. This mostly impacts 
the plots and forecast summary methods. For high frequency data, the user should 
make use of a non-named representation such as \dQuote{matrix} or \dQuote{numeric} 
as the package has yet to implement methods for checking and working with 
frequencies higher than daily (and is unlikely to do so). Finally, the functions
\code{\link{ForwardDates}} and \code{\link{WeekDayDummy}} offer some simple Date 
manipulation methods for working with forecast dates and creating day of the 
week dummy variables for use in GARCH modelling.\cr
Some benchmarks (published and comparison with commercial package), are 
available through the \code{\link{ugarchbench}} function. The \sQuote{inst} 
folder of the source distribution also contains various tests which can be 
sourced and run by the user, also exposing some finer details of the 
functionality of the package. The user should really consult the examples 
supplied in this folder which are quite numerous and instructive with some 
comments.
}
\section{How to cite this package}{
  Whenever using this package, please cite as\cr

 \preformatted{@Manual{Ghalanos_2011,
 author       = {Alexios Ghalanos},
 title        = {{rugarch}: Univariate GARCH models.},
 year         = {2011},
 note 	      = {R package version 1.00.},}}
\emph{}
}
\section{License}{
  The releases of this package is licensed under GPL version 3.
}
\author{Alexios Ghalanos}
\references{
Baillie, R.T. and Bollerslev, T. and Mikkelsen,H.O. 1996, Fractionally 
integrated generalized autoregressive conditional heteroskedasticity, 
\emph{Journal of Econometrics}, 3--30 .\cr
Berkowitz, J. 2001, Testing density forecasts, with applications to risk 
management, \emph{Journal of Business and Economic Statistics}, \bold{19(4)}, 
465--474.\cr
Bollerslev, T. 1986, Generalized Autoregressive Conditional Heteroskedasticity 
1986, \emph{Journal of Econometrics}, \bold{31}, 307--327.\cr
Ding, Z., Granger, C.W.J. and Engle, R.F. 1993, A Long Memory Property of Stock 
Market Returns and a New Model, \emph{Journal of Empirical Finance}, \bold{1}, 
83--106.\cr
Engle, R.F. and Ng, V. K. 1993, Measuring and Testing the Impact of News on 
Volatility, \emph{Journal of Finance}, \bold{48}, 1749--1778.\cr
Glosten, L.R., Jagannathan, R. and Runkle, D.E. 1993, On the Relation between 
the Expected Value and the Volatility of the Nominal Excess Return on Stocks, 
\emph{Journal of Finance}, \bold{48(5)}, 1779--1801.\cr
Hansen, B.E. 1990, Langrange Multiplier Tests for Parameter Instability in 
Non-Linear Models, \emph{mimeo}.\cr
Hentschel, Ludger. 1995, All in the family Nesting symmetric and asymmetric 
GARCH models, \emph{Journal of Financial Economics}, \bold{39(1)}, 71--104.\cr
Nelson, D.B. 1991, Conditional Heteroskedasticity in Asset Returns: A New 
Approach, \emph{Econometrica}, \bold{59}, 347--370.\cr
Pascual, L., Romo, J. and Ruiz, E. 2004, Bootstrap predictive inference for 
ARIMA processes, \emph{Journal of Time Series Analysis}.\cr
Pascual, L., Romo, J. and Ruiz, E. 2006, Bootstrap prediction for returns and 
volatilities in GARCH models, \emph{Computational Statistics and Data Analysis}.\cr
Vlaar, P.J.G. and Palm, F.C. 1993, The Message in Weekly Exchange Rates in the 
European Monetary System: Mean Reversion Conditional Heteroskedasticity and 
Jumps, \emph{Journal of Business and Economic Statistics}, \bold{11}, 351--360.\cr
}