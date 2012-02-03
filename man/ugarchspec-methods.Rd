\name{ugarchspec-methods}
\docType{methods}
\alias{ugarchspec}
\alias{ugarchspec-methods}
\alias{ugarchspec,ANY-method}
\title{function: Univariate GARCH Specification}
\description{
Method for creating a univariate GARCH specification object prior to fitting.
}
\usage{
ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1), 
submodel = NULL, external.regressors = NULL, variance.targeting = FALSE), 
mean.model = list(armaOrder = c(1, 1), include.mean = TRUE, archm = FALSE, 
archpow = 1, arfima = FALSE, external.regressors = NULL, archex = FALSE), 
distribution.model = "norm", start.pars = list(), fixed.pars = list(), ...)
}
\arguments{
\item{variance.model}{
List containing the variance model specification:\cr
\code{model} Valid models (currently implemented) are \dQuote{sGARCH}, 
\dQuote{fGARCH}, \dQuote{eGARCH}, \dQuote{gjrGARCH}, \dQuote{apARCH} and
\dQuote{iGARCH}.\cr
\code{garchOrder} The ARCH (q) and GARCH (p) orders.\cr
\code{submodel} If the model is \dQuote{fGARCH}, valid submodels are 
\dQuote{GARCH}, \dQuote{TGARCH}, \dQuote{AVGARCH}, \dQuote{NGARCH}, 
\dQuote{NAGARCH}, \dQuote{APARCH},\dQuote{GJRGARCH} and \dQuote{ALLGARCH}.\cr
\code{external.regressors} A matrix object containing the external regressors to 
include in the variance equation with as many rows as will be included in the
data (which is passed in the fit function).
\code{variance.targeting} Indicates whether to use variance targeting for the 
sigma intercept \dQuote{omega}.).
}
\item{mean.model}{
List containing the mean model specification:\cr
\code{armaOrder} The autoregressive (ar) and moving average (ma) orders (if any).\cr
\code{include.mean} Whether to include the mean.\cr
\code{archm} Whether to include ARCH volatility in the mean regression.\cr
\code{archpow} Indicates whether to use st.deviation (1) or variance (2) in the 
ARCH in mean regression.\cr
\code{arfima} Whether to fractional differencing in the ARMA regression.\cr
\code{external.regressors} A matrix object containing the external regressors to 
include in the mean equation with as many rows as will be included in the data 
(which is passed in the fit function).\cr
\code{archex} (integer) Whether to multiply the last 'archex' external regressors
by the conditional standard deviation.
}
\item{distribution.model}{
The conditional density to use for the innovations. Valid choices are 
\dQuote{norm} for the normal distibution, \dQuote{snorm} for the skew-normal 
distribution, \dQuote{std} for the student-t, \dQuote{sstd} for the skew-student, 
\dQuote{ged} for the generalized error distribution, \dQuote{sged} for the 
skew-generalized error distribution, \dQuote{nig} for the normal inverse 
gaussian distribution, \dQuote{ghyp} for the Generalized  Hyperbolic, and \dQuote{jsu} 
for Johnson's SU distribution. Note that some of the distributions are taken 
from the fBasics package and implenented locally here for convenience. 
The \dQuote{jsu} distribution is the reparametrized version from the 
\dQuote{gamlss} package.}
\item{start.pars}{
List of staring parameters for the optimization routine. These are not usually 
required unless the optimization has problems converging.}
\item{fixed.pars}{
List of parameters which are to be kept fixed during the optimization. It is 
possible that you designate all parameters as fixed so as to quickly recover 
just the results of some previous work or published work. The optional argument
\dQuote{fixed.se} in the \code{\link{ugarchfit}} function indicates whether to
calculate standard errors for those parameters fixed during the post 
optimization stage.}
\item{...}{ . }
}
\value{
A \code{\linkS4class{uGARCHspec}} object containing details of the GARCH 
specification.
}
\details{
The specification allows for a wide choice in univariate GARCH models, 
distributions, and mean equation modelling. For the \dQuote{fGARCH} model, 
this represents Hentschel's omnibus model which subsumes many others.\cr
For the mean equation, ARFIMAX is fully supported in fitting, forecasting and 
simulation. There is also an option to multiply the external regressors by
the conditional standard deviation, which may be of use for example in 
calculating the  correlation coefficient in a CAPM type setting.\cr
The \dQuote{iGARCH} implements the integrated GARCH model. For the \dQuote{EWMA} 
model just set \dQuote{omega} to zero in the fixed parameters list.\cr
The asymmetry term in the rugarch package, for all implemented models, follows 
the order of the arch parameter \code{alpha}.\cr
Variance targeting, referred to in Engle and Mezrich (1996), replaces the 
intercept \dQuote{omega} in the variance equation by 1 minus the persistence 
multiplied by the unconditional variance which is calculated by its sample 
counterpart in the squared residuals during estimation. In the presence of 
external regressors in the variance equation, the sample average of the external 
regresssors is multiplied by their coefficient and subtracted from the 
variance target.\cr
In order to understand which parameters can be entered in the start.pars and 
fixed.pars optional arguments, the list below exposes the names used for the
parameters across the various models:(note that when a parameter is followed by 
a number, this represents the order of the model. Just increment the number 
for higher orders):\cr
\emph{\bold{Mean Model}}:\cr
   \tabular{rll}{
     \tab constant	\tab \code{mu}\cr
     \tab AR term	\tab \code{ar1}\cr
     \tab MA term	\tab \code{ma1}\cr
     \tab ARCH in mean	\tab \code{archm}\cr
     \tab exogenous regressors	\tab \code{mxreg1}\cr
     \tab arfima	\tab \code{arfima}\cr
   }
\emph{\bold{Distribution Model}}:\cr
   \tabular{rll}{
     \tab ghlambda \tab \code{lambda (for GHYP distribution)}\cr
     \tab skew	\tab \code{skew}\cr
     \tab shape	\tab \code{shape}\cr
   }
\emph{\bold{Variance Model (common specs)}}:\cr
  \tabular{rll}{
     \tab constant	\tab \code{omega}\cr
     \tab ARCH term	\tab \code{alpha1}\cr
     \tab GARCH term \tab \code{beta1}\cr
     \tab exogenous regressors	\tab \code{vxreg1}\cr
   }
\emph{\bold{Variance Model (GJR, EGARCH)}}:\cr
  \tabular{rll}{
     \tab assymetry term \tab \code{gamma1}\cr
   }
\emph{\bold{Variance Model (APARCH)}}:\cr
  \tabular{rll}{
     \tab assymetry term \tab \code{gamma1}\cr
     \tab power term	\tab \code{delta}\cr
   }
\emph{\bold{Variance Model (FGARCH)}}:\cr
  \tabular{rll}{
     \tab assymetry term1 (rotation)	\tab \code{eta11}\cr
     \tab assymetry term2 (shift)	\tab \code{eta21}\cr
     \tab power term1(shock)	\tab \code{delta}\cr
     \tab power term2(variance)\tab \code{lambda}\cr
   }
}
\author{Alexios Ghalanos}
\examples{
# a standard specification
spec1 = ugarchspec()
spec1
# an example which keep the ar1 and ma1 coefficients fixed:
spec2 = ugarchspec(mean.model=list(armaOrder=c(2,2), 
fixed.pars=list(ar1=0.3,ma1=0.3)))
spec2
# an example of the EWMA Model
spec3 = ugarchspec(variance.model=list(model="iGARCH", garchOrder=c(1,1)), 
		mean.model=list(armaOrder=c(0,0), include.mean=TRUE),  
		distribution.model="norm", fixed.pars=list(omega=0))
}
\keyword{methods}