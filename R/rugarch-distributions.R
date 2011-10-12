#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008, 2009, 2010, 2011
##   This file is part of the R package rugarch.
##
##   The R package rugarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rugarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

# Local version START
#################################################################################
## from fBasics library: normal, skew-normal, student, skew-student, ged, skew-ged,
## nig, skew-nig, and sgh/gh distributions locally implemented in rugarch
#################################################################################
## Distributions functions from Rmetrics Libraries
## Copyrights (C)
##   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
##   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
##   info@rmetrics.org
##   www.rmetrics.org
#################################################################################

.dsnorm <-function(x, xi) 
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Compute the density function of the "normalized" skew 
	#   normal distribution
	
	# FUNCTION:
	
	# Standardize:
	m1 = 2/sqrt(2*pi)
	mu = m1 * (xi - 1/xi)
	sigma = sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	z = x*sigma + mu  
	# Compute:
	Xi = xi^sign(z)
	g = 2 / (xi + 1/xi) 
	Density = g * dnorm(x = z/Xi)  
	# Return Value:
	Density * sigma 
}
# ------------------------------------------------------------------------------
dsnorm <-function(x, mean = 0, sd = 1, xi = 1.5)
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Compute the density function of the skew normal distribution
	
	# Arguments:
	#   x - a numeric vector of quantiles.
	#   mean, sd, xi - location parameter, scale parameter, and 
	#       skewness parameter.
	
	# FUNCTION:
	
	# Shift and Scale:
	result = .dsnorm(x = (x-mean)/sd, xi = xi) / sd
	
	# Return Value:
	result
}
# ------------------------------------------------------------------------------
.psnorm <-function(q, xi) 
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Internal Function
	
	# FUNCTION:
	
	# Standardize:
	m1 = 2/sqrt(2*pi)
	mu = m1 * (xi - 1/xi)
	sigma = sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	z = q*sigma + mu
	# Compute:  
	Xi = xi^sign(z)
	g = 2  / (xi + 1/xi)  
	Probability = Heaviside(z) - sign(z) * g * Xi * pnorm(q = -abs(z)/Xi)
	# Return Value:
	Probability 
}
# ------------------------------------------------------------------------------
psnorm <-function(q, mean = 0, sd = 1, xi = 1.5)
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Compute the distribution function of the 
	#   skew normal distribution
	
	# Arguments:
	#   q - a numeric vector of quantiles.
	#   mean, sd, xi - location parameter, scale parameter, and 
	#       skewness parameter.
	
	# FUNCTION:
	
	# Shift and Scale:
	result = .psnorm(q = (q-mean)/sd, xi = xi)
	
	# Return Value:
	result
}
# ------------------------------------------------------------------------------    
.qsnorm <-function(p, xi) 
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Internal Function
	
	# FUNCTION:
	
	# Standardize:
	m1 = 2/sqrt(2*pi)
	mu = m1 * (xi - 1/xi)
	sigma = sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	# Compute:  
	g = 2  / (xi + 1/xi)
	sig = sign(p-1/2) 
	Xi = xi^sig         
	p = (Heaviside(p-1/2)-sig*p) / (g*Xi)
	Quantile = (-sig*qnorm(p = p, sd = Xi) - mu ) / sigma
	# Return Value:
	Quantile 
}
# ------------------------------------------------------------------------------
qsnorm <-function(p, mean = 0, sd = 1, xi = 1.5)
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Compute the quantile function of the 
	#   skew normal distribution
	
	# Arguments:
	#   p - a numeric vector of probabilities.
	#   mean, sd, xi - location parameter, scale parameter, and 
	#       skewness parameter.
	
	# FUNCTION:
	
	# Shift and Scale:
	result = .qsnorm(p = p, xi = xi) * sd + mean
	
	# Return Value:
	result
}
# ------------------------------------------------------------------------------
.rsnorm <-function(n, xi) 
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Internal Function
	
	# FUNCTION:
	
	# Generate Random Deviates:
	weight = xi / (xi + 1/xi)
	z = runif(n, -weight, 1-weight)
	Xi = xi^sign(z)
	Random = -abs(rnorm(n))/Xi * sign(z)  
	# Scale:
	m1 = 2/sqrt(2*pi)
	mu = m1 * (xi - 1/xi)
	sigma = sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	Random = (Random - mu ) / sigma   
	# Return value:
	Random 
}
# ------------------------------------------------------------------------------
rsnorm <-function(n, mean = 0, sd = 1, xi = 1.5)
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Generate random deviates from the 
	#   skew normal distribution
	
	# Arguments:
	#   n - an integer value giving the number of observation.
	#   mean, sd, xi - location parameter, scale parameter, and 
	#       skewness parameter.
	
	# FUNCTION:
	
	# Shift and Scale:
	result = .rsnorm(n = n, xi = xi) * sd + mean
	
	# Return Value:
	result
}

normFit <-function(x, control = list())
{   
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Fit the parameters for a skew Normal distribution
	
	# FUNCTION:
	
	# Start Value:
	ctrl = .solnpctrl(control)
	start = c(mean = mean(x), sd = sqrt(var(x)))
	
	# Log-likelihood Function:
	loglik = function(x, y = x){ 
		f = -sum(log(dnorm(y, x[1], x[2])))
		f }
	
	# Minimization:
	fit = solnp(pars = start, fun = loglik, LB = c(-Inf, 0), UB = c(Inf, Inf), 
			control = ctrl, y = x)
	
	# Add Names to $par
	names(fit$par) = c("mean", "sd")
	
	# Return Value:
	fit
}   
# ------------------------------------------------------------------------------
snormFit = function(x, control = list())
{   
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Fit the parameters for a skew Normal distribution
	
	# FUNCTION:
	
	# Start Value:
	ctrl = .solnpctrl(control)
	start = c(mean = mean(x), sd = sqrt(var(x)), xi = 1)
	
	# Log-likelihood Function:
	loglik = function(x, y = x){ 
		f = -sum(log(dsnorm(y, x[1], x[2], x[3])))
		f }
	
	# Minimization:
	fit = solnp(pars = start, fun = loglik, 
			LB = c(-Inf, 0, 0), upper = c(Inf, Inf, Inf), control = ctrl, y = x)
	
	# Add Names to $par
	names(fit$par) = c("mean", "sd", "xi")
	
	# Return Value:
	fit
}
################################################################################
dged <-function(x, mean = 0, sd = 1, nu = 2)
{   
	# A function imlemented by Diethelm Wuertz
	
	# Description:
	#   Compute the density for the 
	#   generalized error distribution.
	
	# FUNCTION:
	
	# Compute Density:
	z = (x - mean ) / sd
	lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
	g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
	result = g * exp (-0.5*(abs(z/lambda))^nu) / sd
	
	# Return Value
	result
}
# ------------------------------------------------------------------------------
pged <-function(q, mean = 0, sd = 1, nu = 2)
{   
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Compute the probability for the  
	#   generalized error distribution.
	
	# FUNCTION:
	
	# Compute Probability:
	q = (q - mean ) / sd
	lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
	g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
	h = 2^(1/nu) * lambda * g * gamma(1/nu) / nu
	s = 0.5 * ( abs(q) / lambda )^nu
	result = 0.5 + sign(q) * h * pgamma(s, 1/nu)
	
	# Return Value:
	result
}


# ------------------------------------------------------------------------------
qged <-function(p, mean = 0, sd = 1, nu = 2)
{   
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Compute the quantiles for the  
	#   generalized error distribution.
	
	# FUNCTION:
	
	# Compute Quantiles:
	lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
	q = lambda * (2*qgamma((abs(2*p-1)), 1/nu))^(1/nu)
	result = q*sign(2*p-1) * sd + mean
	
	# Return Value:
	result
}


# ------------------------------------------------------------------------------
rged <-function(n, mean = 0, sd = 1, nu = 2)
{   
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Generate GED random deviates. The function uses the 
	#   method based on the transformation of a Gamma random 
	#   variable.
	
	# FUNCTION:
	
	# Generate Random Deviates:
	lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
	# print(lambda)
	r = rgamma(n, 1/nu)
	z =  lambda * (2*r)^(1/nu) * sign(runif(n)-1/2)
	result = z * sd + mean
	
	
	# Return Value:
	result
}


# ------------------------------------------------------------------------------
.dsged <-function(x, nu, xi) 
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Internal Function
	
	# FUNCTION:
	
	# Standardize:
	lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
	g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
	m1 = 2^(1/nu) * lambda * gamma(2/nu) / gamma(1/nu)
	mu = m1*(xi-1/xi)
	sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	z = x*sigma + mu  
	
	# Compute:
	Xi = xi^sign(z)
	g = 2  / (xi + 1/xi)    
	Density = g * dged(x = z/Xi, nu=nu)  
	
	# Return Value:
	Density * sigma 
}

# norm [ nu = 2, xi = 1 ]
# laplace [ nu = 1, xi = 1 ]
dsged <-function(x, mean = 0, sd = 1, nu = 2, xi = 1.5)
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Compute the density function of the 
	#   skewed generalized error distribution
	
	# FUNCTION:
	
	# Shift and Scale:
	result = .dsged(x = (x-mean)/sd, nu = nu, xi = xi) / sd
	
	# Return Value:
	result
}


# ------------------------------------------------------------------------------
.psged <-function(q, nu, xi) 
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Internal Function
	
	# FUNCTION:
	
	# Standardize:
	lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
	g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
	m1 = 2^(1/nu) * lambda * gamma(2/nu) / gamma(1/nu)
	mu = m1*(xi-1/xi)
	sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	z = q*sigma + mu
	
	# Compute:  
	Xi = xi^sign(z)
	g = 2  / (xi + 1/xi)    
	Probability = Heaviside(z) - sign(z) * g * Xi * pged(q = -abs(z)/Xi, nu=nu)
	
	# Return Value:
	Probability 
}

psged <-function(q, mean = 0, sd = 1, nu = 2, xi = 1.5)
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Compute the distribution function of the 
	#   skewed generalized error distribution
	
	# FUNCTION:
	
	# Shift and Scale:
	result = .psged(q = (q-mean)/sd, nu = nu, xi = xi)
	
	# Return Value:
	result
}


# ------------------------------------------------------------------------------    
.qsged <-function(p, nu, xi) 
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Internal Function
	
	# FUNCTION:
	
	# Standardize:
	lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
	g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
	m1 = 2^(1/nu) * lambda * gamma(2/nu) / gamma(1/nu)
	mu = m1*(xi-1/xi)
	sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	
	# Compute:  
	g = 2  / (xi + 1/xi)
	sig = sign(p-1/2) 
	Xi = xi^sig       
	p = (Heaviside(p-1/2)-sig*p) / (g*Xi)
	Quantile = (-sig*qged(p=p, sd=Xi, nu=nu) - mu ) / sigma
	
	# Return Value:
	Quantile 
}

qsged <-function(p, mean = 0, sd = 1, nu = 2, xi = 1.5)
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Compute the quantile function of the 
	#   skewed generalized error distribution
	
	# FUNCTION:
	
	# Shift and Scale:
	result = .qsged(p = p, nu = nu, xi = xi) * sd + mean
	
	# Return Value:
	result
}


# ------------------------------------------------------------------------------
.rsged <-function(n, nu, xi) 
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Internal Function
	
	# FUNCTION:
	
	# Generate Random Deviates:
	weight = xi / (xi + 1/xi)
	z = runif(n, -weight, 1-weight)
	Xi = xi^sign(z)
	Random = -abs(rged(n, nu=nu))/Xi * sign(z)  
	
	# Scale:
	lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
	g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
	m1 = 2^(1/nu) * lambda * gamma(2/nu) / gamma(1/nu)
	mu = m1*(xi-1/xi)
	sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	Random = (Random - mu ) / sigma 
	
	# Return value:
	Random 
}


# ------------------------------------------------------------------------------
rsged <-function(n, mean = 0, sd = 1, nu = 2, xi = 1.5)
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Generate random deviates from the 
	#   skewed generalized error distribution
	
	# FUNCTION:
	
	# Shift and Scale:
	result = .rsged(n = n, nu = nu, xi = xi) * sd + mean
	
	# Return Value:
	result
}

gedFit = function(x, control = list())
{   
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Fit the parameters for a skew Normal distribution
	
	# FUNCTION:
	
	# Start Value:
	ctrl = .solnpctrl(control)
	start = c(mean = mean(x), sd = sqrt(var(x)), nu = 2)
	
	# Log-likelihood Function:
	loglik = function(x, y = x){ 
		f = -sum(log(dged(y, x[1], x[2], x[3])))
		f }
	
	# Minimization:
	fit = solnp(pars = start, fun = loglik, 
			LB = c(-Inf, 0, 0), UB = c(Inf, Inf, Inf), , control  = ctrl, y = x)
	
	# Add Names to $par
	names(fit$par) = c("mean", "sd", "nu")
	
	# Return Value:
	fit
}      


# ------------------------------------------------------------------------------
sgedFit = function(x, control = list())
{   
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Fit the parameters for a skew Normal distribution
	
	# FUNCTION:
	
	# Start Value:
	ctrl = .solnpctrl(control)
	start = c(mean = mean(x), sd = sqrt(var(x)), nu = 2, xi = 1)
	
	# Log-likelihood Function:
	loglik = function(x, y = x){ 
		f = -sum(log(dsged(y, x[1], x[2], x[3], x[4])))
		f }
	
	# Minimization:
	fit = solnp(pars = start, fun = loglik, 
			LB = c(-Inf, 0, 0, 0), UB = c(Inf, Inf, Inf, Inf), , control  = ctrl, y = x)
	
	# Add Names to $par
	names(fit$par) = c("mean", "sd", "nu", "xi")
	
	# Return Value:
	fit
}
################################################################################
Heaviside<-function(x, a = 0) 
{   
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Computes the Heaviside or unit step function.
	
	# Arguments:
	#   x - a numeric vector.
	#   a - the location of the break.
	
	# Details:
	#   The Heaviside step function is 1 for x>a, 1/2 for x=a,
	#   and 0 for x<a.
	
	# Notes:
	#   Heaviside Unit Step Function and Related Functions
	#   See:  http://mathworld.wolfram.com/HeavisideStepFunction.html
	#   Note: sign(x) is part of R's base package
	
	# FUNCTION:
	
	# Compute H:
	result = (sign(x-a) + 1)/2
	
	# Return Value:
	result
}

dstd <- function(x, mean = 0, sd = 1, nu = 5)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Compute the density for the
	#   Student-t distribution.
	
	# FUNCTION:
	
	# Compute Density:
	s = sqrt(nu/(nu-2))
	z = (x - mean) / sd
	# result = .Internal(dnt(x = z*s, df = nu, ncp = 0, log = FALSE)) / (sd/s)
	result = dt(x = z*s, df = nu) * s / sd
	
	# Return Value:
	result
}


# ------------------------------------------------------------------------------
pstd <-function (q, mean = 0, sd = 1, nu = 5)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Compute the probability for the
	#   Student-t distribution.
	
	# FUNCTION:
	
	# Compute Probability:
	s = sqrt(nu/(nu-2))
	z = (q - mean) / sd
	# result = .Internal(pnt(q = z*s, df = nu, ncp = 0, lower.tail = TRUE,
	#   log.p = FALSE))
	result = pt(q = z*s, df = nu)
	
	# Return Value:
	result
}


# ------------------------------------------------------------------------------
qstd <-function (p, mean = 0, sd = 1, nu = 5)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Compute the quantiles for the
	#   Student-t distribution.
	
	# FUNCTION:
	
	# Compute Quantiles:
	s = sqrt(nu/(nu-2))
	# x = .Internal(qt(p = p, df = nu, lower.tail = TRUE, log.p = FALSE)) / s
	# result = x*sd + mean
	result = qt(p = p, df = nu) * sd / s + mean
	
	# Return Value:
	result
}


# ------------------------------------------------------------------------------
rstd <-function(n, mean = 0, sd = 1, nu = 5)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Generate random deviates from the
	#   Student-t distribution.
	
	# FUNCTION:
	
	# Generate Random Deviates:
	s = sqrt(nu/(nu-2))
	# result = .Internal(rt(n = n, df = nu)) * sd / s + mean
	result = rt(n = n, df = nu) * sd / s + mean
	
	# Return Value:
	result
}
# ------------------------------------------------------------------------------
.dsstd <-function(x, nu, xi)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Internal Function
	
	# FUNCTION:
	
	# For SPlus compatibility:
	if (!exists("beta"))
		beta <- function (a, b) exp( lgamma(a) + lgamma(b) -lgamma(a+b) )
	
	# Standardize:
	m1 = 2 * sqrt(nu-2) / (nu-1) / beta(1/2, nu/2)
	mu = m1*(xi-1/xi)
	sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	z = x*sigma + mu
	
	# Compute:
	Xi = xi^sign(z)
	g = 2 / (xi + 1/xi)
	Density = g * dstd(x = z/Xi, nu = nu)
	
	# Return Value:
	Density * sigma
}


# ------------------------------------------------------------------------------
dsstd <-function(x, mean = 0, sd = 1, nu = 5, xi = 1.5)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Compute the density function of the
	#   skewed Student-t distribution
	
	# FUNCTION:
	
	# Shift and Scale:
	result = .dsstd(x = (x-mean)/sd, nu = nu, xi = xi) / sd
	
	# Return Value:
	result
}


# ------------------------------------------------------------------------------
.psstd <-function(q, nu, xi)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Internal Function
	
	# FUNCTION:
	
	# For SPlus compatibility:
	if (!exists("beta"))
		beta <- function (a, b) exp( lgamma(a) + lgamma(b) -lgamma(a+b) )
	
	# Standardize:
	m1 = 2 * sqrt(nu-2) / (nu-1) / beta(1/2, nu/2)
	mu = m1*(xi-1/xi)
	sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	z = q*sigma + mu
	
	# Compute:
	Xi = xi^sign(z)
	g = 2 / (xi + 1/xi)
	Probability = Heaviside(z) - sign(z) * g * Xi * pstd(q = -abs(z)/Xi, nu = nu)
	
	# Return Value:
	Probability
}


# ------------------------------------------------------------------------------
psstd <-function(q, mean = 0, sd = 1, nu = 5, xi = 1.5)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Compute the distribution function of the
	#   skewed Student-t distribution
	
	# FUNCTION:
	
	# Shift and Scale:
	result = .psstd(q = (q-mean)/sd, nu = nu, xi = xi)
	
	# Return Value:
	result
}
# ------------------------------------------------------------------------------
.qsstd <-function(p, nu, xi)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Internal Function
	
	# FUNCTION:
	
	# For SPlus compatibility:
	if (!exists("beta"))
		beta <- function (a, b) exp( lgamma(a) + lgamma(b) -lgamma(a+b) )
	
	# Standardize:
	m1 = 2 * sqrt(nu-2) / (nu-1) / beta(1/2, nu/2)
	mu = m1*(xi-1/xi)
	sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	
	# Compute:
	g = 2  / (xi + 1/xi)
	sig = sign(p-1/2)
	Xi = xi^sig
	p = (Heaviside(p-1/2)-sig*p) / (g*Xi)
	Quantile = (-sig*qstd(p = p, sd = Xi, nu = nu) - mu ) / sigma
	
	# Return Value:
	Quantile
}
# ------------------------------------------------------------------------------
qsstd <-function(p, mean = 0, sd = 1, nu = 5, xi = 1.5)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Compute the quantile function of the
	#   skewed Student-t distribution
	
	# FUNCTION:
	
	# Shift and Scale:
	result = .qsstd(p = p, nu = nu, xi = xi) * sd + mean
	
	# Return Value:
	result
}


# ------------------------------------------------------------------------------
.rsstd <-function(n, nu, xi)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Internal Function
	
	# FUNCTION:
	
	# For SPlus compatibility:
	if (!exists("beta"))
		beta <- function (a, b) exp( lgamma(a) + lgamma(b) -lgamma(a+b) )
	
	# Generate Random Deviates:
	weight = xi / (xi + 1/xi)
	z = runif(n, -weight, 1-weight)
	Xi = xi^sign(z)
	Random = -abs(rstd(n, nu = nu))/Xi * sign(z)
	
	# Scale:
	m1 = 2 * sqrt(nu-2) / (nu-1) / beta(1/2, nu/2)
	mu = m1*(xi-1/xi)
	sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	Random = (Random - mu ) / sigma
	
	# Return value:
	Random
}


# ------------------------------------------------------------------------------
rsstd <-function(n, mean = 0, sd = 1, nu = 5, xi = 1.5)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Generate random deviates from the
	#   skewed Student-t distribution
	
	# FUNCTION:
	
	# Shift and Scale:
	result = .rsstd(n = n, nu = nu, xi = xi) * sd + mean
	
	# Return Value:
	result
}

stdFit = function(x, control = list())
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Fit the parameters for a skew Normal distribution
	
	# FUNCTION:
	
	# Start Value:
	ctrl = .solnpctrl(control)
	start = c(mean = mean(x), sd = sqrt(var(x)), nu = 4)
	
	# Log-likelihood Function:
	loglik = function(x, y = x){
		f = -sum(log(dstd(y, x[1], x[2], x[3])))
		f }
	
	# Minimization:
	fit = solnp(pars = start, fun = loglik,
			LB = c(-Inf, 0, 2.01), UB = c(Inf, Inf, Inf), control = ctrl, y = x)
	
	# Add Names to $par
	names(fit$par) = c("mean", "sd", "nu")
	
	# Return Value:
	fit
}


# ------------------------------------------------------------------------------
sstdFit = function(x, control = list())
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Fit the parameters for a skew Sudent-t distribution
	#   with unit variance
	
	# FUNCTION:
	
	ctrl = .solnpctrl(control)
	
	# Start Value:
	start = c(mean = mean(x), sd = sqrt(var(x)), nu = 4, xi = 1)
	
	# Log-likelihood Function:
	loglik = function(x, y = x){
		f = -sum(log(dsstd(y, x[1], x[2], x[3], x[4])))
		f }
	
	# Minimization:
	#fit = nlm(f = loglik, p = p, y = x, ...)
	fit = solnp(pars = start, fun = loglik, LB = c(-Inf, 0, 2.01, 0.1), 
			UB = c(Inf, Inf, Inf, Inf), control = ctrl, y = x)
	Names = c("mean", "sd", "nu", "xi")
	names(fit$par) = Names
	#names(fit$gradient) = Names
	
	# Return Value:
	fit
}
################################################################################
.unitroot<-function(f, interval, lower = min(interval), upper = max(interval), 
		tol = .Machine$double.eps^0.25, ...)
{   
	if (is.null(args(f))) {
		if (f(lower) * f(upper) >=0) return(NA)
	} else {
		if (f(lower, ...) * f(upper, ...) >= 0) return(NA)
	}
	ans = uniroot(f = f, interval = interval, lower = lower,
			upper = upper, tol = tol, ...)
	ans$root
}

.kappaGH <-function(x, lambda = 1)
{    
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Returns modified Bessel function ratio
	
	# FUNCTION:
	
	# Check:
	stopifnot(x >= 0)
	stopifnot(length(lambda) == 1)
	
	# Ratio:
	if (lambda == -0.5) {
		# NIG:
		kappa = 1/x
	} else {
		# GH:
		kappa = (
					besselK(x, lambda+1, expon.scaled = TRUE) /
					besselK(x, lambda, expon.scaled = TRUE) ) / x
	}
	
	# Return Value:
	kappa
}
# ------------------------------------------------------------------------------
.deltaKappaGH<-function(x, lambda = 1)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Returns difference of Bessel functions ratios
	
	# FUNCTION:
	
	# Difference in Ratios:
	if (lambda == -0.5) {
		# NIG:
		# Replace this with the recursion relation ...
		deltaKappa = .kappaGH(x, lambda+1) - .kappaGH(x, lambda)
	} else {
		# GH:
		deltaKappa = .kappaGH(x, lambda+1) - .kappaGH(x, lambda)
	}
	
	# Return Value:
	deltaKappa
}
# ------------------------------------------------------------------------------
.paramGH <-function(zeta = 1, rho = 0 , lambda = 1)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Change parameterizations to alpha(zeta, rho, lambda)
	
	# FUNCTION:
	
	# Transformation:
	Rho2 = 1 - rho^2
	alpha = zeta^2 * .kappaGH(zeta, lambda) / Rho2 
	alpha = alpha * ( 1 + rho^2 * zeta^2 * .deltaKappaGH(zeta, lambda) / Rho2)
	alpha = sqrt(alpha)  
	beta = alpha * rho
	delta = zeta / ( alpha * sqrt(Rho2) )
	mu = -beta * delta^2 * .kappaGH(zeta, lambda)
	
	# Return Value:
	c(alpha = alpha, beta = beta, delta = delta, mu = mu)  
}

dgh<-function(x, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1, log = FALSE)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Returns density for the generalized hyperbolic distribution
	
	# FUNCTION:
	
	# Checks:
	if (alpha <= 0) stop("alpha must be greater than zero")
	if (delta <= 0) stop("delta must be greater than zero")
	if (abs(beta) >= alpha) stop("abs value of beta must be less than alpha")
	
	# Density:
	arg = delta*sqrt(alpha^2-beta^2)
	a = (lambda/2)*log(alpha^2-beta^2) - (
				log(sqrt(2*pi)) + (lambda-0.5)*log(alpha) + lambda*log(delta) +
				log(besselK(arg, lambda, expon.scaled = TRUE)) - arg )
	f = ((lambda-0.5)/2)*log(delta^2+(x - mu)^2)
	
	# Use exponential scaled form to prevent from overflows:
	arg = alpha * sqrt(delta^2+(x-mu)^2)
	k = log(besselK(arg, lambda-0.5, expon.scaled = TRUE)) - arg
	e = beta*(x-mu)
	
	# Put all together:
	ans = a + f + k + e
	if(!log) ans = exp(ans)
	
	# Return Value:
	ans
}
# ------------------------------------------------------------------------------
pgh<-function(q, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Returns probability for the generalized hyperbolic distribution
	
	# FUNCTION:
	
	# Checks:
	if (alpha <= 0) stop("alpha must be greater than zero")
	if (delta <= 0) stop("delta must be greater than zero")
	if (abs(beta) >= alpha) stop("abs value of beta must be less than alpha")
	
	# Probability:
	ans = NULL
	for (Q in q) {
		Integral = integrate(dgh, -Inf, Q, stop.on.error = FALSE,
				alpha = alpha, beta = beta, delta = delta, mu = mu,
				lambda = lambda)
		ans = c(ans, as.numeric(unlist(Integral)[1]))
	}
	
	# Return Value:
	ans
}
# ------------------------------------------------------------------------------
qgh<-function (p, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Returns quantiles for the generalized hyperbolic distribution
	
	# FUNCTION:
	
	# Checks:
	if (alpha <= 0) stop("alpha must be greater than zero")
	if (delta <= 0) stop("delta must be greater than zero")
	if (abs(beta) >= alpha) stop("abs value of beta must be less than alpha")
	
	# Internal Function:
	.froot <-
			function(x, alpha, beta, delta, mu, lambda, p)
	{
		pgh(q = x, alpha = alpha, beta = beta, delta = delta,
				mu = mu, lambda = lambda) - p
	}
	
	# Quantiles:
	result = NULL
	for (pp in p) {
		lower = -1
		upper = +1
		counter = 0
		iteration = NA
		while (is.na(iteration)) {
			iteration = .unitroot(f = .froot, interval = c(lower,
							upper), alpha = alpha, beta = beta, delta = delta,
					mu = mu, lambda = lambda, p = pp)
			counter = counter + 1
			lower = lower - 2^counter
			upper = upper + 2^counter
		}
		result = c(result, iteration)
	}
	# Return Value:
	result
}

ghFit = function(x, control = list())
{
	ctrl = .solnpctrl(control)
	
	start = c(alpha = 1, beta = 0, delta = sqrt(var(x)), mu = mean(x), lambda = -2)
	eps = .Machine$double.eps^0.5
	BIG = 10000
	# Log-likelihood Function:
	loglik = function(pars, y){
		#.pars<<-pars
		if(abs(pars[2])>pars[1]) return(100000)
		f =  -sum(log(dgh(y, pars[1], pars[2], pars[3], pars[4], pars[5], log = FALSE)))
		f 
	}
	con = function(pars, y){
		ans = abs(pars[2]) - pars[1]
		ans
	}
	# Minimization:
	fit = solnp(pars = start, fun = loglik, LB = c(eps, -BIG, eps, -BIG, -6), 
			UB = c(BIG,  BIG, BIG,  BIG,  6), ineqfun = con, ineqLB = -BIG, 
			ineqUB = -0.1, control = ctrl, y = x)
	Names = c("alpha", "beta", "delta", "mu", "lambda")
	names(fit$par) = Names
	#names(fit$gradient) = Names
	
	# Return Value:
	fit
	
}

.rghyp = function(n, theta)
{	# A function implemented by Diethelm Wuertz
	
	# Author:
	#	Original Version by David Scott
	
	# FUNCTION:
	
	# Settings:
	lambda = theta[1]
	alpha = theta[2]
	beta = theta[3]
	delta = theta[4]
	mu = theta[5]
	chi = delta^2
	psi = alpha^2 - beta^2
	
	# Ckecks:
	if (alpha <= 0) stop("alpha must be greater than zero")  
	if (delta <= 0) stop("delta must be greater than zero")
	if (abs(beta) >= alpha) stop("abs value of beta must be less than alpha")
	
	# Random Numbers:
	if (lambda == 1){
		X = .rgigjd1(n, c(lambda, chi, psi))
	} else{
		X = .rgigjd(n, c(lambda, chi, psi))
	}
	
	# Result:
	sigma = sqrt(X)
	Z = rnorm(n)
	Y = mu + beta*sigma^2 + sigma*Z
	
	# Return Value:
	Y
}
.rgigjd = function(n, theta)
{	# A function implemented by Diethelm Wuertz
	
	# Author:
	#	Original Version by David Scott
	
	# FUNCTION:
	
	# Settings:
	lambda = theta[1]
	chi = theta[2]
	psi = theta[3]
	
	# Checks:
	if (chi < 0) stop("chi can not be negative")
	if (psi < 0) stop("psi can not be negative")
	if ((lambda >= 0)&(psi==0)) stop("When lambda >= 0, psi must be > 0")
	if ((lambda <= 0)&(chi==0)) stop("When lambda <= 0, chi must be > 0")
	if (chi == 0) stop("chi = 0, use rgamma")
	if (psi == 0) stop("algorithm only valid for psi > 0")
	
	alpha = sqrt(psi/chi)
	beta = sqrt(psi*chi)
	
	m = (lambda-1+sqrt((lambda-1)^2+beta^2))/beta
	
	g = function(y){
		0.5*beta*y^3 - y^2*(0.5*beta*m+lambda+1) + y*((lambda-1)*m-0.5*beta) + 0.5*beta*m
	}
	
	upper = m
	while (g(upper) <= 0) upper = 2*upper
	yM = uniroot(g, interval=c(0,m))$root
	yP = uniroot(g, interval=c(m,upper))$root
	
	a = (yP-m)*(yP/m)^(0.5*(lambda-1))*exp(-0.25*beta*(yP+1/yP-m-1/m))
	b = (yM-m)*(yM/m)^(0.5*(lambda-1))*exp(-0.25*beta*(yM+1/yM-m-1/m))
	c = -0.25*beta*(m+1/m) + 0.5*(lambda-1)*log(m)
	
	output = numeric(n)
	
	for(i in 1:n){
		need.value = TRUE
		while(need.value==TRUE){
			R1 = runif (1)
			R2 = runif (1)
			Y = m + a*R2/R1 + b*(1-R2)/R1
			if (Y>0){
				if (-log(R1)>=-0.5*(lambda-1)*log(Y)+0.25*beta*(Y+1/Y)+c){
					need.value = FALSE
				}
			}
		}
		output[i] = Y
	}
	
	# Return Value:
	output/alpha
}

# ------------------------------------------------------------------------------
.rgigjd1 = function(n, theta)
{	# A function implemented by Diethelm Wuertz
	
	# Description:
	# 	Modified version of rgigjd to generate random observations
	# 	from a generalised inverse Gaussian distribution in the
	# 	special case where lambda = 1.
	
	# Author:
	#	Original Version by David Scott
	
	# FUNCTION:
	
	if (length(theta) == 2) theta = c(1, theta)
	
	# Settings:
	lambda = 1
	chi = theta[2]
	psi = theta[3]
	
	# Checks:
	if (chi < 0) stop("chi can not be negative")
	if (psi < 0) stop("psi can not be negative")	
	if (chi == 0) stop("chi = 0, use rgamma")
	if (psi == 0) stop("When lambda >= 0, psi must be > 0")
	
	alpha = sqrt(psi/chi)
	beta = sqrt(psi*chi)
	m = abs(beta)/beta
	g = function(y){
		0.5*beta*y^3 - y^2*(0.5*beta*m+lambda+1) +
				y*(-0.5*beta) + 0.5*beta*m
	}
	
	upper = m
	while (g(upper)<=0) upper = 2*upper
	yM = uniroot(g,interval=c(0,m))$root
	yP = uniroot(g,interval=c(m,upper))$root
	
	a = (yP-m)*exp(-0.25*beta*(yP+1/yP-m-1/m))
	b = (yM-m)*exp(-0.25*beta*(yM+1/yM-m-1/m))
	c = -0.25*beta*(m+1/m)
	
	output = numeric(n)
	
	for(i in 1:n){
		need.value = TRUE
		while(need.value==TRUE){
			R1 = runif (1)
			R2 = runif (1)
			Y = m + a*R2/R1 + b*(1-R2)/R1
			if (Y>0){
				if (-log(R1)>=0.25*beta*(Y+1/Y)+c){
					need.value = FALSE
				}
			}
		}
		output[i] = Y
	}
	
	# Return Value:
	output/alpha
}



# ------------------------------------------------------------------------------
rgh<-function (n, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Returns random variates for the generalized hyperbolic distribution
	
	# FUNCTION:
	
	# Checks:
	if (alpha <= 0) stop("alpha must be greater than zero")
	if (delta <= 0) stop("delta must be greater than zero")
	if (abs(beta) >= alpha) stop("abs value of beta must be less than alpha")
	
	# Settings:
	theta = c(lambda, alpha, beta, delta, mu)
	
	# Random Numbers:
	ans = .rghyp(n, theta)
	# Return Value:
	ans
}

# ------------------------------------------------------------------------------
dsgh <-function(x, zeta = 1, rho = 0, lambda = 1, log = FALSE) 
{
	param = .paramGH(zeta, rho, lambda)
	dgh(x, param[1], param[2], param[3], param[4], lambda, log)
}
# ------------------------------------------------------------------------------
psgh <-function(q, zeta = 1, rho = 0, lambda = 1) 
{
	param = .paramGH(zeta, rho, lambda)
	pgh(q, param[1], param[2], param[3], param[4], lambda)
}
# ------------------------------------------------------------------------------
qsgh <-function(p, zeta = 1, rho = 0, lambda = 1) 
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Returns quantiles of the sgh distribution
	
	# FUNCTION:
	
	# Compute Quantiles:
	param = .paramGH(zeta, rho, lambda)
	qgh(p, param[1], param[2], param[3], param[4], lambda)
}
# ------------------------------------------------------------------------------
rsgh <-function(n, zeta = 1, rho = 0, lambda = 1) 
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Generates sgh distributed random variates
	
	# FUNCTION:
	
	# Generate Random Numbers:
	param = .paramGH(zeta, rho, lambda)
	rgh(n, param[1], param[2], param[3], param[4], lambda)
}
# ------------------------------------------------------------------------------
sghFit = function (x, zeta = 1, rho = 0, lambda = 1, include.lambda = TRUE, 
		scale = TRUE, span = "auto", trace = FALSE, title = NULL, description = NULL, ...) 
{
	x.orig = x
	x = as.vector(x)
	if (scale) x = (x-mean(x)) / sd(x)
	
	eps = .Machine$double.eps^0.5
	BIG = 1000
	
	if (include.lambda) 
	{
		# LLH Function:
		esghmle.include = function(x, y = x) {
			f = -sum(log(dsgh(y, x[1], x[2], x[3], log = FALSE)))
			f
		}
		# LLH Optimization:
		fit = solnp(
				pars = c(zeta, rho, lambda), 
				fun = esghmle.include, 
				LB = c(eps, -0.9999, -2), 
				UB = c(BIG, +0.9999, +5), 
				y = x)
		names(fit$par) <- c("zeta", "rho", "lambda")
		
	} else {
		
		# LLH Function:
		esghmle = function(x, y = x, lambda) {
			f = -sum(log(dsgh(y, x[1], x[2], lambda, log = FALSE)))
			f
		}
		# LLH Optimization:
		fit = solnp(
				pars = c(zeta, rho), 
				Jfun = esghmle, 
				LB = c(eps, -0.9999), 
				UB = c(BIG, +0.9999), 
				y = x, 
				lambda = lambda)
		fit$par = c(fit$par, lambda)
		names(fit$par) <- c("zeta", "rho", "fix.lambda")
		
	}
	
	param = .paramGH(fit$par[1], fit$par[2], fit$par[3])
	
	# Default Title and Description:
	if (is.null(title)) 
		title = "SGH Parameter Estimation"
	# Fit:
	ans = list(
			estimate = fit$par,
			minimum = -fit$values[length(fit$values)], 
			code = fit$convergence,
			param = param, 
			mean = mean(x.orig),
			var = var(x.orig))
	return(ans)
}


nigFit = function(x, control = list())
{
	ctrl = .solnpctrl(control)
	
	start = c(alpha = 1, beta = 0, delta = sqrt(var(x)), mu = mean(x))
	
	# Log-likelihood Function:
	loglik = function(x, y = x){
		f = -sum(log(dnig(y, x[1], x[2], x[3], x[4])))
		f }
	con = function(x, y = x){
		abs(x[2])-x[1]
	}
	# Minimization:
	fit = solnp(pars = start, fun = loglik, LB = c(0.1, -100, eps, -Inf), 
			UB = c(100, 100, Inf, Inf), ineqfun = con, ineqLB = -200, ineqUB = 0,control = ctrl, y = x)
	Names = c("alpha", "beta", "delta", "mu")
	names(fit$par) = Names
	#names(fit$gradient) = Names
	
	# Return Value:
	fit
	
}

dnig = function(x, alpha = 1, beta = 0, delta = 1, mu = 0, log = FALSE)
{   
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Returns density for inverse Gaussian DF
	
	# FUNCTION:
	
	# Density:
	dgh(x = x, alpha = alpha, beta = beta, delta = delta, mu = mu, 
			lambda = -0.5, log = log)
}
# ------------------------------------------------------------------------------
pnig <-function(q, alpha = 1, beta = 0, delta = 1, mu = 0)
{   
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Returns probability for for inverse Gaussian DF
	
	# Function:
	
	# Probability:
	pgh(q = q, alpha = alpha, beta = beta, delta = delta, mu = mu, 
			lambda = -0.5)
}
# ------------------------------------------------------------------------------
qnig <-function(p, alpha = 1, beta = 0, delta = 1, mu = 0)
{   
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Returns quantiles for for inverse Gaussian DF
	
	# FUNCTION:
	
	# Quantiles:
	.qnigC(p = p, alpha = alpha, beta = beta, delta = delta, mu = mu)
}
# ------------------------------------------------------------------------------
rnig <-function(n, alpha = 1, beta = 0, delta = 1, mu = 0)
{   
	rgh(n, alpha = alpha, beta = beta, delta = delta, mu = mu, lambda=-0.5)
}
################################################################################
.qnigC <-function(p, alpha = 1, beta = 0, delta = 1, mu = 0)
{   
	# Description:
	#   Returns quantiles for for inverse Gaussian DF
	
	# FUNCTION:
	
	# Checks:
	if(alpha <= 0) stop("Invalid parameters: alpha <= 0.\n")
	if(alpha^2 <= beta^2) stop("Invalid parameters: alpha^2 <= beta^2.\n")
	if(delta <= 0) stop("Invalid parameters: delta <= 0.\n")
	if((sum(is.na(p)) > 0)) 
		stop("Invalid probabilities:\n",p,"\n")
	else 
	if(sum(p < 0)+sum(p > 1) > 0) stop("Invalid probabilities:\n",p,"\n")
	
	n <- length(p)
	q <- rep(0, n)
	
	# Evaluate NIG cdf by calling C function
	retValues <- .C("qNIG",
			.CArrange(p,1,1,n),
			as.double(mu),
			as.double(delta),
			as.double(alpha),
			as.double(beta),
			as.integer(n),
			.CArrange(q, 1, 1, n), PACKAGE="rugarch")
	quantiles <- retValues[[7]]
	quantiles[quantiles <= -1.78e+308] <- -Inf
	quantiles[quantiles >= 1.78e+308] <- Inf
	
	# Return Value:
	quantiles
}
# ------------------------------------------------------------------------------
.CArrange <-function(obj, i, j, n)
{
	# Description:
	#   Arrange input matrices and vectors in a suitable way for the C program
	#   Matrices are transposed because the C program stores matrices by row 
	#   while R stores matrices by column
	
	# Arguments:
	#   i - length of first dimension
	#   j - length of second dimension
	#   n - length of third dimension
	
	# Value:
	#   out - transformed data set
	
	# Author: 
	#   Daniel Berg <daniel at nr.no> (Kjersti Aas <Kjersti.Aas at nr.no>)
	#   Date: 12 May 2005
	#   Version: 1.0.2
	
	# FUNCTION:
	
	if(is.null(obj)) stop("Missing data")
	
	if(is.vector(obj)) {
		if(i==1 & j==1 & length(obj)==n) out <- as.double(obj)
		else stop("Unexpected length of vector")
	} else if(is.matrix(obj)) {
		if(nrow(obj) == i && ncol(obj) == j) out <- as.double(rep(t(obj), n))
		else stop("Unexpected dimensions of matrix")
	} else {
		stop("Unexpected object")
	}
	
	# Return Value:
	out 
}
################################################################################
# FUNCTION:            DESCRIPTION:
#  dsnig                Returns density of the SNIG distribution
#  psnig                Returns probabilities of the SNIG distribution
#  qsnig                Returns quantiles of the SNIG distribution
#  rsnig                Generates SNIG distributed random variates
# FUNCTION:            DESCRIPTION:
#  .qsnigC              Fast qsnig from C code
################################################################################
dsnig<-function(x, zeta = 1, rho = 0, log = FALSE) 
{
	# Description:
	#   Returns density of the snig distribution
	
	# FUNCTION:
	
	# Compute Density - Quick and Dirty:
	dsgh(x, zeta, rho, lambda = -0.5, log = log)
}
# ------------------------------------------------------------------------------
psnig <-function(q, zeta = 1, rho = 0) 
{
	# Description:
	#   Returns probabilities of the snig distribution
	
	# FUNCTION:
	
	# Compute Probabilities - Quick and Dirty:
	psgh(q, zeta, rho, lambda = -0.5)
}
# ------------------------------------------------------------------------------
qsnig <-function(p, zeta = 1, rho = 0) 
{
	# Description:
	#   Returns quantiles of the snig distribution
	
	# FUNCTION:
	
	# Compute Quantiles:
	qsgh(p, zeta, rho, lambda = -0.5)
}
# ------------------------------------------------------------------------------
rsnig <-function(n, zeta = 1, rho = 0) 
{
	# Description:
	#   Generates snig distributed random variates
	
	# FUNCTION:
	
	# Generate Random Numbers:
	rsgh(n, zeta, rho, lambda = -0.5)
}
################################################################################
.qsnigC <-function(p, zeta = 1, rho = 0) 
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Returns quantiles of the snig distribution
	
	# FUNCTION:
	
	# Compute Quantiles:
	param = .paramGH(zeta, rho, lambda = -0.5)
	.qnigC(p, param[1], param[2], param[3], param[4])
}

snigFit =
		function (x, zeta = 1, rho = 0, scale = TRUE, doplot = TRUE, 
				span = "auto", trace = 0, title = NULL, description = NULL, ...) 
{   
	
	# Update Slots: 
	if (is.null(title)) title = "SNIG Parameter Estimation"
	
	# Quick and dirty ...
	ans = sghFit(x, zeta = zeta, rho = rho, lambda = -0.5, include.lambda = FALSE,
			scale = scale,span = span, trace = trace, 
			title = title, description = description, ...)
	
	# Update Slots:    

	
	# Return Value:
	ans
}



################################################################################
# Local version END

################################################################################
# Johnson SU Distribution (from gamlss package) - reparametrized version so that mu=mean and sigma=sd

jsuFit = function(x, control = list())
{
	# a function implemented by Alexios Ghalanos
	ctrl = .solnpctrl(control)
	start = c(mean(x), sd(x), nu = 1, tau = 0.5)
	
	loglik = function(y, z){
		-sum(djsu(y = y, z[1], z[2], z[3], z[4], log = TRUE))
	}
	
	fit = solnp(pars = start, fun = loglik, LB = c(-Inf, 0, -20, 0.1),
			UB = c(Inf, Inf, 20, 100), control = ctrl, y = x)
	names(fit$par) = c("mean", "sd", "nu","tau")
	
	# Return Value:
	fit
}
# to use in the c-code
djsu = function(y, mu = 0, sigma = 1, nu = 1, tau = 0.5, log = FALSE)
{
	if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
	if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
	rtau <- 1/tau
	if (length(tau)>1)
		w <- ifelse(rtau<0.0000001,1,exp(rtau^2))       
	else w <- if (rtau<0.0000001) 1  else  exp(rtau^2)
	omega <- -nu*rtau 
	c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
	z <- (y-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
	r <- -nu + asinh(z)/rtau
	loglik <- -log(sigma)-log(c)-log(rtau)-.5*log(z*z+1)-.5*log(2*pi)-.5*r*r
	if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
	ft
}

#----------------------------------------------------------------------------------------  
pjsu <- function(q, mu = 0, sigma = 1, nu = 1, tau = .5, lower.tail = TRUE, log.p = FALSE)
{  
	if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
	if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))           
	rtau <- 1/tau
	w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
	omega <- -nu*rtau
	c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
	z <- (q-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
	r <- -nu + asinh(z)/rtau    
	p <- pnorm(r,0,1)
	if(lower.tail==TRUE) p  <- p else  p <- 1-p 
	if(log.p==FALSE) p  <- p else  p <- log(p) 
	p
}

qjsu <-  function(p, mu=0, sigma=1, nu=0, tau=.5, lower.tail = TRUE, log.p = FALSE)
{   
	if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
	if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
	if (log.p==TRUE) p <- exp(p) else p <- p
	if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
	if (lower.tail==TRUE) p <- p else p <- 1-p
	rtau <- 1/tau
	r <- qnorm(p,0,1)
	z <- sinh(rtau*(r+nu))
	w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
	omega <- -nu*rtau 
	c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)     
	q <- (mu+c*sigma*w^(.5)*sinh(omega))+c*sigma*z   
	q
}

rjsu <- function(n, mu=0, sigma=1, nu=0, tau=.5)
{
	if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
	n <- ceiling(n)
	p <- runif(n)
	r <- qjsu(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
	r
}

# Distribution Model functions:
# Distribution Bounds
.DistributionBounds = function(distribution)
{
	ghlambda = 0
	ghlambda.LB = 0
	ghlambda.UB = 0
	if (distribution == "norm"){
		skew 	= 0
		skew.LB = 0
		skew.UB = 0
		shape 	= 0
		shape.LB = 0
		shape.UB = 0}
	if (distribution == "ged"){
		skew 	= 0
		skew.LB = 0
		skew.UB = 10
		shape 	= 2
		shape.LB = 0.01
		shape.UB = 50}
	if (distribution == "std"){
		skew 	= 0
		skew.LB = 0
		skew.UB = 0
		shape 	= 4
		shape.LB = 2.1
		shape.UB = 100}
	if (distribution == "snorm"){
		skew 	= 0.9
		skew.LB	= 0.1
		skew.UB	= 10
		shape 	= 0
		shape.LB = 0
		shape.UB = 0}
	if (distribution == "sged"){
		skew 	= 1
		skew.LB	= 0.01
		skew.UB	= 30
		shape 	= 2
		shape.LB = 0.1
		shape.UB = 60}
	if (distribution == "sstd"){
		skew 	= 1
		skew.LB = 0.01
		skew.UB = 30
		shape 	= 4
		shape.LB = 2.01
		shape.UB = 60}
	if (distribution == "nig"){
		skew 	= 0.2
		skew.LB = -0.99
		skew.UB	= 0.99
		shape 	= 0.4
		shape.LB = 0.01
		shape.UB = 15
		}
	if(distribution == "ghyp"){
		skew 	= 0.2
		skew.LB = -0.99
		skew.UB	= 0.99
		shape 	= 2
		shape.LB = 0.25
		shape.UB = 25
		ghlambda = -0.5
		ghlambda.LB = -6
		ghlambda.UB = 6
	}
	if(distribution == "jsu"){
		skew 	= 0
		skew.LB	= -20
		skew.UB	= 20
		shape 	= 1
		shape.LB = 0.1
		shape.UB = 10
	}
	# johnson has 2 shape parameters. The second one we model with the "skew"
	# representation in rugarch
	skewed.dists = c("snorm", "sged", "sstd", "nig", "ghyp", "jsu")
	shaped.dists = c("ged", "sged", "std", "sstd", "nig", "ghyp", "jsu")
	skew0  = 0
	shape0 = 0
	if(any(skewed.dists == distribution)) include.skew=TRUE else include.skew=FALSE
	if(any(shaped.dists == distribution)) include.shape=TRUE else include.shape=FALSE
	if(distribution == "ghyp") include.ghlambda = TRUE else include.ghlambda = FALSE
	return(list(shape = shape, shape.LB = shape.LB, shape.UB = shape.UB, skew = skew,
					skew.LB = skew.LB, skew.UB = skew.UB, include.skew = include.skew, 
					include.shape = include.shape, skew0 = skew0, shape0 = shape0,
					include.ghlambda = include.ghlambda, ghlambda = ghlambda, ghlambda.LB = ghlambda.LB,
					ghlambda.UB = ghlambda.UB))
}
# ------------------------------------------------------------------------------
.getDensity = function(cond.distribution = "norm")
{
	# Normal Distribution:
	if(cond.distribution == "norm") {
		.garchDensity = function(z, hh, lambda = 0, skew, shape) {
			dnorm(x = z, mean = 0, sd = 1)
		}
	}
	if(cond.distribution == "snorm") {
		.garchDensity = function(z, hh, lambda = 0, skew, shape) {
			dsnorm(x = z, mean = 0, sd = 1, xi = skew)
		}
	}
	
	# Standardized Student-t:
	if(cond.distribution == "std") {
		.garchDensity = function(z, hh, lambda = 0, skew, shape) {
			dstd(x = z, mean = 0, sd = 1, nu = shape)
		}
	}
	if(cond.distribution == "sstd") {
		.garchDensity = function(z, hh, lambda = 0, skew, shape) {
			dsstd(x = z, mean = 0, sd = 1, nu = shape, xi = skew)
		}
	}
	
	# Generalized Error Distribution:
	if(cond.distribution == "ged") {
		.garchDensity = function(z, hh, lambda = 0, skew, shape) {
			dged(x = z, mean = 0, sd = 1, nu = shape)
		}
	}
	if(cond.distribution == "sged") {
		.garchDensity = function(z, hh, lambda = 0, skew, shape) {
			dsged(x = z, mean = 0, sd = 1, nu = shape, xi = skew)
		}
	}
	if(cond.distribution == "nig") {
		.garchDensity = function(z, hh, lambda = 0, skew, shape) {
			dsnig(x = z, zeta = shape, rho = skew)
		}
	}
	
	if(cond.distribution == "ghyp") {
		.garchDensity = function(z, hh, lambda = 0, skew, shape) {
			dsgh(x = z, zeta = shape, rho = skew, lambda = lambda)
		}
	}
	
	if(cond.distribution == "jsu") {
		.garchDensity = function(z, hh, lambda = 0, skew, shape) {
			djsu(z, mu = 0 , sigma = 1, nu = skew, tau = shape)
		}
	}
	# Return Value:
	return(.garchDensity)
}

# ------------------------------------------------------------------------------
.getDistribution = function(cond.distribution = "norm")
{
	# Normal Distribution:
	if(cond.distribution == "norm") {
		.garchDist = function(z, hh, lambda = 0, skew, shape) {
			pnorm(q = z, mean = 0, sd = 1)
		}
	}
	if(cond.distribution == "snorm") {
		.garchDist = function(z, hh, lambda = 0, skew, shape) {
			psnorm(q = z, mean = 0, sd = 1, xi = skew)
		}
	}
	
	# Standardized Student-t:
	if(cond.distribution == "std") {
		.garchDist = function(z, hh, lambda = 0, skew, shape) {
			pstd(q = z, mean = 0, sd = 1, nu = shape)
		}
	}
	if(cond.distribution == "sstd") {
		.garchDist = function(z, hh, lambda = 0, skew, shape) {
			psstd(q = z, mean = 0, sd = 1, nu = shape, xi = skew)
		}
	}
	
	# Generalized Error Distribution:
	if(cond.distribution == "ged") {
		.garchDist = function(z, hh, lambda = 0, skew, shape) {
			pged(q = z, mean = 0, sd = 1, nu = shape)
		}
	}
	if(cond.distribution == "sged") {
		.garchDist = function(z, hh, lambda = 0, skew, shape) {
			psged(q = z, mean = 0, sd = 1, nu = shape, xi = skew)
		}
	}
	if(cond.distribution == "nig") {
		.garchDist = function(z, hh, lambda = 0, skew, shape) {
			psnig(q=z, zeta = shape, rho = skew)
		}
	}

	if(cond.distribution == "ghyp") {
		.garchDist = function(z, hh, lambda = 0, skew, shape) {
			psgh(q = z, zeta = shape, rho = skew, lambda = lambda)
		}
	}
		
	if(cond.distribution == "jsu") {
		.garchDist = function(z, hh, lambda = 0, skew, shape) {
			pjsu(z, mu = 0 , sigma = 1, nu = skew, tau = shape)
		}
	}
	
	# Return Value:
	return(.garchDist)
}

# ------------------------------------------------------------------------------
.getQuantile = function(cond.distribution = "norm")
{
	# Normal Distribution:
	if(cond.distribution == "norm") {
		.garchQuantile = function(z, hh, lambda = 0, skew, shape) {
			qnorm(p = z, mean = 0, sd = 1)
		}
	}
	if(cond.distribution == "snorm") {
		.garchQuantile = function(z, hh, lambda = 0, skew, shape) {
			qsnorm(p = z, mean = 0, sd = 1, xi = skew)
		}
	}
	
	# Standardized Student-t:
	if(cond.distribution == "std") {
		.garchQuantile = function(z, hh, lambda = 0, skew, shape) {
			qstd(p = z, mean = 0, sd = 1, nu = shape)
		}
	}
	if(cond.distribution == "sstd") {
		.garchQuantile = function(z, hh, lambda = 0, skew, shape) {
			qsstd(p = z, mean = 0, sd = 1, nu = shape, xi = skew)
		}
	}
	
	# Generalized Error Distribution:
	if(cond.distribution == "ged") {
		.garchQuantile = function(z, hh, lambda = 0, skew, shape) {
			qged(p = z, mean = 0, sd = 1, nu = shape)
		}
	}
	
	if(cond.distribution == "sged") {
		.garchQuantile = function(z, hh, lambda = 0, skew, shape) {
			qsged(p = z, mean = 0, sd = 1, nu = shape, xi = skew)
		}
	}
	
	if(cond.distribution == "nig") {
		.garchQuantile = function(z, hh, lambda = 0, skew, shape) {
			qsnig(p = z, zeta = shape, rho = skew)
		}
	}
	
	if(cond.distribution == "ghyp") {
		.garchQuantile = function(z, hh, lambda = 0, skew, shape) {
			qsgh(p = z, zeta = shape, rho = skew, lambda = lambda)
		}
	}
	
	if(cond.distribution == "jsu") {
		.garchQuantile = function(z, hh, lambda = 0, skew, shape) {
			qjsu(z, mu = 0 , sigma = 1, nu = skew, tau = shape)
		}
	}
	# Return Value:
	return(.garchQuantile)
}
# ------------------------------------------------------------------------------
.makeSample<-function(distribution, lambda = -0.5, skew, shape, n, seed)
{
	# Normal Distribution:
	set.seed(seed)
	if(distribution == "norm") {
			x = rnorm(n)
		}
	# Skew Normal Distribution
	if(distribution == "snorm") {
			x = rsnorm(n, xi = skew)
	}
	# Student-t:
	if(distribution == "std") {
			x = rstd(n, nu = shape)
	}
	# Skew Student-t:
	if(distribution == "sstd") {
			x = rsstd(n, nu = shape, xi = skew)
	}
	# Generalized Error Distribution:
	if(distribution == "ged") {
			x = rged(n, nu = shape)
	}
	# Skew Generalized Error Distribution:
	if(distribution == "sged") {
			x = rsged(n, nu = shape, xi = skew)
	}
	# Normal Inverse Gaussian Distribution
	if(distribution == "nig") {
			x = rsnig(n, zeta = shape, rho = skew)
	}
	
	if(distribution == "ghyp") {
			x = rsgh(n, zeta = shape, rho = skew, lambda = lambda)
	}
	
	if(distribution == "jsu"){
			x = rjsu(n, mu = 0, sigma = 1, nu = skew, tau = shape)
	}
	# Return Value:
	return(x)
}


# Scaling Transformation
.scaledist = function(dist, mu, sigma, lambda = -0.5, skew, shape)
{
	ans = switch(dist,
			norm  = .normscale(mu, sigma),
			snorm = .snormscale(mu, sigma, skew),
			std   = .stdscale(mu, sigma, shape),
			sstd  = .sstdscale(mu, sigma, skew, shape),
			nig   = .nigscale(mu, sigma, skew, shape),
			ghyp  = .ghypscale(mu, sigma, skew, shape, lambda),
			ged   = .gedscale(mu, sigma, shape),
			sged  = .sgedscale(mu, sigma, skew, shape),
			jsu   = .jsuscale(mu, sigma, skew, shape)
			)
	return(ans)
}

# returns the time varying density parameters of the actual returns
# (rescaled from the (0,1) parametrization)
.nigtransform = function(zeta, rho)
{
	nigpars = t(apply(cbind(zeta, rho), 1, FUN = function(x) .paramGH(zeta = x[1], rho = x[2], lambda = -0.5)))
	colnames(nigpars) = c("alpha", "beta", "delta", "mu")
	return(nigpars)
}

.ghyptransform = function(zeta, rho, lambda)
{
	n = length(zeta)
	ghyppars = t(apply(cbind(zeta, rho), 1, FUN = function(x) .paramGH(zeta = x[1], rho = x[2], lambda = lambda)))
	ghyppars = cbind(ghyppars, rep(lambda, n))
	colnames(ghyppars) = c("alpha", "beta", "delta", "mu", "lambda")
	return(ghyppars)
}

.nigscale2 = function(mu, sigma, nigalpha, nigbeta, nigdelta, nigmu)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	# alpha, beta, delta, mu [alpha = shape, beta = skew]
	xdensity[,1] = nigalpha/sigma
	xdensity[,2] = nigbeta/sigma
	xdensity[,3] = nigdelta*sigma
	xdensity[,4] = nigmu*sigma + mu
	colnames(xdensity) = c("alpha", "beta", "delta", "mu")
	return(xdensity)
}


.nigscale = function(mu, sigma, skew, shape)
{
	nigpars = t(apply(cbind(shape, skew), 1, FUN=function(x) .paramGH(zeta = x[1], rho = x[2], lambda = -0.5)))
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	# alpha, beta, delta, mu
	xdensity[,4] = nigpars[,1]/sigma
	xdensity[,3] = nigpars[,2]/sigma
	xdensity[,2] = nigpars[,3]*sigma
	xdensity[,1] = nigpars[,4]*sigma + mu
	colnames(xdensity) = c("mu", "sigma", "skew", "shape")
	return(xdensity)
}

.ghypscale = function(mu, sigma, skew, shape, lambda)
{
	if(length(lambda)>1){
		ghpars = t(apply(cbind(shape, skew, lambda), 1, FUN=function(x) .paramGH(zeta = x[1], rho = x[2], lambda = x[3])))
	} else{
		ghpars = t(apply(cbind(shape, skew), 1, FUN=function(x) .paramGH(zeta = x[1], rho = x[2], lambda = lambda)))
	}
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	# alpha, beta, delta, mu
	xdensity[,4] = ghpars[,1]/sigma
	xdensity[,3] = ghpars[,2]/sigma
	xdensity[,2] = ghpars[,3]*sigma
	xdensity[,1] = ghpars[,4]*sigma + mu
	colnames(xdensity) =c("mu", "sigma", "skew", "shape")
	return(xdensity)
}

.ghypscale2 = function(mu, sigma, ghalpha, ghbeta, ghdelta, ghmu)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	# alpha, beta, delta, mu
	xdensity[,1] = ghalpha/sigma
	xdensity[,2] = ghbeta/sigma
	xdensity[,3] = ghdelta*sigma
	xdensity[,4] = ghmu*sigma + mu
	colnames(xdensity) = c("alpha", "beta", "delta", "mu")
	return(xdensity)
}


.normscale = function(mu, sigma)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	xdensity[,1] = mu
	xdensity[,2] = sigma
	xdensity[,3] = 0
	xdensity[,4] = 0
	
	colnames(xdensity) = c("mu", "sigma", "skew", "shape")
	return(xdensity)
}

.snormscale = function(mu, sigma, skew)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	xdensity[,1] = mu
	xdensity[,2] = sigma
	xdensity[,3] = skew
	xdensity[,4] = 0
	colnames(xdensity) = c("mu", "sigma", "skew", "shape")
	return(xdensity)
}

.stdscale = function(mu, sigma, shape)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	xdensity[,1] = mu
	xdensity[,2] = sigma
	xdensity[,3] = 0
	xdensity[,4] = shape
	colnames(xdensity) = c("mu", "sigma", "skew", "shape")
	return(xdensity)
}

.sstdscale = function(mu, sigma, skew, shape)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	xdensity[,1] = mu
	xdensity[,2] = sigma
	xdensity[,3] = skew
	xdensity[,4] = shape
	colnames(xdensity) = c("mu", "sigma", "skew", "shape")
	return(xdensity)
}


.gedscale = function(mu, sigma, shape)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	xdensity[,1] = mu
	xdensity[,2] = sigma
	xdensity[,3] = 0
	xdensity[,4] = shape
	colnames(xdensity) = c("mu", "sigma", "skew", "shape")
	return(xdensity)
}

.sgedscale = function(mu, sigma, skew, shape)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	xdensity[,1] = mu
	xdensity[,2] = sigma
	xdensity[,3] = skew
	xdensity[,4] = shape
	colnames(xdensity) = c("mu", "sigma", "skew", "shape")
	return(xdensity)
}



.jsuscale = function(mu, sigma, skew, shape)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	xdensity[,1] = mu
	xdensity[,2] = sigma
	xdensity[,3] = skew
	xdensity[,4] = shape
	colnames(xdensity) = c("mu", "sigma", "skew", "shape")
	return(xdensity)
}


.qdensity = function(z, mu, sigma,  lambda = -0.5, skew, shape, distribution = "norm")
{
	if(distribution == "norm") {
		ans = apply(cbind(z, mu, sigma), 1, FUN=function(x) qnorm(p = x[1], 
							mean = x[2], sd = x[3]))
	}
	
	if(distribution == "snorm") {
		ans = apply(cbind(z, mu, sigma, skew), 1, FUN=function(x) qsnorm(p = x[1], 
							mean = x[2], sd = x[3], xi = x[4]))
	}
	
	if(distribution == "std") {
		ans = apply(cbind(z, mu, sigma, shape), 1, FUN=function(x) qstd(p = x[1], 
							mean = x[2], sd = x[3], nu = x[4]))
	}
	
	if(distribution == "sstd") {
		ans = apply(cbind(z, mu, sigma, skew, shape), 1, FUN=function(x) qsstd(p = x[1], 
							mean = x[2], sd = x[3], nu = x[5], xi = x[4]))
	}
	
	if(distribution == "ged") {
		ans = apply(cbind(z, mu, sigma, shape), 1, FUN=function(x) qged(p = x[1], 
							mean = x[2], sd = x[3], nu = x[4]))
	}
	
	if(distribution == "sged") {
		ans = apply(cbind(z, mu, sigma, skew, shape), 1, FUN=function(x) qsged(p = x[1], 
							mean = x[2], sd = x[3], nu = x[5], xi = x[4]))
	}
	
	if(distribution == "nig") {
		ans = apply(cbind(z, mu, sigma, shape, skew), 1, FUN=function(x) qnig(p = x[1], 
							alpha = x[4], beta = x[5], delta = x[3], mu = x[2]))
	}
	
	if(distribution == "ghyp") {
		if(length(lambda)>1){
			ans = apply(cbind(z, mu, sigma, shape, skew, lambda), 1, FUN=function(x) qgh(p = x[1],
								alpha = x[4], beta = x[5], delta = x[3], mu = x[2], lambda = x[6]))
		} else{
			ans = apply(cbind(z, mu, sigma, shape, skew), 1, FUN=function(x) qgh(p = x[1],
							alpha = x[4], beta = x[5], delta = x[3], mu = x[2], lambda = lambda))
		}
	}
	
	if(distribution == "jsu") {
		ans = apply(cbind(z, mu, sigma, skew, shape), 1, FUN=function(x) qjsu(p = x[1], 
							mu = x[2], sigma = x[3], nu = x[4], tau = x[5]))
	}
	return(ans)
}


.ddensity = function(z, mu, sigma,  lambda = -0.5, skew, shape, distribution = "norm")
{
	if(distribution == "norm") {
		ans = dnorm(x = z, mean = mu, sd = sigma)
	}
	
	if(distribution == "snorm") {
		ans = dsnorm(x = z, mean = mu, sd = sigma, xi = skew)
	}
	
	if(distribution == "std") {
		ans = dstd(x = z, mean = mu, sd = sigma, nu = shape)
	}
	
	if(distribution == "sstd") {
		ans = dsstd(x = z, mean = mu, sd = sigma, nu = shape, xi = skew)
	}
	
	if(distribution == "ged") {
		ans = dged(x = z, mean = mu, sd = sigma, nu = shape)
	}
	
	if(distribution == "sged") {
		ans = dsged(x = z, mean = mu, sd = sigma, nu = shape, xi = skew)
	}
	
	if(distribution == "nig") {
		ans  = dnig(x = z, alpha = shape, beta = skew, delta = sigma, mu = mu)
	}
	
	if(distribution == "ghyp") {
		ans  = dgh(x = z, alpha = shape, beta = skew, delta = sigma, mu = mu, lambda = lambda)
	}
	
	if(distribution == "jsu") {
		ans = djsu(y = z, mu = mu, sigma = sigma, nu = skew, tau = shape)
	}
	return(ans)
}

.dskewness = function(mu, sigma,  lambda = -0.5, skew, shape, distribution)
{
	if(distribution == "sstd"){
		eta = shape
		k2  = skew
		lda = ( k2-1 )/( k2+1 )
		ep1 = ( eta+1 )/2
		lnc = lgamma( ep1 ) - lgamma( eta/2 ) -0.5*log( pi*( eta-2 ) )
		cc  = exp(lnc)
		a   = 4*lda*cc*( eta-2 )/( eta-1 )
		b   = sqrt( 1+3*lda^2-a^2 )		
		my2 = 1+3*lda^2
		my3 = 16*cc*lda*( 1+lda^2 )*( ( eta-2 )^2 )/( ( eta-1 )*( eta-3 ) )
		ans = ( my3-3*a*my2+2*a^3 )/( b^3 )
	}
}
.dexkurtosis = function(mu, sigma, lambda = -0.5,  skew, shape, distribution)
{
	if(distribution == "norm") {
		ans = 0
	}
	if(distribution == "snorm"){

	}
	if(distribution == "std"){
		ans = 6/(shape-4)
	}
	if(distribution == "sstd"){
		eta = shape
		k2  = skew
		lda = ( k2-1 )/( k2+1 )
		ep1 = ( eta+1 )/2
		lnc = lgamma( ep1 ) - lgamma( eta/2 ) -0.5*log( pi*( eta-2 ) )
		cc  = exp(lnc)
		a   = 4*lda*cc*( eta-2 )/( eta-1 )
		b   = sqrt( 1+3*lda^2-a^2 )		
		my2 = 1+3*lda^2
		my3 = 16*cc*lda*( 1+lda^2 )*( (eta-2)^2 )/( ( eta-1 )*( eta-3 ) )
		my4 = 3*( eta-2 )*( 1+10*lda^2+5*lda^4 )/( eta-4 )
		m3  = ( my3-3*a*my2+2*a^3 )/( b^3 )
		ans = ( my4-4*a*my3+6*( a^2 )*my2-3*a^4 )/( b^4 ) - 3
	}
	if(distribution == "ged"){
		ans = gamma(5/shape)*gamma(1/shape)/gamma(3/shape)^2 - 3
	}
	return(ans)
}

#---------------------------------------------------------------------------------
# functions for export:
#---------------------------------------------------------------------------------
nigtransform = function(mu = 0, sigma = 1,  skew = 0, shape = 3)
{
	return(.scaledist(dist = "nig", mu = mu, sigma = sigma, skew = skew,
					shape = shape, lambda = -0.5))
}

ghyptransform = function(mu = 0, sigma = 1,  skew = 0, shape = 3, lambda = -0.5)
{
	return(.scaledist(dist = "ghyp", mu = mu, sigma = sigma, skew = skew,
					shape = shape, lambda = lambda))
}


fitdist = function(distribution = "norm", x, control=list()){
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig",
			"ghyp", "jsu")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distributions\n", call. = FALSE)
	ans = switch(distribution,
			norm = normFit(x, control),
			snorm = snormFit(x, control),
			std = stdFit(x, control),
			sstd = sstdFit(x, control),
			ged = gedFit(x, control),
			sged = sgedFit(x, control),
			nig = nigFit(x, control),
			ghyp = ghFit(x, control),
			jsu = jsuFit(x, control)
	)
	
}
ddist = function(distribution = "norm", y, mu = 0, sigma = 1, lambda = -0.5, skew = 1, shape = 5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig",
			"ghyp", "jsu")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distributions\n", call. = FALSE)
	if(distribution == "nig" | distribution == "ghyp"){
		# transform to alpha-beta
		pars = .paramGH(zeta = shape, rho = skew , lambda = lambda)
		# re-scale
		pars = pars * c(1/sigma, 1/sigma, sigma, sigma)
		pars[4] = pars[4] + mu
		# expand
		alpha = pars[1]
		beta = pars[2]
		delta = pars[3]
		xmu = pars[4]
	}
	ans = switch(distribution,
		norm = dnorm(y, mean = mu, sd = sigma, log = FALSE),
		snorm = dsnorm(y, mean = mu, sd = sigma, xi = skew),
		std = dstd(y, mean = mu, sd = sigma, nu = shape),
		sstd = dsstd(y, mean = mu, sd = sigma, nu = shape, xi = skew),
		ged = dged(y, mean = mu, sd = sigma, nu = shape),
		sged = dsged(y, mean = mu, sd = sigma, nu = shape, xi = skew),
		nig = dnig(y, alpha = alpha, beta = beta, delta = delta, mu = xmu),
		ghyp = dgh(y, alpha = alpha, beta = beta, delta = delta, mu = xmu, lambda = lambda),
		jsu = djsu(y, mu = mu, sigma = sigma, nu = skew, tau = shape)
		)
	return(ans)
}

pdist = function(distribution = "norm", q, mu = 0, sigma = 1, lambda = -0.5, skew = 1, shape = 5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig",
			"ghyp", "jsu")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distributions\n", call. = FALSE)
	if(distribution == "nig" | distribution == "ghyp"){
		# transform to alpha-beta
		pars = .paramGH(zeta = shape, rho = skew , lambda = lambda)
		# re-scale
		pars = pars * c(1/sigma, 1/sigma, sigma, sigma)
		pars[4] = pars[4] + mu
		# expand
		alpha = pars[1]
		beta = pars[2]
		delta = pars[3]
		xmu = pars[4]
	}
	ans = switch(distribution,
			norm = pnorm(q, mean = mu, sd = sigma, log.p = FALSE),
			snorm = psnorm(q, mean = mu, sd = sigma, xi = skew),
			std = pstd(q, mean = mu, sd = sigma, nu = shape),
			sstd = psstd(q, mean = mu, sd = sigma, nu = shape, xi = skew),
			ged = pged(q, mean = mu, sd = sigma, nu = shape),
			sged = psged(q, mean = mu, sd = sigma, nu = shape, xi = skew),
			nig = pnig(q, alpha = alpha, beta = beta, delta = delta, mu = xmu),
			ghyp = pgh(q, alpha = alpha, beta = beta, delta = delta, mu = xmu, lambda = lambda),
			jsu = pjsu(q, mu = mu, sigma = sigma, nu = skew, tau = shape)
	)
	return(ans)
}

qdist = function(distribution = "norm", p, mu = 0, sigma = 1, lambda = -0.5, skew = 1, shape = 5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig",
			"ghyp", "jsu")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distributions\n", call. = FALSE)
	if(distribution == "nig" | distribution == "ghyp"){
		# transform to alpha-beta
		pars = .paramGH(zeta = shape, rho = skew , lambda = lambda)
		# re-scale
		pars = pars * c(1/sigma, 1/sigma, sigma, sigma)
		pars[4] = pars[4] + mu
		# expand
		alpha = pars[1]
		beta = pars[2]
		delta = pars[3]
		xmu = pars[4]
	}
	ans = switch(distribution,
			norm = qnorm(p, mean = mu, sd = sigma, log.p = FALSE),
			snorm = qsnorm(p, mean = mu, sd = sigma, xi = skew),
			std = qstd(p, mean = mu, sd = sigma, nu = shape),
			sstd = qsstd(p, mean = mu, sd = sigma, nu = shape, xi = skew),
			ged = qged(p, mean = mu, sd = sigma, nu = shape),
			sged = qsged(p, mean = mu, sd = sigma, nu = shape, xi = skew),
			nig = qnig(p, alpha = alpha, beta = beta, delta = delta, mu = xmu),
			ghyp = qgh(p, alpha = alpha, beta = beta, delta = delta, mu = xmu, lambda = lambda),
			jsu = qjsu(p, mu = mu, sigma = sigma, nu = skew, tau = shape)
	)
	return(ans)
}

rdist = function(distribution = "norm", n, mu = 0, sigma = 1, lambda = -0.5, skew = 1, shape = 5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig",
			"ghyp", "jsu")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distribution\n", call. = FALSE)
	if(distribution == "nig" | distribution == "ghyp"){
		# transform to alpha-beta
		pars = .paramGH(zeta = shape, rho = skew , lambda = lambda)
		# re-scale
		pars = pars * c(1/sigma, 1/sigma, sigma, sigma)
		pars[4] = pars[4] + mu
		# expand
		alpha = pars[1]
		beta = pars[2]
		delta = pars[3]
		xmu = pars[4]
	}
	ans = switch(distribution,
			norm = rnorm(n, mean = mu, sd = sigma),
			snorm = rsnorm(n, mean = mu, sd = sigma, xi = skew),
			std = rstd(n, mean = mu, sd = sigma, nu = shape),
			sstd = rsstd(n, mean = mu, sd = sigma, nu = shape, xi = skew),
			ged = rged(n, mean = mu, sd = sigma, nu = shape),
			sged = rsged(n, mean = mu, sd = sigma, nu = shape, xi = skew),
			nig = rnig(n, alpha = alpha, beta = beta, delta = delta, mu = xmu),
			ghyp = rgh(n, alpha = alpha, beta = beta, delta = delta, mu = xmu, lambda = lambda),
			jsu = rjsu(n, mu = mu, sigma = sigma, nu = skew, tau = shape)
	)
	return(ans)
}


dskewness = function(distribution = "norm", skew = 1, shape = 5, lambda = -0.5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "jsu")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distribution\n", call. = FALSE)
	if( distribution == "nig" | distribution == "ghyp"){
		pars = .paramGH(zeta = shape, rho = skew, lambda = ifelse(distribution == "nig", -0.5, lambda))
	}
	ans = switch(distribution,
			norm 	= 0,
			snorm 	= .snormskew(skew = skew),
			std 	= 0,
			sstd 	= .sstdskew(skew = skew, shape = shape),
			ged 	= 0,
			sged 	= .sgedskew(skew = skew, shape = shape),
			nig 	= .nigskew(alpha = pars[1], beta = pars[2], delta = pars[3], mu = pars[4]),
			ghyp 	= .ghypskew(lambda = lambda, alpha = pars[1], beta = pars[2],  delta = pars[3], mu = pars[4]),
			jsu 	= .jsuskew(mu = 0, sigma = 1, skew = skew, shape = shape)
	)
	return(as.numeric(ans))
}

dkurtosis = function(distribution = "norm", skew = 1, shape = 5, lambda = -0.5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "jsu")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distribution\n", call. = FALSE)
	if( distribution == "nig" | distribution == "ghyp"){
		pars = .paramGH(zeta = shape, rho = skew, lambda = ifelse(distribution == "nig", -0.5, lambda))
	}
	ans = switch(distribution,
			norm 	= 0,
			snorm 	= 0,
			std 	= .stdexkurt(shape = shape),
			sstd 	= .sstdexkurt(skew = skew, shape = shape),
			ged 	= .gedexkurt(shape = shape),
			sged 	= .sgedexkurt(skew = skew, shape = shape),
			nig 	= .nigexkurt(alpha = pars[1], beta = pars[2], delta = pars[3], mu = pars[4]),
			ghyp 	= .ghypexkurt(lambda = lambda, alpha = pars[1], beta = pars[2],  delta = pars[3], mu = pars[4]),
			jsu 	= .jsuexkurt(mu = 0, sigma = 1, skew = skew, shape = shape)
	)
	return(as.numeric(ans))
}


# NIG Moments
.nigmu = function(alpha, beta, delta, mu){
	gm = sqrt(alpha^2 - beta^2)
	ans = mu + (delta*beta)/gm
	return(ans)
}
.nigsigma = function(alpha, beta, delta, mu){
	gm = sqrt(alpha^2 - beta^2)
	ans = sqrt((delta*alpha^2)/(gm^3))
	return(ans)
}
.nigskew = function(alpha, beta, delta, mu){
	gm = sqrt(alpha^2 - beta^2)
	ans = 3*beta/(alpha*sqrt(delta*gm))
	return(ans)
}

.nigexkurt = function(alpha, beta, delta, mu){
	gm = sqrt(alpha^2 - beta^2)
	ans = 3*(1+4*(beta^2)/(alpha^2))/(delta*gm)
	return(ans)
}

.ghypmu = function(lambda, alpha, beta, delta, mu){
	gm = sqrt(alpha^2 - beta^2)
	ans = mu + (delta * beta * besselK( delta * gm, lambda + 1) )/( gm * besselK( delta * gm, lambda) ) 
	return(ans)
}

.ghypsigma = function(lambda, alpha, beta, delta, mu){
	gm = sqrt(alpha^2 - beta^2)
	x1 = delta * besselK( delta * gm, lambda + 1) / ( gm * besselK( delta * gm, lambda) )
	x2 = ( (beta^2 * delta^2) / (gm^2) ) * ( ( besselK( delta * gm, lambda + 2) / besselK( delta * gm, lambda) ) 
				- ( besselK( delta * gm, lambda + 1)^2 / besselK( delta * gm, lambda)^2 ) )
	ans = sqrt(x1 + x2)
	return(ans)
}

.ghypskew = function(lambda, alpha, beta, delta, mu){
	skew = ghypMom(3, lambda, alpha, beta, delta, mu, momType = "central")/(.ghypsigma(lambda, alpha, beta, delta, mu)^3)
	return(skew)
}

.ghypexkurt = function(lambda, alpha, beta, delta, mu){
	kurt = ghypMom(4, lambda, alpha, beta, delta, mu, momType = "central")/(.ghypsigma(lambda, alpha, beta, delta, mu)^4) - 3
	return(kurt)
}

.norm2snorm1 = function(mu, sigma, skew)
{
	m1 = 2/sqrt(2 * pi)
	mu + m1 * (skew - 1/skew)*sigma
}

.norm2snorm2 = function(mu, sigma, skew)
{
	m1 = 2/sqrt(2 * pi)
	m2 = 1
	sigx = sqrt((1-m1^2)*(skew^2+1/skew^2) + 2*m1^2 - 1)
	
	sigx*sigma
}

.snormskew = function( skew )
{
	m1 = 2/sqrt(2 * pi)
	m2 = 1
	m3 = 4/sqrt(2 * pi)
	(skew - 1/skew) * ( ( m3 + 2 * m1^3 - 3 * m1 * m2 ) * ( skew^2 + (1/skew^2) ) + 3 * m1 * m2 - 4 * m1^3 )/
			( ( (m2 - m1^2) * (skew^2 + 1/skew^2) + 2 * m1^2 - m2) ^ (3/2) )
}

.snormexkurt = function( skew )
{
	0
}

.stdskew = function( shape )
{
	0
}

.stdexkurt = function( shape )
{
	ifelse(shape > 4, 6/(shape - 4), NA)
}

.sstdskew = function( skew, shape )
{
	# Standardize: (these are the absolute moments of dstd NOT dt)
	if(shape > 3 ){
		# Standardize: (these are the absolute moments of dstd NOT dt)
		m1 = 2 * sqrt(shape - 2) / (shape - 1) / beta(1/2, shape/2)
		m2 = 1
		m3 = integrate(f = function(x) 2 * x^3 * dstd(x, 0, 1, shape), 0, Inf, rel.tol=1e-12)$value
		(skew - 1/skew) * ( ( m3 + 2 * m1^3 - 3 * m1 * m2 ) * ( skew^2 + (1/skew^2) ) + 3 * m1 * m2 - 4 * m1^3 )/
				( ( (m2 - m1^2) * (skew^2 + 1/skew^2) + 2 * m1^2 - m2) ^ (3/2) )
	} else{
		NA
	}
}

.sstdexkurt = function( skew, shape )
{
	if(shape > 5 ){
		# Standardize: (these are the absolute moments of dstd NOT dt)
		m1 = 2 * sqrt(shape - 2) / (shape - 1) / beta(1/2, shape/2)
		m2 = 1
		m3 = integrate(f = function(x) 2 * x^3 * dstd(x, 0, 1, shape), 0, Inf, rel.tol=1e-12)$value
		m4 = integrate(f = function(x) 2 * x^4 * dstd(x, 0, 1, shape), 0, Inf, rel.tol=1e-12)$value
		cm4 = (-3 * m1^4 * (skew - 1/skew)^4) + 
				( 6 * m1^2 * (skew - 1/skew)^2 * m2*(skew^3 + 1/skew^3) )/(skew + 1/skew) - 
				( 4 * m1*(skew - 1/skew) * m3 * (skew^4 - 1/skew^4) )/(skew+1/skew) + 
				( m4 * (skew^5 + 1/skew^5) )/(skew + 1/skew)
		( cm4/( ( (m2 - m1^2) * (skew^2 + 1/skew^2) + 2 * m1^2 - m2) ^ 2) ) - 3
	} else{
		NA
	}
}

.gedskew = function( shape )
{
	0
}

.gedexkurt = function( shape )
{
	( ( ( gamma(1/shape)/gamma(3/shape) )^2 ) * ( gamma(5/shape)/gamma(1/shape) ) ) - 3
}

.sgedskew = function( skew, shape )
{
	lambda = sqrt ( 2^(-2/shape) * gamma(1/shape) / gamma(3/shape) )
	m1 = ((2^(1/shape)*lambda)^1 * gamma(2/shape) / gamma(1/shape))
	m2 = 1
	m3 = ((2^(1/shape)*lambda)^3 * gamma(4/shape) / gamma(1/shape))
	(skew - 1/skew) * ( ( m3 + 2 * m1^3 - 3 * m1 * m2 ) * ( skew^2 + (1/skew^2) ) + 3 * m1 * m2 - 4 * m1^3 )/
			( ( (m2 - m1^2) * (skew^2 + 1/skew^2) + 2 * m1^2 - m2) ^ (3/2) )
	
}

.sgedexkurt= function( skew, shape )
{
	lambda = sqrt ( 2^(-2/shape) * gamma(1/shape) / gamma(3/shape) )
	m1 = ((2^(1/shape)*lambda)^1 * gamma(2/shape) / gamma(1/shape))
	m2 = 1
	m3 = ((2^(1/shape)*lambda)^3 * gamma(4/shape) / gamma(1/shape))
	m4 = ((2^(1/shape)*lambda)^4 * gamma(5/shape) / gamma(1/shape))
	cm4 = (-3 * m1^4 * (skew - 1/skew)^4) + 
			( 6 * m1^2 * (skew - 1/skew)^2 * m2*(skew^3 + 1/skew^3) )/(skew + 1/skew) - 
			( 4 * m1*(skew - 1/skew) * m3 * (skew^4 - 1/skew^4) )/(skew+1/skew) + 
			( m4 * (skew^5 + 1/skew^5) )/(skew + 1/skew)
	( cm4/( ( (m2 - m1^2) * (skew^2 + 1/skew^2) + 2 * m1^2 - m2) ^ 2 ) ) - 3
}

.jsuskew = function( mu = 0, sigma = 1, skew, shape )
{	
	f = function(x) ( (x-mu)^3 )*djsu(x, mu = mu, sigma = sigma, nu = skew, tau = shape, log = FALSE)
	integrate(f, -Inf, Inf)$value/sigma^3
}

.jsuexkurt = function( mu = 0, sigma = 1, skew, shape )
{
	f = function(x) ( (x-mu)^4 )*djsu(x, mu = mu, sigma = sigma, nu = skew, tau = shape, log = FALSE)
	(integrate(f, -Inf, Inf)$value/sigma^4)-3
}
