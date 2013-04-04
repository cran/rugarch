#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2013.
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
# ------------------------------------------------------------------------------
# Skew Generalized Hyberolic Student's T 
# alpha = abs(beta)+1e-12, lambda = -nu/2
# ------------------------------------------------------------------------------

# Location-Scale Invariant Parametrization
.paramGHST = function(betabar, nu){
	# Alexios Ghalanos 2012
	# betabar is the skew parameter = beta*delta (parametrization 4 in Prause)
	# nu is the shape parameter
	delta = ( ((2 * betabar^2)/((nu-2)*(nu-2)*(nu-4))) + (1/(nu-2)) )^(-0.5)
	beta = betabar/delta
	mu = -( (beta * (delta^2))/(nu-2))
	return( c(mu, delta, beta, nu) )
}

dsghst = function(x, mean=0, sd=1, skew=1, shape=8, log = FALSE){
	f = Vectorize( .dsghst )
	ans =  f(x, mean, sd, skew, shape, log)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
}

.dsghst = function(x, mean=0, sd=1, skew=1, shape=8, log = FALSE){
	if(abs(skew)<1e-12) skew = 1e-12
	z = (x - mean)/sd
	params = .paramGHST(skew, shape)
	beta = params[3]
	delta = params[2]
	mu = params[1]
	nu = shape
	f = 2^(0.5*(1-nu))*(delta^nu)*abs(beta)^(0.5*(nu+1))*besselK(sqrt(beta^2*(delta^2+(z-mu)^2)), 0.5*(nu+1))*exp(beta*(z-mu))
	f = f/(gamma(nu/2)*sqrt(pi)*(sqrt(delta^2+(z-mu)^2)^(0.5*(nu+1))))
	f = f/sd
	if(log) f = log(f)
	return( f )
}

rsghst = function(n, mean=0, sd=1, skew=1, shape=8){
	if(length(n)>1) stop("\nn cannot be a vector of length greater than 1!")
	f = Vectorize( .rsghst , c("mean","sd","skew","shape"))
	ans = f(n, mean, sd, skew, shape)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
}

.rsghst = function(n, mean=0, sd=1, skew=1, shape=8){
	params = .paramGHST(skew, shape)
	beta = params[3]
	delta = params[2]
	mu = params[1]
	nu = shape
	y <- 1/rgamma(n, shape = nu/2, scale = 2/(delta*sd)^2)
	sigma <- sqrt(y)
	z <- rnorm(n)
	f <- mean+(mu*sd) + beta/sd * sigma^2 + sigma * z
	return( f )
}

####################################################
# functions adapted/imported from the SkewHyperbolic of Scott

.skewhypStepSize = function(dist, delta, beta, nu, side = c("right", "left")) 
{
	side <- match.arg(side)
	if (beta > 0) {
		step = ifelse(side == "left", delta, delta * abs(beta) * 
						(nu * dist)^(-2/nu))
	}
	if (beta < 0) {
		step = ifelse(side == "right", delta, delta * abs(beta) * 
						(nu * dist)^(-2/nu))
	}
	if (isTRUE(all.equal(beta, 0))) {
		step = exp(dist/nu)
	}
	return( step )
}

modeghst = function(mean = 0, sd = 1, skew = 1, shape = 8){
	modeFun <- function(x) {
		dsghst(x, mean, sd, skew, shape, log = TRUE)
	}
	start <- 0
	options(warn = -1)
	opt <- optim(start, modeFun, control = list(fnscale = -1, 
					maxit = 1000, method = "BFGS"))
	ifelse(opt$convergence == 0, distMode <- opt$par, distMode <- NA)
	return( distMode )
}

psghst = function(q, mean=0, sd=1, skew=1, shape=8, lower.tail = TRUE, log = FALSE, cluster = NULL, ...){
	if(!is.null(cluster)){
		if(length(q)>1) stop("\ncluster evaluation cannot have length(q)>1")
		parallel::clusterExport(cluster, c(".paramGHST", "modeghst", 
						"dsghst", ".dsghst", ".skewhypStepSize"), envir = environment())
		ans = parallel::clusterMap(cluster, .psghst, q = q, mean=mean, sd=sd, 
				skew=skew, shape=shape, lower.tail = lower.tail, log = log, 
				SIMPLIFY = TRUE, USE.NAMES = TRUE, 
				.scheduling = c("static", "dynamic")[2])
	} else{
		f = Vectorize( .psghst , c("q", "mean","sd","skew","shape", "lower.tail", "log"))
		ans =  f(q, mean, sd, skew, shape, lower.tail, log, ...)
	}
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
}

.psghst = function(q, mean=0, sd=1, skew=1, shape=8, lower.tail = TRUE, log = FALSE, ...){
	distMode = modeghst(mean, sd, skew, shape)
	qLess = which((q <= distMode) & (is.finite(q)))
	qGreater = which((q > distMode) & (is.finite(q)))
	prob = rep(NA, length(q))
	err = rep(NA, length(q))
	prob[q == -Inf] = 0
	prob[q == Inf]  = 0
	err[q %in% c(-Inf, Inf)] = 0
	dskewhypInt = function(q) {
		dsghst(q, mean, sd, skew, shape)
	}
	for (i in qLess) {
		intRes <- integrate(dskewhypInt, -Inf, q[i], ...)
		prob[i] = intRes$value
		err[i]  = intRes$abs.error
	}
	for (i in qGreater) {
		intRes = integrate(dskewhypInt, q[i], Inf, ...)
		prob[i] = intRes$value
		err[i]  = intRes$abs.error
	}
	if (lower.tail == TRUE) {
		if (length(q > distMode) > 0) {
			prob[q > distMode] = 1 - prob[q > distMode]
		}
	}
	else {
		if (length(q <= distMode) > 0) {
			prob[q <= distMode] = 1 - prob[q <= distMode]
		}
	}
	if(log) prob = log(prob)
	return( prob )
}

qsghst = function(p, mean=0, sd=1, skew=1, shape=8, lower.tail = TRUE, cluster = NULL, ...){
	if(!is.null(cluster)){
		if(length(p)>1) stop("\ncluster evaluation cannot have length(p)>1")
		parallel::clusterExport(cluster, c(".paramGHST", "modeghst", "psghst", 
						".psghst", "dsghst", ".dsghst", ".skewhypStepSize"), envir = environment())
		ans = parallel::clusterMap(cluster, .qsghst, p = p, mean=mean, sd=sd, skew=skew, shape=shape, 
				lower.tail = lower.tail, SIMPLIFY = TRUE, USE.NAMES = TRUE, 
				.scheduling = c("static", "dynamic")[2])
	} else{
		f = Vectorize( .qsghst , c("p", "mean","sd","skew","shape", "lower.tail") )
		ans =  f(p, mean, sd, skew, shape, lower.tail)
	}
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
}


.qsghst = function(p, mean=0, sd=1, skew=1, shape=8, lower.tail = TRUE, ...){
	if (!lower.tail) {
		p <- 1 - p
		lower.tail == TRUE
	}
	params = .paramGHST(skew, shape)
	# scale the parameters
	beta = params[3]/sd
	delta = params[2]*sd
	mu = params[1]*sd + mean
	distMode = modeghst(mean, sd, skew, shape)
	pModeDist <- psghst(distMode, mean, sd, skew, shape)
	ans = rep(NA, length(p))
	invalid = which(p < 0 | p > 1)
	pFinite = which((p > 0) & (p < 1))
	ans[p == 0] = -Inf
	ans[p == 1] =  Inf
	less <- which((p <= pModeDist) & (p > 0))
	if (length(less) > 0) {
		pLow = min(p[less])
		xLow = distMode - .skewhypStepSize(delta, delta, beta, shape, "left")
		while(psghst(xLow, mean, sd, skew, shape, ...) >= pLow){
			xLow = xLow - .skewhypStepSize(distMode - xLow, delta, beta, shape, "left")
		}
		xRange = c(xLow, distMode)
		zeroFn = function(x, mean, sd, skew, shape, p, ...){
			return(psghst(x, mean, sd, skew, shape, ...) - p)
		}
		for(i in less){
			ans[i] = uniroot(zeroFn, mean = mean, sd = sd, skew = skew, 
					shape = shape, p = p[i], interval = xRange, ...)$root
		}
	}
	greater = which((p > pModeDist) & (p < 1))
	p[greater] = 1 - p[greater]
	if(length(greater) > 0) {
		pHigh = min(p[greater])
		xHigh = distMode + .skewhypStepSize(delta, delta, beta, shape, "right")
		while(psghst(xHigh, mean, sd, skew, shape, lower.tail = FALSE, ...) >= pHigh){
			xHigh <- xHigh + .skewhypStepSize(xHigh - distMode, delta, beta, shape, "right")
		}
		xRange = c(distMode, xHigh)
		zeroFn = function(x, mean, sd, skew, shape, p, ...){
			return(psghst(x, mean, sd, skew, shape, lower.tail = FALSE, ...) - p)
		}
		for(i in greater){
			ans[i] = uniroot(zeroFn, mean = mean, sd = sd, skew = skew, 
					shape = shape, p = p[i], interval = xRange, ...)$root
		}
	}
	return( ans )
}
####################################################


ghstFit = function(x, control){
	
	f = function(pars, x){
		return(sum(-dsghst(x, mean=pars[1], sd=pars[2], skew=pars[3], shape=pars[4], log = TRUE)))
	}
	x = as.numeric(x)
	x0 = c(mean(x), sd(x), 0.2, 8)
	fit = try(solnp(x0, fun = f, LB = c(-5, 1e-10, -10, 4.001), UB = c(5, 10, 10, 25), control = control, x = x), 
			silent = TRUE)
	
	# Add Names to $par
	names(fit$par) = c("mean", "sd", "skew", "shape")
	
	# Return Value:
	return( fit )

}
# equivalence: 
# 1. dskewhyp(x, mu = (sigma*params[1]+0.5), delta = params[2]*sigma, beta = params[3]/sigma, nu = shape)
# 2. dskewhyp(z, mu = params[1], delta = params[2], beta = params[3], nu = shape)/sigma

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


# ------------------------------------------------------------------------------

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
	return( Density * sigma )
}
# ------------------------------------------------------------------------------
# vectorized: yes
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
	fun = function(x, mean, sd, xi){
		.dsnorm(x = (x-mean)/sd, xi = xi) / sd
	}
	# Shift and Scale:
	f = Vectorize(fun)
	result = f(x, mean, sd, xi)
	# Return Value:
	return( result )
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
	return( Probability )
}
# ------------------------------------------------------------------------------
# vectorized: yes
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
	fun = function(q, mean, sd, xi){
		.psnorm(q = (q-mean)/sd, xi = xi)
	}
	# Shift and Scale:
	f = Vectorize(fun)
	result = f(q, mean, sd, xi)
	# Return Value:
	return( result )
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
	return( Quantile )
}
# ------------------------------------------------------------------------------
# vectorized: yes
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
	fun = function(p, mean, sd, xi){
		.qsnorm(p = p, xi = xi) * sd + mean
	}
	# Shift and Scale:
	f = Vectorize(fun)
	result = f(p, mean, sd, xi)	
	# Return Value:
	return( result )
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
	return( Random )
}
# ------------------------------------------------------------------------------
# vectorized: yes
rsnorm <-function(n, mean = 0, sd = 1, xi = 1.5)
{   
	# A function implemented by Diethelm Wuertz 
	### modified by Alexios Ghalanos to vectorized form
	# Description:
	#   Generate random deviates from the 
	#   skew normal distribution
	# Arguments:
	#   n - an integer value giving the number of observation.
	#   mean, sd, xi - location parameter, scale parameter, and 
	#       skewness parameter.
	fun = function(n, mean, sd, xi) .rsnorm(n = n, xi = xi) * sd + mean
	if(length(n)>1) stop("\nn cannot be a vector of length greater than 1!")
	f = Vectorize( fun , c("mean","sd","xi"))
	ans = f(n, mean, sd, xi)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
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
	return( fit )
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
	return( fit )
}
################################################################################
# vectorized: yes
dged <-function(x, mean = 0, sd = 1, nu = 2)
{   
	# A function imlemented by Diethelm Wuertz
	
	# Description:
	#   Compute the density for the 
	#   generalized error distribution.
	
	# FUNCTION:
	
	# Compute Density:
	fun = function(x, mean, sd, nu){
		z = (x - mean ) / sd
		lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
		g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
		g * exp (-0.5*(abs(z/lambda))^nu) / sd
	}
	# Shift and Scale:
	f = Vectorize(fun)
	result = f(x, mean, sd, nu)
	# Return Value
	return( result )
}
# ------------------------------------------------------------------------------
# vectorized: yes
pged <-function(q, mean = 0, sd = 1, nu = 2)
{   
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Compute the probability for the  
	#   generalized error distribution.
	
	# FUNCTION:
	
	# Compute Probability:
	fun = function(q, mean, sd, nu){
		q = (q - mean ) / sd
		lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
		g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
		h = 2^(1/nu) * lambda * g * gamma(1/nu) / nu
		s = 0.5 * ( abs(q) / lambda )^nu
		result = 0.5 + sign(q) * h * pgamma(s, 1/nu)
	}
	# Shift and Scale:
	f = Vectorize(fun)
	result = f(q, mean, sd, nu)
	# Return Value:
	return( result )
}


# ------------------------------------------------------------------------------
# vectorized: yes
qged <-function(p, mean = 0, sd = 1, nu = 2)
{   
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Compute the quantiles for the  
	#   generalized error distribution.
	
	# FUNCTION:
	
	# Compute Quantiles:
	fun = function(p, mean, sd, nu){
		lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
		q = lambda * (2*qgamma((abs(2*p-1)), 1/nu))^(1/nu)
		q*sign(2*p-1) * sd + mean
	}
	# Shift and Scale:
	f = Vectorize(fun)
	result = f(p, mean, sd, nu)
	# Return Value:
	return( result )
}


# ------------------------------------------------------------------------------
# vectorized: yes
rged <- function(n, mean = 0, sd = 1, nu = 2)
{   
	# A function implemented by Diethelm Wuertz
	### modified by Alexios Ghalanos to vectorized form
	
	# Description:
	#   Generate GED random deviates. The function uses the 
	#   method based on the transformation of a Gamma random 
	#   variable.
	
	# FUNCTION:
	
	# Generate Random Deviates:
	fun = function(n, mean, sd, nu){
		lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
		r = rgamma(n, 1/nu)
		z =  lambda * (2*r)^(1/nu) * sign(runif(n)-1/2)
		return( z * sd + mean )
	}
	if(length(n)>1) stop("\nn cannot be a vector of length greater than 1!")
	f = Vectorize( fun , c("mean","sd","nu"))
	ans = f(n, mean, sd, nu)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )	
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
	return( Density * sigma )
}

# norm [ nu = 2, xi = 1 ]
# laplace [ nu = 1, xi = 1 ]
# vectorized: yes
dsged <-function(x, mean = 0, sd = 1, nu = 2, xi = 1.5)
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Compute the density function of the 
	#   skewed generalized error distribution
	
	# FUNCTION:
	
	# Shift and Scale:
	fun = function(x, mean, sd, nu, xi){
		.dsged(x = (x-mean)/sd, nu = nu, xi = xi) / sd
	}
	# Shift and Scale:
	f = Vectorize(fun)
	result = f(x, mean, sd, nu, xi)	
	# Return Value:
	return( result )
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
	return( Probability )
}
# vectorized: yes
psged <-function(q, mean = 0, sd = 1, nu = 2, xi = 1.5)
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Compute the distribution function of the 
	#   skewed generalized error distribution
	
	# FUNCTION:
	
	# Shift and Scale:
	fun = function(q, mean, sd, nu, xi){
		.psged(q = (q-mean)/sd, nu = nu, xi = xi)
	}
	# Shift and Scale:
	f = Vectorize(fun)
	result = f(q, mean, sd, nu, xi)
	# Return Value:
	return( result )
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
	return( Quantile )
}
# vectorized: yes
qsged <-function(p, mean = 0, sd = 1, nu = 2, xi = 1.5)
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Compute the quantile function of the 
	#   skewed generalized error distribution
	
	# FUNCTION:
	
	# Shift and Scale:
	fun = function(p, mean, sd, nu, xi){
		.qsged(p = p, nu = nu, xi = xi) * sd + mean
	}
	# Shift and Scale:
	f = Vectorize(fun)
	result = f(p, mean, sd, nu, xi)	
	# Return Value:
	return( result )
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
	return( Random )
}


# ------------------------------------------------------------------------------
# vectorized: yes
rsged <-function(n, mean = 0, sd = 1, nu = 2, xi = 1.5)
{   
	# A function implemented by Diethelm Wuertz 
	
	# Description:
	#   Generate random deviates from the 
	#   skewed generalized error distribution
	
	# FUNCTION:
	
	# Shift and Scale:
	fun = function(n, mean, sd, nu, xi){
		return( .rsged(n = n, nu = nu, xi = xi) * sd + mean )
	}
	if(length(n)>1) stop("\nn cannot be a vector of length greater than 1!")
	f = Vectorize( fun , c("mean","sd","nu","xi"))
	ans = f(n, mean, sd, nu, xi)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
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
	return( fit )
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
	return( fit )
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
	return( result )
}
# ------------------------------------------------------------------------------
# vectorized: yes
dstd <- function(x, mean = 0, sd = 1, nu = 5)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Compute the density for the
	#   Student-t distribution.
	
	# FUNCTION:
	
	# Compute Density:
	fun = function(x, mean, sd, nu){
		s = sqrt(nu/(nu-2))
		z = (x - mean) / sd
		# result = .Internal(dnt(x = z*s, df = nu, ncp = 0, log = FALSE)) / (sd/s)
		dt(x = z*s, df = nu) * s / sd
	}
	# Shift and Scale:
	f = Vectorize(fun)
	result = f(x, mean, sd, nu)		
	# Return Value:
	return( result )
}


# ------------------------------------------------------------------------------
# vectorized: yes
pstd <-function (q, mean = 0, sd = 1, nu = 5)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Compute the probability for the
	#   Student-t distribution.
	
	# FUNCTION:
	
	# Compute Probability:
	fun = function(q, mean, sd, nu){
		s = sqrt(nu/(nu-2))
		z = (q - mean) / sd
		pt(q = z*s, df = nu)
	}
	# Shift and Scale:
	f = Vectorize(fun)
	result = f(q, mean, sd, nu)		
	# Return Value:
	return( result )
}


# ------------------------------------------------------------------------------
# vectorized: yes
qstd <-function (p, mean = 0, sd = 1, nu = 5)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Compute the quantiles for the
	#   Student-t distribution.
	
	# FUNCTION:
	
	# Compute Quantiles:
	fun = function(p, mean, sd, nu){
		s = sqrt(nu/(nu-2))
		qt(p = p, df = nu) * sd / s + mean
	}
	# Shift and Scale:
	f = Vectorize(fun)
	result = f(p, mean, sd, nu)
	# Return Value:
	return( result )
}


# ------------------------------------------------------------------------------
# vectorized: yes
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
	
	# AG: removed SPlus compatibility
	# For SPlus compatibility:
	# if (!exists("beta"))
	#	beta <- function (a, b) exp( lgamma(a) + lgamma(b) -lgamma(a+b) )
	
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
	return( Density * sigma )
}


# ------------------------------------------------------------------------------
# vectorized: yes
dsstd <-function(x, mean = 0, sd = 1, nu = 5, xi = 1.5)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Compute the density function of the
	#   skewed Student-t distribution
	
	# FUNCTION:
	fun = function(x, mean, sd, nu, xi){
		.dsstd(x = (x-mean)/sd, nu = nu, xi = xi) / sd
	}
	# Shift and Scale:
	f = Vectorize(fun)
	result = f(x, mean, sd, nu, xi)
	
	# Return Value:
	return( result )
}


# ------------------------------------------------------------------------------
.psstd <-function(q, nu, xi)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Internal Function
	
	# FUNCTION:
	# AG: removed SPlus compatibility
	# For SPlus compatibility:
	# if (!exists("beta"))
	#	beta <- function (a, b) exp( lgamma(a) + lgamma(b) -lgamma(a+b) )
	
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
	return( Probability )
}


# ------------------------------------------------------------------------------
# vectorized: yes
psstd <-function(q, mean = 0, sd = 1, nu = 5, xi = 1.5)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Compute the distribution function of the
	#   skewed Student-t distribution
	
	# FUNCTION:
	
	# Shift and Scale:
	fun = function(q, mean, sd, nu, xi){
		.psstd(q = (q-mean)/sd, nu = nu, xi = xi)
	}
	# Shift and Scale:
	f = Vectorize(fun)
	result = f(q, mean, sd, nu, xi)
	
	# Return Value:
	return( result )
}
# ------------------------------------------------------------------------------
.qsstd <-function(p, nu, xi)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Internal Function
	
	# FUNCTION:
	# AG: removed SPlus compatibility
	# For SPlus compatibility:
	#if (!exists("beta"))
	#	beta <- function (a, b) exp( lgamma(a) + lgamma(b) -lgamma(a+b) )
	
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
	return( Quantile )
}
# ------------------------------------------------------------------------------
# vectorized: yes
qsstd <-function(p, mean = 0, sd = 1, nu = 5, xi = 1.5)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Compute the quantile function of the
	#   skewed Student-t distribution
	
	# FUNCTION:
	
	# Shift and Scale:
	fun = function(p, mean, sd, nu, xi){
		.qsstd(p = p, nu = nu, xi = xi) * sd + mean
	}
	# Shift and Scale:
	f = Vectorize(fun)
	result = f(p, mean, sd, nu, xi)	
	# Return Value:
	return( result )
}


# ------------------------------------------------------------------------------
.rsstd <-function(n, nu, xi)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Internal Function
	
	# FUNCTION:
	# AG: removed SPlus compatibility
	# For SPlus compatibility:
	#if (!exists("beta"))
	#	beta <- function (a, b) exp( lgamma(a) + lgamma(b) -lgamma(a+b) )
	
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
	return( Random )
}


# ------------------------------------------------------------------------------
# vectorized: yes

rsstd <-function(n, mean = 0, sd = 1, nu = 5, xi = 1.5)
{
	# A function implemented by Diethelm Wuertz
	
	# Description:
	#   Generate random deviates from the
	#   skewed Student-t distribution
	
	# FUNCTION:
	
	# Shift and Scale:
	fun = function(n, mean, sd, nu, xi){
		return(.rsstd(n = n, nu = nu, xi = xi) * sd + mean)
	}
	if(length(n)>1) stop("\nn cannot be a vector of length greater than 1!")
	f = Vectorize( fun , c("mean","sd","nu", "xi"))
	ans = f(n, mean, sd, nu, xi)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )	
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
	return( fit )
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
	return( fit )
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
	return( ans$root )
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
		kappa = (besselK(x, lambda+1, expon.scaled = TRUE)/besselK(x, lambda, expon.scaled = TRUE) ) / x
	}
	# Return Value:
	return( kappa )
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
	return( deltaKappa )
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
	return( c(alpha = alpha, beta = beta, delta = delta, mu = mu) )
}

# vectorized: yes
dgh = function(x, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1, log = FALSE)
{
	f = Vectorize( .dgh  )
	ans =  f(x, alpha, beta, delta, mu, lambda, log)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
}

.dgh<-function(x, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1, log = FALSE)
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
	return( ans )
}
# ------------------------------------------------------------------------------
# vectorized: yes
pgh = function(q, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1, cluster = NULL)
{
	if(!is.null(cluster)){
		if(length(q)>1) stop("\ncluster evaluation cannot have length(q)>1")
		parallel::clusterExport(cluster, c(".dgh"), envir = environment())
		ans = parallel::clusterMap(cluster, .pgh, q = q, alpha=alpha, beta=beta, 
				delta=delta, mu=mu, lambda = lambda, 
				SIMPLIFY = TRUE, USE.NAMES = TRUE, 
				.scheduling = c("static", "dynamic")[2])
	} else{
		f = Vectorize( .pgh , c("q", "alpha","beta","delta","mu", "lambda") )
		ans =  f(q, alpha, beta, delta, mu, lambda)
	}
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )	
}
.pgh<-function(q, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1)
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
		Integral = integrate(.dgh, -Inf, Q, stop.on.error = FALSE,
				alpha = alpha, beta = beta, delta = delta, mu = mu,
				lambda = lambda)
		ans = c(ans, as.numeric(unlist(Integral)[1]))
	}
	
	# Return Value:
	return( ans )
}
# ------------------------------------------------------------------------------
# vectorized: yes

qgh = function(p, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1, cluster = NULL){
	if(!is.null(cluster)){
		if(length(p)>1) stop("\ncluster evaluation cannot have length(p)>1")
		parallel::clusterExport(cluster, c(".pgh", ".dgh", ".unitroot"), envir = environment())
		ans = parallel::clusterMap(cluster, .qgh, p = p, alpha=alpha, beta=beta, 
				delta=delta, mu=mu, lambda = lambda, 
				SIMPLIFY = TRUE, USE.NAMES = TRUE, 
				.scheduling = c("static", "dynamic")[2])
	} else{
		f = Vectorize( .qgh , c("p", "alpha","beta","delta","mu", "lambda") )
		ans =  f(p, alpha, beta, delta, mu, lambda)
	}
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
}

.qgh<-function (p, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1)
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
	.froot <- function(x, alpha, beta, delta, mu, lambda, p)
	{
		.pgh(q = x, alpha = alpha, beta = beta, delta = delta,
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
	return( result )
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
	return( fit )
	
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
	return( Y )
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
	return( output/alpha )
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
	return( output/alpha )
}



# ------------------------------------------------------------------------------
# vectorized: yes
rgh = function (n, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1)
{
	# A function implemented by Diethelm Wuertz
	# modified by Alexios Ghalanos (vectorized)
	
	# Description:
	#   Returns random variates for the generalized hyperbolic distribution
	
	# FUNCTION:
	
	# Checks:
	if (alpha <= 0) stop("alpha must be greater than zero")
	if (delta <= 0) stop("delta must be greater than zero")
	if (abs(beta) >= alpha) stop("abs value of beta must be less than alpha")
	if(length(n)>1) stop("\nn cannot be a vector of length greater than 1!")
	
	# Settings:
	fun = function(n, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1){
		theta = c(lambda, alpha, beta, delta, mu)
		return( .rghyp(n, theta) )
	}
	f = Vectorize( fun , c("alpha","beta","delta", "mu", "lambda"))
	ans = f(n, alpha, beta, delta, mu, lambda)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )	
}

# ------------------------------------------------------------------------------
# vectorized: yes
dsgh = function(x, zeta = 1, rho = 0, lambda = 1, log = FALSE) 
{
	fun = function(x, zeta = 1, rho = 0, lambda = 1, log = FALSE){
		param = .paramGH(zeta, rho, lambda)
		return( .dgh(x, param[1], param[2], param[3], param[4], lambda, log) )
	}
	f = Vectorize( fun )
	ans = f(x, zeta, rho, lambda, log)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )	
}
# ------------------------------------------------------------------------------
# vectorized: yes
psgh = function(q, zeta = 1, rho = 0, lambda = 1) 
{
	fun = function(q, zeta = 1, rho = 0, lambda = 1){
		param = .paramGH(zeta, rho, lambda)
		return( .pgh(q, param[1], param[2], param[3], param[4], lambda) )
	}
	f = Vectorize( fun  )
	ans = f(q, zeta, rho, lambda)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
}
# ------------------------------------------------------------------------------
# vectorized: yes
qsgh = function(p, zeta = 1, rho = 0, lambda = 1) 
{
	# A function implemented by Diethelm Wuertz
	# modified by Alexios Ghalanos (vectorized)
	# Description:
	#   Returns quantiles of the sgh distribution
	
	# FUNCTION:
	
	# Compute Quantiles:
	fun = function(p, zeta = 1, rho = 0, lambda = 1){
		param = .paramGH(zeta, rho, lambda)
		return( .qgh(p, param[1], param[2], param[3], param[4], lambda) )
	}
	f = Vectorize( fun )
	ans = f(p, zeta, rho, lambda)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
}
# ------------------------------------------------------------------------------
# vectorized: yes
rsgh = function(n, zeta = 1, rho = 0, lambda = 1) 
{
	# A function implemented by Diethelm Wuertz
	# modified by Alexios Ghalanos (vectorized)
	# Description:
	#   Generates sgh distributed random variates
	
	# FUNCTION:
	
	# Generate Random Numbers:
	
	if(length(n)>1) stop("\nn cannot be a vector of length greater than 1!")
	
	# Settings:
	fun = function(n,  zeta = 1, rho = 0, lambda = 1){
		param = .paramGH(zeta, rho, lambda)
		theta = c(lambda,  param[1], param[2], param[3], param[4])
		return( .rghyp(n, theta) )		
	}
	f = Vectorize( fun , c("zeta","rho","lambda"))
	ans = f(n, zeta, rho, lambda)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )	
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
		esghmle = function(x, y = x, ghlambda) {
			f = -sum(log(dsgh(y, x[1], x[2], ghlambda, log = FALSE)))
			f
		}
		# LLH Optimization:
		fit = solnp(
				pars = c(zeta, rho), 
				fun = esghmle, 
				LB = c(eps, -0.9999), 
				UB = c(BIG, +0.9999), 
				y = x, 
				ghlambda = lambda)
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
	return( fit )
}

dnig = function(x, alpha = 1, beta = 0, delta = 1, mu = 0, log = FALSE)
{   
	# A function implemented by Diethelm Wuertz
	# modified by Alexios Ghalanos (makes use of vectorized dgh)
	# Description:
	#   Returns density for inverse Gaussian DF
	
	# FUNCTION:
	
	# Density:
	return( dgh(x = x, alpha = alpha, beta = beta, delta = delta, mu = mu, 
			lambda = -0.5, log = log) )
}
# ------------------------------------------------------------------------------
pnig <-function(q, alpha = 1, beta = 0, delta = 1, mu = 0)
{   
	# A function implemented by Diethelm Wuertz
	# modified by Alexios Ghalanos (makes use of vectorized pgh)
	
	# Description:
	#   Returns probability for for inverse Gaussian DF
	
	# Function:
	
	# Probability:
	return( pgh(q = q, alpha = alpha, beta = beta, delta = delta, mu = mu, 
			lambda = -0.5) )
}
# ------------------------------------------------------------------------------
qnig <-function(p, alpha = 1, beta = 0, delta = 1, mu = 0)
{   
	# A function implemented by Diethelm Wuertz
	# modified by Alexios Ghalanos (vectorized)
	
	# Description:
	#   Returns quantiles for for inverse Gaussian DF
	
	# FUNCTION:
	
	# Quantiles:
	fun = function(p, alpha = 1, beta = 0, delta = 1, mu = 0){
		return( .qnigC(p = p, alpha = alpha, beta = beta, delta = delta, mu = mu) )
	}
	f = Vectorize( fun )
	ans = f(p, alpha, beta, delta, mu)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
}
# ------------------------------------------------------------------------------
rnig <-function(n, alpha = 1, beta = 0, delta = 1, mu = 0)
{   
	# (uses vectorized rgh)
	return( rgh(n, alpha = alpha, beta = beta, delta = delta, mu = mu, lambda=-0.5) )
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
			p = as.double(.CArrange(p,1,1,n)),
			i_mu = as.double(mu),
			i_delta = as.double(delta),
			i_alpha = as.double(alpha),
			i_beta = as.double(beta),
			i_n = as.integer(n),
			q = as.double(.CArrange(q, 1, 1, n)), PACKAGE="rugarch")
	quantiles <- retValues[[7]]
	quantiles[quantiles <= -1.78e+308] <- -Inf
	quantiles[quantiles >= 1.78e+308] <- Inf
	
	# Return Value:
	return( quantiles )
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
	return(out) 
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
	# modified by Alexios Ghalanos (makes use of vectorized dsgh)
	
	# FUNCTION:
	
	# Compute Density - Quick and Dirty:
	return( dsgh(x, zeta, rho, lambda = -0.5, log = log) )
}
# ------------------------------------------------------------------------------
psnig <-function(q, zeta = 1, rho = 0) 
{
	# Description:
	#   Returns probabilities of the snig distribution
	# modified by Alexios Ghalanos (makes use of vectorized psgh)
	
	# FUNCTION:
	
	# Compute Probabilities - Quick and Dirty:
	return( psgh(q, zeta, rho, lambda = -0.5) )
}
# ------------------------------------------------------------------------------
qsnig <-function(p, zeta = 1, rho = 0) 
{
	# Description:
	#   Returns quantiles of the snig distribution
	# modified by Alexios Ghalanos (vectorized)
	
	# FUNCTION:
	# 4x faster than calling qsgh with lambda=-0.5
	# Compute Quantiles:
	fun = function(p, zeta = 1, rho = 0){
		return(.qsnigC(p, zeta, rho) )
	}
	f = Vectorize( fun )
	ans = f(p, zeta, rho)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
}
# ------------------------------------------------------------------------------
rsnig <-function(n, zeta = 1, rho = 0) 
{
	# Description:
	#   Generates snig distributed random variates
	# modified by Alexios Ghalanos (makes use of vectorized rsgh)
	
	# FUNCTION:
	
	# Generate Random Numbers:
	return( rsgh(n, zeta, rho, lambda = -0.5) )
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
	return( .qnigC(p, param[1], param[2], param[3], param[4]) )
}

snigFit = function (x, zeta = 1, rho = 0, scale = TRUE, doplot = TRUE, 
				span = "auto", trace = 0, title = NULL, description = NULL, ...) 
{   	
	# Update Slots: 
	if (is.null(title)) title = "SNIG Parameter Estimation"	
	# Quick and dirty ...
	ans = sghFit(x, zeta = zeta, rho = rho, lambda = -0.5, include.lambda = FALSE,
			scale = scale,span = span, trace = trace, 
			title = title, description = description, ...)
	
	return( ans )
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
	return( fit )
}
djsu = function(y, mu = 0, sigma = 1, nu = 1, tau = 0.5, log = FALSE)
{
	if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
	if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
	rtau <- 1/tau
	w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
	omega <- -nu*rtau 
	c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
	z <- (y-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
	r <- -nu + asinh(z)/rtau
	loglik <- -log(sigma)-log(c)-log(rtau)-.5*log(z*z+1)-.5*log(2*pi)-.5*r*r
	if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
	return( ft )
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
	return( p )
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
	return( q )
}

rjsu <- function(n, mu=0, sigma=1, nu=0, tau=.5)
{
	fun = function(n, mu=0, sigma=1, nu=0, tau=.5){
		return(.rjsu(n = n, mu = mu, sigma = sigma, nu = nu, tau = tau) )
	}
	if(length(n)>1) stop("\nn cannot be a vector of length greater than 1!")
	f = Vectorize( fun , c("mu","sigma","nu", "tau"))
	ans = f(n, mu, sigma, nu, tau)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )	
}

.rjsu <- function(n, mu=0, sigma=1, nu=0, tau=.5)
{
	if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
	n <- ceiling(n)
	p <- runif(n)
	r <- qjsu(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
	return( r )
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
		shape.LB = 0.1
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
		shape.UB = 25
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
	if(distribution == "ghst"){
		skew 	= 0
		skew.LB	= -80
		skew.UB	= 80
		shape 	= 8.2
		shape.LB = 4.1
		shape.UB = 25
	}
	# johnson has 2 shape parameters. The second one we model with the "skew"
	# representation in rugarch
	skewed.dists = c("snorm", "sged", "sstd", "nig", "ghyp", "jsu", "ghst")
	shaped.dists = c("ged", "sged", "std", "sstd", "nig", "ghyp", "jsu", "ghst")
	skew0  = 0
	shape0 = 0
	if(any(skewed.dists == distribution)) include.skew=TRUE else include.skew=FALSE
	if(any(shaped.dists == distribution)) include.shape=TRUE else include.shape=FALSE
	if(distribution == "ghyp") include.ghlambda = TRUE else include.ghlambda = FALSE
	return( list(shape = shape, shape.LB = shape.LB, shape.UB = shape.UB, skew = skew,
					skew.LB = skew.LB, skew.UB = skew.UB, include.skew = include.skew, 
					include.shape = include.shape, skew0 = skew0, shape0 = shape0,
					include.ghlambda = include.ghlambda, ghlambda = ghlambda, 
					ghlambda.LB = ghlambda.LB, ghlambda.UB = ghlambda.UB) )
}

# ------------------------------------------------------------------------------
.makeSample = function(distribution, lambda = -0.5, skew, shape, n, seed)
{
	set.seed(seed)
	x = switch(distribution,
			norm = rnorm(n),
			snorm = rsnorm(n, xi = skew),
			std = rstd(n, nu = shape),
			sstd = rsstd(n, nu = shape, xi = skew),
			ged = rged(n, nu = shape),
			sged = rsged(n, nu = shape, xi = skew),
			nig = rsnig(n, zeta = shape, rho = skew),
			ghyp = rsgh(n, zeta = shape, rho = skew, lambda = lambda),
			ghst = rsghst(n, mean=0, sd = 1, skew = skew, shape = shape),
			jsu = rjsu(n, mu = 0, sigma = 1, nu = skew, tau = shape)
			)	
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
			jsu   = .jsuscale(mu, sigma, skew, shape),
			ghst  = .ghstscale(mu, sigma, skew, shape)
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


.nigscale = function(mu, sigma, skew, shape)
{
	nigpars = t(apply(cbind(shape, skew), 1, FUN=function(x) .paramGH(zeta = x[1], rho = x[2], lambda = -0.5)))
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	# alpha, beta, delta, mu
	xdensity[,4] = nigpars[,1]/sigma
	xdensity[,3] = nigpars[,2]/sigma
	xdensity[,2] = nigpars[,3]*sigma
	xdensity[,1] = nigpars[,4]*sigma + mu
	# technically: mu, delta, beta, alpha
	colnames(xdensity) = c("mu", "sigma", "skew", "shape")
	return(xdensity)
}

.ghstscale = function(mu, sigma, skew, shape)
{
	ghpars = t(apply(cbind(skew, shape), 1, FUN=function(x) .paramGHST(betabar = x[1], nu = x[2])))
	xdensity = matrix(0, ncol = 5, nrow = length(sigma))
	# alpha, beta, delta, mu
	xdensity[,5] = -shape/2
	xdensity[,4] = (abs(ghpars[,3])+1e-12)/sigma
	xdensity[,3] = ghpars[,3]/sigma
	xdensity[,2] = ghpars[,2]*sigma
	xdensity[,1] = ghpars[,1]*sigma + mu
	# technically: mu, delta, beta, alpha
	colnames(xdensity) = c("mu", "sigma", "skew", "shape", "lambda")
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
	# technically: mu, delta, beta, alpha
	colnames(xdensity) =c("mu", "sigma", "skew", "shape")
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
			"ghyp", "jsu", "ghst")
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
			jsu = jsuFit(x, control),
			ghst = ghstFit(x, control)
	)
	return(ans)
}
ddist = function(distribution = "norm", y, mu = 0, sigma = 1, lambda = -0.5, skew = 1, shape = 5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig",
			"ghyp", "jsu", "ghst")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distributions\n", call. = FALSE)
	ans = switch(distribution,
		norm = dnorm(y, mean = mu, sd = sigma, log = FALSE),
		snorm = dsnorm(y, mean = mu, sd = sigma, xi = skew),
		std = dstd(y, mean = mu, sd = sigma, nu = shape),
		sstd = dsstd(y, mean = mu, sd = sigma, nu = shape, xi = skew),
		ged = dged(y, mean = mu, sd = sigma, nu = shape),
		sged = dsged(y, mean = mu, sd = sigma, nu = shape, xi = skew),
		nig = dsnig((y-mu)/sigma, rho = skew, zeta = shape)/sigma,
		ghyp = dsgh((y-mu)/sigma, rho = skew, zeta = shape, lambda = lambda)/sigma,
		jsu = djsu(y, mu = mu, sigma = sigma, nu = skew, tau = shape),
		ghst = dsghst(y, mean = mu, sd = sigma, skew = skew, shape = shape)
		)
	return(ans)
}

pdist = function(distribution = "norm", q, mu = 0, sigma = 1, lambda = -0.5, skew = 1, shape = 5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig",
			"ghyp", "jsu","ghst")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distributions\n", call. = FALSE)
	ans = switch(distribution,
			norm = pnorm(q, mean = mu, sd = sigma, log.p = FALSE),
			snorm = psnorm(q, mean = mu, sd = sigma, xi = skew),
			std = pstd(q, mean = mu, sd = sigma, nu = shape),
			sstd = psstd(q, mean = mu, sd = sigma, nu = shape, xi = skew),
			ged = pged(q, mean = mu, sd = sigma, nu = shape),
			sged = psged(q, mean = mu, sd = sigma, nu = shape, xi = skew),
			nig = psnig((q-mu)/sigma, rho = skew, zeta = shape),
			ghyp = psgh((q-mu)/sigma, rho = skew, zeta = shape, lambda = lambda),
			jsu = pjsu(q, mu = mu, sigma = sigma, nu = skew, tau = shape),
			ghst = psghst(q, mean = mu, sd = sigma, skew = skew, shape = shape)
	)
	return(ans)
}

qdist = function(distribution = "norm", p, mu = 0, sigma = 1, lambda = -0.5, skew = 1, shape = 5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig",
			"ghyp", "jsu","ghst")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distributions\n", call. = FALSE)
	ans = switch(distribution,
			norm = qnorm(p, mean = mu, sd = sigma, log.p = FALSE),
			snorm = qsnorm(p, mean = mu, sd = sigma, xi = skew),
			std = qstd(p, mean = mu, sd = sigma, nu = shape),
			sstd = qsstd(p, mean = mu, sd = sigma, nu = shape, xi = skew),
			ged = qged(p, mean = mu, sd = sigma, nu = shape),
			sged = qsged(p, mean = mu, sd = sigma, nu = shape, xi = skew),
			nig = qsnig(p, rho = skew, zeta = shape)*sigma + mu,
			ghyp = qsgh(p, rho = skew, zeta = shape, lambda = lambda)*sigma + mu,
			jsu = qjsu(p, mu = mu, sigma = sigma, nu = skew, tau = shape),
			ghst = qsghst(p, mean = mu, sd = sigma, skew = skew, shape = shape),
	)
	return(ans)
}

# EQUAL:
# set.seed(10)
# rsnig(10, rho = skew, zeta = shape)*sigma + mu
# set.seed(10)
# rdist("nig", 10, mu = mu, sigma = sigma, skew = skew, shape = shape)
rdist = function(distribution = "norm", n, mu = 0, sigma = 1, lambda = -0.5, skew = 1, shape = 5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig",
			"ghyp", "jsu","ghst")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distribution\n", call. = FALSE)
	ans = switch(distribution,
			norm = rnorm(n, mean = mu, sd = sigma),
			snorm = rsnorm(n, mean = mu, sd = sigma, xi = skew),
			std = rstd(n, mean = mu, sd = sigma, nu = shape),
			sstd = rsstd(n, mean = mu, sd = sigma, nu = shape, xi = skew),
			ged = rged(n, mean = mu, sd = sigma, nu = shape),
			sged = rsged(n, mean = mu, sd = sigma, nu = shape, xi = skew),
			nig =  mu + sigma*rsnig(n, rho = skew, zeta = shape),
			ghyp = mu + sigma*rsgh(n, rho = skew, zeta = shape, lambda = lambda),
			jsu = rjsu(n, mu = mu, sigma = sigma, nu = skew, tau = shape),
			ghst = rsghst(n, mean = mu, sd = sigma, skew = skew, shape = shape)
	)
	return(ans)
}


dskewness = function(distribution = "norm", skew = 1, shape = 5, lambda = -0.5)
{
	f = Vectorize(.dskewness)
	ans = f(distribution, skew, shape, lambda)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return(ans)
}
	
.dskewness = function(distribution = "norm", skew = 1, shape = 5, lambda = -0.5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "jsu","ghst")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distribution\n", call. = FALSE)
	ans = switch(distribution,
			norm 	= 0,
			snorm 	= .snormskew(skew = skew),
			std 	= 0,
			sstd 	= .sstdskew(skew = skew, shape = shape),
			ged 	= 0,
			sged 	= .sgedskew(skew = skew, shape = shape),
			nig 	= .snigskew(skew = skew, shape = shape),
			ghyp 	= .sghypskew(skew = skew, shape = shape, lambda = lambda),
			jsu 	= .jsuskew(mu = 0, sigma = 1, skew = skew, shape = shape),
			ghst	= .ghstskew(skew, shape)
	)
	return(as.numeric(ans))
}

dkurtosis = function(distribution = "norm", skew = 1, shape = 5, lambda = -0.5)
{
	f = Vectorize(.dkurtosis)
	ans = f(distribution, skew, shape, lambda)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return(ans)
}

.dkurtosis = function(distribution = "norm", skew = 1, shape = 5, lambda = -0.5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "jsu","ghst")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distribution\n", call. = FALSE)
	ans = switch(distribution,
			norm 	= 0,
			snorm 	= 0,
			std 	= .stdexkurt(shape = shape),
			sstd 	= .sstdexkurt(skew = skew, shape = shape),
			ged 	= .gedexkurt(shape = shape),
			sged 	= .sgedexkurt(skew = skew, shape = shape),
			nig 	= .snigexkurt(skew = skew, shape = shape),
			ghyp 	= .sghypexkurt(skew = skew, shape = shape, lambda = lambda),
			jsu 	= .jsuexkurt(mu = 0, sigma = 1, skew = skew, shape = shape),
			ghst	= .ghstexkurt(skew, shape)
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

.snigskew = function(skew, shape){
	fun = function(skew, shape){
		pars = .paramGH(zeta = shape, rho = skew, lambda = -0.5)
		return(.nigskew(alpha = pars[1], beta = pars[2], delta = pars[3], mu = pars[4]))
	}
	f = Vectorize( fun )
	ans = f(skew, shape)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
}

.nigskew = function(alpha, beta, delta, mu){
	gm = sqrt(alpha^2 - beta^2)
	ans = 3*beta/(alpha*sqrt(delta*gm))
	return(ans)
}

.snigexkurt = function(skew, shape){
	fun = function(skew, shape){
		pars = .paramGH(zeta = shape, rho = skew, lambda = -0.5)
		return(.nigexkurt(alpha = pars[1], beta = pars[2], delta = pars[3], mu = pars[4]))
	}
	f = Vectorize( fun )
	ans = f(skew, shape)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
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

.sghypskew = function(skew, shape, lambda){
	fun = function(skew, shape, lambda){
		pars = .paramGH(zeta = shape, rho = skew, lambda = lambda)
		return(.ghypskew(lambda = lambda, alpha = pars[1], beta = pars[2], delta = pars[3], mu = pars[4]))
	}
	f = Vectorize( fun )
	ans = f(skew, shape, lambda)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
}

.ghypskew = function(lambda, alpha, beta, delta, mu){
	skew = ghypMom(3, lambda, alpha, beta, delta, mu, momType = "central")/(.ghypsigma(lambda, alpha, beta, delta, mu)^3)
	return(skew)
}

.sghypexkurt = function(skew, shape, lambda){
	fun = function(skew, shape, lambda){
		pars = .paramGH(zeta = shape, rho = skew, lambda = lambda)
		return(.ghypexkurt(lambda = lambda, alpha = pars[1], beta = pars[2], delta = pars[3], mu = pars[4]))
	}
	f = Vectorize( fun )
	ans = f(skew, shape, lambda)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
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

.sstdskew = function(skew, shape){
	# Theoretical moments based on bijection betweeen Fernandez and Steel verions
	# and Hansen's Generalized Skew-T (credit due to Michael Rockinger)
	if (shape > 2) {
		eta  = shape
		k2   = skew^2
		lda  = (k2-1)/(k2+1)
		ep1 = (eta+1)/2
		lnc = lgamma(ep1) - lgamma(eta/2) -0.5*log( pi*(eta-2))
		c   = exp(lnc)
		a   = 4*lda*c*(eta-2)/(eta-1)
		b   = sqrt(1+3*lda^2-a^2)
		my2 = 1+3*lda^2
		my3 = 16*c*lda*(1+lda^2)*((eta-2)^2)/((eta-1)*(eta-3))
		my4 = 3*(eta-2)*(1+10*lda^2+5*lda^4)/(eta-4)
		m3  = (my3-3*a*my2+2*a^3)/(b^3)
	} else{
		m3 = NA
	}
	return(m3)
}

.sstdexkurt = function( skew, shape )
{
	# Theoretical moments based on bijection betweeen Fernandez and Steel verions
	# and Hansen's Generalized Skew-T (credit due to Michael Rockinger)
	if(shape > 4 ){
		eta  = shape
		k2   = skew^2
		lda  = (k2-1)/(k2+1)
		ep1 = (eta+1)/2
		lnc = lgamma(ep1) - lgamma(eta/2) -0.5*log( pi*(eta-2))
		c   = exp(lnc)
		a   = 4*lda*c*(eta-2)/(eta-1)
		b   = sqrt(1+3*lda^2-a^2)
		my2 = 1+3*lda^2
		my3 = 16*c*lda*(1+lda^2)*((eta-2)^2)/((eta-1)*(eta-3))
		my4 = 3*(eta-2)*(1+10*lda^2+5*lda^4)/(eta-4)
		m3  = (my3-3*a*my2+2*a^3)/(b^3)
		m4  = -3 + (my4-4*a*my3+6*(a^2)*my2-3*a^4)/(b^4)
	} else{
		m4 = NA
	}
	return(m4)
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
	Omega = -skew/shape
	w = exp(shape^-2)
	s3 = -0.25*sqrt(w)*( (w-1)^2 )*(w*(w+2)*sinh(3*Omega)+3*sinh(Omega))
	s3/(0.5*(w-1)*(w*cosh(2*Omega)+1))^(3/2)
}

.jsuexkurt = function( mu = 0, sigma = 1, skew, shape )
{
	Omega = -skew/shape
	w = exp(shape^-2)
	s4 = 0.125 * (w-1)^2*(w^2*(w^4+2*w^3+3*w^2-3)*cosh(4*Omega)+4*w^2*(w+2)*cosh(2*Omega)+3*(2*w+1))
	ans = s4/(0.5*(w-1)*(w*cosh(2*Omega)+1))^2
	return(ans - 3)
}

.ghstskew = function(skew, shape){
	if(shape<6){
		ans = NA
	} else{
		params = .paramGHST(nu = shape, betabar = skew)
		delta = params[2]
		beta = params[3]
		nu = params[4]
		beta2 = beta*beta
		delta2 = delta*delta
		ans = ( (2 * sqrt(nu - 4)*beta*delta)/( (2*beta2*delta2 + (nu-2)*(nu-4))^(3/2) ) ) * (3*(nu-2) + ((8*beta2*delta2)/(nu-6)))
	}
	return( ans )
}

.ghstexkurt = function(skew, shape){
	if(shape<8){
		ans = NA
	} else{
		params = .paramGHST(nu = shape, betabar = skew)
		delta = params[2]
		beta = params[3]
		nu = params[4]
		beta2 = beta*beta
		delta2 = delta*delta
		k1 = 6/( (2*beta2*delta2+(nu-2)*(nu-4))^2)
		k21 = (nu-2)*(nu-2)*(nu-2)
		k22 = (16*beta2*delta2*(nu-2)*(nu-4))/(nu-6)
		k23 = (8*(beta2^2)*(delta2^2)*(5*nu-22))/((nu-6)*(nu-8))
		ans = k1*(k21+k22+k23)
	}
	return( ans )
}