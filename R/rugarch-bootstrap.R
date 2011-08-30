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

# this function implements a simulation based method for generating standard errors and the
# bootstrapped predictive density for the n.ahead forecasts of a garch model, 
# taking into account parameter uncertainty partly based on the paper by Pascual, Romo and
# Ruiz (2006) : Computational Statistics and Data Analysis 50
# "Bootstrap prediction for returns and volatilities in GARCH models"
# [PRR paper]

# the conditional bootstrap (Partial) does not consider parameter uncertainty
# the Full model does (but is more expensive since we need to simulate the parameter distribution)


.ugarchbootfit = function(fitORspec, data = NULL, method = c("Partial", "Full"), n.ahead = 10, n.bootfit = 100, 
		n.bootpred = 500, out.sample = 0, rseed = NA, solver = "solnp", solver.control = list(), fit.control = list(), 
		external.forecasts =  list(mregfor = NULL, vregfor = NULL), parallel = FALSE, 
		parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2))
{

	method = tolower(method)
	ans = switch(method,
			partial = .ub1p1(fitORspec, data = data, n.ahead = n.ahead, n.bootfit = n.bootfit, 
						n.bootpred = n.bootpred, rseed = rseed, solver.control = solver.control, 
						fit.control = fit.control, external.forecasts =  external.forecasts, 
						parallel = parallel, parallel.control = parallel.control),
			full = .ub1f1(fitORspec, data = data, n.ahead = n.ahead, n.bootfit = n.bootfit, 
					n.bootpred = n.bootpred, rseed = rseed, solver = solver, solver.control = solver.control, 
					fit.control = fit.control, external.forecasts =  external.forecasts, 
					parallel = parallel, parallel.control = parallel.control))
	return(ans)
}

.ugarchbootspec = function(fitORspec, data = NULL, method = c("Partial", "Full"), n.ahead = 10, n.bootfit = 100, 
		n.bootpred = 500, out.sample = 0, rseed = NA, solver = "solnp", solver.control = list(), fit.control = list(), 
		external.forecasts =  list(mregfor = NULL, vregfor = NULL), parallel = FALSE, 
		parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2))
{
	method = tolower(method)
	ans = switch(method,
			partial = .ub1p2(fitORspec, data = data, n.ahead = n.ahead, n.bootfit = n.bootfit, 
					n.bootpred = n.bootpred, out.sample = out.sample, rseed = rseed, solver.control = solver.control, 
					fit.control = fit.control, external.forecasts =  external.forecasts, 
					parallel = parallel, parallel.control = parallel.control),
			full = .ub1f2(fitORspec, data = data, n.ahead = n.ahead, n.bootfit = n.bootfit, 
					n.bootpred = n.bootpred, out.sample = out.sample, rseed = rseed, solver = solver, 
					solver.control = solver.control, fit.control = fit.control, external.forecasts =  external.forecasts, 
					parallel = parallel, parallel.control = parallel.control))
	return(ans)
}

# method from a fit object
.ub1f1 = function(fitORspec, data = NULL, n.ahead = 10, n.bootfit = 100, n.bootpred = 500, rseed = NA, solver = "solnp",
		solver.control = list(), fit.control = list(), external.forecasts =  list(mregfor = NULL, vregfor = NULL), parallel = FALSE, 
		parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2))
{
	
	# inputs:
	# n.bootfit		: the number of simulations for which we will generate the parameter 
	#				uncertainty distribution (i.e. no of refits)
	# n.bootpred	: the number of simulations / h.ahead to use (for averaging to obtain forecast
	#				density)
	# out.sample	: data to keep for out of sample comparison purposes
	# n.ahead		: the horizon of the bootstrap density forecast
	fit = fitORspec
	model = fit@model
	m = model$maxOrder
	data = model$modeldata$data
	N = model$modeldata$T
	ns = model$n.start
	if (is.na(rseed[1])){
		sseed = as.integer(runif(n.bootpred + n.bootfit, 0, 65000))
	} else{
		if(length(rseed) < n.bootpred){
			stop("\nugarchboot-->error: seed length must equal n.bootpred + n.bootfit for full method\n")
		} else {
			sseed = rseed
		}
	}
	sseed1 = sseed[1:n.bootfit]
	sseed2 = sseed[(n.bootfit+1):(n.bootpred + n.bootfit)]
	# generate paths of equal length to data based on empirical re-sampling of z
	# p.2296 equation (5)
	fz = fit@fit$z
	empz = matrix(0, ncol = n.bootfit, nrow = N)
	
	empz = apply(as.data.frame(1:n.bootfit), 1, FUN=function(i) {
				set.seed(sseed1[i]);
				sample(fz, size = N, replace = TRUE);})
		
	if(ns > 0) {
		spec = getspec(fit)
		spec@model$fixed.pars = as.list(coef(fit))
		realized.x = fit@model$modeldata$data[(N+1):(N+ns)]
		filtered.s = tail(sigma(ugarchfilter(data = fit@model$modeldata$data, spec = spec)), fit@model$n.start)
	} else{
		realized.x = NULL
		filtered.s = NULL
	}
	
	# presigma uses the same starting values as the original fit
	# -> in paper they use alternatively the unconditional long run sigma (P.2296 paragraph 2
	# "...marginal variance..."
	paths = ugarchsim(fit, n.sim = N, m.sim = n.bootfit, presigma = tail(fit@fit$sigma, m), 
			prereturns = tail(model$modeldata$data[1:N], m), preresiduals = tail(residuals(fit), m), 
			startMethod = "sample", custom.dist = list(name = "sample", distfit = as.matrix(empz)),
			rseed = sseed1, mexsimdata = external.forecasts$mregfor, vexsimdata = external.forecasts$vregfor)
	fitlist = vector(mode="list", length = n.bootfit)
	path.df = as.data.frame(paths, which = "series")
	spec = getspec(fit)
	rownames(path.df) = as.character(model$modeldata$dates[1:N])
	# help the optimization with good starting parameters
	spec@model$start.pars = as.list(coef(fit))
	nx = dim(path.df)[2]
	
	# get the distribution of the parameters (by fitting to the new path data)
	#-------------------------------------------------------------------------
	cat("\nfitting stage...")
	if( parallel ){
		os = .Platform$OS.type
		if(is.null(parallel.control$pkg)){
			if( os == "windows" ) parallel.control$pkg = "snowfall" else parallel.control$pkg = "multicore"
			if( is.null(parallel.control$cores) ) parallel.control$cores = 2
		} else{
			mtype = match(tolower(parallel.control$pkg[1]), c("multicore", "snowfall"))
			if(is.na(mtype)) stop("\nParallel Package type not recognized in parallel.control\n")
			parallel.control$pkg = tolower(parallel.control$pkg[1])
			if( os == "windows" && parallel.control$pkg == "multicore" ) stop("\nmulticore not supported on windows O/S\n")
			if( is.null(parallel.control$cores) ) parallel.control$cores = 2 else parallel.control$cores = as.integer(parallel.control$cores[1])
		}
		if( parallel.control$pkg == "multicore" ){
			if(!exists("mclapply")){
				require('multicore')
			}
			fitlist = mclapply(1:nx, FUN = function(i) ugarchfit(spec = spec, data = path.df[, i, drop = FALSE], out.sample = 0,
						solver = solver, fit.control = fit.control, solver.control = solver.control), mc.cores = parallel.control$cores)
		} else{
			if(!exists("sfLapply")){
				require('snowfall')
			}
			nx = dim(path.df)[2]
			sfInit(parallel=TRUE, cpus = parallel.control$cores)
			sfExport("path.df", "spec", "solver", "out.sample", "solver.control", local = TRUE)
			fitlist = sfLapply(as.list(1:nx), fun = function(i)
						rugarch::ugarchfit(spec = spec, data = path.df[, i, drop = FALSE], out.sample = 0,
								solver = solver, fit.control = fit.control, solver.control = solver.control))
			sfStop()
		}
	} else{
		for(i in 1:n.bootfit){
			fitlist[[i]] = .safefit(spec = spec, data = path.df[, i, drop = FALSE], out.sample = 0, fit.control = fit.control, 
					solver = solver, solver.control = solver.control)
		}
	}
	# check for any non convergence and remove
	if(any(sapply(fitlist,FUN=function(x) is.null(coef(x))))){
		exc = which(sapply(fitlist,FUN=function(x) is.null(coef(x))))
		fitlist = fitlist[-exc]
		n.bootfit = n.bootfit - length(exc)
		# in case something went very wrong:
		if(n.bootfit == 0) stop("\nugarchboot-->error: the fitting routine failed. No convergene at all!\n", call. = FALSE)
	}
	cat("done!\n")
	# extract the coefficient distribution and generate spec list to feed to path function
	coefdist = matrix(NA, ncol = length(coef(fit)), nrow = n.bootfit)
	speclist = vector(mode = "list", length = n.bootfit)
	for(i in 1:n.bootfit){
		coefdist[i,] = coef(fitlist[[i]])
		speclist[[i]] = spec
		speclist[[i]]@model$fixed.pars = as.list(coef(fitlist[[i]]))
	}
	colnames(coefdist) = names(coef(fit))
	#-------------------------------------------------------------------------

	# generate path based forecast values
	# for each path we generate n.bootpred vectors of resampled data of length n.ahead
	# Equation (6) in the PRR paper (again using z from original fit)
	#-------------------------------------------------------------------------
	empzlist = vector(mode = "list", length = n.bootfit)
	for(i in 1:n.bootfit){
		empz = matrix(0, ncol = n.bootpred, nrow = n.ahead)
		empz = apply(as.data.frame(1:n.bootpred), 1, FUN=function(i) {
					set.seed(sseed2[i]);
					sample(fz, size = n.ahead, replace = TRUE);})
		empzlist[[i]] = matrix(empz, ncol = n.bootpred)
	}
	# we start all forecasts from the last value of sigma based on the original series but
	# the pertrubed parameters as in equation (7) in PRR paper.
	if( parallel ){
		if( parallel.control$pkg == "multicore" ){
			st = mclapply(as.list(1:n.bootfit), FUN = function(x) .sigmat(spec = speclist[[i]], origdata = data[1:N], m), mc.cores = parallel.control$cores)
			st = unlist(st)
		} else{
			nx = length(speclist)
			sfInit(parallel = TRUE, cpus = parallel.control$cores)
			sfExport("speclist", "data", "m", "N", local = TRUE)
			st = sfLapply(as.list(1:n.bootfit), fun = function(i) rugarch:::.sigmat(spec = speclist[[i]], origdata =  data[1:N], m))
			sfStop()
			st = unlist(st)
		}
	} else{
		st = sapply(speclist, FUN=function(x) .sigmat(spec = x, origdata =  data[1:N], m))
	}
	forcseries = NULL
	forcsigma  = NULL
	tmp = vector(mode = "list", length = n.bootfit)
	cat("\nprediction stage...")
	xdat = tail( data[1:N], m )
	if( parallel ){
		if( parallel.control$pkg == "multicore" ){
		tmp = mclapply(as.list(1:n.bootfit), FUN = function(i) .quicksimulate(fitlist[[i]], n.sim = n.ahead, 
							m.sim = n.bootpred, presigma = st[i], prereturns = xdat,  
							n.start = 0, rseed = sseed[(n.bootfit+1):(n.bootfit + n.bootpred)], 
							custom.dist = list(name = "sample", distfit = as.matrix(empzlist[[i]])), 
							mexsimdata = external.forecasts$mregfor, vexsimdata = external.forecasts$vregfor), 
				mc.cores = parallel.control$cores)
		} else{
			sfInit(parallel = TRUE, cpus = parallel.control$cores)
			sfExport("fitlist", "n.ahead", "n.bootpred", "n.bootfit", "st", "xdat", "sseed", "empzlist",
					"external.forecasts", local = TRUE)
			tmp = sfLapply(as.list(1:n.bootfit), fun = function(i) .quicksimulate(fitlist[[i]], n.sim = n.ahead, 
								m.sim = n.bootpred, presigma = st[i], prereturns = xdat,  
								n.start = 0, rseed = sseed[(n.bootfit+1):(n.bootfit + n.bootpred)], 
								custom.dist = list(name = "sample", distfit = as.matrix(empzlist[[i]])), 
								mexsimdata = external.forecasts$mregfor, vexsimdata = external.forecasts$vregfor))
			sfStop()
		}
	} else{
		tmp = vector(mode = "list", length = n.bootfit)
		for(i in 1:n.bootfit){
			tmp[[i]] = .quicksimulate(fit = fitlist[[i]], n.sim = n.ahead, 
					m.sim = n.bootpred, presigma = st[i], prereturns = xdat,  
					n.start = 0, rseed = sseed[(n.bootfit+1):(n.bootfit + n.bootpred)], 
					custom.dist = list(name = "sample", distfit = as.matrix(empzlist[[i]])), 
					mexsimdata = external.forecasts$mregfor, vexsimdata = external.forecasts$vregfor)
			# reduce memory
			empzlist[[i]] = 0
			gc(verbose = FALSE)
		}
	}
	# do some cleanup
	rm(empzlist)
	gc(verbose = FALSE)

	# we will have (n.bootfit x n.bootpred) x n.ahead matrix
	forcseries = lapply(tmp, FUN = function(x) x[1:n.ahead, ,drop = F])
	meanseriesfit = t(sapply(forcseries, FUN = function(x) apply(t(x), 2, "mean")))
	forcseries = matrix(unlist(forcseries), nrow = n.ahead, ncol = n.bootpred * n.bootfit, byrow = FALSE)
	forcseries = t(forcseries)
	forcsigma = lapply(tmp, FUN = function(x) x[(n.ahead+1):(2*n.ahead), ,drop = F])
	meansigmafit = t(sapply(forcsigma, FUN = function(x) apply(t(x), 2, "mean")))
	forcsigma = matrix(unlist(forcsigma), nrow = n.ahead, ncol = n.bootpred * n.bootfit, byrow = FALSE)
	forcsigma = t(forcsigma)
	# do some cleanup
	rm(tmp)
	gc(verbose = FALSE)	
	cat("done!\n")
	#-------------------------------------------------------------------------

	# now we have the bootstrapped distribution of n.ahead forecast values
	# original forecast
	forc = ugarchforecast(fitORspec = fit, n.ahead = n.ahead, n.roll = 0, external.forecasts = external.forecasts) 	
	model$truecoef = coef(fit)
	model$modeldata$realized.x = realized.x
	model$modeldata$filtered.s = filtered.s
	model$n.ahead = n.ahead
	model$n.bootfit = n.bootfit
	model$n.bootpred = n.bootpred
	model$modeldata$meanfit.x = meanseriesfit
	model$modeldata$meanfit.s = meansigmafit
	model$seeds = sseed
	model$type = "full"
	ans = new("uGARCHboot",
			fseries = forcseries,
			fsigma = forcsigma,
			bcoef = as.data.frame(coefdist),
			model = model,
			forc = forc)
	return(ans)
}


# method from a spec object
.ub1f2 = function(fitORspec, data = NULL, n.ahead = 10, n.bootfit = 100, 
		n.bootpred = 500, out.sample = 0, rseed = NA, solver = "solnp", solver.control = list(), fit.control = list(), 
		external.forecasts =  list(mregfor = NULL, vregfor = NULL), parallel = FALSE, 
		parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2))
{

	spec = fitORspec
	model = spec@model
	if(is.null(data))
		stop("\nugarchboot-->error: data must be supplied if SPEC object used.", call. = FALSE)	
	
	flt = ugarchfilter(data = data, spec = spec, out.sample = out.sample)
	xdata = flt@model$modeldata$data
	m = spec@model$maxOrder
	N = flt@model$modeldata$T
	if (is.na(rseed[1])){
		sseed = as.integer(runif(n.bootpred + n.bootfit,0,65000))
	} else{
		if(length(rseed) < n.bootpred){
			stop("\nugarchboot-->error: seed length must equal n.bootpred + n.bootfit for full method\n")
		} else {
			sseed = rseed
		}
	}
	sseed1 = sseed[1:n.bootfit]
	sseed2 = sseed[(n.bootfit+1):(n.bootpred + n.bootfit)]
	
	if(out.sample > 0) {
		flt2 = ugarchfilter(data = data, spec = spec, out.sample = 0)
		realized.x = tail(xdata, out.sample)
		filtered.s = tail(sigma(flt2), out.sample)
	} else{
		realized.x = NULL
		filtered.s = NULL
	}
	
	
	# generate paths of equal length to data based on empirical re-sampling of z
	# p.2296 equation (5)
	# use of filter when using spec/this will also check whether the fixed.pars are correctly
	# specified
	fz = flt@filter$z
	empz = matrix(0, ncol = n.bootfit, nrow = N)
	empz = apply(as.data.frame(1:n.bootfit), 1, FUN=function(i) {
				set.seed(sseed1[i]);
				sample(fz, size = N, replace = TRUE);})
	# presigma uses the same starting values as the original fit
	# -> in paper they use alternatively the unconditional long run sigma (P.2296 paragraph 2
	# "...marginal variance...")
	paths = ugarchpath(spec, n.sim = N, m.sim = n.bootfit, presigma = tail(flt@filter$sigma, m), 
			prereturns = tail(xdata[1:N], m), preresiduals = tail(residuals(flt), m), rseed = sseed1,
			n.start = 0, custom.dist = list(name = "sample", distfit = as.matrix(empz)), 
			mexsimdata = external.forecasts$mregfor, vexsimdata = external.forecasts$vregfor)
	fitlist = vector(mode="list", length = n.bootfit)
	path.df = as.data.frame(paths, which = "series")
	rownames(path.df) = as.character(flt@model$modeldata$dates[1:N])
	nx = dim(path.df)[2]
	
	# help the optimization with good starting parameters
	spex = spec
	spex@model$fixed.pars = NULL
	spex@model$start.pars = spec@model$fixed.pars
	# get the distribution of the parameters (by fitting to the new path data)
	#-------------------------------------------------------------------------
	cat("\nfitting stage...")
	if( parallel ){
		os = .Platform$OS.type
		if(is.null(parallel.control$pkg)){
			if( os == "windows" ) parallel.control$pkg = "snowfall" else parallel.control$pkg = "multicore"
			if( is.null(parallel.control$cores) ) parallel.control$cores = 2
		} else{
			mtype = match(tolower(parallel.control$pkg[1]), c("multicore", "snowfall"))
			if(is.na(mtype)) stop("\nParallel Package type not recognized in parallel.control\n")
			parallel.control$pkg = tolower(parallel.control$pkg[1])
			if( os == "windows" && parallel.control$pkg == "multicore" ) stop("\nmulticore not supported on windows O/S\n")
			if( is.null(parallel.control$cores) ) parallel.control$cores = 2 else parallel.control$cores = as.integer(parallel.control$cores[1])
		}
		if( parallel.control$pkg == "multicore" ){
			if(!exists("mclapply")){
				require('multicore')
			}
			fitlist = mclapply(as.list(1:nx), FUN = function(x) ugarchfit(spec = spex, data = path.df[,i,drop = FALSE], out.sample = 0,
								solver = solver, fit.control = fit.control, solver.control = solver.control), mc.cores = parallel.control$cores)
		} else{
			if(!exists("sfLapply")){
				library('snowfall', pos = "package:base")
			}
			nx = dim(path.df)[2]
			sfInit(parallel=TRUE, cpus = parallel.control$cores)
			sfExport("path.df", "spex", "solver", "out.sample", "solver.control", local = TRUE)
			fitlist = sfLapply(as.list(1:nx), fun = function(i)
						rugarch::ugarchfit(spec = spex, data = path.df[,i,drop = FALSE], out.sample = 0,
								solver = solver, fit.control = fit.control, solver.control = solver.control))
			sfStop()
		}
	} else{
		for(i in 1:n.bootfit){
			fitlist[[i]] = .safefit(spec = spex, data = path.df[,i,drop = FALSE], out.sample = 0, fit.control = fit.control, 
					solver = solver, solver.control = solver.control)
		}
	}
	# check for any non convergence and remove
	if(any(sapply(fitlist,FUN=function(x) is.null(coef(x))))){
		exc = which(sapply(fitlist,FUN=function(x) is.null(coef(x))))
		fitlist = fitlist[-exc]
		n.bootfit = n.bootfit - length(exc)
		# in case something went very wrong:
		if(n.bootfit == 0) stop("\nugarchboot-->error: the fitting routine failed. No convergene at all!\n", call. = FALSE)
	}
	cat("done!\n")
	
	# extract the coefficient distribution and generate spec list to feed to path function
	coefdist = matrix(NA, ncol = length(coef(flt)), nrow = n.bootfit)
	speclist = vector(mode = "list", length = n.bootfit)
	for(i in 1:n.bootfit){
		coefdist[i,] = coef(fitlist[[i]])
		speclist[[i]] = spex
		speclist[[i]]@model$start.pars = NULL
		speclist[[i]]@model$fixed.pars = as.list(coef(fitlist[[i]]))
	}
	colnames(coefdist) = names(coef(flt))
	
	# generate path based forecast values
	# for each path we generate n.bootpred vectors of resampled data of length n.ahead
	# Equation (6) in the PRR paper (again using z from original fit)
	empzlist = vector(mode = "list", length = n.bootfit)
	for(i in 1:n.bootfit){
		empz = matrix(0, ncol = n.bootpred, nrow = n.ahead)
		empz = apply(as.data.frame(1:n.bootpred), 1, FUN=function(i) {
					set.seed(sseed2[i]);
					sample(fz, size = n.ahead, replace = TRUE);})
		empzlist[[i]] = matrix(empz, ncol = n.bootpred)
	}
	
	# we start all forecasts from the last value of sigma based on the original series but
	# the pertrubed parameters as in equation (7) in PRR paper.
	if( parallel ){
		if( parallel.control$pkg == "multicore" ){
			st = mclapply(speclist, FUN=function(x) .sigmat(spec = x, origdata = xdata[1:N], m), mc.cores = parallel.control$cores)
			st = unlist(st)
		} else{
			nx = length(speclist)
			sfInit(parallel = TRUE, cpus = parallel.control$cores)
			sfExport("speclist", "xdata", "m", "N", local = TRUE)
			st = sfLapply(as.list(1:nx), fun = function(i) rugarch:::.sigmat(spec = speclist[[i]], origdata = xdata[1:N], m))
			sfStop()
			st = unlist(st)
		}
	} else{
		st = sapply(speclist, FUN=function(x) .sigmat(spec = x, origdata = xdata[1:N], m))
	}

	forcseries = NULL
	forcsigma  = NULL
	forcseries = NULL
	forcsigma  = NULL
	tmp = vector(mode = "list", length = n.bootfit)
	cat("\nprediction stage...")
	
	xdat =  tail(xdata[1:N], m)
	if( parallel ){
		if( parallel.control$pkg == "multicore" ){			
			tmp = mclapply(as.list(1:n.bootfit), FUN = function(i) .quicksimulate(fitlist[[i]], n.sim = n.ahead, 
								m.sim = n.bootpred, presigma = st[i], prereturns = xdat,  
								n.start = 0, rseed = sseed[(n.bootfit+1):(n.bootfit + n.bootpred)], 
								custom.dist = list(name = "sample", distfit = as.matrix(empzlist[[i]])), 
								mexsimdata = external.forecasts$mregfor, vexsimdata = external.forecasts$vregfor), 
					mc.cores = parallel.control$cores)
		} else{
			sfInit(parallel = TRUE, cpus = parallel.control$cores)
			sfExport("fitlist", "n.ahead", "n.bootpred", "n.bootfit", "st", "xdat", "sseed", "empzlist",
					"external.forecasts", local = TRUE)
			tmp = sfLapply(as.list(1:n.bootfit), fun = function(i) .quicksimulate(fitlist[[i]], n.sim = n.ahead, 
								m.sim = n.bootpred, presigma = st[i], prereturns = xdat,  
								n.start = 0, rseed = sseed[(n.bootfit+1):(n.bootfit + n.bootpred)], 
								custom.dist = list(name = "sample", distfit = as.matrix(empzlist[[i]])), 
								mexsimdata = external.forecasts$mregfor, vexsimdata = external.forecasts$vregfor))
			sfStop()
		}
	} else{
		tmp = vector(mode = "list", length = n.bootfit)
		for(i in 1:n.bootfit){
			tmp[[i]] = .quicksimulate(fit = fitlist[[i]], n.sim = n.ahead, 
					m.sim = n.bootpred, presigma = st[i], prereturns = xdat,  
					n.start = 0, rseed = sseed[(n.bootfit+1):(n.bootfit + n.bootpred)], 
					custom.dist = list(name = "sample", distfit = as.matrix(empzlist[[i]])), 
					mexsimdata = external.forecasts$mregfor, vexsimdata = external.forecasts$vregfor)
			empzlist[[i]] = 0
			gc(verbose = FALSE)
		}
	}
	# we will have (n.bootfit x n.bootpred) x n.ahead matrix
	forcseries = lapply(tmp, FUN = function(x) x[1:n.ahead, ,drop = F])
	meanseriesfit = t(sapply(forcseries, FUN = function(x) apply(t(x), 2, "mean")))
	forcseries = matrix(unlist(forcseries), nrow = n.ahead, ncol = n.bootpred * n.bootfit, byrow = FALSE)
	forcseries = t(forcseries)
	forcsigma = lapply(tmp, FUN = function(x) x[(n.ahead+1):(2*n.ahead), ,drop = F])
	meansigmafit = t(sapply(forcsigma, FUN = function(x) apply(t(x), 2, "mean")))
	forcsigma = matrix(unlist(forcsigma), nrow = n.ahead, ncol = n.bootpred * n.bootfit, byrow = FALSE)
	forcsigma = t(forcsigma)
	cat("done!\n")
	rm(empzlist)
	rm(tmp)
	rm(fitlist)
	gc(verbose = FALSE)
	#-------------------------------------------------------------------------
	# now we have the bootstrapped distribution of n.ahead forecast values
	# original forecast
	forc = ugarchforecast(fitORspec = spec, data = head(data, N), n.ahead = n.ahead, n.roll = 0, external.forecasts = external.forecasts)	
	model$truecoef = coef(flt)
	model$modeldata$realized.x = realized.x
	model$modeldata$filtered.s = filtered.s
	model$n.ahead = n.ahead
	model$n.bootfit = n.bootfit
	model$n.bootpred = n.bootpred
	model$modeldata$meanfit.x = meanseriesfit
	model$modeldata$meanfit.s = meansigmafit
	model$seeds = sseed
	model$type = "full"
	ans = new("uGARCHboot",
			fseries = forcseries,
			fsigma = forcsigma,
			bcoef = as.data.frame(coefdist),
			model = model,
			forc = forc)
	return(ans)
}

# Partial method (very fast but does not take into account the parameter uncertainty)
.ub1p1 = function(fitORspec, data = NULL, n.ahead = 10, n.bootfit = 100, 
		n.bootpred = 500, rseed = NA, solver.control = list(), fit.control = list(),
		external.forecasts =  list(mregfor = NULL, vregfor = NULL), parallel = FALSE, 
		parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2))
{
	fit = fitORspec
	model = fit@model
	m = fit@model$maxOrder
	ns = fit@model$n.start
	data = data.frame(fit@model$modeldata$data)
	rownames(data) = as.character( fit@model$modeldata$dates )
	xdata = fit@model$modeldata$data
	N = length(xdata) - ns
	spec = getspec(fit)
	if(is.na(rseed[1])){
		sseed = as.integer(runif(n.bootpred,0,65000))
	} else{
		if(length(rseed) < n.bootpred){
			stop("\nugarchboot-->error: seed length must equal n.bootpred for partial method\n")
		} else {
			sseed = rseed
		}
	}
	# generate path based forecast values
	# for each path we generate n.bootpred vectors of resampled data of length n.ahead
	# Equation (6) in the PRR paper
	fz = fit@fit$z
	empz = matrix(apply(as.data.frame(1:n.bootpred), 1, FUN = function(i) {
														set.seed(sseed[i]);
														sample(fz, size = n.ahead, replace = TRUE);}), ncol = n.bootpred)
	
	#empz[1:m,] = matrix(rep(tail(fit@fit$z, m), n.bootpred), nrow = m)
	
	# we start all forecasts from the last value of sigma based on the original series
	spec@model$fixed.pars = as.list(coef(fit))
	st = .sigmat(spec, xdata[1:N], m)
	
	if(ns > 0){
		realized.x = xdata[(N+1):(N+ns)]
		filtered.s = tail(sigma(ugarchfilter(data = xdata, spec = spec)), ns)
	} else{
		realized.x = NULL
		filtered.s = NULL
	}
	
	forcseries = matrix(NA, ncol = n.ahead, nrow = n.bootpred)
	forcsigma  = matrix(NA, ncol = n.ahead, nrow = n.bootpred)
	sim = vector(mode = "list", length = 1)
	sim = ugarchsim(fit, n.sim = n.ahead, m.sim = n.bootpred, 
				presigma = st, prereturns = tail(xdata[1:N], m), preresiduals = tail(residuals(fit), m),
				rseed = sseed, n.start = 0, startMethod = "sample", custom.dist = list(name = "sample", distfit = empz),
				mexsimdata = external.forecasts$mregfor, vexsimdata = external.forecasts$vregfor)
	# we transpose to get n.boot x n.ahead
	forcseries = t(sim@simulation$seriesSim)
	forcsigma =  t(sim@simulation$sigmaSim)
	
	# now we have the bootstrapped distribution of n.ahead forecast values
	# original forecast
	forc = ugarchforecast(fitORspec = fit, n.ahead = n.ahead, n.roll = 0, external.forecasts = external.forecasts) 
	coefdist = as.data.frame(coef(fit))
	model$truecoef = coef(fit)
	model$modeldata$realized.x = realized.x
	model$modeldata$filtered.s = filtered.s
	model$n.ahead = n.ahead
	model$n.bootfit = n.bootfit
	model$n.bootpred = n.bootpred
	model$type = "partial"
	ans = new("uGARCHboot",
			fseries = (forcseries),
			fsigma = (forcsigma),
			bcoef = coefdist,
			model = model,
			forc = forc)
	return(ans)
}

# using spec method
.ub1p2 = function(fitORspec, data = NULL, n.ahead = 10, n.bootfit = 100, 
		n.bootpred = 500, out.sample = 0, rseed = NA, solver.control = list(), fit.control = list(), 
		external.forecasts =  list(mregfor = NULL, vregfor = NULL), parallel = FALSE, 
		parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2))
{
	spec = fitORspec
	if(is.null(data))
		stop("\nugarchboot-->error: data must be supplied if SPEC object used.", call. = FALSE)
	model = spec@model
	flt = ugarchfilter(data = data, spec = spec, out.sample = out.sample)
	tmp = .extractdata(data)
	dates = tmp$dates
	xdata = tmp$data
	dformat = tmp$dformat
	ns = out.sample
	N = length(xdata) - out.sample
	sigma = flt@filter$sigma
	m = spec@model$maxOrder
	
	if (is.na(rseed[1])){
		sseed = as.integer(runif(n.bootpred,0,65000))
	} else{
		if(length(rseed) < n.bootpred){
			stop("\nugarchboot-->error: seed length must equal n.bootpred for partial method\n", call. = FALSE)
		} else {
			sseed = rseed
		}
	}
	# generate path based forecast values
	empz = matrix(0, ncol = n.bootpred, nrow = n.ahead)
	fz = flt@filter$z
	empz = matrix(apply(as.data.frame(1:n.bootpred), 1, FUN = function(i) {
						set.seed(sseed[i]);
						sample(fz, size = n.ahead, replace = TRUE);}), ncol = n.bootpred)
	
	if(ns > 0) {
		realized.x = xdata[(N+1):(N+ns)]
		filtered.s = tail(sigma(ugarchfilter(data = data, spec = spec)), ns)
	} else{
		realized.x = NULL
		filtered.s = NULL
	}
	
	# we start all forecasts from the last value of sigma based on the original series
	#st = .sigmat(spec, data, m)
	
	
	forcseries = matrix(NA, ncol = n.ahead, nrow = n.bootpred)
	forcsigma  = matrix(NA, ncol = n.ahead, nrow = n.bootpred)
	sim = ugarchpath(spec, n.sim = n.ahead, m.sim = n.bootpred, presigma = tail(sigma, m),
				prereturns = tail(xdata[1:N], m), preresiduals = tail(flt@filter$residuals, m),
				rseed = sseed, n.start = 0, custom.dist = list(name = "sample", 
						distfit = as.matrix(empz)), mexsimdata = external.forecasts$mregfor, 
				vexsimdata = external.forecasts$vregfor)
	# we transpose to get n.boot x n.ahead
	forcseries = t(sim@path$seriesSim)
	forcsigma =  t(sim@path$sigmaSim)
	
	# now we have the bootstrapped distribution of n.ahead forecast values
	# original forecast
	forc = ugarchforecast(fitORspec = spec, data = head(data, N), n.ahead = n.ahead, n.roll = 0, external.forecasts = external.forecasts) 
	
	model$truecoef = coef(flt)
	model$modeldata$realized.x = realized.x
	model$modeldata$filtered.s = filtered.s
	model$n.ahead = n.ahead
	model$n.bootfit = n.bootfit
	model$n.bootpred = n.bootpred
	coefdist = data.frame(NULL)
	model$type = "partial"
	ans = new("uGARCHboot",
			fseries = (forcseries),
			fsigma = (forcsigma),
			bcoef = coefdist,
			model = model,
			forc = forc)
	return(ans)
}



.sigmat = function(spec, origdata, m)
{
	flt = ugarchfilter(data = origdata, spec = spec)
	st = tail(flt@filter$sigma, m)
	return(st)
}

.quicksimulate = function(fit, n.sim, m.sim, presigma = NA, prereturns = NA,  n.start = 0, 
		rseed = NA, custom.dist = list(name = "sample", distfit = NULL), mexsimdata = NULL, vexsimdata = NULL)
{
	ans = ugarchsim(fit = fit, n.sim = n.sim, m.sim = m.sim, presigma = presigma, prereturns = prereturns,  
			n.start = n.start, rseed = rseed, custom.dist = custom.dist, mexsimdata = mexsimdata, 
			vexsimdata = vexsimdata)
	ret = rbind(ans@simulation$seriesSim, ans@simulation$sigmaSim)
	return(ret)
}