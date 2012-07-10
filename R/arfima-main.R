#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008, 2009, 2010, 2011, 
##	 2012
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

#---------------------------------------------------------------------------------
# SECTION sGARCH fit
#---------------------------------------------------------------------------------
.arfimafit = function(spec, data, out.sample = 0, solver = "solnp", solver.control = list(), 
		fit.control = list(fixed.se = 0, scale = 0))
{
	tic = Sys.time()
	if(is.null(solver.control$trace)) trace = 0 else trace = solver.control$trace
	
	if(is.null(fit.control$stationarity)) fit.control$stationarity = TRUE
	if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
	if(is.null(fit.control$scale)) fit.control$scale = FALSE
	# if we have external regressors in variance turn off scaling
	if(spec@model$modelinc[15] > 0) fit.control$scale = FALSE
	# if we have arch-in-mena turn off scaling
	if(spec@model$modelinc[5] > 0) fit.control$scale = FALSE
	
	# if there are fixed pars we do no allow scaling as there would be no way of mixing scaled
	# amd non scaled parameters	
	if(sum(spec@model$pars[,2]) > 0) fit.control$scale = FALSE
	xdata = .extractdata(data)
	if(!is.numeric(out.sample)) stop("\narfimafit-->error: out.sample must be numeric\n")
	if(as.numeric(out.sample)<0) stop("\narfimafit-->error: out.sample must be positive\n")
	n.start = round(out.sample,0)
	n = length(xdata$data)
	#if((n-n.start)<100) stop("\narfimafit-->error: function requires at least 100 data\n points to run\n")
	data = xdata$data[1:(n-n.start)]
	dates = xdata$pos[1:(n-n.start)]
	origdata = xdata$data
	origdates = xdata$pos
	dformat = xdata$dformat
	# create a temporary environment to store values (deleted at end of function)
	tempenvir = new.env(hash = TRUE, parent = .GlobalEnv)
	model = spec@model
	modelinc = model$modelinc
	pidx = model$pidx
	# expand the spec object and assign spec lists
	if(modelinc[6] > 0){
		mexdata = model$modeldata$mexdata[1:(n-n.start), , drop = FALSE]
	} else{
		mexdata = NULL
	}
	
	assign("dates", dates, envir = tempenvir)
	assign("trace", trace, envir = tempenvir)
	assign("fit.control", fit.control, envir = tempenvir)
	
	m =  model$maxOrder
	model$modeldata$T = T = length(as.numeric(data))
	dist = model$modeldesc$distribution
	if(fit.control$scale) dscale = sd(data) else dscale = 1
	zdata = data/dscale
	assign("dscale", dscale, envir = tempenvir)
	assign("model", model, envir = tempenvir)
	ipars = model$pars
	# Optimization Starting Parameters Vector & Bounds
	ipars = .arfimastart(ipars, data = zdata, garchenv = tempenvir)
	assign("ipars", ipars, envir = tempenvir)
	# we now split out any fixed parameters
	estidx = as.logical( ipars[,4] )
	assign("estidx", estidx, envir = tempenvir)
	
	npars = sum(estidx)
	if(any(ipars[,2]==1)){
		if(npars == 0){
			if(fit.control$fixed.se==0) {
				# if all parameters are fixed an no standard erros are to
				# be calculated then we return a ugarchfilter object
				warning("\narfimafit-->warning: all parameters fixed...returning ugarchfilter 
								object instead\n")
				return(arfimafilter(data = data, spec = spec, out.sample = out.sample))
			} else{
				# if all parameters are fixed but we require standard errors, we
				# skip the solver
				use.solver = 0
				ipars[ipars[,2]==1, 4] = 1
				ipars[ipars[,2]==1, 2] = 0
				assign("ipars", ipars, envir = tempenvir)
				estidx = as.logical( ipars[,4] )
				assign("estidx", estidx, envir = tempenvir)
			}
		} else{
			# with some parameters fixed we extract them (to be rejoined at end)
			# so that they do not enter the solver
			use.solver = 1
		}
	} else{
		use.solver = 1
	}
	
	# start counter
	assign(".llh", 1, envir = tempenvir)
	
	if(fit.control$stationarity==1 && modelinc[2]>0){
		Ifn = function(pars, data, returnType, garchenv){
			arx = 0.5
			modelinc = get("model", garchenv)$modelinc
			idx = get("model", garchenv)$pidx
			ipars = get("ipars", garchenv)
			estidx = get("estidx", garchenv)
			ipars[estidx, 1] = pars		
			if(modelinc[2] > 0 && modelinc[4] == 0) arx = ipars[idx["ar",1]:idx["ar",2],1]
			min(Mod(polyroot(c(1, -arx))))
		}
		ILB = 1.01
		IUB = 100
		if(solver == "solnp" | solver == "gosolnp") fit.control$stationarity = 0
	} else{
		Ifn = ILB = IUB = NULL
	}
	assign("fit.control", fit.control, , envir = tempenvir)
	
	if(use.solver){
		parscale = rep(1, length = npars)
		names(parscale) = rownames(ipars[estidx,])
		if(modelinc[1] > 0) parscale["mu"] = abs(mean(zdata))
		if(modelinc[7] > 0) parscale["sigma"] = sd(zdata)
		solution = .garchsolver(solver, pars = ipars[estidx, 1], fun = .arfimaLLH, Ifn, ILB, IUB, 
				gr = NULL, hessian = NULL, parscale = parscale, 
				control = solver.control, LB = ipars[estidx, 5], UB = ipars[estidx, 6], 
				ux = NULL, ci = NULL, mu = NULL, data = zdata, returnType = "llh", garchenv = tempenvir)
		sol = solution$sol
		hess = solution$hess
		timer = Sys.time()-tic
		ipars = get("ipars", tempenvir)
		if(!is.null(sol$par)) ipars[estidx, 1] = sol$par else ipars[estidx, 1] = NA
		if(sum(ipars[,2]) == 0){
			if(modelinc[1] > 0) ipars[pidx["mu",1]:pidx["mu",2], 1] = ipars[pidx["mu",1]:pidx["mu",2], 1] * dscale
			if(modelinc[5] > 0){
				ipars[pidx["mxreg", 1]:pidx["mxreg", 2], 1] = ipars[pidx["mxreg", 1]:pidx["mxreg", 2], 1] * dscale
			}
			ipars[pidx["sigma",1], 1] = ipars[pidx["sigma",1],1] * dscale
		}
		assign("ipars", ipars, envir = tempenvir)
		convergence = sol$convergence
	} else{
		hess = NULL
		timer = Sys.time()-tic
		convergence = 0
		sol = list()
		sol$message = "all parameters fixed"
	}
	fit = list()
	# check convergence else write message/return
	ipars2 = ipars
	if(convergence == 0){
		if(sum(ipars[,2]) > 0 && fit.control$fixed.se == 1){
			ipars[ipars[,2]==1, 4] = 1
			ipars[ipars[,2]==1, 2] = 0
			assign("ipars", ipars, envir = tempenvir)
			estidx = as.logical( ipars[,4] )
			assign("estidx", estidx, envir = tempenvir)
		}
		fit = .makearfimafitmodel(f = .arfimaLLH, data = data,  T = T, m = m, timer = timer, 
				convergence = convergence, message = sol$message, hess, garchenv = tempenvir)
		model$modelinc[7] = modelinc[7]
		model$modeldata$data = origdata
		model$modeldata$dates = origdates
		model$modeldata$date.format = dformat
		
		model$pars = ipars
		model$pars[, 1] = fit$ipars[,1]
		fit$ipars[, 4] = ipars2[, 4]
		fit$ipars[, 2] = ipars2[, 2]
	} else{
		fit$message = sol$message
		fit$convergence = 1
		model$modeldata$data = origdata
		model$modeldata$dates = origdates
		model$modeldata$date.format = dformat
	}
	
	# make model list to return some usefule information which
	# will be called by other functions (show, plot, sim etc)
	model$n.start = n.start

	ans = new("ARFIMAfit",
			fit = fit,
			model = model)
	rm(tempenvir)
	return(ans)
}

#---------------------------------------------------------------------------------
# SECTION sGARCH LLH
#---------------------------------------------------------------------------------
.arfimaLLH = function(pars, data, returnType = "llh", garchenv)
{
	# prepare inputs
	# rejoin fixed and pars
	model = get("model", garchenv)
	estidx = get("estidx", garchenv)
	idx = model$pidx
	ipars = get("ipars", garchenv)
	ipars[estidx, 1] = pars
	trace = get("trace", garchenv)
	T = length(data)
	
	fit.control = get("fit.control", garchenv)
	m = model$maxOrder
	N = c(m,T)
	mexdata = model$modeldata$mexdata[1:T,, drop = FALSE]
	distribution = model$modeldesc$distribution
	modelinc = model$modelinc
	# distribution number
	# 1 = norm, 2=snorm, 3=std, 4=sstd, 5=ged, 6=sged, 7=nig
	dist = model$modeldesc$distno
	hm = 0
	if(modelinc[4]>0){
		rx = .arfimaxfilter(modelinc, ipars[,1], idx, mexdata = mexdata, h = 0, data = data, N = N, garchenv)
		res = rx$res
		zrf = rx$zrf
		res[is.na(res) | !is.finite(res) | is.nan(res)] = 0
	} else{
		res = rep(0, T)
		zrf = 0
	}
	# unconditional sigma value
	assign("ipars", ipars, envir = garchenv)
	if(fit.control$stationarity == 1 && modelinc[4] == 0 && modelinc[2] > 0){
		arx = ipars[idx["ar",1]:idx["ar",2],1]
		kappa = min(Mod(polyroot(c(1, -arx))))
		if(!is.na(kappa) && (kappa < 1 | kappa > 100) ) return(llh = get(".llh", garchenv) + 0.1*(abs(get(".llh", garchenv))))
	}
	if(modelinc[6]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)

	ans = try( .C("arfimafitC", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					x = as.double(data), res = as.double(res), mexdata = mexdata, zrf = as.double(zrf),
					constm = double(T), condm = double(T), m = as.integer(m), T = as.integer(T),
					z = double(T), llh = double(1), LHT = double(T), PACKAGE = "rugarch"), silent = TRUE )
	
	if(inherits(ans, "try-error")){
		assign(".csol", 1, envir = garchenv)
		assign(".filtermessage", ans, envir = garchenv)
		if( trace > 0 ) cat(paste("\narfimafit-->warning: ", get(".filtermessage",garchenv),"\n", sep=""))
		return(llh = (get(".llh", garchenv) + 0.1*(abs(get(".llh", garchenv)))))
	} else{
		assign(".csol", 0, envir = garchenv)
	}
	z = ans$z
	epsx = ans$res
	llh = ans$llh
	
	if(is.finite(llh) && !is.na(llh) && !is.nan(llh)){
		assign(".llh", llh, envir = garchenv)
	} else {
		llh = (get(".llh", garchenv) + 0.1*(abs(get(".llh", garchenv))))
	}
	
	# LHT = raw scores
	#LHT = -ans$LHT[(m):T]
	LHT = -ans$LHT
	ans = switch(returnType,
			llh = llh,
			LHT = LHT,
			all = list(llh = llh, res = epsx, z = z, LHT = LHT))
	return(ans)
}

#---------------------------------------------------------------------------------
# SECTION sGARCH filter
#---------------------------------------------------------------------------------
.arfimafilter = function(spec, data, out.sample = 0, n.old = NULL)
{
	# n.old is optional and indicates the length of the original dataseries (in
	# cases when this represents a dataseries augmented by newer data). The reason
	# for using this is so that the old and new datasets agree since the original
	# recursion uses the sum of the residuals to start the recursion and therefore
	# is influenced by new data. For a small augmentation the values converge after
	# x periods, but it is sometimes preferable to have this option so that there is
	# no forward looking information contaminating the study.
	.garchenv = environment()
	xdata = .extractdata(data)
	data = xdata$data
	dates = xdata$pos
	dformat = xdata$dformat
	origdata = data
	origdates = dates
	T = length(origdata)  - out.sample
	data = origdata[1:T]
	dates = origdates[1:T]
	if(!is.null(n.old)) Nx = n.old else Nx = length(data)
	
	model = spec@model
	ipars = model$pars
	pars = unlist(model$fixed.pars)
	parnames = names(pars)
	modelnames = .checkallfixed(spec)
	if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
		cat("\narfimafilter-->error: parameters names do not match specification\n")
		cat("Expected Parameters are: ")
		cat(paste(modelnames))
		cat("\n")
		stop("Exiting", call. = FALSE)
	}
	
	# once more into the spec
	setfixed(spec)<-as.list(pars)
	model = spec@model
	model$modeldata$T = T
	ipars = model$pars
	idx = model$pidx
	modelinc = model$modelinc
	m = model$maxOrder
	N = c(m,T)
	mexdata = model$modeldata$mexdata[1:T, , drop = FALSE]
	distribution = model$modeldesc$distribution
	dist = model$modeldesc$distno
	
	
	if(modelinc[4]>0){
		rx = .arfimaxfilter(modelinc, ipars[,1], idx, mexdata = mexdata, h = 0, data = data, N = N, .garchenv)
		res = rx$res
		zrf = rx$zrf
		res[is.na(res) | !is.finite(res) | is.nan(res)] = 0
	} else{
		res = rep(0, T)
		zrf = 0
	}
	
	if(modelinc[6]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	
	ans = try( .C("arfimafitC", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					x = as.double(data), res = as.double(res), mexdata = mexdata, zrf = as.double(zrf),
					constm = double(T), condm = double(T), m = as.integer(m), T = as.integer(T),
					z = double(T), llh = double(1), LHT = double(T), PACKAGE = "rugarch"), silent = TRUE )
	
	if(inherits(ans, "try-error")) stop("\nugarchfilter-->error: problem in C code...exiting.")
	filter = list()
	filter$z = ans$z
	filter$residuals = ans$res
	filter$LLH = -ans$llh
	filter$log.likelihoods = ans$LHT
	filter$distribution = distribution
	filter$ipars = ipars
	model$modeldata$data = origdata
	model$modeldata$dates = origdates
	model$modeldata$date.format = dformat
	model$n.start = out.sample
	
	rm(.garchenv)
	sol = new("ARFIMAfilter",
			filter = filter,
			model = model)
	return(sol)
}

#---------------------------------------------------------------------------------
# SECTION sGARCH forecast
#---------------------------------------------------------------------------------
.arfimaforecast = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0, 
		external.forecasts = list(mregfor = NULL), ...)
{
	fit = fitORspec
	data = fit@model$modeldata$data
	Nor = length(as.numeric(data))
	dates = fit@model$modeldata$dates
	dformat = fit@model$modeldata$date.format
	ns = fit@model$n.start
	N = Nor - ns
	model = fit@model
	ipars = fit@fit$ipars
	modelinc = model$modelinc
	idx = model$pidx
	if( n.roll > ns ) stop("\nugarchforecast-->error: n.roll must not be greater than out.sample!")
	pars = fit@fit$coef
	ipars = fit@fit$ipars
	
	# check if necessary the external regressor forecasts provided first
	xreg = .forcregressors(model, external.forecasts$mregfor, NULL, ipars, n.ahead, Nor, out.sample = ns, n.roll)
	mxf = xreg$mxf
	
	# filter data (check external regressor data - must equal length of origData)
	fcreq = ifelse(ns >= (n.ahead+n.roll), n.ahead+n.roll, ns)
	fspec = arfimaspec(
			mean.model = list(armaOrder = c(modelinc[2], modelinc[3]),
					include.mean = modelinc[1], arfima = modelinc[4], 
					external.regressors = mxf[1:(N + fcreq), , drop = FALSE]), 
			distribution.model = model$modeldesc$distribution, fixed.pars = as.list(pars))
	tmp =  data.frame(data[1:(N + fcreq)])
	rownames(tmp) = as.character(dates[1:(N + fcreq)])
	flt = .arfimafilter(data = tmp[,1,drop=FALSE], spec = fspec, n.old = N)
	resfilter = flt@filter$residuals
	zfilter = flt@filter$z
	
	
	# forecast ARFIMA process
	forecasts = vector(mode = "list", length = n.roll+1)
	fwdd = vector(mode="list", length = n.roll+1)
	for(i in 1:(n.roll+1)){
		np = N + i - 1
		if(modelinc[1] > 0){
			mu = rep(ipars[idx["mu",1]:idx["mu",2], 1], N+i+n.ahead-1)
		} else{
			mu = rep(0, N+i+n.ahead-1)
		}
		epsx = c(resfilter[1:(N+i-1)], rep(0, n.ahead))
		x = c(data[1:(N+i-1)], rep(0, n.ahead))
		z = c(zfilter[1:(N+i-1)], rep(0, n.ahead))
		# forecast of externals is provided outside the system
		mxfi = mxf[1:(N+i-1+n.ahead), , drop = FALSE]
		
		if(modelinc[4]>0){
			res = arfimaf(ipars, modelinc, idx, mu, mxfi, h=0, epsx, z, data = x, N = np, n.ahead)
		} else{
			res = armaf(ipars, modelinc, idx, mu, mxfi, h=0, epsx, z, data = x, N = np, n.ahead)
		}
		ans = res[(np + 1):(np + n.ahead)]
		fdf = data.frame(series = ans)
		fwdd[[i]] = .forcdates( dates, n.ahead, N, i, ns, dformat )
		rownames(fdf) = as.character(fwdd[[i]])
		forecasts[[i]] = fdf
	}
	fcst = list()
	fcst$n.ahead = n.ahead
	fcst$N = N+ns
	fcst$n.start = ns
	fcst$n.roll = n.roll
	fcst$forecasts = forecasts
	fcst$fdates = fwdd
	model$modeldata$residuals = flt@filter$residuals
	ans = new("ARFIMAforecast",
			forecast = fcst,
			model = model)
	return(ans)
}

#---------------------------------------------------------------------------------
# 2nd dispatch method for forecast
.arfimaforecast2 = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0, 
		external.forecasts = list(mregfor = NULL), ...)
{
	# first we filter the data to get the results:
	spec = fitORspec
	if(is.null(data)) stop("\nugarchforecast-->error: data must not be NULL when using a specification!")
	# we then extract the data/coefs etc
	xdata = .extractdata(data)
	Nor = length(as.numeric(xdata$data))
	data = xdata$data
	N = length(as.numeric(data))
	dates = xdata$pos
	dformat = xdata$dformat	
	ns = out.sample
	if( n.roll > ns ) stop("\nugarchforecast-->error: n.roll must not be greater than out.sample!")
	N = Nor - ns
	
	model = spec@model
	ipars = model$pars
	pars = unlist(model$fixed.pars)
	parnames = names(pars)
	modelnames = .checkallfixed(spec)
	if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
		cat("\nugarchforecast-->error: parameters names do not match specification\n")
		cat("Expected Parameters are: ")
		cat(paste(modelnames))
		cat("\n")
		stop("Exiting", call. = FALSE)
	}
	# once more into the spec
	setfixed(spec)<-as.list(pars)
	model = spec@model
	idx = model$pidx
	ipars = model$pars
	modelinc = model$modelinc
	model$modeldata$data = data
	model$modeldata$dates = dates
	model$modeldata$dformat = dformat
	
	# check if necessary the external regressor forecasts provided first
	xreg = .forcregressors(model, external.forecasts$mregfor, NULL, ipars, n.ahead, Nor, out.sample = ns, n.roll)
	mxf = xreg$mxf
	vxf = xreg$vxf
	
	# filter data (check external regressor data - must equal length of origData)
	fcreq = ifelse(ns >= (n.ahead+n.roll), n.ahead+n.roll, ns)
	fspec = arfimaspec(
			mean.model = list(armaOrder = c(modelinc[2], modelinc[3]),
					include.mean = modelinc[1], arfima = modelinc[4], 
					external.regressors = mxf[1:(N + fcreq), , drop = FALSE]), 
			distribution.model = model$modeldesc$distribution, fixed.pars = as.list(pars))
	tmp =  data.frame(data[1:(N + fcreq)])
	rownames(tmp) = as.character(dates[1:(N + fcreq)])
	flt = .arfimafilter(data = tmp[,1,drop=FALSE], spec = fspec, n.old = N)
	resfilter = flt@filter$residuals
	zfilter = flt@filter$z
	# forecast GARCH process
	forecasts = vector(mode="list", length = n.roll+1)
	fwdd = vector(mode="list", length = n.roll+1)
	
	for(i in 1:(n.roll+1)){
		np = N + i - 1
		if(modelinc[1] > 0){
			mu = rep(ipars[idx["mu",1]:idx["mu",2], 1], N+i+n.ahead-1)
		} else{
			mu = rep(0, N+i+n.ahead-1)
		}
		epsx = c(resfilter[1:(N+i-1)], rep(0, n.ahead))
		x = c(data[1:(N+i-1)], rep(0, n.ahead))
		z = c(zfilter[1:(N+i-1)], rep(0, n.ahead))
		# forecast of externals is provided outside the system
		mxfi = mxf[1:(N+i-1+n.ahead), , drop = FALSE]
		
		if(modelinc[4]>0){
			res = arfimaf(ipars, modelinc, idx, mu, mxfi, h=0, epsx, z, data = x, N = np, n.ahead)
		} else{
			res = armaf(ipars, modelinc, idx, mu, mxfi, h=0, epsx, z, data = x, N = np, n.ahead)
		}
		ans = res[(np + 1):(np + n.ahead)]
		fdf = data.frame(series = ans)
		fwdd[[i]] = .forcdates( dates, n.ahead, N, i, ns, dformat )
		rownames(fdf) = as.character(fwdd[[i]])
		forecasts[[i]] = fdf
	}
	fcst = list()
	fcst$n.ahead = n.ahead
	fcst$N = N+ns
	fcst$n.start = ns
	fcst$n.roll = n.roll
	fcst$forecasts = forecasts
	fcst$fdates = fwdd
	model$modeldata$residuals = flt@filter$residuals
	ans = new("ARFIMAforecast",
			forecast = fcst,
			model = model)
	return(ans)
	
}
#---------------------------------------------------------------------------------
# SECTION sGARCH simulate
#---------------------------------------------------------------------------------

.arfimasim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, startMethod = 
				c("unconditional","sample"), prereturns = NA, 
		preresiduals = NA, rseed = NA, custom.dist = list(name = NA, distfit = NA, type = "z"), 
		mexsimdata = NULL)
{
	# some checks
	if(is.na(rseed[1])){
		sseed = as.integer(runif(1,0,as.integer(Sys.time())))
	} else{
		if(length(rseed) != m.sim) sseed = as.integer(rseed[1]) else sseed = rseed[1:m.sim]
	}
	
	model = fit@model
	modelinc = model$modelinc
	# Enlarge Series:
	# need to allow for arfima case:
	if(modelinc[4]>0) n.start = max(modelinc[2]+modelinc[3], n.start)
	n = n.sim + n.start
	startMethod = startMethod[1]
	data = fit@model$modeldata$data
	N = length(as.numeric(data))
	data = data[1:(N - fit@model$n.start)]
	N = length(as.numeric(data))
	m = fit@model$maxOrder
	resids = fit@fit$residuals

	idx = model$pidx
	ipars = fit@fit$ipars
	# check if necessary the external regressor forecasts provided first
	xreg = .simregressors(model, mexsimdata, NULL, fit@model$ipars, N, n, m.sim, m)	
	mexsim = xreg$mexsimlist
	Sigma = ipars[idx["sigma", 1], 1]
	if(N < n.start){
		startmethod[1] = "unconditional"
		warning("\narfimasim-->warning: n.start greater than length of data...using unconditional start method...\n")
	}
	
	# Random Samples from the Distribution
	# For the arfima method this is the residuals, NOT the standardized residuals
	if(length(sseed) == 1){
		zmatrix = data.frame(dist = model$modeldesc$distribution, lambda = ipars[idx["ghlambda",1], 1], 
				skew = ipars[idx["skew",1], 1], shape = ipars[idx["shape",1], 1], n = n * m.sim, seed = sseed[1])
		z = .custzdist(custom.dist, zmatrix, m.sim, n)
	} else{
		zmatrix = data.frame(dist = rep(model$modeldesc$distribution, m.sim), lambda = rep(ipars[idx["ghlambda",1], 1], m.sim), 
				skew = rep(ipars[idx["skew",1], 1], m.sim), shape = rep(ipars[idx["shape",1], 1], m.sim), 
				n = rep(n, m.sim), seed = sseed)
		z = .custzdist(custom.dist, zmatrix, m.sim, n)
	}
	if( is.matrix( custom.dist$distfit ) && custom.dist$type != "z") cres = TRUE else cres = FALSE
	
	if(startMethod == "unconditional"){
		z = rbind(matrix(0, nrow = m, ncol = m.sim), z)
	} else{
		z = rbind(matrix(tail(fit@fit$z, m), nrow = m, ncol = m.sim), z) 
	}
	
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<m) stop(paste("\narfimasim-->error: prereturns must be of length ", m, sep=""))
	}
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<m) stop(paste("\narfimasim-->error: preresiduals must be of length ", m, sep=""))
		preres = matrix(preresiduals, nrow = m)
	}
	if(is.na(prereturns[1])){
		if(startMethod[1] == "unconditional"){
			prereturns = as.numeric(rep(uncmean(fit), m))
		}
		else{
			prereturns = tail(data, m)
		}
	}
	
	# input vectors/matrices
	x = c(prereturns, rep(0, n))
	constm = matrix(ipars[idx["mu",1]:idx["mu",2], 1], ncol = m.sim, nrow = n + m)
	
	# outpus matrices
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	residSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	
	for(i in 1:m.sim){
		if(any(is.na(preresiduals))){
			if(startMethod[1] == "sample"){
				preres 	= tail(resids, m)
				res 	= c(preres, if(cres) tail(z[(m+1):(n+m), i], n) else tail(z[(m+1):(n+m), i], n) * Sigma)
				residSim[,i] = tail(res, n.sim)
			} else{
				res 	= if(cres) as.numeric(z[,i]) else as.numeric(z[,i])*Sigma
				residSim[,i] = tail(res, n.sim)
			}
		}
		if(modelinc[6]>0){
			mxreg = matrix( ipars[idx["mxreg",1]:idx["mxreg",2], 1], ncol = modelinc[6] )
			constm[,i] = constm[,i] + mxreg %*%t( matrix( mexsim[[i]], ncol = modelinc[6] ) )
		}
		
		if(modelinc[4]>0){
			x = c(prereturns, rep(0, n +  modelinc[3]))
			res = c(if(cres) as.numeric(z[(m+1):(n+m),i]) else as.numeric(z[(m+1):(n+m),i]) * Sigma, if(modelinc[3]>0) rep(0,  modelinc[3]) else NULL)
			ans2 = .arfimaxsim(modelinc, ipars, idx, constm[1:n, i], res[1:(n+modelinc[3])], T = n)	
			seriesSim[,i] = tail(ans2$series, n.sim)
		} else{
			ans2 = .armaxsim(modelinc, ipars, idx, constm[,i],  x, res, T = n + m, m)
			seriesSim[,i] = tail(ans2$x, n.sim)
		}
	}
	sim = list(seriesSim = seriesSim, residSim = residSim)
	sol = new("ARFIMAsim",
			simulation = sim,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}


#---------------------------------------------------------------------------------
# SECTION sGARCH path
#---------------------------------------------------------------------------------

.arfimapath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, prereturns = NA, preresiduals = NA, 
		rseed = NA, custom.dist = list(name = NA, distfit = NA, type = "z"), mexsimdata = NULL)
{
	# some checks
	if(is.na(rseed[1])){
		sseed = as.integer(runif(1,0,as.integer(Sys.time())))
	} else{
		if(length(rseed) != m.sim) sseed = as.integer(rseed[1]) else sseed = rseed[1:m.sim]
	}
	model = spec@model
	ipars = model$pars
	pars = unlist(model$fixed.pars)
	parnames = names(pars)
	modelnames = .checkallfixed(spec)
	if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
		cat("\narfimapath-->error: parameters names do not match specification\n")
		cat("Expected Parameters are: ")
		cat(paste(modelnames))
		cat("\n")
		stop("Exiting", call. = FALSE)
	}
	# once more into the spec
	setfixed(spec)<-as.list(pars)
	model = spec@model
	ipars = model$pars
	idx = model$pidx
	modelinc = model$modelinc
	Sigma = ipars[idx["sigma", 1], 1]
	# Enlarge Series:
	n = n.sim + n.start
	m = model$maxOrder
	N = 0
	if(modelinc[6]>0) {
		mexdata = matrix(model$modeldata$mexdata, ncol = modelinc[6])
		N = dim(mexdata)[1]
	} else { mexdata = NULL }
	distribution = model$modeldesc$distribution	
	# check if necessary the external regressor forecasts provided first
	xreg = .simregressors(model, mexsimdata, NULL, ipars, N, n, m.sim, m)	
	mexsim = xreg$mexsimlist
	
	if(length(sseed) == 1){
		zmatrix = data.frame(dist = model$modeldesc$distribution, lambda = ipars[idx["ghlambda",1], 1], 
				skew = ipars[idx["skew",1], 1], shape = ipars[idx["shape",1], 1], n = n * m.sim, seed = sseed[1])
		z = .custzdist(custom.dist, zmatrix, m.sim, n)
	} else{
		zmatrix = data.frame(dist = rep(model$modeldesc$distribution, m.sim), lambda = rep(ipars[idx["ghlambda",1], 1], m.sim), 
				skew = rep(ipars[idx["skew",1], 1], m.sim), shape = rep(ipars[idx["shape",1], 1], m.sim), 
				n = rep(n, m.sim), seed = sseed)
		z = .custzdist(custom.dist, zmatrix, m.sim, n)
	}
	z = rbind(matrix(0, nrow = m, ncol = m.sim), z)
	
	if( is.matrix( custom.dist$distfit ) && custom.dist$type != "z") cres = TRUE else cres = FALSE
	
		
	# create the presample information
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<m) stop(paste("\narfimapath-->error: prereturns must be of length ", m, sep=""))
	}
	
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<m) stop(paste("\narfimapath-->error: preresiduals must be of length ", m, sep=""))
		preres = matrix(preresiduals, nrow = m)
	}
	if(is.na(prereturns[1])){
		prereturns = as.numeric(rep(uncmean(spec), times = m))
	}
	
	
	# input vectors/matrices
	x = c(prereturns, rep(0, n))
	constm = matrix(ipars[idx["mu",1]:idx["mu",2],1], ncol = m.sim, nrow = n + m)

	# outpus matrices
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	residSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	
	for(i in 1:m.sim){
		if(any(is.na(preresiduals))){
			preres = if(cres) as.numeric(z[1:m,i]) else as.numeric(z[1:m,i])*Sigma
		} 
		res = c(preres, if(cres) as.numeric(z[(m+1):(n+m),i]) else as.numeric(z[(m+1):(n+m),i]) * Sigma)
		residSim[,i] = tail(res, n.sim)
		
		if(modelinc[6]>0){
			mxreg = matrix( ipars[idx["mxreg",1]:idx["mxreg",2], 1], ncol = modelinc[6] )
			constm[,i] = constm[,i] + mxreg %*%t( matrix( mexsim[[i]], ncol = modelinc[6] ) )
		}
		if(modelinc[4]>0){
			x = c(prereturns, rep(0, n + modelinc[3]))
			res = c(if(cres) as.numeric(z[(m+1):(n+m),i]) else as.numeric(z[(m+1):(n+m),i]) * Sigma, if(modelinc[3]>0) rep(0, modelinc[3]) else NULL)
			ans2 = .arfimaxsim(modelinc, ipars, idx, constm[1:n, i], res[1:(n+modelinc[3])], T = n)
			seriesSim[,i] = tail(ans2$series, n.sim)
		} else{
			ans2 = .armaxsim(modelinc, ipars, idx, constm[,i],  x, res, T = n + m, m)
			seriesSim[,i] = tail(ans2$x, n.sim)
		}
	}
	path = list(seriesSim = seriesSim, residSim = residSim)
	sol = new("ARFIMApath",
			path = path,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}

.arfimaroll = function(spec,  data, n.ahead = 1, forecast.length = 500, refit.every = 25, 
		refit.window = c("recursive", "moving"), parallel = FALSE, 
		parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2), solver = "solnp", 
		fit.control = list(), solver.control = list(), calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), ...)
{

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
	}
	datanames = names(data)
	xdata = .extractdata(data)
	data = xdata$data
	dates = xdata$pos
	xdata = data.frame(data)
	rownames(xdata) = as.character(dates)
	if(is.null(fit.control$stationarity)) fit.control$stationarity=1
	if(is.null(fit.control$fixed.se)) fit.control$fixed.se=0
	T = length(data)
	startT = T-forecast.length
	# forecast window index
	fwindex = t( .embed((T - forecast.length - n.ahead + 2):T, refit.every, by = refit.every, TRUE ) )
	#forecast.length = dim(fwindex)[1] * dim(fwindex)[2]
	fitindx.start = c(fwindex[1,]-1)
	fwindex.end = fwindex[refit.every,]
	nf = length(fwindex.end)
	fitdata = vector(mode="list",length = nf)
	mexdata = vector(mode="list",length = nf)
	fitlist = vector(mode="list",length = nf)
	outdata = vector(mode="list",length = nf)
	coefv = vector(mode="list", length = nf)
	distribution = spec@model$modeldesc$distribution
	model = spec@model
	modelinc = model$modelinc
	forecastlist = vector(mode="list", length = nf)
	forecast.length = floor(forecast.length/refit.every) * refit.every
	VaR.list = NULL
	
	for(i in 1:nf){
		if(refit.window[1]=="recursive"){
			fitdata[[i]] = xdata[1:(fitindx.start[i]+refit.every), , drop = FALSE]
			outdata[[i]] = xdata[(fitindx.start[i]+1):(fitindx.start[i]+refit.every), , drop = FALSE]
			if( modelinc[6]>0 ){
				mexdata[[i]] = model$modeldata$mexdata[1:(fitindx.start[i]+refit.every), , drop = FALSE]
			} else{
				mexdata[[i]] = NULL
			}
		} else{
			fitdata[[i]] = xdata[(i*refit.every):(fitindx.start[i]+refit.every), , drop = FALSE]
			outdata[[i]] = xdata[(fitindx.start[i]+1):(fitindx.start[i]+refit.every), , drop = FALSE]
			if( modelinc[6]>0 ){
				mexdata[[i]] = model$modeldata$mexdata[(i*refit.every):(fitindx.start[i]+refit.every), , drop = FALSE]
			} else{
				mexdata[[i]] = NULL
			}
		}
	}
	
	cat("\n...estimating refit windows...\n")
	specx = vector(mode = "list", length = nf)
	for(i in 1:nf){
		specx[[i]] = spec
		if(modelinc[6]>0)  specx[[i]]@model$modeldata$mexdata = mexdata[[i]] else specx[[i]]@model$modeldata$mexdata = NULL
	}
	if( parallel ){
		if( parallel.control$pkg == "multicore" ){
			if(!exists("mclapply")){
				require('multicore')
			}
			fitlist = multicore::mclapply(1:nf, FUN = function(i) arfimafit(spec = specx[[i]], data = fitdata[[i]], 
								solver = solver, out.sample = refit.every, solver.control = solver.control,
								fit.control = fit.control), mc.cores = parallel.control$cores)
		} else{
			if(!exists("sfLapply")){
				require('snowfall')
			}
			sfInit(parallel = TRUE, cpus = parallel.control$cores)
			sfExport("specx", "fitdata", "refit.every", "solver", "solver.control", "fit.control",local = TRUE)
			fitlist = sfLapply(as.list(1:nf), fun = function(i) rugarch::arfimafit(spec = specx[[i]], 
								data = fitdata[[i]], solver = solver, out.sample = refit.every, 
								solver.control = solver.control, fit.control = fit.control))
			sfStop()
		}
	} else{
		for(i in 1:nf){
			fitlist[[i]] = arfimafit(data = fitdata[[i]], spec = specx[[i]], solver = solver, 
					out.sample = refit.every, solver.control = solver.control,
					fit.control = fit.control)
			if(fitlist[[i]]@fit$convergence!=0){
				if(i>1){
					specx[[i]]@model$start.pars = as.list(coef(fitlist[[i-1]]))
					fitlist[[i]] = arfimafit(data = fitdata[[i]], spec = specx[[i]], solver = solver, 
							out.sample = refit.every, solver.control = solver.control,
							fit.control = fit.control)
					specx[[i]]@model$start.pars = NULL
				}
			}
		}
	}
	converge = sapply(fitlist, FUN=function(x) x@fit$convergence)
	if(any(converge!=0)){
		ncon = which(converge!=0)
		stop(paste("\nno convergence in the following fits: ", ncon, "...exiting", sep=""), call.=FALSE)
	}
	cat("\nDone!...all converged.\n")
	coefv = t(sapply(fitlist,FUN=function(x) x@fit$coef))
	coefmat = lapply(fitlist,FUN=function(x) x@fit$robust.matcoef)
	LLH = sapply(fitlist, FUN=function(x) x@fit$LLH)
	# filtered sigma/series(actual) estimates
	filter.sigma = vector(mode = "numeric", length = forecast.length)
	filter.series = unlist(outdata)
	
	if( parallel ){
		if( parallel.control$pkg == "multicore" ){
			forecastlist = multicore::mclapply(fitlist, FUN = function(x) arfimaforecast(x, 
								n.ahead = n.ahead, n.roll = refit.every-1), mc.cores = parallel.control$cores)
		} else{
			nx = length(fitlist)
			sfInit(parallel = TRUE, cpus = parallel.control$cores)
			sfExport("fitlist", "n.ahead", "refit.every", local = TRUE)
			forecastlist = sfLapply(as.list(1:nx), fun = function(i) rugarch::arfimaforecast(fitlist[[i]], 
								n.ahead = n.ahead, n.roll = refit.every-1))
			sfStop()
		}
	} else{
		forecastlist = lapply(fitlist, FUN = function(x) arfimaforecast(x, n.ahead = n.ahead, n.roll = refit.every-1))
	}
	
	eindex = t(.embed(1:forecast.length, refit.every, refit.every, TRUE))
	# collect forecasts [mu sigma (skew shape)]
	# 2 cases : n.ahead = 1 we use a matrix, else
	# n.ahead > 1 we use a list object
	if(n.ahead == 1){
		fser = as.numeric(sapply(forecastlist,FUN=function(x) sapply(x@forecast$forecasts, FUN=function(x) x)))
		fdates = dates[as.numeric(fwindex)]
		rdat = data[(fitindx.start[1]+1):(fitindx.start[1] + forecast.length + 1 - 1)]
		eindex = t(.embed(1:forecast.length, refit.every, refit.every, TRUE))
		if(modelinc[16]>0) skew.f = sapply(fitlist,FUN=function(x) coef(x)["skew"])
		if(modelinc[17]>0) shape.f = sapply(fitlist,FUN=function(x) coef(x)["shape"])
		if(modelinc[18]>0) ghlambda.f = sapply(fitlist,FUN=function(x) coef(x)["ghlambda"])
		sigma.f = sapply(fitlist,FUN=function(x) coef(x)["sigma"])
		fmatrix = matrix(NA, ncol = 5, nrow = forecast.length)
		colnames(fmatrix) = c("muf", "sigmaf", "ghlambdaf", "skewf", "shapef")
		fmatrix[,1] = fser
		for(i in 1:dim(eindex)[2]){
			fmatrix[eindex[,i], 2] = rep(sigma.f[i], refit.every)
			fmatrix[eindex[,i], 3] = rep(if(modelinc[18]>0) ghlambda.f[i] else 0, refit.every)
			fmatrix[eindex[,i], 4] = rep(if(modelinc[16]>0) skew.f[i] else 0, refit.every)
			fmatrix[eindex[,i], 5] = rep(if(modelinc[17]>0) shape.f[i] else 0, refit.every)
		}
		#re - scale the fmatrix to returns-based density
		smatrix = .scaledist(dist = distribution, fmatrix[,1], fmatrix[,2], fmatrix[,3], fmatrix[,4], fmatrix[,5])
		#fdates = dates[as.numeric(fwindex)]
		fdensity = vector(mode = "list", length = 1)
		f01density = vector(mode = "list", length = 1)
		fdensity[[1]] = data.frame(fdate = fdates, fmu=smatrix[,1], fsigma=smatrix[,2], fdlambda = fmatrix[,3], fskew=smatrix[,3], fshape=smatrix[,4])
		f01density[[1]] = data.frame(f01date=fdates, f01mu=fmatrix[,1], f01sigma=fmatrix[,2], f01ghlambda = fmatrix[,3], f01skew=fmatrix[,3], f01shape=fmatrix[,4])
		if(calculate.VaR){
			VaR.list = vector(mode="list", length = 1)
			cat("\n...calculating VaR...\n")
			if(is.null(VaR.alpha)) VaR.alpha = c(0.01, 0.05)
			n.v = dim(fdensity[[1]])[1]
			m.v = length(VaR.alpha)
			VaR.matrix = matrix(NA, ncol = m.v + 1, nrow = n.v)
			for(i in 1:m.v){
				VaR.matrix[,i] = .qdensity(rep(VaR.alpha[i], n.v) , mu = smatrix[,1], sigma = smatrix[,2], lambda = fmatrix[,3], 
						skew = smatrix[,3], shape = smatrix[,4], distribution = distribution)
			}
			VaR.matrix[,m.v+1] = unlist(outdata)
			colnames(VaR.matrix) = c(paste("alpha(", round(VaR.alpha,2)*100, "%)",sep=""), "actual")
			rownames(VaR.matrix) = as.character(fdates)
			VaR.list[[1]] = VaR.matrix
			cat("\nDone!\n")
		} else{
			VaR.matrix = NULL
			VaR.alpha = NULL
		}
		rolllist = list()
		rolllist$n.refit = nf
		rolllist$refit.every = refit.every
		rolllist$fdensity = fdensity
		rolllist$f01density = f01density
		rolllist$coefs = coefv
		rolllist$coefmat = coefmat
		rolllist$LLH = LLH
		rolllist$VaR.out = VaR.list
		rolllist$n.ahead = n.ahead
		rolllist$forecast.length = forecast.length
		rolllist$VaR.alpha =VaR.alpha
		model$modeldata$filterseries = filter.series
		model$modeldata$data = data
		model$modeldata$dates = dates
		model$modeldata$dates.format = xdata$dformat
		model$modeldata$datanames = datanames
		ans = new("ARFIMAroll",
				roll = rolllist,
				forecast = forecastlist,
				model = model)
	} else{
		eindex = t(.embed(1:forecast.length, refit.every, refit.every, TRUE))
		sigma.f = sapply(fitlist,FUN=function(x) coef(x)["sigma"])
		if(modelinc[16]>0) skew.f = sapply(fitlist,FUN=function(x) coef(x)["skew"])
		if(modelinc[17]>0) shape.f = sapply(fitlist,FUN=function(x) coef(x)["shape"])
		if(modelinc[18]>0) ghlambda.f = sapply(fitlist,FUN=function(x) coef(x)["ghlambda"])
		# split the array into 1:n-day ahead forecasts
		flist = vector(mode="list", length = n.ahead)
		f01density = vector(mode="list", length = n.ahead)
		fdensity = vector(mode="list", length = n.ahead)
		for(i in 1:n.ahead){
			fser =  unlist(lapply(forecastlist, FUN=function(x) sapply(x@forecast$forecast, FUN=function(x) x[i,])))
			rdat = data[(fitindx.start[1]+i):(fitindx.start[1] + forecast.length + i - 1)]
			fdates = as.character(dates[(fitindx.start[1]+i):(fitindx.start[1] + forecast.length + i - 1)])
			flist[[i]] = fdensity[[i]] = f01density[[i]] = matrix(NA, ncol = 5, nrow = forecast.length)	
			f01density[[i]][,1] = fser
			for(j in 1:dim(eindex)[2]){
				f01density[[i]][eindex[,j],2] = rep(sigma.f[j], refit.every)
				f01density[[i]][eindex[,j],3] = rep(if(modelinc[18]>0) ghlambda.f[j] else 0, refit.every)
				f01density[[i]][eindex[,j],4] = rep(if(modelinc[16]>0) skew.f[j] else 0, refit.every)
				f01density[[i]][eindex[,j],5] = rep(if(modelinc[17]>0) shape.f[j] else 0, refit.every)
			}
			fdensity[[i]] = .scaledist(dist = distribution, f01density[[i]][,1], f01density[[i]][,2], f01density[[i]][,3], f01density[[i]][,4], f01density[[i]][,5])
			fdensity[[i]] = as.data.frame(fdensity[[i]])
			f01density[[i]] = as.data.frame(f01density[[i]])
			
			fdensity[[i]] = cbind(fdates, fdensity[[i]])
			
			f01density[[i]] = cbind(fdates, f01density[[i]])
			
			colnames(fdensity[[i]])  = c("fdate", "fmu", "fsigma", "fskew", "fshape")
			colnames(f01density[[i]]) = c("f01date", "f01mu", "f01sigma", "f01ghlambda", "f01skew", "f01shape")
		}
		if(calculate.VaR){
			# should consider implementing mclapply here as well
			cat("\n...calculating VaR...\n")
			if(is.null(VaR.alpha)) VaR.alpha = c(0.01, 0.05)
			n.v = forecast.length
			m.v = length(VaR.alpha)
			VaR.list = vector(mode="list", length = n.ahead)
			for(j in 1:n.ahead){
				VaR.list[[j]] = matrix(NA, ncol = m.v+1, nrow = n.v)
				for(i in 1:m.v){
					VaR.list[[j]][,i] = .qdensity(rep(VaR.alpha[i], n.v) , mu = fdensity[[j]][,2],
							sigma = fdensity[[j]][,3], lambda = f01density[[i]][,3], skew = fdensity[[j]][,4], 
							shape = fdensity[[j]][,5],distribution = distribution)
				}
				VaR.list[[j]][,m.v+1] = data[(fitindx.start[1]+j):(fitindx.start[1] + forecast.length + j - 1)]
				colnames(VaR.list[[j]]) = c(paste("alpha(", round(VaR.alpha,2)*100, "%)",sep=""), "actual")
				rownames(VaR.list[[j]]) = as.character(dates[(fitindx.start[1]+j):(fitindx.start[1] + forecast.length + j - 1)])
			}
			cat("\nDone!\n")
		} else{
			VaR.list = NULL
			VaR.alpha = NULL
		}
			rolllist = list()
			rolllist$n.refit = nf
			rolllist$refit.every = refit.every
			rolllist$garchmodel = model
			rolllist$fdensity = fdensity
			rolllist$f01density = f01density
			rolllist$coefs = coefv
			rolllist$coefmat = coefmat
			rolllist$LLH = LLH
			rolllist$VaR.out = VaR.list
			rolllist$n.ahead = n.ahead
			rolllist$forecast.length = forecast.length
			rolllist$VaR.alpha =VaR.alpha
			model$modeldata$filterseries = filter.series
			model$modeldata$data = data
			model$modeldata$dates = dates
			model$modeldata$dates.format = xdata$dformat
			model$modeldata$datanames = datanames
			
			ans = new("ARFIMAroll",
					roll = rolllist,
					model = model,
					forecast = forecastlist)
	}
	return(ans)
}

.arfimadistribution = function(fitORspec, n.sim = 2000, n.start = 1, 
		m.sim = 100, recursive = FALSE, recursive.length = 6000, recursive.window = 1000, 
		prereturns = NA, preresiduals = NA, rseed = NA, custom.dist = list(name = NA, distfit = NA, type = "z"), 
		mexsimdata = NULL, fit.control = list(), solver = "solnp", 
		solver.control = list(), parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2), ...)
{
	if(recursive){
		nwindows = 1 + round( (recursive.length - n.sim) / recursive.window )
		swindow = vector(mode = "list", length = nwindows)
		rwindow = vector(mode = "list", length = nwindows)
	} else{
		nwindows = 1
		swindow = vector(mode = "list", length = nwindows)
		rwindow = vector(mode = "list", length = nwindows)
		recursive.window = 0
	}
	if(is(fitORspec, "ARFIMAfit")){
		for(i in 1:nwindows){
			sim = arfimasim(fitORspec, n.sim = n.sim + (i-1)*recursive.window, 
					n.start = n.start, m.sim = m.sim, prereturns = prereturns, 
					preresiduals = preresiduals, rseed = rseed, custom.dist = custom.dist, 
					mexsimdata = mexsimdata)
			swindow[[i]]$path.df = as.data.frame(sim)
			swindow[[i]]$seed = sim@seed
		}
		xmodel = fitORspec@model
		fixpars = as.list(coef(fitORspec))
		truecoef = fitORspec@fit$robust.matcoef
		spec = getspec(fitORspec)
	}
	# simulate series paths
	if(is(fitORspec, "ARFIMAspec")){
		
		for(i in 1:nwindows){
			sim  = arfimapath(fitORspec, n.sim =  n.sim + (i-1)*recursive.window, 
					n.start = n.start, m.sim = m.sim, prereturns = prereturns, 
					preresiduals = preresiduals, rseed = rseed, custom.dist = custom.dist, 
					mexsimdata = mexsimdata)
			swindow[[i]]$path.df = as.data.frame(sim)
			swindow[[i]]$seed = sim@seed
		}
		spec = fitORspec
		xmodel = spec@model
		setfixed(spec) <- list(NA)
		spec@model$pars[,4] = spec@model$pars[,2]
		spec@model$pars[,2] = 0
		fixpars = fitORspec@model$fixed.pars
		truecoef = as.matrix(cbind(unlist(fitORspec@model$fixed.pars),rep(0,length(fixpars)),
						rep(10, length(fixpars)),rep(0,length(fixpars))))
	}
	fitlist = vector( mode = "list", length = m.sim )
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
			for(i in 1:nwindows){
				rwindow[[i]]$fitlist = multicore::mclapply(swindow[[i]]$path.df, FUN = function(x) .fitandextractarfima(spec, x, out.sample = 0, 
									solver = solver, fit.control = fit.control, solver.control = solver.control), mc.cores = parallel.control$cores)
			}
		} else{
			for(i in 1:nwindows){
				nx = dim(swindow[[i]]$path.df)[2]
				sfInit(parallel = TRUE, cpus = parallel.control$cores)
				sfExport("spec", "swindow", "solver", "fit.control", "solver.control", local = TRUE)
				rwindow[[i]]$fitlist = sfLapply(as.list(1:nx), fun = function(j) rugarch:::.fitandextractarfima(spec, swindow[[i]]$path.df[,j], 
									out.sample = 0, solver = solver, fit.control = fit.control, solver.control = solver.control))
				sfStop()
			}
		}
	} else{
		for(i in 1:nwindows){
			rwindow[[i]]$fitlist = lapply(swindow[[i]]$path.df, FUN = function(x) .fitandextractarfima(spec, x, out.sample = 0, 
								solver = solver, fit.control = fit.control, solver.control = solver.control))
		}
	}
	
	
	reslist = vector(mode = "list", length = nwindows)
	for(j in 1:nwindows){
		reslist[[j]]$simcoef = 	matrix(NA, ncol = length(fixpars), nrow = m.sim)
		reslist[[j]]$rmse = 	rep(NA, length = length(fixpars))
		reslist[[j]]$simcoefse = matrix(NA, ncol = length(fixpars), nrow = m.sim)
		reslist[[j]]$likelist = rep(NA, length = m.sim)
		reslist[[j]]$mlongrun = rep(NA, length = m.sim)
		reslist[[j]]$simmaxdata  = matrix(NA, ncol = 2, nrow = m.sim)
		reslist[[j]]$simmindata  = matrix(NA, ncol = 2, nrow = m.sim)
		reslist[[j]]$simmeandata  = matrix(NA, ncol = 2, nrow = m.sim)
		reslist[[j]]$simmomdata = matrix(NA, ncol = 2, nrow = m.sim)
		reslist[[j]]$convergence = 	rep(1, length = m.sim)
		reslist[[j]]$seeds = rep(1, length = m.sim)
		for(i in 1:m.sim){
			if(rwindow[[j]]$fitlist[[i]]$convergence!=0) next()
			reslist[[j]]$simcoef[i, ] = rwindow[[j]]$fitlist[[i]]$simcoef
			reslist[[j]]$simcoefse[i,] = rwindow[[j]]$fitlist[[i]]$simcoefse
			reslist[[j]]$likelist[i] = 	rwindow[[j]]$fitlist[[i]]$llh
			reslist[[j]]$mlongrun[i] = 	rwindow[[j]]$fitlist[[i]]$mlongrun
			reslist[[j]]$simmaxdata[i, ] = rwindow[[j]]$fitlist[[i]]$maxdata
			reslist[[j]]$simmindata[i, ] = rwindow[[j]]$fitlist[[i]]$mindata
			reslist[[j]]$simmeandata[i, ] = rwindow[[j]]$fitlist[[i]]$meandata
			reslist[[j]]$simmomdata[i, ] = rwindow[[j]]$fitlist[[i]]$momdata
			reslist[[j]]$convergence[i] = rwindow[[j]]$fitlist[[i]]$convergence
		}
		reslist[[j]]$seed = swindow[[j]]$seed
		reslist[[j]]$rmse = .rmse(reslist[[j]]$simcoef, unlist(fixpars))
	}
	reslist$details = list(n.sim = n.sim, n.start = n.start, m.sim = m.sim,  
			recursive = recursive, recursive.length = recursive.length, recursive.window = recursive.window,
			nwindows = nwindows)
	ans = new("ARFIMAdistribution",
			dist = reslist,
			truecoef = truecoef,
			model = xmodel)
	
	return(ans)
}

.fitandextractarfima = function(spec, x, out.sample = 0,  solver = "solnp", fit.control = list(), solver.control = list())
{
	dist = list()
	fit = .safefitarfima(spec, x, out.sample = 0, solver = solver, fit.control = fit.control, solver.control = solver.control)
	if( is.null(fit) || fit@fit$convergence == 1 || !is( fit, "ARFIMAfit" ) || any( is.na( coef( fit ) ) ) ){
		dist$convergence = 1
		return(dist)
	}
	dist$simcoef = coef(fit)
	dist$simcoefse = fit@fit$robust.matcoef[, 2]
	dist$llh = likelihood(fit)
	dist$mlongrun = uncmean(fit)
	tmp = as.data.frame(fit)
	dist$maxdata = apply(tmp[, -2], 2, "max")
	dist$mindata = apply(tmp[, -2], 2, "min")
	dist$meandata = apply(tmp[, -2], 2, "mean")
	dist$momdata = c(.kurtosis(tmp[,1]), .skewness(tmp[,1]))
	dist$convergence = fit@fit$convergence
	return(dist)
}

# autoarfima tests for best penalized insample fit based on information criteria
autoarfima = function(data, ar.max = 2, ma.max = 2, criterion = c("AIC", "BIC", "SIC", "HQIC"),
		method = c("partial", "full"), arfima = FALSE, include.mean = NULL, 
		distribution.model = "norm",
		parallel = FALSE, parallel.control = list(pkg = "snowfall", cores = 2), 
		external.regressors = NULL, solver = "solnp", solver.control=list(), 
		fit.control=list(), return.all = FALSE){
	m = tolower(method)
	ans = switch(m, 
			partial = .autoarfima1(Data = data, ar.max = ar.max, ma.max = ma.max, 
					criterion = criterion[1], arfima = arfima, include.mean = include.mean, 
					distribution.model = distribution.model,
					parallel = parallel, parallel.control = parallel.control, 
					external.regressors = external.regressors, solver = solver, 
					solver.control=solver.control, fit.control=fit.control, 
					return.all = return.all),
			full = .autoarfima2(Data = data, ar.max = ar.max, ma.max = ma.max, 
					criterion = criterion[1], arfima = arfima, include.mean = include.mean, 
					distribution.model = distribution.model,
					parallel = parallel, parallel.control = parallel.control, 
					external.regressors = external.regressors, solver = solver, 
					solver.control=solver.control, fit.control=fit.control, 
					return.all = return.all))
	return(ans)
}


.autoarfima1 = function(Data, ar.max = 2, ma.max = 2, criterion = c("AIC", "BIC", "SIC", "HQIC"),
		arfima = FALSE, include.mean = NULL, distribution.model = "norm", 
		parallel = FALSE, parallel.control = list(pkg = "snowfall", cores = 2), 
		external.regressors = NULL, solver = "solnp", solver.control=list(), fit.control=list(), return.all = FALSE){
	# combinations
	ar = 0:ar.max
	ma = 0:ma.max
	if(is.null(include.mean)) im = c(0,1) else im = as.integer(include.mean)
	if(arfima) arf = c(0, 1) else arf = 0
	d = expand.grid(ar=ar, ma=ma, im = im, arf = arf)
	# eliminate the zero row
	check = apply(d, 1, "sum")
	if(any(check == 0)){
		idx=which(check==0)
		d = d[-idx,]
	}
	n = dim(d)[1]
	IC = match(criterion[1], c("AIC", "BIC", "SIC", "HQIC"))
	fitlist = vector(mode = "list", length = n)
	if( parallel ){
		if( parallel.control$pkg == "multicore" ){
			if(!exists("mclapply")){
				require('multicore')
			}
			fitlist = multicore::mclapply(1:n, FUN = function(i){
						spec = arfimaspec(mean.model = list(armaOrder = c(d[i,1], d[i,2]),
										include.mean =  as.logical(d[i,3]), arfima = as.logical(d[i,4])),
								external.regressors = external.regressors, 
								distribution.model = distribution.model)
						fit = try(arfimafit(spec = spec, data = Data, 
								solver = solver, solver.control = solver.control, 
								fit.control = fit.control), silent = TRUE)
						if(!is(fit, "try-error") && fit$convergence!=0 && solver == "solnp"){
							fit = try(arfimafit(spec = spec, data = Data, 
											solver = "nlminb", solver.control = list(trace=0), 
											fit.control = fit.control), silent = TRUE)
						}
						return(fit)}, mc.cores = parallel.control$cores)
		} else{
			if(!exists("sfLapply")){
				require('snowfall')
			}
			sfInit(parallel = TRUE, cpus = parallel.control$cores)
			sfExport("d", "Data", "n", "solver", "solver.control", "fit.control",local = TRUE)
			fitlist = sfLapply(as.list(1:n), fun = function(i){
						spec = rugarch::arfimaspec(mean.model = list(armaOrder = c(d[i,1], d[i,2]),
										include.mean =  as.logical(d[i,3]), arfima = as.logical(d[i,4])),
								external.regressors = external.regressors, 
								distribution.model = distribution.model)
						fit = try(rugarch::arfimafit(spec = spec, data = Data, 
								solver = solver, solver.control = solver.control, 
								fit.control = fit.control), silent = TRUE)
						if(!is(fit, "try-error") && fit$convergence!=0 && solver == "solnp"){
							fit = try(arfimafit(spec = spec, data = Data, 
											solver = "nlminb", solver.control = list(trace=0), 
											fit.control = fit.control), silent = TRUE)
						}
						return(fit)})
			sfStop()
		}
	} else{
		for(i in 1:n){
			spec = arfimaspec(mean.model = list(armaOrder = c(d[i,1], d[i,2]),
							include.mean =  as.logical(d[i,3]), arfima = as.logical(d[i,4])),
					external.regressors = external.regressors, 
					distribution.model = distribution.model)
			fitlist[[i]] = try(arfimafit(spec = spec, data = Data, 
					solver = solver, solver.control = solver.control, 
					fit.control = fit.control), silent = TRUE)
			if(!is(fitlist[[i]], "try-error") && fitlist[[i]]$convergence!=0 && solver == "solnp"){
				fitlist[[i]] = try(arfimafit(spec = spec, data = Data, 
								solver = "nlminb", solver.control = list(trace=0), 
								fit.control = fit.control), silent = TRUE)
			}
		}
	}
	rankmat = matrix(NA, ncol = 6, nrow = n)
	colnames(rankmat) = c("AR", "MA", "Mean", "ARFIMA", criterion[1], "converged")
	rankmat[,1:4] = as.matrix(d[1:4])
	for(i in 1:n){
		if(inherits(fitlist[[i]], 'try-error') || fitlist[[i]]@fit$convergence!=0){
			rankmat[i,6] = FALSE
		} else{
			rankmat[i,5] = infocriteria(fitlist[[i]])[IC]
			rankmat[i,6] = TRUE
		}
	}
	rk = rankmat[order(rankmat[,5]),]
	rownames(rk) = 1:dim(rk)[1]
	i = which(rankmat[,5] == min(rankmat[,5], na.rm = TRUE))
	if(return.all){
		ans = list(fit = fitlist, rank.matrix = rankmat)
	} else{
		ans = list(fit = fitlist[[i]], rank.matrix = rk)	
	}
	return(ans)
}




.autoarfima2 = function(Data, ar.max = 2, ma.max = 2, criterion = c("AIC", "BIC", "SIC", "HQIC"),
		arfima = FALSE, include.mean = NULL, distribution.model = "norm", 
		parallel = FALSE, parallel.control = list(pkg = "snowfall", cores = 2), 
		external.regressors = NULL, solver = "solnp", solver.control=list(), fit.control=list(),
		return.all = FALSE)
{
	# combinations
	arnames = paste("ar", 1:ar.max, sep = "")
	manames = paste("ma", 1:ma.max, sep = "")
	.str = NULL
	if(ar.max>0){
		for(i in 1:ar.max){
			.str = c(.str, paste(arnames[i],"=c(0,1),",sep=""))
		}
	}
	if(ma.max>0){
		for(i in 1:ma.max){
			.str = c(.str, paste(manames[i],"=c(0,1),",sep=""))
		}
	}
	if(is.null(include.mean)){
		.str = c(.str, "im = c(0,1)")
	} else{
		.str = c(.str, "im = as.integer(include.mean)")
	}
	if(is.null(arfima)){
		.str = c(.str, ",arf = c(0,1)")
	} else{
		.str = c(.str, ",arf = as.integer(arfima)")
	}
	str = c("d = expand.grid(", paste(.str), ')')
	xstr = paste(str, sep="", collapse="")
	eval(parse(text=xstr))
	# eliminate the zero row
	check = apply(d, 1, "sum")
	if(any(check == 0)){
		idx=which(check==0)
		d = d[-idx,]
	}
	sumar = apply(d[,1:ar.max,], 1, "sum")
	summa = apply(d[,(ar.max+1):(ma.max+ar.max),], 1, "sum")
	n = dim(d)[1]
	IC = match(criterion[1], c("AIC", "BIC", "SIC", "HQIC"))
	fitlist = vector(mode = "list", length = n)
	if( parallel ){
		if( parallel.control$pkg == "multicore" ){
			if(!exists("mclapply")){
				require('multicore')
			}
			fitlist = multicore::mclapply(1:n, FUN = function(i){
						if(ar.max>0){
							arr = d[i,1:ar.max]
							if(ma.max>0){
								mar = d[i,(ar.max+1):(ma.max+ar.max)]
							} else{
								mar = 0
							}
						} else{
							arr = 0
							if(ma.max>0){
								mar = d[i,1:ma.max]
							} else{
								mar = 0
							}
						}
						spec = .zarfimaspec( arOrder = arr, maOrder = mar, 
								include.mean = d[i,'im'], arfima = d[i,'arf'], 
								external.regressors = external.regressors, 
								distribution.model = distribution.model)
						fit = try(arfimafit(spec = spec, data = Data, 
										solver = solver, solver.control = solver.control, 
										fit.control = fit.control), silent = TRUE)
						return(fit)}, mc.cores = parallel.control$cores)
		} else{
			if(!exists("sfLapply")){
				require('snowfall')
			}
			sfInit(parallel = TRUE, cpus = parallel.control$cores)
			sfExport("d", "Data", "n", "solver", "external.regressors",  "distribution.model",
					"ar.max", "ma.max", "solver.control", "fit.control",local = TRUE)
			fitlist = sfLapply(as.list(1:n), fun = function(i){
						if(ar.max>0){
							arr = d[i,1:ar.max]
							if(ma.max>0){
								mar = d[i,(ar.max+1):(ma.max+ar.max)]
							} else{
								mar = 0
							}
						} else{
							arr = 0
							if(ma.max>0){
								mar = d[i,1:ma.max]
							} else{
								mar = 0
							}
						}
						spec = rugarch:::.zarfimaspec( arOrder = arr, maOrder = mar, 
								include.mean = d[i,'im'], arfima = d[i,'arf'], 
								external.regressors = external.regressors, 
								distribution.model = distribution.model)
						fit = try(rugarch::arfimafit(spec = spec, data = Data, 
										solver = solver, solver.control = solver.control, 
										fit.control = fit.control), silent = TRUE)
						return(fit)})
			sfStop()
		}
	} else{
		for(i in 1:n){
			if(ar.max>0){
				arr = d[i,1:ar.max]
				if(ma.max>0){
					mar = d[i,(ar.max+1):(ma.max+ar.max)]
				} else{
					mar = 0
				}
			} else{
				arr = 0
				if(ma.max>0){
					mar = d[i,1:ma.max]
				} else{
					mar = 0
				}
			}
			spec = .zarfimaspec( arOrder = arr, maOrder = mar, 
					include.mean = d[i,'im'], arfima = d[i,'arf'], 
					external.regressors = external.regressors, 
					distribution.model = distribution.model)
			fitlist[[i]] = try(arfimafit(spec = spec, data = Data, 
							solver = solver, solver.control = solver.control, 
							fit.control = fit.control), silent = TRUE)
		}
	}
	m = dim(d)[2]
	rankmat = matrix(NA, ncol = m+2, nrow = n)
	colnames(rankmat) = c(colnames(d), criterion[1], "converged")
	rankmat[,1:m] = as.matrix(d)
	for(i in 1:n){
		if(inherits(fitlist[[i]], 'try-error') || fitlist[[i]]@fit$convergence!=0){
			rankmat[i,m+2] = 0
		} else{
			rankmat[i,m+1] = infocriteria(fitlist[[i]])[IC]
			rankmat[i,m+2] = 1
		}
	}
	rk = rankmat[order(rankmat[,m+1]),]
	rownames(rk) = 1:dim(rk)[1]
	i = which(rankmat[,m+1] == min(rankmat[,m+1], na.rm = TRUE))
	if(return.all){
		ans = list(fit = fitlist, rank.matrix = rankmat)
	} else{
	 	ans = list(fit = fitlist[[i]], rank.matrix = rk)	
	}
	return(ans)
}