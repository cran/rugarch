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


.rollfdensity = function(spec,  data, n.ahead = 1, forecast.length = 500, 
		refit.every = 25, refit.window = c("recursive", "moving"), parallel = FALSE, 
		parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2), 
		solver = "solnp", fit.control = list(), solver.control = list() ,
		calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05))
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
	vexdata = vector(mode="list",length = nf)
	fitlist = vector(mode="list",length = nf)
	outdata = vector(mode="list",length = nf)
	coefv = vector(mode="list", length = nf)
	distribution = spec@model$modeldesc$distribution
	model = spec@model
	modelinc = model$modelinc
	forecastlist = vector(mode="list", length = nf)
	forecast.length = floor(forecast.length/refit.every) * refit.every
	VaR.list = NULL
	NL = vector(mode = "list", length = nf)
	for(i in 1:nf){
		if(refit.window[1]=="recursive"){
			fitdata[[i]] = xdata[1:(fitindx.start[i]+refit.every), , drop = FALSE]
			outdata[[i]] = xdata[(fitindx.start[i]+1):(fitindx.start[i]+refit.every), , drop = FALSE]
			if( modelinc[6]>0 ){
				mexdata[[i]] = model$modeldata$mexdata[1:(fitindx.start[i]+refit.every), , drop = FALSE]
			} else{
				mexdata[[i]] = NULL
			}
			if( modelinc[15]>0 ){
				vexdata[[i]] = model$modeldata$vexdata[1:(fitindx.start[i]+refit.every), , drop = FALSE]
			} else{
				vexdata[[i]] = NULL
			}
			NL[[i]] = 1:(fitindx.start[i]+refit.every)
		} else{
			fitdata[[i]] = xdata[(i*refit.every):(fitindx.start[i]+refit.every), , drop = FALSE]
			outdata[[i]] = xdata[(fitindx.start[i]+1):(fitindx.start[i]+refit.every), , drop = FALSE]
			if( modelinc[6]>0 ){
				mexdata[[i]] = model$modeldata$mexdata[(i*refit.every):(fitindx.start[i]+refit.every), , drop = FALSE]
			} else{
				mexdata[[i]] = NULL
			}
			if( modelinc[15]>0 ){
				vexdata[[i]] = model$modeldata$vexdata[(i*refit.every):(fitindx.start[i]+refit.every), , drop = FALSE]
			} else{
				vexdata[[i]] = NULL
			}
			NL[[i]] = (i*refit.every):(fitindx.start[i]+refit.every)
		}
	}
	
	cat("\n...estimating refit windows...\n")
	specx = vector(mode = "list", length = nf)
	for(i in 1:nf){
		specx[[i]] = spec
		if(modelinc[6]>0)  specx[[i]]@model$modeldata$mexdata = mexdata[[i]] else specx[[i]]@model$modeldata$mexdata = NULL
		if(modelinc[15]>0) specx[[i]]@model$modeldata$vexdata = vexdata[[i]] else specx[[i]]@model$modeldata$vexdata = NULL
	}
	if( parallel ){
		if( parallel.control$pkg == "multicore" ){
			if(!exists("mclapply")){
				require('multicore')
			}
			fitlist = multicore::mclapply(1:nf, FUN = function(i){
						xdat = data.frame(fitdata[[i]])
						rownames(xdat) = dates[NL[[i]]]
					ans = ugarchfit(spec = specx[[i]], data = xdat, solver = solver,
							out.sample = refit.every, solver.control = solver.control,
							fit.control = fit.control)
					if(ans@fit$convergence!=0){
						if(solver == "solnp") solverx = "nlminb" else solverx = "gosolnp"
						ans = ugarchfit(spec = specx[[i]], data = fitdata[[i]], solver = solverx,
							out.sample = refit.every, solver.control = list(),
							fit.control = fit.control)
					}
					return(ans)
				}, mc.cores = parallel.control$cores)
		} else{
			if(!exists("sfLapply")){
				require('snowfall')
			}
			sfInit(parallel = TRUE, cpus = parallel.control$cores)
			sfExport("NL", "specx", "fitdata", "dates", "refit.every", "solver", "solver.control", "fit.control",local = TRUE)
			
			fitlist = sfLapply(as.list(1:nf), fun = function(i){
						xdat = data.frame(fitdata[[i]])
						rownames(xdat) = dates[NL[[i]]]
						ans = rugarch::ugarchfit(spec = specx[[i]], data = xdat, solver = solver,
								out.sample = refit.every, solver.control = solver.control,
								fit.control = fit.control)
						if(ans@fit$convergence!=0){
							if(solver == "solnp") solverx = "nlminb" else solverx = "gosolnp"
							ans = rugarch::ugarchfit(spec = specx[[i]], data = fitdata[[i]], solver = solverx,
									out.sample = refit.every, solver.control = list(),
									fit.control = fit.control)
						}
						return(ans)
					})
			sfStop()
		}
	} else{
		for(i in 1:nf){
			xdat = data.frame(fitdata[[i]])
			rownames(xdat) = dates[NL[[i]]]
			fitlist[[i]] = ugarchfit(data = xdat, spec = specx[[i]], solver = solver, 
					out.sample = refit.every, solver.control = solver.control,
					fit.control = fit.control)
			if(fitlist[[i]]@fit$convergence!=0){
				if(i>1){
					specx[[i]]@model$start.pars = as.list(coef(fitlist[[i-1]]))
					fitlist[[i]] = ugarchfit(data = xdat, spec = specx[[i]], solver = solver, 
						out.sample = refit.every, solver.control = solver.control,
						fit.control = fit.control)
					specx[[i]]@model$start.pars = NULL
				}
			}
		}
	}
	# TODO: return fitlist and allow to resume
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
			forecastlist = multicore::mclapply(fitlist, FUN = function(x) ugarchforecast(x, 
							n.ahead = n.ahead, n.roll = refit.every-1), mc.cores = parallel.control$cores)
		} else{
			nx = length(fitlist)
			sfInit(parallel = TRUE, cpus = parallel.control$cores)
			sfExport("fitlist", "n.ahead", "refit.every", local = TRUE)
			forecastlist = sfLapply(as.list(1:nx), fun = function(i) rugarch::ugarchforecast(fitlist[[i]], 
								n.ahead = n.ahead, n.roll = refit.every-1))
			sfStop()
		}
	} else{
		forecastlist = lapply(fitlist, FUN = function(x) ugarchforecast(x, n.ahead = n.ahead, n.roll = refit.every-1))
	}

	eindex = t(.embed(1:forecast.length, refit.every, refit.every, TRUE))
	for(i in 1:nf){
		filter.sigma[eindex[,i]] = forecastlist[[i]]@model$modeldata$sigma[(fitindx.start[i]
							+1):(fitindx.start[i]+refit.every)]
	}
	# collect forecasts [mu sigma (skew shape)]
	# 2 cases : n.ahead = 1 we use a matrix, else
	# n.ahead > 1 we use a list object
	if(n.ahead == 1){
		fser = as.numeric(sapply(forecastlist,FUN=function(x) sapply(x@forecast$forecasts, FUN=function(x) x[,2])))
		fsig = as.numeric(sapply(forecastlist,FUN=function(x) sapply(x@forecast$forecasts, FUN=function(x) x[,1])))
		fdates = dates[as.numeric(fwindex)]
		rdat = data[(fitindx.start[1]+1):(fitindx.start[1] + forecast.length + 1 - 1)]
		eindex = t(.embed(1:forecast.length, refit.every, refit.every, TRUE))
		if(modelinc[16]>0) skew.f = sapply(fitlist,FUN=function(x) coef(x)["skew"])
		if(modelinc[17]>0) shape.f = sapply(fitlist,FUN=function(x) coef(x)["shape"])
		if(modelinc[18]>0) ghlambda.f = sapply(fitlist,FUN=function(x) coef(x)["ghlambda"])
		fmatrix = matrix(NA, ncol = 5, nrow = forecast.length)
		colnames(fmatrix) = c("muf", "sigmaf", "ghlambdaf", "skewf", "shapef")
		fmatrix[,1] = fser
		fmatrix[,2] = fsig
		for(i in 1:dim(eindex)[2]){
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
		model$modeldata$filtersigma = filter.sigma
		model$modeldata$data = data
		model$modeldata$dates = dates
		model$modeldata$dates.format = xdata$dformat
		model$modeldata$datanames = datanames
		ans = new("uGARCHroll",
				roll = rolllist,
				forecast = forecastlist,
				model = model)
	} else{
		eindex = t(.embed(1:forecast.length, refit.every, refit.every, TRUE))
		if(modelinc[16]>0) skew.f = sapply(fitlist,FUN=function(x) coef(x)["skew"])
		if(modelinc[17]>0) shape.f = sapply(fitlist,FUN=function(x) coef(x)["shape"])
		if(modelinc[18]>0) ghlambda.f = sapply(fitlist,FUN=function(x) coef(x)["ghlambda"])
		# split the array into 1:n-day ahead forecasts
		flist = vector(mode="list", length = n.ahead)
		f01density = vector(mode="list", length = n.ahead)
		fdensity = vector(mode="list", length = n.ahead)
		for(i in 1:n.ahead){
			fsig  = unlist(lapply(forecastlist, FUN=function(x) sapply(x@forecast$forecast, FUN=function(x) x[i,1])))
			fser =  unlist(lapply(forecastlist, FUN=function(x) sapply(x@forecast$forecast, FUN=function(x) x[i,2])))
			rdat = data[(fitindx.start[1]+i):(fitindx.start[1] + forecast.length + i - 1)]
			fdates = as.character(dates[(fitindx.start[1]+i):(fitindx.start[1] + forecast.length + i - 1)])
			flist[[i]] = fdensity[[i]] = f01density[[i]] = matrix(NA, ncol = 5, nrow = forecast.length)	
			f01density[[i]][,1] = fser
			f01density[[i]][,2] = fsig
			for(j in 1:dim(eindex)[2]){
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
		model$modeldata$filtersigma = filter.sigma
		model$modeldata$data = data
		model$modeldata$dates = dates
		model$modeldata$dates.format = xdata$dformat
		model$modeldata$datanames = datanames
		
		ans = new("uGARCHroll",
				roll = rolllist,
				forecast = forecastlist,
				model = model)
	}
	return(ans)
}

# forecast performance measures
.ugarchrollreport = function(object, type = "VaR", n.ahead = 1, VaR.alpha = 0.01, conf.level = 0.95)
{
	switch(type,
			VaR = .rollVaRreport(object, n.ahead, VaR.alpha, conf.level),
			fpm = .rollfpmreport(object))
	invisible(object)
}

.rollfpmreport = function(object)
{
	vmodel = object@model$modeldesc$vmodel
	vsubmodel = object@model$modeldesc$vsubmodel
	cat(paste("\nGARCH Roll Mean Forecast Performance Measures", sep = ""))
	cat(paste("\n---------------------------------------------", sep = ""))
	cat(paste("\nModel : ", vmodel, sep = ""))
	if(vmodel == "fGARCH"){
		cat(paste("\nSubModel : ", vsubmodel, sep = ""))
	}
	cat(paste("\nno.refits : ", object@roll$n.refit, sep = ""))
	cat(paste("\nn.ahead   : ", object@roll$n.ahead, sep = ""))
	cat(paste("\nn.rolls   : ", object@roll$n.refit*object@roll$refit.every, sep = ""))
	cat("\n\n")
	tmp = fpm(object)
	print(signif(tmp, 4))
	cat("\n\n")
}

.rollVaRreport = function(object, n.ahead = 1, VaR.alpha = 0.01, conf.level = 0.95)
{
	n = object@roll$n.ahead
	v.a = object@roll$VaR.alpha
	
	if(is.null(object@roll$VaR.out)) stop("\nplot-->error: VaR was not calculated for this object\n", call.=FALSE)
	if(n.ahead > n) stop("\nplot-->error: n.ahead chosen is not valid for object\n", call.=FALSE)
	if(!is.null(v.a) && !any(v.a==VaR.alpha[1])) stop("\nplot-->error: VaR.alpha chosen is invalid for the object\n", call.=FALSE)
	if(is.list(object@roll$VaR.out)){
		dvar = object@roll$VaR.out[[n.ahead]]
		m = dim(dvar)[2]
		idx = which(colnames(dvar) == paste("alpha(", round(VaR.alpha,2)*100, "%)",sep=""))
		.VaRreport(object@model$modeldata$datanames, object@model$modeldesc$vmodel, object@model$modeldesc$distribution, p = VaR.alpha, actual=dvar[,m], VaR = dvar[, idx], 
				conf.level = conf.level)
	} else{
		dvar = object@roll$VaR.out
		m = dim(dvar)[2]
		idx = which(colnames(dvar) == paste("alpha(", round(VaR.alpha,2)*100, "%)",sep=""))
		.VaRreport(object@model$modeldata$datanames, object@model$modeldesc$vmodel, object@model$modeldesc$distribution, p = VaR.alpha, actual=dvar[,m], VaR = dvar[, idx], 
				conf.level = conf.level)
	}
	invisible(object)
}

.embed = function(data, k, by = 1, ascending = FALSE) 
{
	# n = no of time points, k = number of columns
	# by = increment. normally =1 but if =b calc every b-th point 
	# ascending If TRUE, points passed in ascending order else descending.
	# Note that embed(1:n,k) corresponds to embedX(n,k,by=1,rev=TRUE)
	# e.g. embedX(10,3)
	if(is.null(dim(data)[1])) n<-length(data) else n<-dim(data)[1]
	s <- seq(1,n-k+1,by)
	lens <- length(s)
	cols <- if (ascending) 1:k else k:1
	return(matrix(data[s + rep(cols,rep(lens,k))-1],lens))
}
