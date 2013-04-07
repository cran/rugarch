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

.multispecall = function( speclist ){
	model = unlist( strsplit(class(speclist[[1]]), "spec") )
	if( model == "ARFIMA" ){
		ans = .multispecarfima( speclist )
	} else{
		ans = .multispecgarch( speclist )
	}
	return( ans )
}

.multispecgarch = function( speclist )
{
	# first create a spec which goes through validation process
	tp = 1
	if( !all(unlist(lapply(speclist, FUN = function(x) is(x, "uGARCHspec"))) ) ){
		stop("\nNot a valid list of univariate GARCH specs.")
	}
	# then check type
	n = length(speclist)
	for(i in 2:n){
		modelnames1 = rownames( speclist[[i]]@model$pars[speclist[[i]]@model$pars[,3]==1, ] )
		modelnames2 = rownames( speclist[[i-1]]@model$pars[speclist[[i-1]]@model$pars[,3]==1, ] )
		if(length(modelnames1) != length(modelnames2))
		{
			tp  = 0
			break()
		} else{
			if(!all(modelnames1 == modelnames2))
			{
				tp  = 0
				break()
			}
		}
	}
	if(tp) type = "equal" else type = "unequal"
	ans = new("uGARCHmultispec",
			spec = speclist,
			type = type)
	return(ans)
}

# a multifit function possible utilizing parallel execution returning a fitlist
# object
.multifitgarch = function(multispec, data, out.sample = 0, solver = "solnp", 
		solver.control = list(), fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, 
				rec.init = "all"), cluster = NULL, ...)
{
	n = length(multispec@spec)
	if(is.null(data)) stop("\nmultifit GARCH-->error: multifit requires a data object", call. = FALSE)
	if(!is.matrix(data) & !is.data.frame(data)) stop("\nmultifit GARCH-->error: multifit only supports matrix or data.frame objects for the data", call. = FALSE)
	if(is.matrix(data)) data = as.data.frame(data)
	asset.names = colnames(data)
	if(dim(data)[2] != n)
		stop("\nmultifit GARCH-->error: speclist length not equal to data length", call. = FALSE)
	fitlist = vector(mode = "list", length = n)
	if(length(out.sample) == 1 | length(out.sample) < n) out.sample = rep(out.sample, n)
	
	##################
	# Parallel Execution Prelim Checks
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("multispec", "data", "out.sample", "solver", 
							"solver.control", "fit.control"), envir = environment())
			fitlist = parLapply(cluster, as.list(1:n), fun = function(i){
						ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
								out.sample = out.sample[i], solver = solver, 
								solver.control = solver.control, fit.control = fit.control)
					})
		} else{
			fitlist = lapply(as.list(1:n), FUN = function(i){
						ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
								out.sample = out.sample[i], solver = solver, 
								solver.control = solver.control, fit.control = fit.control)
					})
	}
	# converged: print
	desc = list()
	desc$type = multispec@type
	desc$asset.names = asset.names
	ans = new("uGARCHmultifit",
			fit  = fitlist,
			desc = desc)
	return(ans)
}

.multifiltergarch1 = function(multifitORspec, data = NULL, out.sample = 0, 
		n.old = NULL, rec.init = 'all', cluster = NULL, ...)
{
	fitlist = multifitORspec
	n = length(fitlist@fit)
	if(is.null(data)) 
		stop("\nmultifilter GARCH-->error: multifilter requires a data object", call. = FALSE)
	if(!is.matrix(data) & !is.data.frame(data)) 
		stop("\nmultifilter GARCH-->error: multifilter only supports matrix or data.frame objects for the data", call. = FALSE)
	if(dim(data)[2] != n)
		stop("\nmultifilter GARCH-->error: fitlist length not equal to data length", call. = FALSE)
	if(is.matrix(data)) data = as.data.frame(data)
	asset.names = colnames(data)
	filterlist = vector(mode = "list", length = n)
	if(length(out.sample) == 1 | length(out.sample) < n) out.sample = rep(out.sample, n)
	if(length(rec.init) == 1 | length(rec.init) < n) rec.init = rep(rec.init, n)
	
	specx = vector(mode = "list", length = n)
	for(i in 1:n){
		specx[[i]] = getspec(fitlist@fit[[i]])
		specx[[i]]@model$fixed.pars = as.list(coef(fitlist@fit[[i]]))
	}
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("specx", "data", "out.sample", "n.old", "rec.init"), envir = environment())
		filterlist = parLapply(cluster, as.list(1:n), fun = function(i){
					ugarchfilter(data = data[, i, drop = FALSE], spec = specx[[i]], 
							out.sample =  out.sample[i], n.old = n.old, rec.init = rec.init[i])
				})
	} else{
		filterlist = lapply(as.list(1:n), FUN = function(i){
					ugarchfilter(data = data[, i, drop = FALSE], spec = specx[[i]], 
							out.sample =  out.sample[i], n.old = n.old, rec.init = rec.init[i])
				})
	}
	desc = list()
	desc$type = "equal"
	desc$asset.names = asset.names
	
	ans = new("uGARCHmultifilter",
			filter = filterlist,
			desc = desc)
	return(ans)
}

.multifiltergarch2 = function(multifitORspec, data = NULL, out.sample = 0, 
		n.old = NULL, rec.init = 'all', cluster = NULL, ...)
{
	speclist = multifitORspec
	n = length(speclist@spec)
	if(is.null(data)) 
		stop("\nmultifilter GARCH-->error: multifilter requires a data object", call. = FALSE)
	if(!is.matrix(data) & !is.data.frame(data)) 
		stop("\nmultifilter GARCH-->error: multifilter only supports matrix or data.frame objects for the data", call. = FALSE)
	if(dim(data)[2] != n)
		stop("\nmultifilter GARCH-->error: multispec length not equal to data length", call. = FALSE)
	if(is.matrix(data)) data = as.data.frame(data)
	asset.names = colnames(data)
	filterlist = vector(mode = "list", length = n)
	if(length(out.sample) == 1 | length(out.sample) < n) out.sample = rep(out.sample, n)
	if(length(rec.init) == 1 | length(rec.init) < n) rec.init = rep(rec.init, n)
	
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("speclist", "data", "out.sample", "n.old", "rec.init"), envir = environment())
		filterlist = parLapply(cluster, as.list(1:n), fun = function(i){
					ugarchfilter(data = data[, i, drop = FALSE], spec = speclist@spec[[i]], 
							out.sample =  out.sample[i], n.old = n.old, rec.init = rec.init[i])
				})
	} else{
		filterlist = lapply(as.list(1:n), FUN = function(i){
					ugarchfilter(data = data[, i, drop = FALSE], spec = speclist@spec[[i]], 
							out.sample =  out.sample[i], n.old = n.old, rec.init = rec.init[i])
				})	
	}
	# converged: print
	desc = list()
	desc$type = speclist@type
	desc$asset.names = asset.names
	
	ans = new("uGARCHmultifilter",
			filter  = filterlist,
			desc = desc)
	return(ans)
}

.multiforecastgarch1 = function(multifitORspec, data = NULL, n.ahead = 1, 
		n.roll = 0, out.sample = 0, external.forecasts = list(mregfor = NULL, vregfor = NULL), 
		cluster = NULL, ...)
{
	multifit = multifitORspec
	n = length(multifit@fit)
	asset.names = multifit@desc$asset.names
	forecastlist = vector(mode = "list", length = n)
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("multifit", "n.ahead", "n.roll", "external.forecasts"), envir = environment())
		forecastlist = parLapply(cluster, as.list(1:n), fun = function(i){
					ugarchforecast(fitORspec = multifit@fit[[i]], data = NULL, n.ahead = n.ahead, 
							n.roll = n.roll, external.forecasts = external.forecasts)
				})
	} else{
		forecastlist = lapply(as.list(1:n), FUN = function(i){
					ugarchforecast(fitORspec = multifit@fit[[i]], data = NULL, n.ahead = n.ahead, 
							n.roll = n.roll, external.forecasts = external.forecasts)
				})
	}
	desc = list()
	desc$type = "equal"
	desc$asset.names = asset.names
	ans = new("uGARCHmultiforecast",
			forecast  = forecastlist,
			desc = desc)
	return(ans)
}

.multiforecastgarch2 = function(multifitORspec, data = NULL, n.ahead = 1, n.roll = 0, out.sample = 0,
		external.forecasts = list(mregfor = NULL, vregfor = NULL), cluster =NULL, ...)
{
	multispec = multifitORspec
	n = length(multispec@spec)
	if(is.null(data)) 
		stop("\nmultiforecast GARCH-->error: multiforecast with multiple spec requires a data object", call. = FALSE)
	if(!is.matrix(data) & !is.data.frame(data)) 
		stop("\nmultiforecast GARCH-->error: multiforecast only supports matrix or data.frame objects for the data", call. = FALSE)
	# data.frame allows to call data with x[i] notation
	if(is.matrix(data)) data = as.data.frame(data)
	if(dim(data)[2] != n)
		stop("\nmultiforecast GARCH-->error: multispec length not equal to data length", call. = FALSE)
	asset.names = colnames(data)
	forecastlist = vector(mode = "list", length = n)
	if(is.null(out.sample)) out.sample = 0
	if(length(out.sample) == 1) out.sample = rep(out.sample, n)
	if(length(out.sample) !=n ) 
		stop("\nmultiforecast GARCH-->error: out.sample length not equal to data length", call. = FALSE)
	
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("multispec", "data", "n.ahead", "n.roll", 
						"out.sample", "external.forecasts"), envir = environment())
		forecastlist = parLapply(cluster, as.list(1:n), fun = function(i){
					ugarchforecast(fitORspec = multispec@spec[[i]], 
							data = data[, i, drop = FALSE], n.ahead = n.ahead, n.roll = n.roll, 
							out.sample = out.sample[i], external.forecasts = external.forecasts)
				})
	} else{
		forecastlist = lapply(as.list(1:n), FUN = function(i){
					ugarchforecast(fitORspec = multispec@spec[[i]], 
							data = data[, i, drop = FALSE], n.ahead = n.ahead, n.roll = n.roll, 
							out.sample = out.sample[i], external.forecasts = external.forecasts, ...)
				})
	}
	desc = list()
	desc$type = multispec@type
	desc$asset.names = asset.names
	
	ans = new("uGARCHmultiforecast",
			forecast  = forecastlist,
			desc = desc)
	return(ans)
}