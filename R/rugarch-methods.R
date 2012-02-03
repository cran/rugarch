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

#----------------------------------------------------------------------------------
# univariate spec method
#----------------------------------------------------------------------------------
ugarchspec = function(variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
				submodel = NULL, external.regressors = NULL, variance.targeting = FALSE), 
		mean.model = list(armaOrder = c(1,1), include.mean = TRUE, archm = FALSE, 
				archpow = 1, arfima = FALSE, external.regressors = NULL, archex = FALSE), 
		distribution.model = "norm", start.pars = list(), fixed.pars = list(), ...)
{
	UseMethod("ugarchspec")
}

# [mu ar ma arfima im mxreg omega alpha beta gamma gamma11 gamma21 delta lambda vxreg skew shape dlamda aux aux aux aux]

.expand.model = function(model){
	modelnames = NULL
	for(i in 1:21){
		if(model[i]>0){
			if(any(c(2,3,6,8,9,10,11,12,15) == i)){
				modelnames = c(modelnames, paste(names(model)[i], 1:model[i], sep = ""))
			} else{
				modelnames = c(modelnames, names(model)[i])
			}
		}
	}
	return( modelnames )
}

# Changelog:
# 06-12-2011 Added "archex" option in mean.model so that external regressor might
# be multiplied by conditional variance. If used must be integer and represents the
# number of series from the end of the supplied external regressors.

.ugarchspec = function(variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
				submodel = NULL, external.regressors = NULL, variance.targeting = FALSE), 
		mean.model = list(armaOrder = c(1,1), include.mean = TRUE, archm = FALSE, 
				archpow = 1, arfima = FALSE, external.regressors = NULL, archex = FALSE), 
		distribution.model = "norm", start.pars = list(), fixed.pars = list())
{
	# some checks and preparation to be passed on to specific models by switch
	modelinc = rep(0, 21)
	names(modelinc) = c("mu", "ar", "ma", "arfima", "archm", "mxreg", "omega", "alpha",
			"beta", "gamma", "eta1", "eta2", "delta", "lambda", "vxreg", "skew", "shape",
			"ghlambda", "aux", "aux", "aux")
	
	modeldesc = list()
	modeldata = list()
	# distribution model
	if(is.null(distribution.model)) modeldesc$distribution = "norm"
	valid.distribution = c("norm", "snorm", "std", "sstd","ged", "sged", "nig", "ghyp", "jsu", "ghst")
	distribution = distribution.model
	
	if(!is.character(distribution[1]))
		stop("\nugarchspec-->error: cond.distribution argument must be a character")
	
	if(!any(distribution==valid.distribution)) 
		stop("\nugarchspec-->error: the cond.distribution does not appear to be a valid choice.")
	
	if(length(distribution)!=1) distribution = distribution[1]
	
	modeldesc$distribution = distribution
	modeldesc$distno = which(distribution == valid.distribution)
	di = .DistributionBounds(distribution)
	modelinc[16] = di$include.skew
	modelinc[17] = di$include.shape
	modelinc[18] = di$include.ghlambda
	# the last aux value is the distribution number
	modelinc[21] = modeldesc$distno
	# variance model:
	vmodel = list()
		
	valid.model = c("sGARCH", "eGARCH", "gjrGARCH", "tGARCH", "fGARCH", "iGARCH", "apARCH")
	if(is.null(variance.model$model)){
		modeldesc$vmodel = "sGARCH"
	} else{
		modeldesc$vmodel = variance.model$model[1]
		
		if(!is.character(modeldesc$vmodel)) 
			stop("\nugarchspec-->error: garch model argument must be a character.\n", call. = FALSE)
		
		if(!any(modeldesc$vmodel == valid.model)) 
			stop("\nugarchspec-->error: the garch model does not appear to be a valid choice.\n", call. = FALSE)
		
		if(modeldesc$vmodel == "fGARCH"){
			modeldesc$vsubmodel = variance.model$submodel
			valid.submodel = c("GARCH","TGARCH","AVGARCH","NGARCH","NAGARCH","APARCH","ALLGARCH","GJRGARCH")
			if(is.null(modeldesc$vsubmodel)) 
				stop("\nugarchspec-->error: NULL not allowed for the submodel when model is of type fGARCH.\n", call. = FALSE)
			if(!any(modeldesc$vsubmodel == valid.submodel)) 
				stop("\nugarchspec-->error: the fGARCH submodel does not appear to be a valid choice. See documentation for valid choices.\n", call. = FALSE)
		}
	}
	
	# depending on model include the additional parameters
	if(is.null(variance.model$garchOrder)){
		modelinc[8] = 1
		modelinc[9] = 1
	} else{
		modelinc[8] = variance.model$garchOrder[1]
		modelinc[9] = variance.model$garchOrder[2]
	}
	
	if( modeldesc$vmodel == "gjrGARCH" ) modelinc[10] = modelinc[8]
	if( modeldesc$vmodel == "eGARCH" ) modelinc[10] = modelinc[8]
	if( modeldesc$vmodel == "apARCH" ){
		modelinc[10] = modelinc[8]
		modelinc[13] = 1
	}
	if( modeldesc$vmodel == "fGARCH" ){
		if(modeldesc$vsubmodel == "AVGARCH"){
			modelinc[12] = modelinc[11] = modelinc[8]
		}
		if(modeldesc$vsubmodel == "GJRGARCH") modelinc[11] = modelinc[8]
		if(modeldesc$vsubmodel == "TGARCH") modelinc[11] = modelinc[8]
		if(modeldesc$vsubmodel == "NGARCH") modelinc[14] = 1
		if(modeldesc$vsubmodel == "NAGARCH") modelinc[12] = modelinc[8]
		if(modeldesc$vsubmodel == "APARCH"){
			modelinc[14] = 1
			modelinc[11] = modelinc[8]
		}
		if(modeldesc$vsubmodel == "ALLGARCH"){
			modelinc[12] = modelinc[11] = modelinc[8]
			modelinc[14] = 1
		}
	}
	if( modeldesc$vmodel == "iGARCH" && modelinc[9] == 0 ) 	stop("\nugarchspec-->error: the iGARCH model requires the GARCH beta parameter.\n", call. = FALSE)
	
	modeldata$vexdata = variance.model$external.regressors
	if( !is.null(variance.model$external.regressors) ) modelinc[15] = dim( variance.model$external.regressors )[2]
	
	if(is.null(variance.model$variance.targeting)) modelinc[7] = 1 else modelinc[7] = as.integer( 1-variance.model$variance.targeting )
	
	# mean model:
	if(is.null(mean.model$armaOrder)){
		modelinc[2] = modelinc[3] = 1
	} else{
		modelinc[2] = mean.model$armaOrder[1]
		modelinc[3] = mean.model$armaOrder[2]
	}
	if(is.null(mean.model$include.mean)) modelinc[1] = 1 else modelinc[1] = as.integer( mean.model$include.mean )
	
	if(is.null(mean.model$archm) || !mean.model$archm){
		modelinc[5] = 0
	} else{
		if(is.null(mean.model$archpow)) mean.model$archpow = 1
		modelinc[5] = as.integer( mean.model$archpow )
	}
	

	if(is.null(mean.model$arfima)) modelinc[4] = 0 else modelinc[4] = as.integer( mean.model$arfima )
	
	modeldata$mexdata = mean.model$external.regressors
	if( !is.null(mean.model$external.regressors) ) modelinc[6] = dim( mean.model$external.regressors )[2]
	
	if(is.null(mean.model$archex) || !mean.model$archex){
		modelinc[20] = 0
	} else{
		modelinc[20] = as.integer( mean.model$archex )
		if(modelinc[6] == 0) stop("\narchex cannot be used without external.regressors!!\n", call. = FALSE)
		if(modelinc[6] < modelinc[20]) stop("\narchex cannot be greater than number of external.regressors!!\n", call. = FALSE)
	}
	
	maxOrder = max(modelinc[c(2,3,8,9)])
	modelnames = .expand.model(modelinc)
	
	
	fmodel = NULL
	if(modeldesc$vmodel == "fGARCH"){
		valid.submodels = c("GARCH","AVGARCH","NGARCH","NAGARCH","TGARCH","GJRGARCH","APARCH","ALLGARCH")
		submodel = modeldesc$vsubmodel
		if(!any(submodel == valid.submodels)) stop("not a valid fmodel name for fGARCH specification. See documentation.")
		fspec = .fgarchModel(submodel)
		fmodel$fpars = fspec$parameters
		# lambda, delta, fb, fc, fk
		fmodel$finclude = fspec$indicator
		fmodel$fgarchtype = fspec$garchtype
		fmodel$fbounds = .fmodelBounds(submodel)
	}
	
	
	pos = 1
	pos.matrix = matrix(0, ncol = 3, nrow = 21)
	colnames(pos.matrix) = c("start", "stop", "include")
	rownames(pos.matrix) = c("mu", "ar", "ma", "arfima", "archm", "mxreg", "omega", "alpha",
			"beta", "gamma", "eta1", "eta2", "delta", "lambda", "vxreg", "skew", "shape",
			"ghlambda", "aux", "aux", "aux")
	
	# check if there are starting or fixed
	# check that it is included in the optimization
	# check that it is needed in the model


	if(modeldesc$vmodel == "fGARCH"){
		for(i in 1:21){
			if( modelinc[i] > 0 ){
				if(i == 11 && fmodel$finclude[3]==1){
					pos.matrix[i,1:3] = c(pos, pos+modelinc[i]-1, 1)
					pos = max(pos.matrix[1:i,2]+1)
				} else if(i == 12 && fmodel$finclude[4]==1){
					pos.matrix[i,1:3] = c(pos, pos+modelinc[i]-1, 1)
					pos = max(pos.matrix[1:i,2]+1)
				} else if(i == 13 && fmodel$finclude[2]==1){
					pos.matrix[i,1:3] = c(pos, pos+modelinc[i]-1, 1)
					pos = max(pos.matrix[1:i,2]+1)
				} else if(i == 14 && fmodel$finclude[1]==1){
					pos.matrix[i,1:3] = c(pos, pos+modelinc[i]-1, 1)
					pos = max(pos.matrix[1:i,2]+1)
				} else{
					pos.matrix[i,1:3] = c(pos, pos+modelinc[i]-1, 1)
					pos = max(pos.matrix[1:i,2]+1)
				}
			}
		}
	} else{
		for(i in 1:21){
			if( modelinc[i] > 0 ){
				pos.matrix[i,1:3] = c(pos, pos+modelinc[i]-1, 1)
				pos = max(pos.matrix[1:i,2]+1)
			}
		}
	}
	nn = length(modelnames)
	modelmatrix = matrix(0, ncol = 3, nrow = nn)
	rownames(modelmatrix) = modelnames
	colnames(modelmatrix) = c("opt", "fixed", "start")
	fixed.names = names(fixed.pars)
	fp = charmatch(fixed.names, modelnames)
	
	if(!is.null(fixed.names) && any(!is.na(fp))){
		fixed = fp[!is.na(fp)]
		modelmatrix[fixed,2] = 1
		fz = charmatch(modelnames, fixed.names)
		fz = fz[!is.na(fz)]
		fixed.pars = fixed.pars[fz]
		names(fixed.pars) = fixed.names[fz]
	} else{
		fixed.pars = NULL
	}
	modelmatrix[,1] = 1 - modelmatrix[,2]
	start.names = names(start.pars)
	sp = charmatch(start.names, modelnames)
	if(!is.null(start.names) && any(!is.na(sp))){
		start = sp[!is.na(sp)]
		modelmatrix[start,3] = 1
		sz = charmatch(modelnames, start.names)
		sz = sz[!is.na(sz)]
		start.pars = start.pars[sz]
	} else{
		start.pars = NULL
	}
	
	
	##################################################################
	# Parameter Matrix
	mm = sum(modelinc[c(2,3,6,8,9,10,11,12,15)])
	mm = mm - length( which(modelinc[c(2,3,6,8,9,10,11,12,15)]>0) )
	pars = matrix(0, ncol = 6, nrow = 18 + mm)
	colnames(pars) = c("Level", "Fixed", "Include", "Estimate", "LB", "UB")
	pidx = matrix(NA, nrow = 18, ncol = 2)
	colnames(pidx) = c("begin", "end")
	rownames(pidx) =  c("mu", "ar", "ma", "arfima", "archm", "mxreg", "omega", "alpha",
			"beta", "gamma", "eta1", "eta2", "delta", "lambda", "vxreg", "skew", "shape",
			"ghlambda")
	fixed.names = names(fixed.pars)
	pnames = NULL
	nx = 0
	if(pos.matrix[1,3]==1){
		pars[1, 3] = 1
		pars[1, 1] = 0
		if(any(substr(fixed.names, 1, 2)=="mu")) pars[1,2] = 1 else pars[1,4] = 1

	}
	pidx[1,1] = 1
	pidx[1,2] = 1
	pnames = c(pnames, "mu")
	nx = 1
	pn = 1
	pidx[2,1] = 2	
	if(pos.matrix[2,3] == 1){
		pn = length( seq(pos.matrix[2,1], pos.matrix[2,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("ar", i, sep="")
			# TODO: FIX "ar" now might conflict with "arfima" when using substr
			if(any(substr(fixed.names, 1, nchar(nnx))==nnx)) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "ar")
	}
	pidx[2,2] = 1+pn
	
	nx = nx + pn
	pn = 1
	pidx[3,1] = nx+1
	if(pos.matrix[3,3] == 1){
		pn = length( seq(pos.matrix[3,1], pos.matrix[3,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("ma", i, sep="")
			if(any(substr(fixed.names, 1, nchar(nnx))==nnx)) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "ma")
	}
	pidx[3,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[4,1] = nx+1
	if(pos.matrix[4,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(substr(fixed.names, 1, 6)=="arfima")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "arfima")
	pidx[4,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[5,1] = nx+1
	if(pos.matrix[5,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(substr(fixed.names, 1, 5)=="archm")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "archm")
	pidx[5,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[6,1] = nx+1
	if(pos.matrix[6,3]==1){
		pn = length( seq(pos.matrix[6,1], pos.matrix[6,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("mxreg", i, sep="")
			if(any(substr(fixed.names, 1, nchar(nnx))==nnx)) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "mxreg")
	}
	pidx[6,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[7,1] = nx+1
	if(pos.matrix[7,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(substr(fixed.names, 1, 5)=="omega")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "omega")
	pidx[7,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[8,1] = nx+1
	if(pos.matrix[8,3]==1){
		pn = length( seq(pos.matrix[8,1], pos.matrix[8,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("alpha", i, sep="")
			if(any(substr(fixed.names, 1, nchar(nnx))==nnx)) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "alpha")
	}
	pidx[8,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[9,1] = nx+1
	if(pos.matrix[9,3]==1){
		pn = length( seq(pos.matrix[9,1], pos.matrix[9,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("beta", i, sep="")
			if(any(substr(fixed.names, 1, nchar(nnx))==nnx)) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
		#-------------------------------------------
		# special consideration for the iGARCH model
		#-------------------------------------------
		if(modeldesc$vmodel == "iGARCH"){
			# last beta not estimated
			pars[nx+pn, 4] = 0
			nnx = paste("beta", pn, sep="")
			# do not allow the last beta to be fixed
			if(any(substr(fixed.names, 1, nchar(nnx))==nnx)) pars[(nx+pn), 2] = 0
		}
	} else{
		pnames = c(pnames, "beta")
	}
	pidx[9,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[10,1] = nx+1
	
	if(pos.matrix[10,3]==1){
		pn = length( seq(pos.matrix[10,1], pos.matrix[10,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("gamma", i, sep="")
			if(any(substr(fixed.names, 1, nchar(nnx))==nnx)) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "gamma")
	}
	pidx[10,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[11,1] = nx+1
	
	if(pos.matrix[11,3]==1){
		pn = length( seq(pos.matrix[11,1], pos.matrix[11,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("eta1", i, sep="")
			if(any(substr(fixed.names, 1, nchar(nnx))==nnx)) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "eta1")
	}
	pidx[11,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[12,1] = nx+1
	
	if(pos.matrix[12,3]==1){
		pn = length( seq(pos.matrix[12,1], pos.matrix[12,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("eta2", i, sep="")
			if(any(substr(fixed.names, 1, nchar(nnx))==nnx)) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "eta2")
	}
	pidx[12,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[13,1] = nx+1
	
	if(pos.matrix[13,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(substr(fixed.names, 1, 5)=="delta")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	} else{
		#-------------------------------------------
		# special consideration for the fGARCH model
		#-------------------------------------------
		if(modeldesc$vmodel == "fGARCH")
		{
			pars[nx+pn, 3] = 1
			pars[nx+pn, 1] = fmodel$fpars$delta
			pars[nx+pn, 4] = pars[nx+pn, 2] = 0
		}
	}
	pidx[13,2] = nx+pn
	
	pnames = c(pnames, "delta")
	
	nx = nx + pn
	pn = 1
	pidx[14,1] = nx+1
	
	if(pos.matrix[14,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(substr(fixed.names, 1, 6)=="lambda")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	} else{
		#-------------------------------------------
		# special consideration for the fGARCH model
		#-------------------------------------------
		if(modeldesc$vmodel == "fGARCH")
		{
			pars[nx+pn, 3] = 1
			pars[nx+pn, 1] = fmodel$fpars$lambda
			pars[nx+pn, 4] = pars[nx+pn, 2] = 0
		}
	}
	pidx[14,2] = nx+pn
	
	pnames = c(pnames, "lambda")
	
	nx = nx + pn
	pn = 1
	pidx[15,1] = nx+1
	
	if(pos.matrix[15,3]==1){
		pn = length( seq(pos.matrix[15,1], pos.matrix[15,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("vxreg", i, sep="")
			if(any(substr(fixed.names, 1, nchar(nnx))==nnx)) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "vxreg")
		
	}
	pidx[15,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[16,1] = nx+1
	
	if(pos.matrix[16,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(substr(fixed.names, 1, 4)=="skew")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pidx[16,2] = nx+pn
	
	pnames = c(pnames, "skew")
	
	nx = nx + pn
	pn = 1
	pidx[17,1] = nx+1
	
	if(pos.matrix[17,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(substr(fixed.names, 1, 5)=="shape")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "shape")
	pidx[17,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[18,1] = nx+1
	
	if(pos.matrix[18,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(substr(fixed.names, 1, 8)=="ghlambda")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pidx[18,2] = nx+pn
	
	
	# Once more for fgarch model pars
	pnames = c(pnames, "ghlambda")
	rownames(pars) = pnames
	
	zf = match(fixed.names, rownames(pars))
	if( length(zf)>0 ) pars[zf, 1] = unlist(fixed.pars)
		
	model = list(modelinc = modelinc, modeldesc = modeldesc, modeldata = modeldata, pars = pars, 
			start.pars = start.pars, fixed.pars = fixed.pars, maxOrder = maxOrder, 
			pos.matrix = pos.matrix, fmodel = fmodel, pidx = pidx)
	ans = new("uGARCHspec", model = model)
	
	return(ans)
}

setMethod(f = "ugarchspec", definition = .ugarchspec)

# extract spec from fit object
getspec = function(object)
{
	UseMethod("getspec")
}

.getspec = function(object)
{
	spec = ugarchspec(variance.model = list(model = object@model$modeldesc$vmodel, garchOrder = c(object@model$modelinc[8],object@model$modelinc[9]), 
				submodel = object@model$modeldesc$vsubmodel, external.regressors = object@model$modeldata$vexdata), 
		mean.model = list(armaOrder = c(object@model$modelinc[2],object@model$modelinc[3]), 
				include.mean = object@model$modelinc[1], 
				archm = ifelse(object@model$modelinc[5]>0,TRUE,FALSE), archpow = object@model$modelinc[5], 
				arfima = object@model$modelinc[4], external.regressors = object@model$modeldata$mexdata,
				archex = object@model$modelinc[20]), 
		distribution.model = object@model$modeldesc$distribution, start.pars  = object@model$start.pars, 
		fixed.pars = object@model$fixed.pars)
	return(spec)
}

setMethod(f = "getspec", signature(object = "uGARCHfit"), definition = .getspec)


# internal function:
.model2spec = function(pars, model, type = "GARCH"){
	if(type == "GARCH"){
		ans = ugarchspec(variance.model = list(model = model$modeldesc$vmodel, garchOrder = c(model$modelinc[8],model$modelinc[9]), 
						submodel = model$modeldesc$vsubmodel, external.regressors = model$modeldata$vexdata), 
				mean.model = list(armaOrder = c(model$modelinc[2],model$modelinc[3]), 
						include.mean = model$modelinc[1], 
						archm = ifelse(model$modelinc[5]>0,TRUE,FALSE), archpow = model$modelinc[5], 
						arfima = model$modelinc[4], external.regressors = model$modeldata$mexdata,
						archex = model$modelinc[20]), 
				distribution.model = model$modeldesc$distribution)
		setfixed(ans)<-pars
	} else{
		ans = arfimaspec( 
				mean.model = list(armaOrder = c(model$modelinc[2],model$modelinc[3]), 
						include.mean = model$modelinc[1], 
						arfima = model$modelinc[4], external.regressors = model$modeldata$mexdata), 
				distribution.model = model$modeldesc$distribution)
		setfixed(ans)<-pars
	}
	return(ans)
}


# Set methods replace any existing fixed or starting parameters already present
setGeneric("setfixed<-", function(object, value){standardGeneric("setfixed<-")})

.setfixed = function(object, value){
	# get parameter values
	model = object@model
	ipars = model$pars
	pars = unlist(value)
	names(pars) = parnames = tolower(names(pars))
	# included parameters in model
	modelnames = rownames(ipars[which(ipars[,4]==1 | ipars[,2]==1), ])
	inc = NULL
	for(i in seq_along(parnames)){
		if(is.na(match(parnames[i], modelnames))){
			warning( (paste("Unrecognized Parameter in Fixed Values: ", parnames[i], "...Ignored", sep = "")))
		} else{
			inc = c(inc, i)
		}
	}
	fixed.pars = pars[inc]
	names(fixed.pars) = tolower(names(pars[inc]))
	# set parameter values
	tmp = ugarchspec(variance.model = list(model = model$modeldesc$vmodel, garchOrder = c(model$modelinc[8], model$modelinc[9]), 
					submodel = model$modeldesc$vsubmodel, external.regressors = model$modeldata$vexdata,
					variance.targeting = ifelse(model$modelinc[7]==0, TRUE, FALSE)), 
			mean.model = list(armaOrder = c(model$modelinc[2], model$modelinc[3]), 
					include.mean = model$modelinc[1], 
					archm = ifelse(model$modelinc[5]>0,TRUE,FALSE), archpow = model$modelinc[5], 
					arfima = model$modelinc[4], external.regressors = model$modeldata$mexdata,
					archex = model$modelinc[20]), 
			distribution.model = model$modeldesc$distribution, start.pars  = model$start.pars, 
			fixed.pars = as.list(fixed.pars))
	return(tmp)
}
setReplaceMethod(f="setfixed", signature= c(object = "uGARCHspec", value = "vector"), definition = .setfixed)


setGeneric("setstart<-", function(object, value){standardGeneric("setstart<-")})

.setstart = function(object, value){
	# get parameter values
	model = object@model
	ipars = model$pars
	pars = unlist(value)
	names(pars) = parnames = tolower(names(pars))
	# included parameters in model
	modelnames = rownames(ipars[which(ipars[,4]==1 | ipars[,2]==1), ])
	inc = NULL
	for(i in seq_along(parnames)){
		if(is.na(match(parnames[i], modelnames))){
			warning( (paste("Unrecognized Parameter in Start Values: ", parnames[i], "...Ignored", sep = "")))
		} else{
			inc = c(inc, i)
		}
	}
	start.pars = pars[inc]
	names(start.pars) = tolower(names(pars[inc]))
	# set parameter values
	tmp = ugarchspec(variance.model = list(model = model$modeldesc$vmodel, garchOrder = c(model$modelinc[8], model$modelinc[9]), 
					submodel = model$modeldesc$vsubmodel, external.regressors = model$modeldata$vexdata,
					variance.targeting = ifelse(model$modelinc[7]==0, TRUE, FALSE)), 
			mean.model = list(armaOrder = c(model$modelinc[2], model$modelinc[3]), 
					include.mean = model$modelinc[1], 
					archm = ifelse(model$modelinc[5]>0,TRUE,FALSE), archpow = model$modelinc[5], 
					arfima = model$modelinc[4], external.regressors = model$modeldata$mexdata,
					archex = model$modelinc[20]), 
			distribution.model = model$modeldesc$distribution, fixed.pars  = model$fixed.pars, 
			start.pars = as.list(start.pars))
	return(tmp)
}

setReplaceMethod(f="setstart", signature= c(object = "uGARCHspec", value = "vector"), definition = .setstart)

.checkallfixed = function( spec ){
	# check that a given spec with fixed parameters
	model = spec@model
	pars = model$pars
	pnames = rownames(pars)
	estpars = pnames[as.logical(pars[,2] * pars[,3] + pars[,3] * pars[,4])]
	return( estpars )
	
}
#----------------------------------------------------------------------------------
# univariate model dispatch methods
#----------------------------------------------------------------------------------

.ugarchfit = function(spec, data, out.sample = 0, solver = "solnp", solver.control = list(), 
		fit.control = list(stationarity = 1, fixed.se = 0, scale = 0), ...)
{
	return( switch(spec@model$modeldesc$vmodel,
					sGARCH = .sgarchfit(spec = spec, data = data, out.sample = out.sample, solver = solver, 
							solver.control = solver.control, fit.control = fit.control, ...),
					iGARCH = .igarchfit(spec = spec, data = data, out.sample = out.sample, solver = solver, 
							solver.control = solver.control, fit.control = fit.control, ...),
					eGARCH = .egarchfit(spec = spec, data = data, out.sample = out.sample, solver = solver, 
							solver.control = solver.control, fit.control = fit.control, ...),
					gjrGARCH = .gjrgarchfit(spec = spec, data = data, out.sample = out.sample, solver = solver, 
							solver.control = solver.control, fit.control = fit.control, ...),
					apARCH = .aparchfit(spec = spec, data = data, out.sample = out.sample, solver = solver, 
							solver.control = solver.control, fit.control = fit.control, ...),
					fGARCH = .fgarchfit(spec = spec, data = data, out.sample = out.sample, solver = solver, 
							solver.control = solver.control, fit.control = fit.control, ...)) )
}

.ugarchfilter = function(spec, data, out.sample = 0, n.old = NULL, ...)
{
	return( switch(spec@model$modeldesc$vmodel,
					sGARCH = .sgarchfilter(spec = spec, data = data, out.sample = out.sample, n.old = n.old, ...),
					iGARCH = .igarchfilter(spec = spec, data = data, out.sample = out.sample, n.old = n.old, ...),
					eGARCH = .egarchfilter(spec = spec, data = data, out.sample = out.sample, n.old = n.old, ...),
					gjrGARCH = .gjrgarchfilter(spec = spec, data = data, out.sample = out.sample, n.old = n.old, ...),
					apARCH = .aparchfilter(spec = spec, data = data, out.sample = out.sample, n.old = n.old, ...),
					fGARCH = .fgarchfilter(spec = spec, data = data, out.sample = out.sample, n.old = n.old, ...)) )
}

.ugarchforecast1 = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0, 
		external.forecasts = list(mregfor = NULL, vregfor = NULL), ...)
{
	return( switch(fitORspec@model$modeldesc$vmodel,
					sGARCH = .sgarchforecast(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll, 
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					iGARCH = .igarchforecast(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll, 
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					eGARCH = .egarchforecast(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll, 
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					gjrGARCH = .gjrgarchforecast(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll, 
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					apARCH = .aparchforecast(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll, 
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					fGARCH = .fgarchforecast(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll, 
							out.sample = out.sample, external.forecasts = external.forecasts, ...)) )
}

.ugarchforecast2 = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0, 
		external.forecasts = list(mregfor = NULL, vregfor = NULL), ...)
{
	return( switch(fitORspec@model$modeldesc$vmodel,
					sGARCH = .sgarchforecast2(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll, 
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					iGARCH = .igarchforecast2(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll, 
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					eGARCH = .egarchforecast2(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll, 
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					gjrGARCH = .gjrgarchforecast2(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll, 
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					apARCH = .aparchforecast2(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll, 
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					fGARCH = .fgarchforecast2(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll, 
							out.sample = out.sample, external.forecasts = external.forecasts, ...)) )
}

.ugarchsim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, startMethod = c("unconditional","sample"), 
		presigma = NA, prereturns = NA, preresiduals = NA, rseed = NA,  custom.dist = list(name = NA, distfit = NA), 
		mexsimdata = NULL, vexsimdata = NULL, ...)
{
	return( switch(fit@model$modeldesc$vmodel,
					sGARCH = .sgarchsim(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
							startMethod = startMethod, presigma = presigma, prereturns = prereturns, 
							preresiduals = preresiduals, rseed = rseed,  custom.dist = custom.dist, 
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, ...),
					iGARCH = .igarchsim(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
							startMethod = startMethod, presigma = presigma, prereturns = prereturns, 
							preresiduals = preresiduals, rseed = rseed,  custom.dist = custom.dist, 
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, ...),
					eGARCH = .egarchsim(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
							startMethod = startMethod, presigma = presigma, prereturns = prereturns, 
							preresiduals = preresiduals, rseed = rseed,  custom.dist = custom.dist, 
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, ...),
					gjrGARCH = .gjrgarchsim(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
							startMethod = startMethod, presigma = presigma, prereturns = prereturns, 
							preresiduals = preresiduals, rseed = rseed,  custom.dist = custom.dist, 
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, ...),
					apARCH = .aparchsim(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
							startMethod = startMethod, presigma = presigma, prereturns = prereturns, 
							preresiduals = preresiduals, rseed = rseed,  custom.dist = custom.dist, 
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, ...),
					fGARCH = .fgarchsim(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
							startMethod = startMethod, presigma = presigma, prereturns = prereturns, 
							preresiduals = preresiduals, rseed = rseed,  custom.dist = custom.dist, 
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, ...) ) )
}

.ugarchpath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA, prereturns = NA, preresiduals = NA, 
		rseed = NA, custom.dist = list(name = NA, distfit = NA), mexsimdata = NULL,  vexsimdata = NULL, ...)
{
	return( switch(spec@model$modeldesc$vmodel,
					sGARCH = .sgarchpath(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
							presigma = presigma, prereturns = prereturns, preresiduals = preresiduals, 
							rseed = rseed,  custom.dist = custom.dist, mexsimdata = mexsimdata, 
							vexsimdata = vexsimdata, ...),
					iGARCH = .igarchpath(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
							presigma = presigma, prereturns = prereturns, preresiduals = preresiduals, 
							rseed = rseed,  custom.dist = custom.dist, mexsimdata = mexsimdata, 
							vexsimdata = vexsimdata, ...),
					eGARCH = .egarchpath(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
							presigma = presigma, prereturns = prereturns, preresiduals = preresiduals, 
							rseed = rseed,  custom.dist = custom.dist, mexsimdata = mexsimdata, 
							vexsimdata = vexsimdata, ...),
					gjrGARCH = .gjrgarchpath(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
							presigma = presigma, prereturns = prereturns, preresiduals = preresiduals, 
							rseed = rseed,  custom.dist = custom.dist, mexsimdata = mexsimdata, 
							vexsimdata = vexsimdata, ...),
					apARCH = .aparchpath(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
							presigma = presigma, prereturns = prereturns, preresiduals = preresiduals, 
							rseed = rseed,  custom.dist = custom.dist, mexsimdata = mexsimdata, 
							vexsimdata = vexsimdata, ...),
					fGARCH = .fgarchpath(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
							presigma = presigma, prereturns = prereturns, preresiduals = preresiduals, 
							rseed = rseed,  custom.dist = custom.dist, mexsimdata = mexsimdata, 
							vexsimdata = vexsimdata, ...) ) )
}


#----------------------------------------------------------------------------------
# univariate filter method
#----------------------------------------------------------------------------------
ugarchfilter = function(spec, data, out.sample = 0, n.old = NULL, ...)
{
	UseMethod("ugarchfilter")
}

setMethod("ugarchfilter", signature(spec = "uGARCHspec"), .ugarchfilter)
#----------------------------------------------------------------------------------
# univariate fit method
#----------------------------------------------------------------------------------
ugarchfit = function(spec, data, out.sample = 0, solver = "solnp", solver.control = list(), 
		fit.control = list(stationarity = 1, fixed.se = 0, scale = 0), ...)
{
	UseMethod("ugarchfit")
}

setMethod("ugarchfit", signature(spec = "uGARCHspec"), .ugarchfit)
#----------------------------------------------------------------------------------
# univariate forecast method
#----------------------------------------------------------------------------------
ugarchforecast = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0, 
		external.forecasts = list(mregfor = NULL, vregfor = NULL), ...)
{
	UseMethod("ugarchforecast")
}

setMethod("ugarchforecast", signature(fitORspec = "uGARCHfit"), .ugarchforecast1)

#setMethod("ugarchforecast", signature(fitORspec = "anstGARCHfit"), .anstgarchforecast)
#alternative dispath method:
# we use the fitORspec rather than a method with fit, spec and data with missing
# methods as this requires implicit declaration of arguments

setMethod("ugarchforecast", signature(fitORspec = "uGARCHspec"), .ugarchforecast2)

#----------------------------------------------------------------------------------
# univariate simulation method
#----------------------------------------------------------------------------------

ugarchsim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1,
		startMethod = c("unconditional","sample"), 
		presigma = NA, prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA), mexsimdata = NULL, 
		vexsimdata = NULL, ...)
{
	UseMethod("ugarchsim")
}


setMethod("ugarchsim", signature(fit = "uGARCHfit"), .ugarchsim)
#----------------------------------------------------------------------------------
# univariate path simulation method
#----------------------------------------------------------------------------------

ugarchpath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1,
		presigma = NA, prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA), mexsimdata = NULL, 
		vexsimdata = NULL, ...)
{
	UseMethod("ugarchpath")
}


setMethod("ugarchpath", signature(spec = "uGARCHspec"), .ugarchpath)
#----------------------------------------------------------------------------------
# univariate garch roll
#----------------------------------------------------------------------------------
# methods to recursively predict/filter/compare with refitting at every N points.
ugarchroll = function(spec,  data, n.ahead = 1, forecast.length = 500, 
		refit.every = 25, refit.window = c("recursive", "moving"), parallel = FALSE, 
		parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2), 
		solver = "solnp", fit.control = list(), solver.control = list() ,
		calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), ...)
{
	setMethod("ugarchroll")
}
setMethod("ugarchroll", signature(spec = "uGARCHspec"),  definition = .rollfdensity)
#----------------------------------------------------------------------------------
# univariate garch parameter distribution
#----------------------------------------------------------------------------------
ugarchdistribution = function(fitORspec, n.sim = 2000, n.start = 1, 
		m.sim = 100, recursive = FALSE, recursive.length = 6000, recursive.window = 1000, 
		presigma = NA, prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA), mexsimdata = NULL, 
		vexsimdata = NULL, fit.control = list(), solver = "solnp", 
		solver.control = list(), parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2), ...)
{
	setMethod("ugarchdistribution")
}
setMethod("ugarchdistribution", signature(fitORspec = "uGARCHfit"), .ugarchdistribution)
setMethod("ugarchdistribution", signature(fitORspec = "uGARCHspec"), .ugarchdistribution)


#----------------------------------------------------------------------------------
# univariate garch bootstrap based forecast distribution
#----------------------------------------------------------------------------------
ugarchboot = function(fitORspec, data = NULL, method = c("Partial", "Full"), n.ahead = 10, 
		n.bootfit = 100, n.bootpred = 500, out.sample = 0, rseed = NA, solver = "solnp", 
		solver.control = list(), fit.control = list(), external.forecasts =  list(mregfor = NULL, 
				vregfor = NULL), parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2))
{
	setMethod("ugarchboot")
}

setMethod("ugarchboot", signature(fitORspec = "uGARCHfit"), .ugarchbootfit)
setMethod("ugarchboot", signature(fitORspec = "uGARCHspec"), .ugarchbootspec)

#----------------------------------------------------------------------------------
# univariate plot method / seperate for fit,sim and forecast
#----------------------------------------------------------------------------------
setMethod(f = "plot", signature(x = "uGARCHfit", y = "missing"), .plotgarchfit)

setMethod(f = "plot", signature(x = "uGARCHfilter", y = "missing"), .plotgarchfilter)

setMethod(f = "plot", signature(x = "uGARCHsim", y = "missing"), .plotgarchsim)

setMethod(f = "plot", signature(x = "uGARCHforecast", y = "missing"), .plotgarchforecast)

setMethod(f = "plot", signature(x = "uGARCHpath", y = "missing"), .plotgarchpath)

setMethod(f = "plot", signature(x = "uGARCHroll", y = "missing"), .plotgarchroll)

setMethod(f = "plot", signature(x = "uGARCHdistribution", y = "missing"), .plotgarchdist)

setMethod(f = "plot", signature(x = "uGARCHboot", y = "missing"), .plotgarchboot)

#----------------------------------------------------------------------------------
# univariate show method / seperate for fit,sim and forecast
#----------------------------------------------------------------------------------

# spec show
setMethod("show",
		signature(object = "uGARCHspec"),
		function(object){
			model = object@model
			vmodel = model$modeldesc$vmodel
			modelinc = model$modelinc
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*       GARCH Model Spec          *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			cat("\n\nConditional Variance Dynamics \t")
			cat(paste("\n------------------------------------\n",sep=""))
			cat(paste("GARCH Model\t\t: ", vmodel, "(",modelinc[8],",",modelinc[9],")\n", sep = ""))
			if(vmodel == "fGARCH"){
				cat(paste("fGARCH Sub-Model\t: ", model$modeldesc$vsubmodel, "\n", sep = ""))
			}
			cat("Variance Targeting\t:", as.logical(1-modelinc[7]), "\n")
			if(modelinc[15]>0) cat(paste("Exogenous Regressor Dimension: ", modelinc[15], "\n",sep = ""))
			cat("\nConditional Mean Dynamics")
			cat(paste("\n------------------------------------\n",sep=""))
			cat("Mean Model\t\t: ARFIMA(", modelinc[2],",", ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Include Mean\t\t:", as.logical(modelinc[1]),"\n")
			cat("GARCH-in-Mean\t\t:", as.logical(modelinc[5]),"\n")
			if(modelinc[6]>0) cat(paste("Exogenous Regressor Dimension: ", modelinc[6],"\n",sep=""))
			cat("\nConditional Distribution")
			cat(paste("\n------------------------------------\n",sep=""))
			cat("Distribution\t: ", model$modeldesc$distribution,"\n")
			cat("Includes Skew\t: ", as.logical(modelinc[16]),"\n")
			cat("Includes Shape\t: ", as.logical(modelinc[17]),"\n")
			cat("Includes Lambda\t: ", as.logical(modelinc[18]),"\n\n")
			invisible(object)
		})

# fit show
setMethod("show",
		signature(object = "uGARCHfit"),
		function(object){
			vmodel = object@model$modeldesc$vmodel
			model = object@model
			modelinc = object@model$modelinc
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*          GARCH Model Fit        *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			cat("\n\nConditional Variance Dynamics \t")
			cat(paste("\n-----------------------------------", sep = ""))
			cat(paste("\nGARCH Model\t: ", vmodel,"(", modelinc[8], ",", modelinc[9], ")\n", sep=""))
			if(vmodel == "fGARCH"){
				cat(paste("fGARCH Sub-Model\t: ", model$modeldesc$vsubmodel, "\n", sep = ""))
			}
			cat("Mean Model\t: ARFIMA(", modelinc[2],",",ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Distribution\t:", model$modeldesc$distribution,"\n")
			if(object@fit$convergence == 0){
				cat("\nOptimal Parameters")
				cat(paste("\n------------------------------------\n",sep=""))
				print(round(object@fit$matcoef,6), digits = 5)
				cat("\nRobust Standard Errors:\n")
				print(round(object@fit$robust.matcoef,6), digits = 5)
				if(!is.null(object@fit$hessian.message)){
					cat(paste("\n", object@fit$hessian.message))
				}
				cat("\nLogLikelihood :", object@fit$LLH, "\n")
				stdresid = object@fit$residuals/object@fit$sigma
				itest = .information.test(object@fit$LLH, nObs = object@model$modeldata$T, nPars = length(object@fit$coef))
				itestm = matrix(0, ncol = 1, nrow = 4)
				itestm[1,1] = itest$AIC
				itestm[2,1] = itest$BIC
				itestm[3,1] = itest$SIC
				itestm[4,1] = itest$HQIC
				colnames(itestm) = ""
				rownames(itestm) = c("Akaike", "Bayes", "Shibata", "Hannan-Quinn")
				cat("\nInformation Criteria")
				cat(paste("\n------------------------------------\n",sep=""))
				print(itestm,digits=5)
				cat("\nQ-Statistics on Standardized Residuals")
				cat(paste("\n------------------------------------\n",sep=""))
				tmp1 = .box.test(stdresid, p = 1, df = sum(modelinc[2:3]))
				print(tmp1, digits = 4)
				cat("\nH0 : No serial correlation\n")
				cat("\nQ-Statistics on Standardized Squared Residuals")
				cat(paste("\n------------------------------------\n",sep=""))
				tmp2 = .box.test(stdresid, p = 2, df = sum(modelinc[2:3]))
				print(tmp2, digits = 4)
				cat("\nARCH LM Tests")
				cat(paste("\n------------------------------------\n",sep=""))
				L2 = .archlmtest(stdresid, lags = 2)
				L5 = .archlmtest(stdresid, lags = 5)
				L10 = .archlmtest(stdresid, lags = 10)
				alm = matrix(0,ncol = 3,nrow = 3)
				alm[1,1:3] = c(L2$statistic, L2$parameter, L2$p.value)
				alm[2,1:3] = c(L5$statistic, L5$parameter, L5$p.value)
				alm[3,1:3] = c(L10$statistic, L10$parameter, L10$p.value)
				colnames(alm) = c("Statistic", "DoF", "P-Value")
				rownames(alm) = c("ARCH Lag[2]", "ARCH Lag[5]", "ARCH Lag[10]")
				print(alm,digits = 4)
				nyb = .nyblomTest(object)
				colnames(nyb$IndividualStat)<-""
				cat("\nNyblom stability test")
				cat(paste("\n------------------------------------\n",sep=""))
				cat("Joint Statistic: ",round(nyb$JointStat,4))
				cat("\nIndividual Statistics:")
				print(nyb$IndividualStat, digits = 4)
				cat("\nAsymptotic Critical Values (10% 5% 1%)")
				cat("\nJoint Statistic:     \t", round(nyb$JointCritical, 3))
				cat("\nIndividual Statistic:\t", round(nyb$IndividualCritical, 2))
				cat("\n\n")
				cat("Sign Bias Test")
				cat(paste("\n------------------------------------\n",sep=""))
				sgtest = signbias(object)
				print(sgtest, digits = 4)
				cat("\n")
				cat("\nAdjusted Pearson Goodness-of-Fit Test:")
				cat(paste("\n------------------------------------\n",sep=""))
				gofm = gof(object,c(20, 30, 40, 50))
				print(gofm, digits = 4)
				cat("\n")
				cat("\nElapsed time :", object@fit$timer,"\n\n")
			} else{
				cat("\nConvergence Problem:")
				cat("\nSolver Message:", object@fit$message,"\n\n")
				
			}
			invisible(object)
		})

		
		
# filter show
setMethod("show",
		signature(object = "uGARCHfilter"),
		function(object){
			vmodel = object@model$modeldesc$vmodel
			model = object@model
			modelinc = object@model$modelinc
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\n*          GARCH Model Filter        *", sep = ""))
			cat(paste("\n*------------------------------------*", sep = ""))
			cat("\n\nConditional Variance Dynamics \t")
			cat(paste("\n--------------------------------------", sep = ""))
			cat(paste("\nGARCH Model\t: ", vmodel,"(", modelinc[8], ",", modelinc[9], ")\n", sep=""))
			if(vmodel == "fGARCH"){
				cat(paste("fGARCH Sub-Model\t: ", model$modeldesc$vsubmodel, "\n", sep = ""))
			}
			cat("Mean Model\t: ARFIMA(", modelinc[2],",",ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Distribution\t:", model$modeldesc$distribution,"\n")
			cat("\nFilter Parameters")
			cat(paste("\n---------------------------------------\n",sep=""))
			print(matrix(coef(object), ncol=1, dimnames = list(names(coef(object)), "")), digits = 5)
			cat("\nLogLikelihood :", object@filter$LLH, "\n")
			stdresid = object@filter$residuals/object@filter$sigma
			itest = .information.test(object@filter$LLH, nObs = object@model$modeldata$T, nPars = length(coef(object)))
			itestm = matrix(0, ncol = 1, nrow = 4)
			itestm[1,1] = itest$AIC
			itestm[2,1] = itest$BIC
			itestm[3,1] = itest$SIC
			itestm[4,1] = itest$HQIC
			colnames(itestm) = ""
			rownames(itestm) = c("Akaike", "Bayes", "Shibata", "Hannan-Quinn")
			cat("\nInformation Criteria")
			cat(paste("\n---------------------------------------\n",sep=""))
			print(itestm,digits=5)
			cat("\nQ-Statistics on Standardized Residuals")
			cat(paste("\n---------------------------------------\n",sep=""))
			tmp1 = .box.test(stdresid, p = 1, df = sum(modelinc[2:3]))
			print(tmp1, digits = 4)
			cat("\nH0 : No serial correlation\n")
			cat("\nQ-Statistics on Standardized Squared Residuals")
			cat(paste("\n---------------------------------------\n",sep=""))
			tmp2 = .box.test(stdresid, p = 2, df = sum(modelinc[2:3]))
			print(tmp2, digits = 4)
			cat("\nARCH LM Tests")
			cat(paste("\n---------------------------------------\n",sep=""))
			L2 = .archlmtest(stdresid, lags = 2)
			L5 = .archlmtest(stdresid, lags = 5)
			L10 = .archlmtest(stdresid, lags = 10)
			alm = matrix(0,ncol = 3,nrow = 3)
			alm[1,1:3] = c(L2$statistic, L2$parameter, L2$p.value)
			alm[2,1:3] = c(L5$statistic, L5$parameter, L5$p.value)
			alm[3,1:3] = c(L10$statistic, L10$parameter, L10$p.value)
			colnames(alm) = c("Statistic", "DoF", "P-Value")
			rownames(alm) = c("ARCH Lag[2]", "ARCH Lag[5]", "ARCH Lag[10]")
			print(alm,digits = 4)
			cat("\n\n")
			cat("Sign Bias Test")
			cat(paste("\n---------------------------------------\n",sep=""))
			sgtest = signbias(object)
			print(sgtest, digits = 4)
			cat("\n")
			cat("\nAdjusted Pearson Goodness-of-Fit Test:")
			cat(paste("\n---------------------------------------\n",sep=""))
			gofm = gof(object,c(20, 30, 40, 50))
			print(gofm, digits = 4)
			cat("\n")
			invisible(object)
		})
# sim show
setMethod("show",
			signature(object = "uGARCHsim"),
			function(object){
			vmodel = object@model$modeldesc$vmodel
			model = object@model
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\n*       GARCH Model Simulation       *", sep = ""))
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\nModel : ", vmodel,sep=""))
			if(vmodel == "fGARCH"){
				cat(paste("fGARCH Sub-Model\t: ", model$modeldesc$vsubmodel, "\n", sep = ""))
			}
			sim = object@simulation
			dates = object@model$dates
			sigma = sim$sigmaSim
			series = sim$seriesSim
			resids = sim$residSim
			m = dim(sigma)[2]
			N = dim(sigma)[1]
			cat(paste("\nHorizon: ",N))
			cat(paste("\nSimulations: ",m,"\n",sep=""))
			sd1 = apply(sigma^2, 2, FUN=function(x) mean(x))
			sd2 = apply(sigma^2, 2, FUN=function(x) range(x))	
			rx1 = apply(series, 2, FUN=function(x) mean(x))
			rx2 = apply(series, 2, FUN=function(x) range(x))
			T = object@model$modeldata$T
			actual = c(0,mean(object@model$modeldata$sigma^2), min(object@model$modeldata$sigma^2),
					max(object@model$modeldata$sigma^2), mean(object@model$modeldata$data[1:T]), 
					min(object@model$modeldata$data[1:T]), max(object@model$modeldata$data[1:T]))
			xspec = .model2spec(as.list(object@model$pars[object@model$pars[,3]==1,1]), object@model, type = "GARCH")
			setfixed(xspec)<-as.list(object@model$pars[which(object@model$pars[,3]==1),1])
			uv = uncvariance(xspec)
			um = uncmean(xspec)
			uncond = c(0, uv, NA, NA, um, NA, NA)
			dd = data.frame(Seed = object@seed, Sigma2.Mean = sd1, Sigma2.Min = sd2[1,],
					Sigma2.Max = sd2[2,], Series.Mean = rx1, Series.Min = rx2[1,],
					Series.Max = rx2[2,])
			meansim = apply(dd, 2, FUN = function(x) mean(x))
			meansim[1] = 0
			dd = rbind(dd, meansim, actual, uncond)
			rownames(dd) = c(paste("sim", 1:m, sep = ""), "Mean(All)", "Actual", "Unconditional")
			print(dd,digits = 3)
			cat("\n\n")
			})
			
# forecast show
setMethod("show",
		signature(object = "uGARCHforecast"),
		function(object){
			vmodel = object@model$modeldesc$vmodel
			model = object@model
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\n*       GARCH Model Forecast         *", sep = ""))
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\nModel: ", vmodel, sep = ""))
			if(vmodel == "fGARCH"){
				cat(paste("fGARCH Sub-Model: ", model$modeldesc$vsubmodel, "\n", sep = ""))
			}
			n.ahead = object@forecast$n.ahead
			cat(paste("\nHorizon: ", n.ahead, sep = ""))
			cat(paste("\nRoll Steps: ",object@forecast$n.roll, sep = ""))
			n.start = object@forecast$n.start
			if(n.start>0) infor = ifelse(n.ahead>n.start, n.start, n.ahead) else infor = 0
			cat(paste("\nOut of Sample: ", infor, "\n", sep = ""))
			cat("\n0-roll forecast: \n")
			zz = object@forecast$forecast[[1]]
			print(zz, digits = 4)
			cat("\n\n")
		})
		
# path show
setMethod("show",
			signature(object = "uGARCHpath"),
			function(object){
			vmodel = object@model$modeldesc$vmodel
			model = object@model
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\n*     GARCH Model Path Simulation    *", sep = ""))
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\nModel: ", vmodel, sep = ""))
			if(vmodel == "fGARCH"){
				cat(paste("fGARCH Sub-Model: ", model$modeldesc$vsubmodel, "\n", sep = ""))
			}
			sim = object@path
			sigma = sim$sigmaSim
			series = sim$seriesSim
			resids = sim$residSim
			m = dim(sigma)[2]
			N = dim(sigma)[1]
			cat(paste("\nHorizon: ", N))
			cat(paste("\nSimulations: ", m, "\n", sep = ""))
			sd1 = apply(sigma^2, 2, FUN = function(x) mean(x))
			sd2 = apply(sigma^2, 2, FUN = function(x) range(x))			
			rx1 = apply(series, 2, FUN = function(x) mean(x))
			rx2 = apply(series, 2, FUN = function(x) range(x))
			xspec = .model2spec(as.list(object@model$pars[object@model$pars[,3]==1,1]), object@model, type = "GARCH")
			setfixed(xspec)<-as.list(object@model$pars[which(object@model$pars[,3]==1),1])
			uv = uncvariance(xspec)
			um = uncmean(xspec)
			uncond = c(NA, uv, NA, NA, um, NA, NA)
			dd = data.frame(Seed = object@seed, Sigma2.Mean = sd1, Sigma2.Min = sd2[1,],
					Sigma2.Max = sd2[2,], Series.Mean = rx1, Series.Min = rx2[1,], 
					Series.Max = rx2[2,])
			meansim = apply(dd, 2, FUN = function(x) mean(x))
			meansim[1] = 0
			dd = rbind(dd, meansim, uncond)			
			rownames(dd) = c(paste("sim", 1:m, sep = ""), "Mean(All)", "Unconditional")
			print(dd, digits = 3)
			cat("\n\n")
			})

# distribution show
setMethod("show",
		signature(object = "uGARCHdistribution"),
		function(object){
			vmodel = object@model$modeldesc$vmodel
			vsubmodel = object@model$modeldesc$vsubmodel
			model = object@model
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\n*    GARCH Parameter Distribution    *", sep = ""))
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\nModel : ", vmodel, sep = ""))
			if(vmodel == "fGARCH"){
				cat(paste("\nSubModel : ", vsubmodel, sep = ""))
			}
			cat(paste("\nNo. Paths (m.sim) : ", object@dist$details$m.sim, sep = ""))
			cat(paste("\nLength of Paths (n.sim) : ", object@dist$details$n.sim, sep = ""))
			cat(paste("\nRecursive : ", object@dist$details$recursive, sep = ""))
			if(object@dist$details$recursive){
				cat(paste("\nRecursive Length : ", object@dist$details$recursive.length, sep = ""))
				cat(paste("\nRecursive Window : ", object@dist$details$recursive.window, sep = ""))
			}
			cat("\n\n")
			cat("Coefficients: True vs Simulation Mean (Window-n)\n")
			nwindows = object@dist$details$nwindows
			nm = object@dist$details$n.sim + (0:(nwindows-1))*object@dist$details$recursive.window
			ns = matrix(0, ncol = dim(object@truecoef)[1], nrow = nwindows)
			for(i in 1:nwindows){
				ns[i,] = apply(as.data.frame(object, window = i), 2, FUN = function(x) mean(x, na.rm = T))
			}
			ns = rbind(object@truecoef[,1], ns)
			colnames(ns) = rownames(object@truecoef)
			rownames(ns) = c("true-coef",paste("window-", nm, sep=""))
			print(as.data.frame(ns), digits=5)
			for(i in 1:nwindows){
				if(any(object@dist[[i]]$convergence==1)) n = length(which(object@dist[[i]]$convergence==1)) else n = 0
				if(n>0) cat(paste("\nwindow-",nm[i]," no. of non-converged fits: ", n, "\n",sep=""))
			}
			cat("\n\n")
		})
			

# boot show
setMethod("show",
		signature(object = "uGARCHboot"),
		function(object){
			vmodel = object@model$modeldesc$vmodel
			vsubmodel = object@model$modeldesc$vsubmodel
			cat(paste("\n*-----------------------------------*", sep = ""))
			cat(paste("\n*     GARCH Bootstrap Forecast      *", sep = ""))
			cat(paste("\n*-----------------------------------*", sep = ""))
			cat(paste("\nModel : ", vmodel, sep = ""))
			if(vmodel == "fGARCH"){
				cat(paste("\nSubModel : ", vsubmodel, sep = ""))
			}
			cat(paste("\nn.ahead : ", object@model$n.ahead, sep = ""))
			cat(paste("\nBootstrap method: ",object@model$type))
			forc = object@forc@forecast$forecasts[[1]]
			zs = rbind(as.data.frame(object, which = "sigma", type = "summary"),  forc[,"sigma"])
			zr = rbind(as.data.frame(object, which = "series", type = "summary"), forc[,"series"])
			rownames(zr)[6] = rownames(zs)[6] = "forecast"
			hh = min(object@model$n.ahead, 10)
			cat("\n\nSeries (summary):\n")
			print(head(round(t(zr), 6), hh), digits = 5)
			cat(".....................\n")
			cat("\nSigma (summary):\n")
			print(head(round(t(zs), 6),hh), digits = 5)
			cat(".....................")
			cat("\n\n")
		})

setMethod("show",
		signature(object = "uGARCHroll"),
		function(object){
			cat(paste("\n*-------------------------------------*", sep = ""))
			cat(paste("\n*              GARCH Roll             *", sep = ""))
			cat(paste("\n*-------------------------------------*", sep = ""))
			N = object@roll$n.refit
			model = object@model
			modelinc = model$modelinc
			vmodel = object@model$modeldesc$vmodel
			cat("\nNo.Refits\t\t:", N)
			cat("\nRefit Horizon\t:", object@roll$refit.every)
			cat("\nForecast Horizon:", object@roll$n.ahead)
			cat(paste("\nGARCH Model\t\t: ", vmodel, "(",modelinc[8],",",modelinc[9],")\n", sep = ""))
			if(vmodel == "fGARCH"){
				cat(paste("fGARCH Sub-Model\t: ", model$modeldesc$vsubmodel, "\n", sep = ""))
			}
			cat("Distribution\t:", model$modeldesc$distribution,"\n")
			mp = sum(model$modelinc[1:18])
			if(N>4){
				cmat = matrix(NA, ncol = 4, nrow = mp)
				colnames(cmat) = c(paste("refit-", 1:2, sep = ""), paste("refit", (N-1):N, sep = ""))
				rownames(cmat) = rownames(object@roll$coefmat[[1]])
				cmat[,1:2] = cbind(object@roll$coefmat[[1]][,1], object@roll$coefmat[[2]][,1])
				cmat[,3:4] = cbind(object@roll$coefmat[[(N-1)]][,1], object@roll$coefmat[[N]][,1])
			} else{
				cmat = matrix(NA, ncol = N, nrow = mp)
				colnames(cmat) = paste("refit-", 1:N, sep = "")
				rownames(cmat) = rownames(object@roll$coefmat[[1]])
				for(i in 1:N) cmat[,i] = object@roll$coefmat[[i]][,1]
			}

			cat("\nCoefficients (first 2 , last 2):\n")
			print(round(cmat,5), digit = 5)
			cat("\n")
			invisible(object)
		})
#-------------------------------------------------------------------------
# multi-methods
setMethod("show",
		signature(object = "uGARCHmultispec"),
		function(object){
			cat(paste("\n*-----------------------------*", sep = ""))
			cat(paste("\n*     GARCH Multi-Spec        *", sep = ""))
			cat(paste("\n*-----------------------------*", sep = ""))
			N = length(object@spec)
			cat(paste("\nMultiple Specifications\t: ", N, sep=""))
			cat(paste("\nMulti-Spec Type\t\t\t: ", object@type, sep=""))
			cat("\n")
			invisible(object)
		})
		
setMethod("show",
		signature(object = "uGARCHmultifit"),
		function(object){
			cat(paste("\n*----------------------------*", sep = ""))
			cat(paste("\n*     GARCH Multi-Fit        *", sep = ""))
			cat(paste("\n*----------------------------*", sep = ""))
			cat(paste("\nNo. Assets :", length(object@fit), sep=""))
			asset.names = object@desc$asset.names
			
			if(object@desc$type == "equal"){
				vmodel = object@fit[[1]]@model$modeldesc$vmodel
				cat(paste("\nGARCH Multi-Spec Type : Equal",sep=""))
				cat(paste("\nGARCH Model Spec",sep=""))
				cat(paste("\n--------------------------",sep=""))
				cat(paste("\nModel : ", vmodel,sep=""))
				if(vmodel == "fGARCH" ){
					cat(paste(" Sub-Model : ", object@fit[[1]]@model$modeldesc$vsubmodel, "\n", sep = ""))
				}
				if(object@fit[[1]]@model$modelinc[15]>0){
					cat("\nExogenous Regressors in variance equation: ", object@fit[[1]]@model$modelinc[15], "\n")
				} else{
					cat("\nExogenous Regressors in variance equation: none\n")
				}
				cat("\nMean Equation :")
				cat("\nInclude Mean : ", object@fit[[1]]@model$modelinc[1])
				cat(paste("\nAR(FI)MA Model : (",object@fit[[1]]@model$modelinc[2],",",
								ifelse(object@fit[[1]]@model$modelinc[4]>0, 1, "d"),
								",",object@fit[[1]]@model$modelinc[3],")",sep=""))
				cat("\nGARCH-in-Mean : ", as.logical(object@fit[[1]]@model$modelinc[5]))
				if(object@fit[[1]]@model$modelinc[6]>0){
					cat("\nExogenous Regressors in mean equation: ", object@fit[[1]]@model$modelinc[6])
				} else{
					cat("\nExogenous Regressors in mean equation: none")
				}
				cat("\nConditional Distribution: ",object@fit[[1]]@model$modeldesc$distribution,"\n")
				cv = sapply(object@fit, FUN = function(x) x@fit$convergence)
				if(any(cv != 0)){
					ncv = which(cv != 0)
					nncv = length(ncv)
					cat("\nNo. of non converged fits: ", ncv,"\n")
					if(ncv>0) cat("\nNon converged fits: ", nncv,"\n")
					
				} else{
				cat(paste("\nGARCH Model Fit", sep = ""))
				cat(paste("\n--------------------------", sep = ""))
				cat("\nOptimal Parameters:\n")
				ll = t(likelihood(object))
				rownames(ll) = "Log-Lik"
				cf = coef(object)
				colnames(cf) = asset.names
				print(round(rbind(cf, ll), digits = 5))
				cat("\n")
			}
		} else{
			cat(paste("\nGARCH Model Fit", sep = ""))
			cat(paste("\n--------------------------", sep = ""))
			cat("\nOptimal Parameters:\n")
			print(coef(object), digits = 5)
		}
		invisible(object)
	})
	
setMethod("show",
		signature(object = "uGARCHmultifilter"),
		function(object){
			asset.names = object@desc$asset.names
			cat(paste("\n*-------------------------------*", sep = ""))
			cat(paste("\n*     GARCH Multi-Filter        *", sep = ""))
			cat(paste("\n*-------------------------------*", sep = ""))
			cat(paste("\nNo. Assets :", length(object@filter), sep=""))
			if(object@model$type == "equal"){
					cat(paste("\nGARCH Model Filter", sep = ""))
					cat(paste("\n--------------------------", sep = ""))
					cat("\nParameters:\n")
					cf = coef(object)
					colnames(cf) = asset.names
					print(round(cf, digits = 5))
			} else{
				cat(paste("\nGARCH Model Filter", sep = ""))
				cat(paste("\n--------------------------", sep = ""))
				cat("\nOptimal Parameters:\n")
				print(coef(object), digits = 5)
			}
			invisible(object)
		})

setMethod("show",
		signature(object = "uGARCHmultiforecast"),
		function(object){
			asset.names = object@desc$asset.names
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*       GARCH Multi-Forecast      *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\nNo. Assets :", length(object@forecast), sep=""))
			cat(paste("\n--------------------------\n",sep=""))
			fc = as.array(object)
			print(fc, digits = 5)
			invisible(object)
		})
#----------------------------------------------------------------------------------
# report method
#----------------------------------------------------------------------------------
report = function(object, ...)
{
	UseMethod("report")
}

setMethod("report", signature(object = "uGARCHroll"), .ugarchrollreport)
#----------------------------------------------------------------------------------
# univariate fit extractors
#----------------------------------------------------------------------------------
# coef methods
.ugarchfitcoef = function(object)
{
	object@fit$coef
}

setMethod("coef", signature(object = "uGARCHfit"), .ugarchfitcoef)

.ugarchfiltercoef = function(object)
{
	object@model$pars[object@model$pars[,2]==1, 1]
}

setMethod("coef", signature(object = "uGARCHfilter"), .ugarchfiltercoef)

# multi-fit and multi-filter coefficients:
.ugarchmultifitcoef = function(object)
{
	if(object@desc$type == "equal"){
		ans = sapply(object@fit, FUN = function(x) coef(x), simplify = TRUE)
	} else{
		ans = lapply(object@fit, FUN = function(x) coef(x))
	}
	return(ans)
}
setMethod("coef", signature(object = "uGARCHmultifit"), .ugarchmultifitcoef)

.ugarchmultifiltercoef = function(object)
{
	
	ans = sapply(object@filter, FUN = function(x) coef(x), simplify = TRUE)
	return(ans)
}

setMethod("coef", signature(object = "uGARCHmultifilter"), .ugarchmultifiltercoef)

#----------------------------------------------------------------------------------
# as.data.frame method for fitted object
.ugarchfitdf = function(x, row.names = NULL, optional = FALSE, ...)
{
	fit = x@fit
	xdata = x@model$modeldata$data
	T = x@model$modeldata$T
	ans = data.frame(data = xdata[1:T], fitted = fit$fitted.values, 
			residuals = fit$residuals, sigma = fit$sigma)
	rownames(ans) = as.character(x@model$modeldata$dates[1:T])
	return( ans )
}

setMethod("as.data.frame", signature(x = "uGARCHfit"), .ugarchfitdf)
#----------------------------------------------------------------------------------
# as.data.frame method for filter object
.ugarchfilterdf = function(x, row.names = NULL, optional = FALSE, ...)
{
	flt = x@filter
	T = x@model$modeldata$T
	ans = data.frame(data = flt$model$modeldata$data[1:T], 
			fitted = flt$model$modeldata$data[1:T] - flt$residuals, 
			residuals = flt$residuals, sigma = flt$sigma)
	rownames(ans) = as.character(x@model$modeldata$dates[1:T])
	return( ans )
}

setMethod("as.data.frame", signature(x = "uGARCHfilter"), .ugarchfilterdf)
#----------------------------------------------------------------------------------
# as.data.frame method for forecast object
.ugarchfordfall = function(x, which = "sigma", aligned = TRUE, prepad = TRUE, type = 0)
{
	if(aligned){
		ans = .ugarchfordf1(x = x, which = which, prepad = prepad)
	} else{
		ans = .ugarchfordf2(x = x, which = which, type = type)
	}
	return(ans)
}

.ugarchfordf1 = function(x, which = "sigma", prepad = TRUE)
{
	# return full frame with columns the rolls and rows the unique dates
	# padded with NA for forecasts beyond n.ahead and actual filtered values
	# for values before n.roll
	forc = x@forecast
	n.start = forc$n.start
	N = x@forecast$N - n.start
	n.ahead = forc$n.ahead
	n.roll = forc$n.roll
	tp = ifelse(which == "sigma", 1, 2)
	fdates = sort(unique(as.vector(sapply(forc$forecast, FUN=function(x) rownames(x)))))
	tmp = matrix(NA, ncol = n.roll+1, nrow = length(fdates))
	tmp[1:n.ahead, 1] = forc$forecast[[1]][,tp]
	if( n.roll > 0 ){
		for(i in 1:(n.roll)){
			if(which == "sigma"){
				if(prepad) {
					tmp[1:i, i+1] = x@model$modeldata$sigma[N:(N+i-1)]
				}
			} else {
				if(prepad) {
					tmp[1:i, i+1] = x@model$modeldata$data[N:(N+i-1)]
				} 
			}
			tmp[(i+1):(i+n.ahead), i+1] = forc$forecast[[i+1]][,tp]
		}
	}
	tmp = as.data.frame(tmp)
	rownames(tmp) = fdates
	colnames(tmp) = paste("roll-", seq(0, n.roll,by = 1), sep="")
	tmp
}

# retval: 0 is the standard, returns all values
# retval: 1 returns only those values which have in sample equivalent data (for testing purposes)
# retval: 2 returns only those values which are truly forecasts without in sample data
.ugarchfordf2 = function(x, which = "sigma", type = 0, ...)
{
	n.ahead = x@forecast$n.ahead
	# n.start == out.sample
	n.start = x@forecast$n.start
	n.roll =  x@forecast$n.roll
	tlength = n.ahead + n.roll
	mat = matrix(NA, ncol = n.roll + 1, nrow = n.ahead)
	mat = sapply(x@forecast$forecast, FUN = function(x) x[, which])
	if(is.vector(mat)) mat = matrix(mat, ncol = n.roll+1)
	colnames(mat) = paste("roll-", seq(0, n.roll, by = 1), sep = "")
	rownames(mat) = paste("t+", 1:n.ahead, sep = "")
	if(n.roll==0 | type==0) return(as.data.frame(mat))
	indices = apply(.embed(1:tlength, n.ahead, by = 1), 1, FUN = function(x) rev(x))
	exc = which(indices>n.start)
	if(type == 1){
		if(length(exc)>0) mat[exc] = NA
	} else{
		if(length(exc)>0) mat[-exc] = NA
	}
	return(as.data.frame(mat))
}

# pad applies to when rollframe = "all", whether to pre-pad forecast with actual values else NA's
# post forecast values are always NA's (this relates to the out.sample option and n.roll)
.ugarchfordf = function(x, row.names = NULL, optional = FALSE, rollframe = 0, 
		aligned = TRUE, which = "sigma", prepad = TRUE, type = 0)
{
	forc = x@forecast
	n.start = forc$n.start
	n.ahead = forc$n.ahead
	n.roll = forc$n.roll
	if(!any(type==(0:2)))
		stop("invalid argument for type in as.data.frame", call. = FALSE)
	#if(rollframe == "all" && n.roll<=0) rollframe = 0
	if(rollframe == "all") return(.ugarchfordfall(x, which, aligned, prepad, type))
	# allow type to be used here as well
	rollframe = as.integer(rollframe)
	if(rollframe>n.roll | rollframe<0) stop("rollframe out of bounds", call. = FALSE)
	indx = (rollframe+1):(rollframe+n.ahead)
	exc = which(indx>n.start)
	ans = forc$forecast[[rollframe+1]]
	if(type == 1){
		if(length(exc)>0) ans[exc,] = NA
	}
	if(type == 2){
		if(length(exc)>0) ans[-exc,] = NA
	}
	return(ans)
}

setMethod("as.data.frame", signature(x = "uGARCHforecast"), .ugarchfordf)

# as.array method for forecast object
.ugarchforar = function(x, ...)
{
	forc = x@forecast
	n.start = forc$n.start
	n.ahead = forc$n.ahead
	n.roll = forc$n.roll
	forar = array(NA, dim = c(n.ahead, 2, n.roll+1))
	for(i in 1:(n.roll+1)){
		forar[,,i] = as.matrix(forc$forecast[[i]])
	}
	return(forar)
}

setMethod("as.array", signature(x = "uGARCHforecast"), .ugarchforar)

# as.list method for forecast object
.ugarchforlist = function(x, ...)
{
	x@forecast$forecast
}

setMethod("as.list", signature(x = "uGARCHforecast"), .ugarchforlist)



# as.list method for forecast object
.ugarchmultiforlist = function(x, ...)
{
	ans = lapply(x@forecast, FUN = function(x) x@forecast$forecast)
	ans
}

setMethod("as.list", signature(x = "uGARCHmultiforecast"), .ugarchmultiforlist)

.ugarchmultiforarray = function(x, which = "sigma", ...)
{	
	n = length(x@forecast)
	asset.names = x@desc$asset.names
	ns = x@forecast[[1]]@forecast$n.roll
	nah = x@forecast[[1]]@forecast$n.ahead
	forar = array(NA, dim = c(nah, n, ns+1))
	rnames = paste("n.ahead-", 1:nah, sep = "")
	dnames = paste("n.roll-", 0:ns, sep = "")
	dimnames(forar) = list(rnames, asset.names, dnames)
	for(i in 1:(ns+1)) forar[,,i] = sapply(x@forecast, FUN = function(x) x@forecast$forecasts[[i]][,which])
	return(forar)
}

setMethod("as.array", signature(x = "uGARCHmultiforecast"), .ugarchmultiforarray)

#----------------------------------------------------------------------------------
# as.data.frame method for distribution object
.ugarchdistdf = function(x, row.names = NULL, optional = FALSE, which = "coef", window = 1, ...)
{
	n = x@dist$details$nwindows
	if(window > n) stop("window size greater than actual available", call. = FALSE)
	
	if(which == "rmse"){
		ans = as.data.frame(t(x@dist[[window]]$rmse))
		colnames(ans) = rownames(x@truecoef)
	}
	
	if(which == "stats"){	
		llh = x@dist[[window]]$likelist
		persist = x@dist[[window]]$persist
		uncvar = x@dist[[window]]$vlongrun
		uncmean = x@dist[[window]]$mlongrun
		maxret = x@dist[[window]]$simmaxdata[,1]
		minret = x@dist[[window]]$simmindata[,1]
		meanret = x@dist[[window]]$simmeandata[,1]
		kurtosis = x@dist[[window]]$simmomdata[,1]
		skewness = x@dist[[window]]$simmomdata[,2]
		maxsigma = x@dist[[window]]$simmaxdata[,3]
		minsigma = x@dist[[window]]$simmindata[,3]
		meansigma = x@dist[[window]]$simmeandata[,3]
		ans = data.frame(llh = llh, persist = persist, uncvar = uncvar, uncmean = uncmean,
				maxret = maxret, minret = minret, meanret = meanret, kurtosis = kurtosis,
				skewness  = skewness, maxsigma = maxsigma, minsigma = minsigma, meansigma = meansigma)
	}
	
	if(which == "coef"){
		cf = x@dist[[window]]$simcoef
		ans = data.frame(coef = cf)
		colnames(ans) = rownames(x@truecoef)
	}
	
	if(which == "coefse"){
		cfe = x@dist[[window]]$simcoefse
		ans = data.frame(coefse = cfe)
		colnames(ans) = rownames(x@truecoef)
	}
	
	ans
}

setMethod("as.data.frame", signature(x = "uGARCHdistribution"), .ugarchdistdf)

#----------------------------------------------------------------------------------

# as.data.frame method for simulation object
.ugarchsimdf = function(x, row.names = NULL, optional = FALSE, which = "sigma")
{
	# simframe: sigma, series, residuals
	#sigmaSim=sigmaSim, seriesSim=seriesSim, residSim=residSim
	sim = x@simulation
	T = x@model$modeldata$T
	dates = x@model$modeldata$dates[1:T]
	sigma = sim$sigmaSim
	m = dim(sigma)[2]
	N = dim(sigma)[1]
	fwdd = .generatefwd(dates, N = N, dformat = "%Y-%m-%d", periodicity = "days")
	if(which == "sigma"){
		ans = data.frame(sigmasim = sigma)
	}
	if(which == "series"){
		series = sim$seriesSim
		ans = data.frame(seriessim = series)
	}
	if(which == "residuals"){
		resids = sim$residSim
		ans = data.frame(residsim = resids)
	}
	rownames(ans) = as.character(fwdd)
	ans
}

setMethod("as.data.frame", signature(x = "uGARCHsim"), .ugarchsimdf)
#----------------------------------------------------------------------------------
# as.data.frame method for path simulation object
.ugarchpathdf = function(x, row.names = NULL, optional = FALSE, which = "sigma")
{
	# simframe: sigma, series, residuals
	#sigmaSim=sigmaSim, seriesSim=seriesSim, residSim=residSim
	sim = x@path
	sigma = sim$sigmaSim
	m = dim(sigma)[2]
	N = dim(sigma)[1]
	if(which == "sigma"){
		ans = data.frame(sigmasim = sigma)
	}
	if(which == "series"){
		series = sim$seriesSim
		ans = data.frame(seriessim = series)
	}
	if(which == "residuals"){
		resids = sim$residSim
		ans = data.frame(residsim = resids)
	}
	ans
}

setMethod("as.data.frame", signature(x = "uGARCHpath"), .ugarchpathdf)
#----------------------------------------------------------------------------------
# as.data.frame method for bootstrap object
.ugarchbootdf = function(x, row.names = NULL, optional = FALSE, which = "sigma", type = "raw", qtile = c(0.01, 0.099))
{
	n.ahead = x@model$n.ahead
	if(which == "sigma")
	{
		if(type == "raw"){
			sigma = x@fsigma
			ans = data.frame(bootsigma = sigma)
			colnames(ans) = paste("t+", 1:n.ahead, sep="")
		}
		if(type == "q"){
			if(all(is.numeric(qtile)) && (all(qtile<1.0) && all(qtile >0.0))){
				sigma = x@fsigma
				ans = apply(sigma, 2, FUN = function(x) quantile(x, qtile))
				ans = as.data.frame(ans)
				colnames(ans) = paste("t+", 1:n.ahead, sep="")
				rownames(ans) = paste("q", qtile, sep = "")
			} else{
				stop("\nfor type q, the qtile value must be numeric and between (>)0 and 1(<)\n", call.  = FALSE)
			}
		} 
		if(type == "summary"){
			sigma = x@fsigma
			ans = apply(sigma, 2, FUN = function(x) c(min(x), quantile(x, 0.25), mean(x), quantile(x, 0.75), max(x) ))
			ans = as.data.frame(ans)
			colnames(ans) = paste("t+", 1:n.ahead, sep="")
			rownames(ans) = c("min", "q0.25", "mean", "q0.75", "max")
		}
	}
	if(which == "series")
	{
		if(type == "raw"){
			series = x@fseries
			ans = data.frame(bootseries = series)
			colnames(ans) = paste("t+", 1:n.ahead, sep="")
		}
		if(type == "q"){
			if(all(is.numeric(qtile)) && (all(qtile<1.0) && all(qtile >0.0))){
						series = x@fseries
						ans = apply(series, 2, FUN = function(x) quantile(x, qtile))
						ans = as.data.frame(ans)
						colnames(ans) = paste("t+", 1:n.ahead, sep="")
						rownames(ans) = paste("q", qtile, sep = "")
					} else{
						stop("\nfor type q, the qtile value must be numeric and between (>)0 and 1(<)\n", call.  = FALSE)
					} 
		}
		if(type == "summary"){
			series = x@fseries
			ans = apply(series, 2, FUN = function(x) c(min(x), quantile(x, 0.25), mean(x), quantile(x, 0.75), max(x) ))
			ans = as.data.frame(ans)
			colnames(ans) = paste("t+", 1:n.ahead, sep="")
			rownames(ans) = c("min", "q.25", "mean", "q.75", "max")
		}
	}
	ans
}
setMethod("as.data.frame", signature(x = "uGARCHboot"), .ugarchbootdf)


#----------------------------------------------------------------------------------
# as.data.frame method for roll object
# valid which = density, fpm, coefs
.ugarchrolldf = function(x, row.names = NULL, optional = FALSE, which = "density", n.ahead = 1, refit = 1, aligned = FALSE,
		prepad = FALSE)
{
	n = x@roll$n.ahead
	if(n.ahead>n)
		stop("n.ahead chosen exceeds roll object specification", call. = FALSE)
	if(which == "sigma"){
		if(refit == "all"){
			nr = x@roll$n.refit
			fn = NULL
			dt = NULL
			for(i in 1:nr){
				fn = c(fn, sapply(x@forecast[[i]]@forecast$forecasts, FUN = function(y) y[n.ahead, 1]))
				dt = c(dt, sapply(x@forecast[[i]]@forecast$forecasts, FUN = function(y) rownames(y[n.ahead, ])))
			}
			ans = data.frame(sigma = fn)
			rownames(ans) = dt
		} else{
			ans = as.data.frame(x@forecast[[refit]], which = "sigma", rollframe = "all", aligned = aligned, prepad = prepad)
		}
	}
	if(which == "series"){
		if(refit == "all"){
			nr = x@roll$n.refit
			fn = NULL
			dt = NULL
			for(i in 1:nr){
				fn = c(fn, sapply(x@forecast[[i]]@forecast$forecasts, FUN = function(y) y[n.ahead, 2]))
				dt = c(dt, sapply(x@forecast[[i]]@forecast$forecasts, FUN = function(y) rownames(y[n.ahead, ])))
			}
			ans = data.frame(series = fn)
			rownames(ans) = dt
		} else{
			ans = as.data.frame(x@forecast[[refit]], which = "series", rollframe = "all", aligned = aligned, prepad = prepad)
		}
	}
	if(which == "coefs"){
		ans = as.data.frame(x@roll$coefs)
		rownames(ans) = paste("refit-", 1:dim(ans)[1], sep = "")
	}
	if(which == "density"){
		ans =  as.data.frame(x@roll$fdensity[[n.ahead]])
		rownames(ans) = paste("roll-", 1:dim(ans)[1], sep = "")
	}
	if(which == "coefmat"){
		ans = as.data.frame(x@roll$coefmat[[refit]])
	}
	if(which == "LLH"){
		ans = as.data.frame(x@roll$LLH)
		rownames(ans) = paste("refit-", 1:dim(ans)[1], sep = "")
		colnames(ans) = "LLH"
	}
	if(which == "VaR"){
		ans = as.data.frame(x@roll$VaR.out[[n.ahead]])
	}
	return(ans)
}
setMethod("as.data.frame", signature(x = "uGARCHroll"), .ugarchrolldf)

as.uGARCHforecast = function(object, ...)
{
	setMethod("as.uGARCHforecast")
}

.roll2forc = function(object, refit = 1)
{
	n = object@roll$n.refit
	if(refit>n)
		stop("refit chosen exceeds roll object specification", call. = FALSE)
	object@forecast[[refit]]
}

setMethod("as.uGARCHforecast", signature(object = "uGARCHroll"), .roll2forc)

#----------------------------------------------------------------------------------
# residuals method
.ugarchfitresids = function(object)
{
	object@fit$residuals
}

setMethod("residuals", signature(object = "uGARCHfit"), .ugarchfitresids)

.ugarchfilterresids = function(object)
{
	object@filter$residuals
}

setMethod("residuals", signature(object = "uGARCHfilter"), .ugarchfilterresids)

.ugarchmultifitresids = function(object)
{
	sapply(object@fit, FUN = function(x) residuals(x), simplify = TRUE)
}

setMethod("residuals", signature(object = "uGARCHmultifit"), .ugarchmultifitresids)

.ugarchmultifilterresids = function(object)
{
	sapply(object@filter, FUN = function(x) residuals(x), simplify = TRUE)
}

setMethod("residuals", signature(object = "uGARCHmultifilter"), .ugarchmultifilterresids)




#----------------------------------------------------------------------------------
# sigma method
sigma = function(object, ...)
{
	UseMethod("sigma")
}

.ugarchfitsigma = function(object)
{
	object@fit$sigma
}

setMethod("sigma", signature(object = "uGARCHfit"), .ugarchfitsigma)

.ugarchfiltersigma = function(object)
{
	object@filter$sigma
}

setMethod("sigma", signature(object = "uGARCHfilter"), .ugarchfiltersigma)

.ugarchmultifitsigma = function(object)
{
	sapply(object@fit, FUN = function(x) sigma(x), simplify = TRUE)
}

setMethod("sigma", signature(object = "uGARCHmultifit"), .ugarchmultifitsigma)


.ugarchmultifiltersigma = function(object)
{
	sapply(object@filter, FUN = function(x) sigma(x), simplify = TRUE)
}

setMethod("sigma", signature(object = "uGARCHmultifilter"), .ugarchmultifiltersigma)

#----------------------------------------------------------------------------------
# nyblom method
nyblom = function(object)
{
	UseMethod("nyblom")
}

setMethod("nyblom", signature(object = "uGARCHfit"), .nyblomTest)
#----------------------------------------------------------------------------------
# signbias method
signbias = function(object)
{
	UseMethod("signbias")
}

setMethod("signbias", signature(object = "uGARCHfit"), .signbiasTest)
setMethod("signbias", signature(object = "uGARCHfilter"), .signbiasTest)
#----------------------------------------------------------------------------------
# goodness of fit method
gof = function(object,groups)
{
	UseMethod("gof")
}

setMethod("gof", signature(object = "uGARCHfit", groups = "numeric"), .gofTest)
setMethod("gof", signature(object = "uGARCHfilter", groups = "numeric"), .gofTest)

#----------------------------------------------------------------------------------
# Info Criteria method
infocriteria = function(object)
{
	UseMethod("infocriteria")
}
.ugarchinfocriteria = function(object)
{
	# indicator object@fit$ipars[,4] denotes the estimated parameters
	if(is(object, "uGARCHfilter")){
		np = sum(object@filter$ipars[,2]) 
	} else{
		np = sum(object@fit$ipars[,4])
	}
	itest = .information.test(likelihood(object), nObs = length(fitted(object)), 
			nPars = np)
	itestm = matrix(0, ncol = 1, nrow = 4)
	itestm[1,1] = itest$AIC
	itestm[2,1] = itest$BIC
	itestm[3,1] = itest$SIC
	itestm[4,1] = itest$HQIC
	colnames(itestm) = ""
	rownames(itestm) = c("Akaike", "Bayes", "Shibata", "Hannan-Quinn")
	return(itestm)
}

setMethod("infocriteria", signature(object = "uGARCHfit"), .ugarchinfocriteria)
setMethod("infocriteria", signature(object = "uGARCHfilter"), .ugarchinfocriteria)

#----------------------------------------------------------------------------------
# Likelihood method
likelihood = function(object)
{
	UseMethod("likelihood")
}
.ugarchfitLikelihood = function(object)
{
	return(object@fit$LLH)
}

setMethod("likelihood", signature(object = "uGARCHfit"), .ugarchfitLikelihood)

.ugarchfilterLikelihood = function(object)
{
	return(object@filter$LLH)
}

setMethod("likelihood", signature(object = "uGARCHfilter"), .ugarchfilterLikelihood)


.ugarchmultifilterLikelihood = function(object)
{
	sapply(object@filter, FUN = function(x) likelihood(x), simplify = TRUE)
}

setMethod("likelihood", signature(object = "uGARCHmultifilter"), .ugarchmultifilterLikelihood)

.ugarchmultifitLikelihood = function(object)
{
	sapply(object@fit, FUN = function(x) likelihood(x), simplify = TRUE)
}

setMethod("likelihood", signature(object = "uGARCHmultifit"), .ugarchmultifitLikelihood)

#----------------------------------------------------------------------------------
# Fitted method
.ugarchfitted = function(object)
{
	return(object@fit$fitted.values)
}

setMethod("fitted", signature(object = "uGARCHfit"), .ugarchfitted)

.ugarchfiltered = function(object)
{
	T = object@model$modeldata$T
	return(object@model$modeldata$data[1:T] - object@filter$residuals)
}

setMethod("fitted", signature(object = "uGARCHfilter"), .ugarchfiltered)

.ugarchmultifiltered = function(object)
{
	sapply(object@filter, FUN = function(x) fitted(x), simplify = TRUE)
}

setMethod("fitted", signature(object = "uGARCHmultifilter"), .ugarchmultifiltered)


.ugarchmultifitted = function(object)
{
	sapply(object@fit, FUN = function(x) fitted(x), simplify = TRUE)
}

setMethod("fitted", signature(object = "uGARCHmultifit"), .ugarchmultifitted)
#----------------------------------------------------------------------------------

# newsimpact curve method (not multiplied by unconditional sigma)
newsimpact = function(object, z = seq(-0.3, 0.3, length.out = 100))
{
	UseMethod("newsimpact")
}

.newsimpact = function(object, z = seq(-0.3, 0.3, length.out = 100))
{
	vmodel = object@model$modeldesc$vmodel
	ans = switch(vmodel,
			sGARCH = .sgarchni(z, object),
			fGARCH = .fgarchni(z, object),
			gjrGARCH = .gjrgarchni(z, object),
			eGARCH = .egarchni(z, object),
			apARCH = .aparchni(z, object),
			iGARCH = .sgarchni(z, object))
	return(ans)
}

setMethod("newsimpact", signature(object = "uGARCHfit"), .newsimpact)
setMethod("newsimpact", signature(object = "uGARCHfilter"), .newsimpact)

# the underlying news impact methods
.sgarchni = function(z, object)
{
	if(is.null(z)){
		if(is(object, "uGARCHfit")) mz = min(object@fit$residuals) else mz = min(object@filter$residuals)
		if(abs(mz) <1){
			zz = seq(round(mz, 2) , round(abs(mz), 2), length.out = 101)
		} else{
			zz = seq(round(mz, 0) , round(abs(mz), 0), length.out = 101)
		}
	} else{
		zz = z
	}
	if(is(object, "uGARCHfit")) ipars = object@fit$ipars else ipars = object@filter$ipars
	modelinc = object@model$modelinc
	idx = object@model$pidx
	omega = ipars[idx["omega",1],1]
	alpha = beta = 0
	if(modelinc[8]>0) alpha = ipars[idx["alpha",1],1]
	if(modelinc[9]>0) beta  = ipars[idx["beta",1],1]
	lrvar = rep(as.numeric(uncvariance(object)), length(zz))
	ans = omega + beta*lrvar + alpha*zz^2
	yexpr = expression(sigma[t]^2)
	xexpr = expression(epsilon[t-1])
	return(list(zy = as.numeric(ans), zx = zz, yexpr = yexpr, xexpr = xexpr))
}

.fgarchni = function(z, object)
{
	if(is.null(z)){
		if(is(object, "uGARCHfit")) mz = min(object@fit$residuals) else mz = min(object@filter$residuals)
		if(abs(mz) <1){
			zz = seq(round(mz, 2) , round(abs(mz), 2), length.out = 101)
		} else{
			zz = seq(round(mz, 0) , round(abs(mz), 0), length.out = 101)
		}
	} else{
		zz = z
	}
	if(is(object, "uGARCHfit")) ipars = object@fit$ipars else ipars = object@filter$ipars
	modelinc = object@model$modelinc
	idx = object@model$pidx
	omega = ipars[idx["omega",1],1]
	alpha = beta = lambda = delta = eta1 = eta2 = 0
	alpha = ipars[idx["alpha",1],1]
	beta = ipars[idx["beta",1],1]
	lambda = ipars[idx["lambda",1],1]
	delta = ipars[idx["delta",1],1]
	eta1 = ipars[idx["eta1",1],1]
	eta2 = ipars[idx["eta2",1],1]
	fk = object@model$fmodel$fpars$fk
	kdelta = delta + fk*lambda
	lrvar = rep(as.numeric(uncvariance(object))^(1/2), length(zz))
	ans = omega + alpha[1]*(lrvar^lambda)*(abs(zz/lrvar - eta2[1]) - eta1[1]*(zz/lrvar - eta2[1]))^kdelta + beta[1]*(lrvar^lambda)
	yexpr = expression(sigma[t]^2)
	xexpr = expression(z[t-1])
	return(list(zy = ans^(2/lambda), zx = zz, yexpr = yexpr, xexpr = xexpr))
}

.egarchni = function(z, object)
{
	if(is.null(z)){
		if(is(object, "uGARCHfit")) mz = min(object@fit$residuals) else mz = min(object@filter$residuals)
		if(abs(mz) <1){
			zz = seq(round(mz, 2) , round(abs(mz), 2), length.out = 101)
		} else{
			zz = seq(round(mz, 0) , round(abs(mz), 0), length.out = 101)
		}
	} else{
		zz = z
	}
	if(is(object, "uGARCHfit")) ipars = object@fit$ipars else ipars = object@filter$ipars
	modelinc = object@model$modelinc
	idx = object@model$pidx
	omega = ipars[idx["omega",1],1]
	alpha = beta = gamma = skew = shape = ghlambda = 0
	if(modelinc[8]>0) alpha = ipars[idx["alpha",1],1]
	if(modelinc[9]>0) beta = ipars[idx["beta",1],1]
	if(modelinc[10]>0) gamma = ipars[idx["gamma",1],1]
	if(modelinc[16]>0) skew = ipars[idx["skew",1],1]
	if(modelinc[17]>0) shape = ipars[idx["shape",1],1]
	if(modelinc[18]>0) ghlambda = ipars[idx["ghlambda",1],1]
	k = egarchKappa(ghlambda, shape, skew, object@model$modeldesc$distribution)
	lrvar = rep(as.numeric(uncvariance(object)), length(zz))
	sqlr = sqrt(lrvar)
	ans=exp(omega + alpha[1]*zz/sqlr + gamma[1]*(abs(zz/sqlr)-k) + beta[1]*log(lrvar))
	yexpr = expression(sigma[t]^2)
	xexpr = expression(epsilon[t-1])
	return(list(zy = ans, zx = zz, yexpr = yexpr, xexpr = xexpr))
}

.gjrgarchni = function(z, object)
{
	if(is.null(z)){
		if(is(object, "uGARCHfit")) mz = min(object@fit$residuals) else mz = min(object@filter$residuals)
		if(abs(mz) <1){
			zz = seq(round(mz, 2) , round(abs(mz), 2), length.out = 101)
		} else{
			zz = seq(round(mz, 0) , round(abs(mz), 0), length.out = 101)
		}
	} else{
		zz = z
	}
	if(is(object, "uGARCHfit")) ipars = object@fit$ipars else ipars = object@filter$ipars
	modelinc = object@model$modelinc
	idx = object@model$pidx
	omega = ipars[idx["omega",1],1]
	alpha = beta = gamma = 0
	if(modelinc[8]>0) alpha = ipars[idx["alpha",1],1]
	if(modelinc[9]>0) beta = ipars[idx["beta",1],1]
	if(modelinc[10]>0) gamma = ipars[idx["gamma",1],1]
	lrvar = rep(as.numeric(uncvariance(object)), length(zz))
	ans = omega + alpha[1]*zz^2 + gamma[1]*(zz^2)*(zz<0) + beta[1]*lrvar
	yexpr = expression(sigma[t]^2)
	xexpr = expression(epsilon[t-1])
	return(list(zy = ans, zx = zz, yexpr = yexpr, xexpr = xexpr))
}

.aparchni = function(z, object)
{
	
	if(is.null(z)){
		if(is(object, "uGARCHfit")) mz = min(object@fit$residuals) else mz = min(object@filter$residuals)
		if(abs(mz) <1){
			zz = seq(round(mz, 2) , round(abs(mz), 2), length.out = 101)
		} else{
			zz = seq(round(mz, 0) , round(abs(mz), 0), length.out = 101)
		}
	} else{
		zz = z
	}
	if(is(object, "uGARCHfit")) ipars = object@fit$ipars else ipars = object@filter$ipars
	modelinc = object@model$modelinc
	idx = object@model$pidx
	omega = ipars[idx["omega",1],1]
	alpha = beta = gamma = 0
	delta = 2
	if(modelinc[8]>0) alpha = ipars[idx["alpha",1],1]
	if(modelinc[9]>0) beta = ipars[idx["beta",1],1]
	if(modelinc[10]>0) gamma = ipars[idx["gamma",1],1]
	if(modelinc[13]>0) delta = ipars[idx["delta",1],1]
	lrvar = rep(as.numeric(uncvariance(object))^(1/2), length(zz))
	ans = omega + alpha[1]*(abs(zz) - gamma[1]*(zz))^delta + beta[1]*(lrvar^delta)
	yexpr = expression(sigma[t]^2)
	xexpr = expression(epsilon[t-1])
	return(list(zy = ans^(2/delta), zx = zz, yexpr = yexpr, xexpr = xexpr))
}
#----------------------------------------------------------------------------------
# Half-Life Method for the various garch models
# ln(0.5)/log(persistence)
halflife = function(object, pars, distribution = "norm", model = "sGARCH", 
		submodel = "GARCH")
{
	UseMethod("halflife")
}

.halflife1<-function(object)
{
	ps = persistence(object)
	hlf = log(0.5)/log(ps)
	names(hlf) = "Half-Life"
	return(hlf)
}

.halflife2 = function(pars, distribution = "norm", model = "sGARCH", 
		submodel = "GARCH")
{
	ps = persistence(pars = pars, distribution = distribution, model = model, submodel = submodel)
	hlf = log(0.5)/log(ps)
	names(hlf) = "Half-Life"
	return(hlf)
}
setMethod("halflife",signature(object = "uGARCHfilter", pars = "missing", 
				distribution = "missing", model = "missing", submodel = "missing"), 
		definition = .halflife1)
		
setMethod("halflife",signature(object = "uGARCHfit", pars = "missing", 
				distribution = "missing", model = "missing", submodel = "missing"), 
		definition = .halflife1)

setMethod("halflife",signature(object = "uGARCHspec", pars = "missing", 
				distribution = "missing", model = "missing", submodel = "missing"), 
		definition = .halflife1)

setMethod("halflife",signature(object = "missing", pars = "numeric", 
				distribution = "character", model = "character", submodel = "ANY"), 
		definition = .halflife2)
#----------------------------------------------------------------------------------
# Persistence
persistence = function(object, pars, distribution = "norm", model = "sGARCH", 
		submodel = "GARCH")
{
	UseMethod("persistence")
}

# filter method
.filterpersistence = function(object)
{
	ans = object@filter$persistence
	names(ans) = "persistence"
	return(ans)
}
setMethod("persistence", signature(object = "uGARCHfilter", pars = "missing", 
				distribution = "missing", model = "missing", submodel = "missing"), 
		definition = .filterpersistence)


# fit method
.persistence1 = function(object)
{
		pars = object@fit$ipars
		idx = object@model$pidx
		distribution = object@model$modeldesc$distribution
		vmodel = object@model$modeldesc$vmodel
		vsubmodel = object@model$modeldesc$vsubmodel
	ans = switch(vmodel,
			sGARCH = .persistsgarch1(pars, idx, distribution),
			eGARCH = .persistegarch1(pars, idx, distribution),
			gjrGARCH = .persistgjrgarch1(pars, idx, distribution),
			apARCH = .persistaparch1(pars, idx, distribution),
			fGARCH = .persistfgarch1(pars, idx, distribution, vsubmodel),
			iGARCH = 1)
	names(ans) = "persistence"
	return(ans)
}

.persistence2 = function(object)
{
	if(is.null(object@model$fixed.pars))
		stop("\nuncvariance with spec required fixed.pars list\n", call. = FALSE)
	# no other checks for now.
	model = object@model
	pars = unlist(model$fixed.pars)
	parnames = names(pars)
	modelnames = .checkallfixed(object)
	if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
		cat("\npersistence-->error: parameters names do not match specification\n")
		cat("Expected Parameters are: ")
		cat(paste(modelnames))
		cat("\n")
		stop("Exiting", call. = FALSE)
	}
	# once more into the spec
	setfixed(object)<-as.list(pars)
	model = object@model
	idx = model$pidx
	modelinc = model$modelinc
	pars = object@model$pars
	distribution = object@model$modeldesc$distribution
	vmodel = object@model$modeldesc$vmodel
	vsubmodel = object@model$modeldesc$vsubmodel
	ans = switch(vmodel,
			sGARCH = .persistsgarch1(pars, idx, distribution),
			eGARCH = .persistegarch1(pars, idx, distribution),
			gjrGARCH = .persistgjrgarch1(pars, idx, distribution),
			apARCH = .persistaparch1(pars, idx, distribution),
			fGARCH = .persistfgarch1(pars, idx, distribution, vsubmodel),
			iGARCH = 1)
	names(ans) = "persistence"
	return(ans)
}


.persistence3 = function(pars, distribution = "norm", model = "sGARCH", 
		submodel = "GARCH")
{
	ans = switch(model,
			sGARCH = .persistsgarch2(pars, distribution),
			eGARCH = .persistegarch2(pars, distribution),
			gjrGARCH = .persistgjrgarch2(pars, distribution),
			apARCH = .persistaparch2(pars, distribution),
			fGARCH = .persistfgarch2(pars, distribution, submodel),
			iGARCH = 1)
	names(ans) = "persistence"
	return(ans)
}


setMethod("persistence",signature(object = "uGARCHfit", pars = "missing", 
				distribution = "missing", model = "missing", submodel = "missing"), 
		definition = .persistence1)
setMethod("persistence",signature(object = "uGARCHspec", pars = "missing", 
				distribution = "missing", model = "missing", submodel = "missing"), 
		definition = .persistence2)
setMethod("persistence",signature(object = "missing", pars = "numeric", 
				distribution = "character", model = "character", submodel = "ANY"), 
		definition = .persistence3)

.persistsgarch1 = function(pars, idx, distribution = "norm"){	
	ps = sum(pars[idx["alpha",1]:idx["alpha",2]]) + sum(pars[idx["beta",1]:idx["beta",2]])
	return(ps)
}

.persistsgarch2 = function(pars, idx, distribution = "norm"){	
	Names = names(pars)
	if(any(substr(Names, 1, 5)=="alpha")){
		i = which(substr(Names, 1, 5)=="alpha")
		alpha = pars[i]
	} else{
		alpha = 0
	}	
	if(any(substr(Names, 1, 4)=="beta")){
		i = which(substr(Names, 1, 4)=="beta")
		beta = pars[i]
	} else{
		beta = 0
	}
	ps = sum(alpha) + sum(beta)	
	return(ps)
}

.persistgjrgarch1 = function(pars, idx, distribution = "norm"){	
	alpha = pars[idx["alpha",1]:idx["alpha",2]]
	beta  = pars[idx["beta",1]:idx["beta",2]]
	gamma = pars[idx["gamma",1]:idx["gamma",2]]
	skew = pars[idx["skew",1]:idx["skew",2]]
	shape = pars[idx["shape",1]:idx["shape",2]]
	ghlambda = pars[idx["ghlambda",1]:idx["ghlambda",2]]

	ps = sum(beta)+ sum(alpha)+sum(apply(as.data.frame(gamma),1,FUN=function(x) 
						x*pneg(ghlambda, shape, skew, distribution)))
	return(ps)
}

.persistgjrgarch2 = function(pars, distribution = "norm"){
	Names = names(pars)
	if(any(substr(Names, 1, 5)=="alpha")){
		i = which(substr(Names, 1, 5)=="alpha")
		alpha = pars[i]
	} else{
		alpha = 0
	}
	
	if(any(substr(Names, 1, 5)=="gamma")){
		i = which(substr(Names, 1, 5)=="gamma")
		gamma = pars[i]
	} else{
		gamma = 0
	}
	
	if(any(substr(Names, 1, 4)=="beta")){
		i = which(substr(Names, 1, 4)=="beta")
		beta = pars[i]
	} else{
		beta = 0
	}
	
	
	if(any(substr(Names, 1, 8)=="ghlambda")){
		i = which(substr(Names, 1, 8)=="ghlambda")
		ghlambda = pars[i]
	} else{
		ghlambda = 0
	}
	
	if(any(substr(Names, 1, 4)=="skew")){
		i = which(substr(Names, 1, 4)=="skew")
		skew = pars[i]
	} else{
		skew = 0
	}
	
	if(any(substr(Names, 1, 5)=="shape")){
		i = which(substr(Names, 1, 5)=="shape")
		shape = pars[i]
	} else{
		shape = 0
	}
	
	ps = sum(beta)+ sum(alpha)+sum(apply(as.data.frame(gamma),1,FUN=function(x) 
						x*pneg(ghlambda, shape, skew, distribution)))
	return(ps)
}

.persistegarch1 = function(pars, idx, distribution = "norm"){	
	beta  = pars[idx["beta",1]:idx["beta",2]]
	ps = sum(beta)
	return(ps)
}

.persistegarch2 = function(pars, distribution = "norm"){
	Names=names(pars)
	if(any(substr(Names, 1, 4)=="beta")){
		i=which(substr(Names, 1, 4)=="beta")
		beta=pars[i]
	} else{
		beta=0
	}
	ps=sum(beta)
	return(ps)
}

.persistaparch1 = function(pars, idx, distribution = "norm"){	
	alpha = pars[idx["alpha",1]:idx["alpha",2]]
	beta  = pars[idx["beta",1]:idx["beta",2]]
	gamma = pars[idx["gamma",1]:idx["gamma",2]]
	delta = pars[idx["delta",1]:idx["delta",2]]
	skew = pars[idx["skew",1]:idx["skew",2]]
	shape = pars[idx["shape",1]:idx["shape",2]]
	ghlambda = pars[idx["ghlambda",1]:idx["ghlambda",2]]
	ps = sum(beta) + sum(apply(cbind(gamma,alpha), 1, FUN=function(x) 
						x[2]*aparchKappa(x[1], delta, ghlambda, shape, skew,distribution)))
	return(ps)
}

.persistaparch2 = function(pars, distribution = "norm"){
	Names = names(pars)
	if(any(substr(Names, 1, 5)=="alpha")){
		i = which(substr(Names, 1, 5)=="alpha")
		alpha = pars[i]
	} else{
		alpha = 0
	}
	
	if(any(substr(Names, 1, 5)=="gamma")){
		i = which(substr(Names, 1, 5)=="gamma")
		gamma = pars[i]
	} else{
		gamma = rep(0, length(alpha))
	}
	
	if(any(substr(Names, 1, 5)=="delta")){
		i = which(substr(Names, 1, 5)=="delta")
		delta = pars[i]
	} else{
		delta=  2
	}
	
	if(any(substr(Names, 1, 4)=="beta")){
		i = which(substr(Names, 1, 4)=="beta")
		beta = pars[i]
	} else{
		beta = 0
	}
	if(any(substr(Names, 1, 8)=="ghlambda")){
		i = which(substr(Names, 1, 8)=="ghlambda")
		ghlambda = pars[i]
	} else{
		ghlambda = 0
	}
	
	if(any(substr(Names, 1, 4)=="skew")){
		i = which(substr(Names, 1, 4)=="skew")
		skew = pars[i]
	} else{
		skew = 0
	}
	
	if(any(substr(Names, 1, 5)=="shape")){
		i = which(substr(Names, 1, 5)=="shape")
		shape = pars[i]
	} else{
		shape = 0
	}
	ps = sum(beta) + sum(apply(cbind(gamma,alpha), 1, FUN=function(x) 
						x[2]*aparchKappa(x[1], delta, ghlambda, shape, skew,distribution)))
	return(ps)
}

.persistfgarch1 = function(pars, idx, distribution = "norm", submodel){
	fm = .fgarchModel(submodel)
	alpha = pars[idx["alpha",1]:idx["alpha",2]]
	beta  = pars[idx["beta",1]:idx["beta",2]]
	eta1 = pars[idx["eta1",1]:idx["eta1",2]]
	eta2 = pars[idx["eta2",1]:idx["eta2",2]]
	lambda = pars[idx["lambda",1]:idx["lambda",2]]
	delta = pars[idx["delta",1]:idx["delta",2]]
	skew = pars[idx["skew",1]:idx["skew",2]]
	shape = pars[idx["shape",1]:idx["shape",2]]
	ghlambda = pars[idx["ghlambda",1]:idx["ghlambda",2]]
	fk = fm$parameters$fk
	ps = sum(beta) + sum(apply(cbind(alpha, eta1, eta2), 1, FUN=function(x) 
						x[1]*fgarchKappa(lambda, delta, x[2], x[3], fk, ghlambda, shape, skew, distribution)))
	return(ps)
}

.persistfgarch2 = function(pars, distribution = "norm", submodel){
	fm = .fgarchModel(submodel)
	Names = names(pars)
	if(any(substr(Names, 1, 5)=="alpha")){
		i = which(substr(Names, 1, 5)=="alpha")
		alpha = pars[i]
	} else{
		alpha = 0
	}
	
	if(any(substr(Names, 1, 6)=="lambda")){
		i = which(substr(Names, 1, 6)=="lambda")
		lambda = pars[i]
	} else{
		lambda = fm$parameters$lambda
	}
	
	if(any(substr(Names, 1, 5)=="delta")){
		i = which(substr(Names, 1, 5)=="delta")
		delta = pars[i]
	} else{
		delta = fm$parameters$delta
	}
	
	if(any(substr(Names, 1, 4)=="eta1")){
		i = which(substr(Names, 1, 4)=="eta1")
		eta1 = pars[i]
	} else{
		eta1 = rep(0, length(alpha))
	}
	
	if(any(substr(Names, 1, 4)=="eta2")){
		i = which(substr(Names, 1, 4)=="eta2")
		eta2 = pars[i]
	} else{
		eta2 = rep(0, length(alpha))
	}
	
	if(any(substr(Names, 1, 4)=="beta")){
		i = which(substr(Names, 1, 4)=="beta")
		beta = pars[i]
	} else{
		beta = 0
	}
	
	if(any(substr(Names, 1, 8)=="ghlambda")){
		i = which(substr(Names, 1, 8)=="ghlambda")
		ghlambda = pars[i]
	} else{
		ghlambda = 0
	}
	
	if(any(substr(Names, 1, 4)=="skew")){
		i = which(substr(Names, 1, 4)=="skew")
		skew = pars[i]
	} else{
		skew = 0
	}
	
	if(any(substr(Names, 1, 5)=="shape")){
		i = which(substr(Names, 1, 5)=="shape")
		shape = pars[i]
	} else{
		shape = 0
	}
	fk = fm$parameters$fk
	ps = sum(beta) + sum(apply(cbind(alpha, eta1, eta2), 1, FUN=function(x) 
						x[1]*fgarchKappa(lambda, delta, x[2], x[3], fk, ghlambda, shape, skew, distribution)))
	return(ps)
}

.persistanstgarch = function(pars, distribution = "norm", submodel){
	Names = names(pars)
		
	if(any(substr(Names, 1, 6)=="alpha1")){
		i = which(substr(Names, 1, 6)=="alpha1")
		alpha1 = pars[i]
	} else{
		alpha1 = 0
	}
	
	if(any(substr(Names, 1, 6)=="alpha2")){
		i = which(substr(Names, 1, 6)=="alpha2")
		alpha2 = pars[i]
	} else{
		alpha2 = 0
	}
	
		if(any(substr(Names, 1, 5)=="beta1")){
		i = which(substr(Names, 1, 5)=="beta1")
		beta1 = pars[i]
	} else{
		beta1 = 0
	}
	
	if(any(substr(Names, 1, 5)=="beta2")){
		i = which(substr(Names, 1, 5)=="beta2")
		beta2 = pars[i]
	} else{
		beta2 = 0
	}
	ps = 0.5 * (sum(alpha1) + sum(alpha2) + sum(beta1) + sum(beta2))
	
	return(ps)
}

#----------------------------------------------------------------------------------
# Unconditional Variance
uncvariance = function(object, pars, distribution = "norm", model = "sGARCH", 
		submodel = "GARCH", vexdata = NULL)
{
	UseMethod("uncvariance")
}

.unconditional1 = function(object)
{
	pars = object@fit$ipars[,1]
	idx = object@model$pidx
	distribution = object@model$modeldesc$distribution
	vmodel = object@model$modeldesc$vmodel
	vsubmodel = object@model$modeldesc$vsubmodel
	T = object@model$modeldata$T
	vexdata = object@model$modeldata$vexdata[1:T, ,drop = FALSE]
	ans = switch(vmodel,
			sGARCH = .uncsgarch1(pars, idx, distribution, vexdata),
			eGARCH = .uncegarch1(pars, idx, distribution, vexdata),
			gjrGARCH = .uncgjrgarch1(pars, idx, distribution, vexdata),
			apARCH = .uncaparch1(pars, idx, distribution, vexdata),
			fGARCH = .uncfgarch1(pars, idx, distribution, vsubmodel, vexdata),
			iGARCH = Inf)
	names(ans) = "unconditional"
	return(ans)
}

.unconditional2 = function(object)
{
	if(is.null(object@model$fixed.pars))
		stop("\nuncvariance with spec required fixed.pars list\n", call. = FALSE)
	# no other checks for now.
	model = object@model
	pars = unlist(model$fixed.pars)
	parnames = names(pars)
	modelnames = .checkallfixed(object)
	if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
		cat("\nuncmean-->error: parameters names do not match specification\n")
		cat("Expected Parameters are: ")
		cat(paste(modelnames))
		cat("\n")
		stop("Exiting", call. = FALSE)
	}
	# once more into the spec
	setfixed(object)<-as.list(pars)
	model = object@model
	idx = model$pidx
	modelinc = model$modelinc
	pars = object@model$pars[,1]
	distribution = object@model$modeldesc$distribution
	vmodel = object@model$modeldesc$vmodel
	vsubmodel = object@model$modeldesc$vsubmodel
	#T = object@model$modeldata$T
	# assume no out of sample since this cannot be provided in the spec
	vexdata = object@model$modeldata$vexdata
	ans=switch(vmodel,
			sGARCH = 	.uncsgarch1(pars, idx, distribution, vexdata),
			eGARCH = 	.uncegarch1(pars, idx, distribution, vexdata),
			gjrGARCH = 	.uncgjrgarch1(pars, idx, distribution, vexdata),
			apARCH = 	.uncaparch1(pars, idx, distribution, vexdata),
			fGARCH = 	.uncfgarch1(pars, idx, distribution, vsubmodel, vexdata),
			iGARCH = 	Inf)
	names(ans) = "unconditional"
	return(ans)
}

.unconditional3 = function(object)
{
	pars = object@filter$ipars[,1]
	distribution = object@model$modeldesc$distribution
	vmodel = object@model$modeldesc$vmodel
	vsubmodel = object@model$modeldesc$vsubmodel
	idx = object@model$pidx
	T = object@model$modeldata$T
	vexdata = object@model$modeldata$vexdata[1:T, ,drop = FALSE]
	ans=switch(vmodel,
			sGARCH = 	.uncsgarch1(pars, idx, distribution, vexdata),
			eGARCH = 	.uncegarch1(pars, idx, distribution, vexdata),
			gjrGARCH = 	.uncgjrgarch1(pars, idx, distribution, vexdata),
			apARCH = 	.uncaparch1(pars, idx, distribution, vexdata),
			fGARCH = 	.uncfgarch1(pars, idx, distribution, vsubmodel, vexdata),
			iGARCH = 	Inf)
	names(ans) = "unconditional"
	return(ans)
}

.unconditional4 = function(pars, distribution = "norm", model = "sGARCH", 
		submodel = "GARCH", vexdata = NULL)
{
	ans = switch(model,
			sGARCH = 	.uncsgarch2(pars, distribution, vexdata),
			eGARCH = 	.uncegarch2(pars, distribution, vexdata),
			gjrGARCH = 	.uncgjrgarch2(pars, distribution, vexdata),
			apARCH = 	.uncaparch2(pars, distribution, vexdata),
			fGARCH = 	.uncfgarch2(pars, distribution, submodel, vexdata),
			iGARCH = 	Inf)
	names(ans) = "unconditional"
	return(ans)
}


setMethod("uncvariance", signature(object = "uGARCHfit", pars = "missing", 
				distribution = "missing", model = "missing", submodel = "missing", vexdata = "missing"), 
		definition = .unconditional1)

setMethod("uncvariance", signature(object = "uGARCHspec", pars = "missing", 
				distribution = "missing", model = "missing", submodel = "missing", vexdata = "missing"), 
		definition = .unconditional2)

setMethod("uncvariance", signature(object = "uGARCHfilter", pars = "missing", 
				distribution = "missing", model = "missing", submodel = "missing", vexdata = "missing"), 
		definition = .unconditional3)

setMethod("uncvariance", signature(object = "missing", pars = "numeric", 
				distribution = "character", model = "character", submodel = "ANY", vexdata = "ANY"),
		definition = .unconditional4)


.uncsgarch1 = function(pars, idx, distribution = "norm", vexdata = NULL){
	if(!is.null(vexdata)){
		vxreg = pars[idx["vxreg",1]:idx["vxreg",2]]
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	omega = pars[idx["omega",1]:idx["omega",2]]
	ps = sum(pars[idx["alpha",1]:idx["alpha",2]]) + sum(pars[idx["beta",1]:idx["beta",2]])
	uvol = (umeanvex+omega)/(1-ps)
	return(uvol)
}

.uncsgarch2 = function(pars, distribution = "norm", vexdata = NULL){
	Names = names(pars)
	if(!is.null(vexdata)){
		if(any(substr(Names, 1, 5) == "vxreg")){
			i = which(substr(Names, 1, 5) == "vxreg")
			vxreg = pars[i]
		} else{
			vxreg = 0
		}
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	if(any(substr(Names, 1, 5) == "omega")){
		i = which(substr(Names, 1, 5) == "omega")
		omega = pars[i]
	} else{
		omega = 0
	}
	ps = persistence(pars = pars, distribution = distribution, model = "sGARCH")
	uvol = (umeanvex+omega)/(1-ps)
	return(uvol)
}
.uncgjrgarch1 = function(pars, idx, distribution="norm", vexdata = NULL){
	if(!is.null(vexdata)){
		vxreg = pars[idx["vxreg",1]:idx["vxreg",2]]
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	omega = pars[idx["omega",1]:idx["omega",2]]
	ps = .persistgjrgarch1(pars, idx, distribution)
	uvol = (umeanvex+omega)/(1-ps)
	return(uvol)
}

.uncgjrgarch2 = function(pars, distribution = "norm", vexdata = NULL){
	Names = names(pars)
	if(!is.null(vexdata)){
		if(any(substr(Names, 1, 5) == "vxreg")){
			i = which(substr(Names, 1, 5) == "vxreg")
			vxreg = pars[i]
		} else{
			vxreg = 0
		}
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	if(any(substr(Names, 1, 5) == "omega")){
		i = which(substr(Names, 1, 5) == "omega")
		omega = pars[i]
	} else{
		omega = 0
	}
	ps = persistence(pars = pars, distribution = distribution, model = "gjrGARCH")
	uvol = (umeanvex+omega)/(1-ps)
	return(uvol)
}

.uncegarch1 = function(pars, idx, distribution="norm", vexdata = NULL){
	if(!is.null(vexdata)){
		vxreg = pars[idx["vxreg",1]:idx["vxreg",2]]
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	omega = pars[idx["omega",1]:idx["omega",2]]
	beta = pars[idx["beta",1]:idx["beta",2]]
	#gamma = pars[idx["gamma",1]:idx["gamma",2]]
	#kappa = egarchKappa(pars[idx["ghlambda",1]], pars[idx["shape",1]], ipars[idx["skew",1]], distribution)
	uvol = exp( (umeanvex+omega)/(1-sum(beta)) )
}

.uncegarch2 = function(pars,distribution="norm", vexdata = NULL){
	Names = names(pars)
	if(!is.null(vexdata)){
		if(any(substr(Names, 1, 5) == "vxreg")){
			i = which(substr(Names, 1, 5) == "vxreg")
			vxreg = pars[i]
		} else{
			vxreg = 0
		}
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	if(any(substr(Names, 1, 5)=="omega")){
		i=which(substr(Names, 1, 5)=="omega")
		omega=pars[i]
	} else{
		omega=0
	}
	if(any(substr(Names, 1, 4)=="beta")){
		i=which(substr(Names, 1, 4)=="beta")
		beta=pars[i]
	} else{
		beta=0
	}
	
	#if(any(substr(Names, 1, 5)=="gamma")){
	#	i=which(substr(Names, 1, 5)=="gamma")
	#	gamma=pars[i]
	#} else{
	#	gamma=0
	#}
	
	#if(any(substr(Names, 1, 8)=="ghlambda")){
	#	i=which(substr(Names, 1, 8)=="ghlambda")
	#	ghlambda=pars[i]
	#} else{
	#	ghlambda=0
	#}
	
	#if(any(substr(Names, 1, 5)=="shape")){
	#	i=which(substr(Names, 1, 5)=="shape")
	#	shape=pars[i]
	#} else{
	#	shape=0
	#}
	
	#if(any(substr(Names, 1, 4)=="skew")){
	#	i=which(substr(Names, 1, 4)=="skew")
	#	skew=pars[i]
	#} else{
	#	skew=0
	#}
	#kappa = egarchKappa(ghlambda, shape, skew, distribution)
	uvol = exp( (umeanvex+omega )/(1-sum(beta)) )
}

.uncaparch1 = function(pars, idx, distribution="norm", vexdata = NULL){
	if(!is.null(vexdata)){
		vxreg = pars[idx["vxreg",1]:idx["vxreg",2]]
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	omega = pars[idx["omega",1]:idx["omega",2]]
	delta = pars[idx["delta",1]:idx["delta",2]]
	ps = .persistaparch1(pars, idx, distribution)
	uvol = ((omega+umeanvex)/(1-ps))^(2/delta)
	return(uvol)
}


.uncaparch2 = function(pars,distribution="norm", vexdata = NULL){
	Names = names(pars)
	if(!is.null(vexdata)){
		if(any(substr(Names, 1, 5) == "vxreg")){
			i = which(substr(Names, 1, 5) == "vxreg")
			vxreg = pars[i]
		} else{
			vxreg = 0
		}
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	if(any(substr(Names, 1, 5)=="omega")){
		i=which(substr(Names, 1, 5)=="omega")
		omega=pars[i]
	} else{
		omega=0
	}
	if(any(substr(Names, 1, 5)=="delta")){
		i=which(substr(Names, 1, 5)=="delta")
		delta=pars[i]
	} else{
		delta=2
	}
	ps=persistence(pars=pars,distribution=distribution,model="apARCH")
	uvol = ((omega+umeanvex)/(1-ps))^(2/delta)
	return(uvol)
}

.uncfgarch1 = function(pars, idx, distribution="norm", submodel, vexdata = NULL){
	if(!is.null(vexdata)){
		vxreg = pars[idx["vxreg",1]:idx["vxreg",2]]
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	omega = pars[idx["omega",1]:idx["omega",2]]
	lambda = pars[idx["lambda",1]:idx["lambda",2]]
	ps = .persistfgarch1(pars, idx, distribution , submodel)
	uvol = ((omega+umeanvex)/(1-ps))^(2/lambda)
	return(uvol)
}

.uncfgarch2 = function(pars, distribution="norm", submodel, vexdata = NULL){
	Names = names(pars)
	if(!is.null(vexdata)){
		if(any(substr(Names, 1, 5) == "vxreg")){
			i = which(substr(Names, 1, 5) == "vxreg")
			vxreg = pars[i]
		} else{
			vxreg = 0
		}
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	fpars = .fgarchModel(submodel)$parameters
	if(any(substr(Names, 1, 6)=="lambda")){
		i=which(substr(Names, 1, 6)=="lambda")
		lambda=pars[i]
	} else{
		lambda = fpars$lambda
	}
	
	if(any(substr(Names, 1, 5)=="omega")){
		i=which(substr(Names, 1, 5)=="omega")
		omega=pars[i]
	} else{
		omega=0
	}
	
	ps = persistence(pars = pars, distribution = distribution, model = "fGARCH", submodel = submodel)
	uvol = ((omega+umeanvex)/(1-ps))^(2/lambda)
	return(uvol)
}

# Unconditional Mean
uncmean = function(object, method = c("analytical", "simulation"), n.sim = 20000, rseed = NA)
{
	UseMethod("uncmean")
}

.unconditionalmean1 = function(object, method = c("analytical", "simulation"), n.sim = 20000, rseed = NA)
{
	method = method[1]
	if(method == "analytical"){
		modelinc = object@model$modelinc
		idx = object@model$pidx
		if(is(object, "uGARCHfilter")){
			h = object@filter$sigma
			pars = object@filter$ipars[,1]
		} else{
			h = object@fit$sigma
			pars = object@fit$ipars[,1]
		}
		N = length(h)
		
		if(modelinc[6]>0){
			mxreg = matrix( pars[idx["mxreg",1]:idx["mxreg",2]], ncol = modelinc[6] )
			if(modelinc[20]==0){
				mexdata = matrix(object@model$modeldata$mexdata[1:N, ], ncol = modelinc[6])
				meanmex = apply(mexdata, 2, "mean")
				umeanmex = sum(mxreg*meanmex)			
			} else{
				if(modelinc[20] == modelinc[6]){
					mexdata = matrix(object@model$modeldata$mexdata[1:N, ], ncol = modelinc[6])
					meanmex = apply(mexdata, 2, "mean")*(uncvariance(object)^(1/2))
					umeanmex = sum(mxreg*meanmex)	
				} else{
					mexdata = matrix(object@model$modeldata$mexdata[1:N, ], ncol = modelinc[6])
					meanmex1 = apply(mexdata[,1:(modelinc[6]-modelinc[20]),drop=FALSE], 2, "mean")
					meanmex2 = apply(mexdata[,(modelinc[6]-modelinc[20]+1):modelinc[6],drop=FALSE], 2, "mean")*(uncvariance(object)^(1/2))
					umeanmex = sum(mxreg[,1:(modelinc[6]-modelinc[20])]*meanmex1)+sum(mxreg[,(modelinc[6]-modelinc[20]+1):modelinc[6]]*meanmex2)
				}
			}
		} else{
			umeanmex = 0	
		}
		if(modelinc[5]>0){
			# this is obviously an approximation....
			if(modelinc[5] == 2){
				mh = uncvariance(object) 
			} else{
				mh = uncvariance(object)^(1/2)
			}
			umeangim = mh * pars[idx["archm",1]]
		} else{
			umeangim = 0
		}
		if(modelinc[1]>0) mu = pars[idx["mu",1]] else mu=0
		umean = (mu + umeangim + umeanmex)
		return(umean)
	 } else{
		 sim = ugarchsim(object, n.sim = n.sim, n.start = 1000, startMethod = "sample", rseed = rseed)
		 umean = mean(sim@simulation$seriesSim[,1])
		 return(umean)
	 }
}

.unconditionalmean2 = function(object, method = c("analytical", "simulation"), n.sim = 20000, rseed = NA)
{
	method = method[1]
	if(is.null(object@model$fixed.pars)) stop("uncmean with uGARCHspec requires fixed.pars list", call. = FALSE)
	if(method == "analytical"){
		model = object@model
		pars = unlist(model$fixed.pars)
		parnames = names(pars)
		modelnames = .checkallfixed(object)
		if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
			cat("\nuncmean-->error: parameters names do not match specification\n")
			cat("Expected Parameters are: ")
			cat(paste(modelnames))
			cat("\n")
			stop("Exiting", call. = FALSE)
		}
		# once more into the spec
		setfixed(object)<-as.list(pars)
		model = object@model
		idx = model$pidx
		modelinc = model$modelinc
		pars = object@model$pars[,1]
		
		
		if(modelinc[6]>0){
			mxreg = matrix( pars[idx["mxreg",1]:idx["mxreg",2]], ncol = modelinc[6] )
			if(modelinc[20]==0){
				mexdata = matrix(object@model$modeldata$mexdata, ncol = modelinc[6])
				meanmex = apply(mexdata, 2, "mean")
				umeanmex = sum(mxreg*meanmex)			
			} else{
				if(modelinc[20] == modelinc[6]){
					mexdata = matrix(object@model$modeldata$mexdata, ncol = modelinc[6])
					meanmex = apply(mexdata, 2, "mean")*(uncvariance(object)^(1/2))
					umeanmex = sum(mxreg*meanmex)	
				} else{
					mexdata = matrix(object@model$modeldata$mexdata, ncol = modelinc[6])
					meanmex1 = apply(mexdata[,1:(modelinc[6]-modelinc[20]),drop=FALSE], 2, "mean")
					meanmex2 = apply(mexdata[,(modelinc[6]-modelinc[20]+1):modelinc[6],drop=FALSE], 2, "mean")*(uncvariance(object)^(1/2))
					umeanmex = sum(mxreg[,1:(modelinc[6]-modelinc[20])]*meanmex1)+sum(mxreg[,(modelinc[6]-modelinc[20]+1):modelinc[6]]*meanmex2)
				}
			}
		} else{
			umeanmex = 0	
		}
		
		if(modelinc[5]>0){
			# this is obviously an approximation....
			if(modelinc[5] == 2){
				mh = uncvariance(object) 
			} else{
				mh = uncvariance(object)^(1/2)
			}
			umeangim = mh * pars[idx["archm",1]]
		} else{
			umeangim = 0
		}
		if(modelinc[1]>0) mu = pars[idx["mu",1]] else mu=0
		umean = (mu + umeangim + umeanmex)
		return(umean)
	}  else{
		sim = ugarchpath(object, n.sim = n.sim, n.start = 1000, rseed = rseed)
		umean = mean(sim@path$seriesSim[,1])
		return(umean)
	}
}
setMethod("uncmean", signature(object = "uGARCHfit"),    definition = .unconditionalmean1)
setMethod("uncmean", signature(object = "uGARCHfilter"), definition = .unconditionalmean1)
setMethod("uncmean", signature(object = "uGARCHspec"),   definition = .unconditionalmean2)

#----------------------------------------------------------------------------------
# The mult- methods
#----------------------------------------------------------------------------------
multispec = function( speclist )
{
	UseMethod("multispec")
}

setMethod("multispec", signature(speclist = "vector"),  definition = .multispecall)


multifit = function(multispec, data, out.sample = 0, solver = "solnp", solver.control = list(), fit.control = list(stationarity = 1, 
				fixed.se = 0, scale = 0), parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2), ...)
{
	UseMethod("multifit")
}

setMethod("multifit", signature(multispec = "uGARCHmultispec"),  definition = .multifitgarch)


multifilter = function(multifitORspec, data = NULL, out.sample = 0, n.old = NULL, 
		parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2), ...)
{
	UseMethod("multifilter")
}

setMethod("multifilter", signature(multifitORspec = "uGARCHmultifit"),  definition = .multifiltergarch1)
setMethod("multifilter", signature(multifitORspec = "uGARCHmultispec"),  definition = .multifiltergarch2)


multiforecast = function(multifitORspec, data = NULL, n.ahead = 1, n.roll = 0, out.sample = 0, 
		external.forecasts = list(mregfor = NULL, vregfor = NULL), 
		parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2), ...)
{
	UseMethod("multiforecast")
}

setMethod("multiforecast", signature(multifitORspec = "uGARCHmultifit"),  definition = .multiforecastgarch1)
setMethod("multiforecast", signature(multifitORspec = "uGARCHmultispec"),  definition = .multiforecastgarch2)

#----------------------------------------------------------------------------------
fpm = function( object, summary = TRUE, ...)
{
	UseMethod("fpm")
}


.fpm1 = function(object, summary = TRUE)
{
	n.ahead = object@forecast$n.ahead
	if(n.ahead == 1){
		n.roll = object@forecast$n.roll
		N = object@forecast$N
		ns = object@forecast$n.start
		if(n.roll == 0 | n.roll<4) stop("\nfpm-->error: Forecast Performance Measures require at least 5 out of sample data points (n.roll>3).")
		forecast = sapply(object@forecast$forecasts, FUN = function(x) x[1,2])
		# get only the forecasts for which out.of.sample data is available
		forecast = forecast[1:ns]
		#dates = sapply(object@forecast$forecasts, FUN = function(x) rownames(x[1,]))
		#dates = dates[1:ns]
		actual = object@model$modeldata$data[(N - ns + 1):(N - ns + n.roll)]
		DAC = apply(cbind(actual, forecast), 1, FUN = function(x) as.integer(sign(x[1]) == sign(x[2])))
		if(summary){
			ans = data.frame(MSE = mean((forecast - actual)^2), MAE = mean(abs(forecast - actual)), DAC = mean(DAC))
		} else{
			ans = data.frame(SE = (forecast - actual)^2, AE = abs(forecast - actual), DAC = DAC)
			rownames(ans) = object@model$modeldata$dates[(N - ns + 1):(N - ns + n.roll)]
		}
	} else{
		if(n.ahead<4) stop("\nfpm-->error: Forecast Performance Measures require at least 5 out of sample data points (n.ahead>4).")
		T = object@model$modeldata$T
		if(summary){
			n.roll = object@forecast$n.roll
			actual = object@model$modeldata$data
			dates = as.character( object@model$modeldata$dates )
			ans = matrix(NA, ncol = n.roll+1, nrow = 4)
			colnames(ans) = paste("roll-", seq(0, n.roll))
			rownames(ans) = c("MSE", "MAE", "DAC", "N")
			for(i in 1:(n.roll+1)){
				tmp = as.data.frame(object, which = "series", aligned = FALSE, rollframe = i-1, type = 1)
				tmp = tmp[,2, drop = FALSE]
				tmp = tmp[which(!is.na(tmp)), , drop = FALSE]
				dt = rownames(tmp)
				actd = actual[match(dt, dates)]
				if(length(actd)>0 && length(actd)>4){
					ans[1,i] = mean((tmp[,1] - actd)^2)
					ans[2,i] = mean(abs(tmp[,1] - actd))
					DAC = apply(cbind(tmp[,1], actd), 1, FUN = function(x) as.integer(sign(x[1]) == sign(x[2])))
					ans[3,i] = mean(DAC)
				}
				ans[4,i] = length(actd)
			}
		} else{
			n.roll = object@forecast$n.roll
			actual = object@model$modeldata$data
			dates = as.character( object@model$modeldata$dates )
			sol = vector(mode = "list", length = n.roll+1)
			names(sol) = paste("roll-", seq(0, n.roll))
			for(i in 1:(n.roll+1)){
				tmp = as.data.frame(object, which = "series", aligned = FALSE, rollframe = i-1, type = 1)
				tmp = tmp[,2, drop = FALSE]
				tmp = tmp[which(!is.na(tmp)), , drop = FALSE]
				ans = matrix(NA, nrow = dim(tmp)[1], ncol = 3)
				colnames(ans) = c("SE", "AE", "DAC")
				dt = rownames(tmp)
				actd = actual[match(dt, dates)]
				if(length(actd)>0 && length(actd)>4){
					ans[,1] = (tmp[,1] - actd)^2
					ans[,2] = abs(tmp[,1] - actd)
					ans[,3] = apply(cbind(tmp[,1], actd), 1, FUN = function(x) as.integer(sign(x[1]) == sign(x[2])))
					ans[3,i] = mean(ans[,3])
				}
				sol[[i]] = ans
			}
			ans = sol
		}
	}
	return( ans )
}


.fpm2 = function(object, summary = TRUE)
{
	if(summary){
		n.ahead = object@roll$n.ahead
		Data = object@model$modeldata$data
		Dates = object@model$modeldata$dates
		ans = matrix(NA, ncol = n.ahead, nrow = 4)
		rownames(ans) = c("MSE", "MAE", "DAC", "N")
		colnames(ans) = paste("n.ahead-", 1:n.ahead, sep = "")
		dtx = vector(mode = "character", length = n.ahead)
		for(i in 1:n.ahead){
			forecast = as.data.frame(object, which = "series", refit = "all", n.ahead = i)
			dt = rownames(forecast)
			actual = Data[match(dt, as.character(Dates))]
			DAC = apply(cbind(actual, forecast), 1, FUN = function(x) as.integer(sign(x[1]) == sign(x[2])))
			tmp = c(mean( (forecast - actual)^2 ), mean(abs(forecast - actual)), mean(DAC), length(forecast[,1]))
			ans[, i] = tmp
			dtx[i] = paste(dt[1], " - ", dt[length(dt)], sep = "")
		}
		ans = data.frame(ans)
		attr(ans, "dates") = dtx
	} else{
		n.ahead = object@roll$n.ahead
		Data = object@model$modeldata$data
		Dates = object@model$modeldata$dates
		ans = NULL
		#rownames(ans) = c("MSE", "MAE", "DAC", "N")
		#colnames(ans) = paste("n.ahead-", 1:n.ahead, sep = "")
		dtx = vector(mode = "character", length = n.ahead)
		for(i in 1:n.ahead){
			forecast = as.data.frame(object, which = "series", refit = "all", n.ahead = i)
			dt = rownames(forecast)
			actual = Data[match(dt, as.character(Dates))]
			DAC = apply(cbind(actual, forecast), 1, FUN = function(x) as.integer(sign(x[1]) == sign(x[2])))
			tmp = cbind( (forecast - actual)^2 , abs(forecast - actual), DAC )
			ans[[i]] = tmp
			rownames(ans[[i]]) = dt
			colnames(ans[[i]]) = c("SE", "AE", "DAC")
		}
		names(ans) = paste("n.ahead-", 1:n.ahead, sep = "")
	}
	return( ans )
}
setMethod("fpm", signature(object = "uGARCHforecast"),  definition = .fpm1)
setMethod("fpm", signature(object = "uGARCHroll"),  definition = .fpm2)

convergence = function(object){
	UseMethod("convergence")
}
.convergence = function(object){
	return( object@fit$convergence )
}
setMethod("convergence", signature(object = "uGARCHfit"),  definition = .convergence)

.vcov = function(object, robust = FALSE){
	if(robust){
		return( object@fit$robust.cvar )
	} else{
		return( object@fit$cvar)
	}
}

setMethod("vcov", signature(object = "uGARCHfit"),  definition = .vcov)
