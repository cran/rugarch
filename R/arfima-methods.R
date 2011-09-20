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


#----------------------------------------------------------------------------------
# univariate spec method
#----------------------------------------------------------------------------------
arfimaspec = function( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, arfima = FALSE, 
				external.regressors = NULL), 
		distribution.model = "norm", start.pars = list(), fixed.pars = list(), ...)
{
	UseMethod("arfimaspec")
}


.xarfimaspec = function( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, arfima = FALSE, 
				external.regressors = NULL), 
		distribution.model = "norm", start.pars = list(), fixed.pars = list(), ...)
{
	mmodel = mean.model
	dmodel = distribution.model
	# make the temporary subsitution so that it will be accepted by ugarchspec
	if(!is.null(fixed.pars$sigma)){
		fixed.pars$omega = fixed.pars$sigma
		fixed.pars$sigma = NULL
	}
	if(!is.null(start.pars$sigma)){
		start.pars$omega = start.pars$sigma
		start.pars$sigma = NULL
	}
	ans = ugarchspec(mean.model = list(armaOrder = mmodel$armaOrder, include.mean = mmodel$include.mean,
					arfima = mmodel$arfima, external.regressors = mmodel$external.regressors, archm = FALSE,
	archpow = 1), variance.model = list(garchOrder = c(0,0), model = "sGARCH"), distribution.model = dmodel,
	start.pars = start.pars, fixed.pars = fixed.pars)

	if(!is.null(ans@model$fixed.pars$omega)){
		ans@model$fixed.pars$sigma = ans@model$fixed.pars$omega
		ans@model$fixed.pars$omega = NULL
	}
	if(!is.null(ans@model$start.pars$omega)){
		ans@model$start.pars$sigma = ans@model$start.pars$omega
		ans@model$start.pars$omega = NULL
	}
	model = ans@model
	# change the omega name to sigma
	names(model$modelinc)[7] = "sigma"
	model$modeldesc$vmodel = "constant"
    idx = which(rownames(model$pars) == "omega")
	rownames(model$pars)[idx] = "sigma"
	rownames(model$pos.matrix)[7] = "sigma"
	rownames(model$pidx)[7] = "sigma"
	sol = new("ARFIMAspec", model = model)
	return(sol)
}

setMethod(f = "arfimaspec", definition = .xarfimaspec)


.getarfimaspec = function(object)
{
	spec = arfimaspec(mean.model = list(armaOrder = c(object@model$modelinc[2], object@model$modelinc[3]),
					include.mean = object@model$modelinc[1], 
					arfima = object@model$modelinc[4], external.regressors = object@model$modeldata$mexdata), 
			distribution.model = object@model$modeldesc$distribution, start.pars  = object@model$start.pars, 
			fixed.pars = object@model$fixed.pars)
	return(spec)
}

setMethod(f = "getspec", signature(object = "ARFIMAfit"), definition = .getarfimaspec)



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
	tmp = arfimaspec(mean.model = list(armaOrder = c(model$modelinc[2], model$modelinc[3]), 
					include.mean = model$modelinc[1], 
					arfima = model$modelinc[4], external.regressors = model$modeldata$mexdata), 
			distribution.model = model$modeldesc$distribution, start.pars  = model$start.pars, 
			fixed.pars = as.list(fixed.pars))
	return(tmp)
}
setReplaceMethod(f="setfixed", signature= c(object = "ARFIMAspec", value = "vector"), definition = .setfixed)

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
			warning( (paste("Unrecognized Parameter in Fixed Values: ", parnames[i], "...Ignored", sep = "")))
		} else{
			inc = c(inc, i)
		}
	}
	start.pars = pars[inc]
	names(start.pars) = tolower(names(pars[inc]))
	# set parameter values
	tmp = arfimaspec(mean.model = list(armaOrder = c(model$modelinc[2], model$modelinc[3]), 
					include.mean = model$modelinc[1], 
					arfima = model$modelinc[4], external.regressors = model$modeldata$mexdata), 
			distribution.model = model$modeldesc$distribution, fixed.pars  = model$fixed.pars, 
			start.pars = as.list(start.pars))
	return(tmp)
}

setReplaceMethod(f="setstart", signature= c(object = "ARFIMAspec", value = "vector"), definition = .setstart)

arfimafit = function(spec, data, out.sample = 0, solver = "solnp", solver.control = list(), 
		fit.control = list(fixed.se = 0, scale = 0), ...)
{
	UseMethod("arfimafit")
}

setMethod("arfimafit", signature(spec = "ARFIMAspec"), .arfimafit)

arfimafilter = function(spec, data, out.sample = 0, n.old = NULL, ...)
{
	UseMethod("arfimafilter")
}

setMethod("arfimafilter", signature(spec = "ARFIMAspec"), .arfimafilter)


arfimaforecast = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0, 
		external.forecasts = list(mregfor = NULL), ...)
{
	UseMethod("arfimaforecast")
}

setMethod("arfimaforecast", signature(fitORspec = "ARFIMAfit"), .arfimaforecast)

setMethod("arfimaforecast", signature(fitORspec = "ARFIMAspec"), .arfimaforecast2)

arfimasim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, startMethod = c("unconditional","sample"), 
		prereturns = NA, preresiduals = NA, rseed = NA, custom.dist = list(name = NA, distfit = NA, type = "z"), 
		mexsimdata = NULL, ...)
{
	UseMethod("arfimasim")
}


setMethod("arfimasim", signature(fit = "ARFIMAfit"), .arfimasim)

arfimapath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, prereturns = NA, preresiduals = NA, 
		rseed = NA,  custom.dist = list(name = NA, distfit = NA, type = "z"), mexsimdata = NULL, ...)
{
	UseMethod("arfimapath")
}


setMethod("arfimapath", signature(spec = "ARFIMAspec"), .arfimapath)


arfimaroll = function(spec,  data, n.ahead = 1, forecast.length = 500, 
		refit.every = 25, refit.window = c("recursive", "moving"), parallel = FALSE, 
		parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2), 
		solver = "solnp", fit.control = list(), solver.control = list() ,
		calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), ...)
{
	setMethod("arfimaroll")
}

setMethod("arfimaroll", signature(spec = "ARFIMAspec"),  definition = .arfimaroll)


arfimadistribution = function(fitORspec, n.sim = 2000, n.start = 1, 
		m.sim = 100, recursive = FALSE, recursive.length = 6000, recursive.window = 1000, 
		prereturns = NA, preresiduals = NA, rseed = NA, custom.dist = list(name = NA, distfit = NA, type = "z"), 
		mexsimdata = NULL, fit.control = list(), solver = "solnp", 
		solver.control = list(), parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2), ...)
{
	setMethod("arfimadistribution")
}
setMethod("arfimadistribution", signature(fitORspec = "ARFIMAfit"), .arfimadistribution)
setMethod("arfimadistribution", signature(fitORspec = "ARFIMAspec"), .arfimadistribution)

setMethod("show",
		signature(object = "ARFIMAspec"),
		function(object){
			model = object@model
			modelinc = model$modelinc
			cat(paste("\n*----------------------------------*", sep = ""))
			cat(paste("\n*       ARFIMA Model Spec          *", sep = ""))
			cat(paste("\n*----------------------------------*", sep = ""))
			cat("\nConditional Mean Dynamics")
			cat(paste("\n------------------------------------\n",sep=""))
			cat("Mean Model\t\t\t: ARFIMA(", modelinc[2],",", ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Include Mean\t\t:", as.logical(modelinc[1]),"\n")
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
# fit show
setMethod("show",
		signature(object = "ARFIMAfit"),
		function(object){
			model = object@model
			modelinc = model$modelinc
			cat(paste("\n*----------------------------------*", sep = ""))
			cat(paste("\n*          ARFIMA Model Fit        *", sep = ""))
			cat(paste("\n*----------------------------------*", sep = ""))
			cat("\nMean Model\t\t\t: ARFIMA(", modelinc[2],",",ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Distribution\t\t:", model$modeldesc$distribution,"\n")
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
				stdresid = object@fit$residuals/coef(object)["sigma"]
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
				cat("\nElapsed time :", object@fit$timer,"\n\n")
			} else{
				cat("\nConvergence Problem:")
				cat("\nSolver Message:", object@fit$message,"\n\n")
				
			}
			invisible(object)
		})

# filter show
setMethod("show",
		signature(object = "ARFIMAfilter"),
		function(object){
			model = object@model
			modelinc = object@model$modelinc
			cat(paste("\n*-------------------------------------*", sep = ""))
			cat(paste("\n*          ARFIMA Model Filter        *", sep = ""))
			cat(paste("\n*-------------------------------------*", sep = ""))
			cat("\nMean Model\t\t\t: ARFIMA(", modelinc[2],",",ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Distribution\t\t:", model$modeldesc$distribution,"\n")
			cat("\nFilter Parameters")
			cat(paste("\n---------------------------------------\n",sep=""))
			print(matrix(coef(object), ncol=1, dimnames = list(names(coef(object)), "")), digits = 5)
			cat("\nLogLikelihood :", object@filter$LLH, "\n")
			stdresid = object@filter$residuals/object@model$pars["sigma", 1]
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
			invisible(object)
		})
# sim show
setMethod("show",
		signature(object = "ARFIMAsim"),
		function(object){
			model = object@model
			modelinc = model$modelinc
			cat(paste("\n*-------------------------------------*", sep = ""))
			cat(paste("\n*       ARFIMA Model Simulation       *", sep = ""))
			cat(paste("\n*-------------------------------------*", sep = ""))
			sim = object@simulation
			dates = object@model$dates
			series = sim$seriesSim
			resids = sim$residSim
			m = dim(sigma)[2]
			N = dim(sigma)[1]
			cat(paste("\nHorizon: ",N))
			cat(paste("\nSimulations: ",m,"\n",sep=""))
			rx1 = apply(series, 2, FUN=function(x) mean(x))
			rx2 = apply(series, 2, FUN=function(x) range(x))
			T = object@model$modeldata$T
			xspec = .model2spec(as.list(object@model$pars[object@model$pars[,3]==1,1]), object@model, type = "ARFIMA")
			actual = c(0, mean(object@model$modeldata$data[1:T]), 
					min(object@model$modeldata$data[1:T]), max(object@model$modeldata$data[1:T]))
			uncond = c(0, uncmean(xspec), NA, NA)
			dd = data.frame(Seed = object@seed, Series.Mean = rx1, Series.Min = rx2[1,],
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
		signature(object = "ARFIMAforecast"),
		function(object){
			cat(paste("\n*----------------------------------*", sep = ""))
			cat(paste("\n*        ARFIMA Model Forecast     *", sep = ""))
			cat(paste("\n*----------------------------------*", sep = ""))
			n.ahead = object@forecast$n.ahead
			cat(paste("\n\nHorizon: ", n.ahead, sep = ""))
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
		signature(object = "ARFIMApath"),
		function(object){
			model = object@model
			modelinc = model$modelinc
			cat(paste("\n*--------.........------------------------*", sep = ""))
			cat(paste("\n*        ARFIMA Model Path Simulation     *", sep = ""))
			cat(paste("\n*-----------------------------------------*", sep = ""))
			sim = object@path
			series = sim$seriesSim
			resids = sim$residSim
			m = dim(series)[2]
			N = dim(series)[1]
			cat(paste("\n\nHorizon: ", N))
			cat(paste("\nSimulations: ", m, "\n", sep = ""))
			T = object@model$modeldata$T
			xspec = .model2spec(as.list(object@model$pars[object@model$pars[,3]==1,1]), object@model, type = "ARFIMA")
			uncond = c(0, uncmean(xspec), NA, NA)
			rx1 = apply(series, 2, FUN = function(x) mean(x))
			rx2 = apply(series, 2, FUN = function(x) range(x))
			dd = data.frame(Seed = object@seed, Series.Mean = rx1, Series.Min = rx2[1,], 
					Series.Max = rx2[2,])
			meansim = apply(dd, 2, FUN = function(x) mean(x))
			meansim[1] = 0			
			dd = rbind(dd, meansim, uncond)
			rownames(dd) = c(paste("sim", 1:m, sep = ""), "Mean(All)", "Unconditional")			
			print(dd, digits = 3)
			cat("\n\n")
		})

# distribution show
# distribution show
setMethod("show",
		signature(object = "ARFIMAdistribution"),
		function(object){
			cat(paste("\n*-------------------------------------*", sep = ""))
			cat(paste("\n*    ARFIMA Parameter Distribution    *", sep = ""))
			cat(paste("\n*-------------------------------------*", sep = ""))
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

#-------------------------------------------------------------------------
# multi-methods
setMethod("show",
		signature(object = "ARFIMAmultispec"),
		function(object){
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*        ARFIMA Multi-Spec        *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			N = length(object@spec)
			cat(paste("\n\nMultiple Specifications\t: ", N, sep=""))
			cat(paste("\nMulti-Spec Type\t\t\t: ", object@type, sep=""))
			cat("\n")
			invisible(object)
		})		

setMethod("show",
		signature(object = "ARFIMAmultifit"),
		function(object){
			cat(paste("\n*--------------------------------*", sep = ""))
			cat(paste("\n*       ARFIMA Multi-Fit         *", sep = ""))
			cat(paste("\n*--------------------------------*", sep = ""))
			cat(paste("\n\nNo. Assets :", length(object@fit), sep=""))
			asset.names = object@desc$asset.names
			if(object@desc$type == "equal"){
				cat(paste("\nMulti-Spec Type : Equal",sep=""))
				cat(paste("\n\nModel Spec",sep=""))
				cat(paste("\n-------------------------------\n",sep=""))
				cat("\nInclude Mean\t: ", as.logical( object@fit[[1]]@model$modelinc[1] ) )
				cat(paste("\nAR(FI)MA Model : (",object@fit[[1]]@model$modelinc[2],",",
								ifelse(object@fit[[1]]@model$modelinc[4]>0, 1, "d"),
								",",object@fit[[1]]@model$modelinc[3],")",sep=""))
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
					cat(paste("\nModel Fit", sep = ""))
					cat(paste("\n-------------------------------\n",sep=""))
					cat("\n")
					ll = t(likelihood(object))
					rownames(ll) = "Log-Lik"
					cf = coef(object)
					colnames(cf) = asset.names
					print(round(rbind(cf, ll), digits = 5))
					cat("\n")
				}
			} else{
				cat(paste("\nARFIMA Model Fit", sep = ""))
				cat(paste("\n-------------------------------\n",sep=""))
				cat(paste("\nModel Fit", sep = ""))
				cat(paste("\n-------------------------------\n",sep=""))
				cat("\n")
				print(coef(object), digits = 5)
			}
			invisible(object)
		})

setMethod("show",
		signature(object = "ARFIMAmultifilter"),
		function(object){
			cat(paste("\n*--------------------------------*", sep = ""))
			cat(paste("\n*       ARFIMA Multi-Filter         *", sep = ""))
			cat(paste("\n*--------------------------------*", sep = ""))
			cat(paste("\n\nNo. Assets :", length(object@filter), sep=""))
			asset.names = object@desc$asset.names
			if(object@desc$type == "equal"){
				cat(paste("\nMulti-Spec Type : Equal",sep=""))
				cat(paste("\n\nModel Spec",sep=""))
				cat(paste("\n-------------------------------\n",sep=""))
				cat("\nInclude Mean\t: ", as.logical( object@filter[[1]]@model$modelinc[1] ) )
				cat(paste("\nAR(FI)MA Model : (",object@filter[[1]]@model$modelinc[2],",",
								ifelse(object@filter[[1]]@model$modelinc[4]>0, 1, "d"),
								",",object@filter[[1]]@model$modelinc[3],")",sep=""))
				if(object@filter[[1]]@model$modelinc[6]>0){
					cat("\nExogenous Regressors in mean equation: ", object@filter[[1]]@model$modelinc[6])
				} else{
					cat("\nExogenous Regressors in mean equation: none")
				}
				cat("\nConditional Distribution: ",object@filter[[1]]@model$modeldesc$distribution,"\n")			
				cat(paste("\nModel Filter", sep = ""))
				cat(paste("\n-------------------------------\n",sep=""))
				cat("\n")
				ll = t(likelihood(object))
				rownames(ll) = "Log-Lik"
				cf = coef(object)
				colnames(cf) = asset.names
				print(round(rbind(cf, ll), digits = 5))
				cat("\n")
			} else{
				cat(paste("\nARFIMA Model Filter", sep = ""))
				cat(paste("\n-------------------------------\n",sep=""))
				cat(paste("\nModel Fit", sep = ""))
				cat(paste("\n-------------------------------\n",sep=""))
				cat("\n")
				print(coef(object), digits = 5)
			}
			invisible(object)
		})

setMethod("show",
		signature(object = "ARFIMAmultiforecast"),
		function(object){
			asset.names = object@desc$asset.names
			cat(paste("\n*----------------------------------*", sep = ""))
			cat(paste("\n*       ARFIMA Multi-Forecast      *", sep = ""))
			cat(paste("\n*----------------------------------*", sep = ""))
			cat(paste("\n\nNo. Assets :", length(object@forecast), sep=""))
			cat(paste("\n--------------------------\n",sep=""))
			fc = as.array(object)
			print(fc, digits = 5)
			invisible(object)
		})

#----------------------------------------------------------------------------------
# univariate fit extractors
#----------------------------------------------------------------------------------
# coef methods
.arfimafitcoef = function(object)
{
	object@fit$coef
}

setMethod("coef", signature(object = "ARFIMAfit"), .arfimafitcoef)

.arfimafiltercoef = function(object)
{
	object@model$pars[object@model$pars[,2]==1, 1]
}

setMethod("coef", signature(object = "ARFIMAfilter"), .arfimafiltercoef)

# multi-fit and multi-filter coefficients:
.arfimamultifitcoef = function(object)
{
	if(object@desc$type == "equal"){
		ans = sapply(object@fit, FUN = function(x) coef(x), simplify = TRUE)
	} else{
		ans = lapply(object@fit, FUN = function(x) coef(x))
	}
	return(ans)
}
setMethod("coef", signature(object = "ARFIMAmultifit"), .arfimamultifitcoef)

.arfimamultifiltercoef = function(object)
{
	
	ans = sapply(object@filter, FUN = function(x) coef(x), simplify = TRUE)
	return(ans)
}

setMethod("coef", signature(object = "ARFIMAmultifilter"), .arfimamultifiltercoef)


# as.data.frame method for fitted object
.arfimafitdf = function(x, row.names = NULL, optional = FALSE, ...)
{
	T = x@model$modeldata$T
	ans = data.frame(data = x@model$modeldata$data[1:T], fitted = x@fit$fitted.values, 
			residuals = x@fit$residuals)
	rownames(ans) = as.character(x@model$modeldata$dates[1:T])
	ans
}

setMethod("as.data.frame", signature(x = "ARFIMAfit"), .arfimafitdf)
#----------------------------------------------------------------------------------
# as.data.frame method for filter object
.arfimafilterdf = function(x, row.names = NULL, optional = FALSE, ...)
{
	xdata = x@model$modeldata$data	
	T = x@model$modeldata$T
	res = x@filter$residuals
	ans = data.frame(data = xdata[1:T], fitted = xdata[1:T] - res, 
			residuals = res)
	rownames(ans) = as.character(x@model$modeldata$dates[1:T])
	ans
}

setMethod("as.data.frame", signature(x = "ARFIMAfilter"), .arfimafilterdf)



#----------------------------------------------------------------------------------
# as.data.frame method for forecast object
.arfimafordfall = function(x, aligned = TRUE, prepad = TRUE, type = 0)
{
	if(aligned){
		ans = .arfimafordf1(x = x, prepad = prepad)
	} else{
		ans = .arfimafordf2(x = x, type = type)
	}
	return(ans)
}

.arfimafordf1 = function(x, prepad = TRUE)
{
	# return full frame with columns the rolls and rows the unique dates
	# padded with NA for forecasts beyond n.ahead and actual filtered values
	# for values before n.roll
	forc = x@forecast
	n.start = forc$n.start
	N = x@forecast$N - n.start
	n.ahead = forc$n.ahead
	n.roll = forc$n.roll
	fdates = sort(unique(as.vector(sapply(forc$forecast, FUN=function(x) rownames(x)))))
	tmp = matrix(NA, ncol = n.roll+1, nrow = length(fdates))
	tmp[1:n.ahead, 1] = forc$forecast[[1]][,1]
	for(i in 1:n.roll){
		if(prepad) {
			tmp[1:i, i+1] = x@model$modeldata$data[(N):(N+i-1)]
		} 
		tmp[(i+1):(i+n.ahead), i+1] = forc$forecast[[i+1]][,1]
	}
	tmp = as.data.frame(tmp)
	rownames(tmp) = fdates
	colnames(tmp) = paste("roll-", seq(0, n.roll,by = 1), sep="")
	tmp
}

# retval: 0 is the standard, returns all values
# retval: 1 returns only those values which have in sample equivalent data (for testing purposes)
# retval: 2 returns only those values which are truly forecasts without in sample data
.arfimafordf2 = function(x, type = 0, ...)
{
	n.ahead = x@forecast$n.ahead
	# n.start == out.sample
	n.start = x@forecast$n.start
	n.roll =  x@forecast$n.roll
	tlength = n.ahead + n.roll
	mat = matrix(NA, ncol = n.roll + 1, nrow = n.ahead)
	mat = sapply(x@forecast$forecast, FUN = function(x) x[,1])
	if(is.vector(mat)) mat = matrix(mat, ncol = n.roll+1)
	colnames(mat) = paste("roll-", seq(0, n.roll, by = 1), sep = "")
	rownames(mat) = paste("t+", 1:n.ahead, sep = "")
	if(n.roll == 0 | type == 0) return(as.data.frame(mat))
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
.arfimafordf = function(x, row.names = NULL, optional = FALSE, rollframe = 0, aligned = TRUE, 
		prepad = TRUE, type = 0)
{
	forc = x@forecast
	n.start = forc$n.start
	n.ahead = forc$n.ahead
	n.roll = forc$n.roll
	if(!any(type==(0:2)))
		stop("invalid argument for type in as.data.frame", call. = FALSE)
	#if(rollframe == "all" && n.roll<=0) rollframe = 0
	if(rollframe == "all") return(.arfimafordfall(x, aligned, prepad, type))
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

setMethod("as.data.frame", signature(x = "ARFIMAforecast"), .arfimafordf)

# as.array method for forecast object
.arfimaforar = function(x, ...)
{
	forc = x@forecast
	n.start = forc$n.start
	n.ahead = forc$n.ahead
	n.roll = forc$n.roll
	forar = array(NA, dim = c(n.ahead, 1, n.roll+1))
	for(i in 1:(n.roll+1)){
		forar[,,i] = as.matrix(forc$forecast[[i]])
	}
	return(forar)
}

setMethod("as.array", signature(x = "ARFIMAforecast"), .arfimaforar)

# as.list method for forecast object
.arfimaforlist = function(x, ...)
{
	x@forecast$forecast
}

setMethod("as.list", signature(x = "ARFIMAforecast"), .arfimaforlist)

# as.list method for forecast object
.arfimamultiforlist = function(x, ...)
{
	ans = lapply(x@forecast, FUN = function(x) x@forecast$forecast)
	names(ans) = x@desc$asset.names
	ans
}

setMethod("as.list", signature(x = "ARFIMAmultiforecast"), .arfimamultiforlist)

.arfimamultiforarray = function(x, ...)
{
	
	n = length(x@forecast)
	ns = x@forecast[[1]]@forecast$n.roll
	nah = x@forecast[[1]]@forecast$n.ahead
	forar = array(NA, dim = c(nah, n, ns+1))
	cnames = x@desc$asset.names
	rnames = paste("n.ahead-", 1:nah, sep = "")
	dnames = paste("n.roll-", 0:ns, sep = "")
	dimnames(forar) = list(rnames, cnames, dnames)
	for(i in 1:(ns+1)) forar[,,i] = as.matrix(sapply(x@forecast, FUN = function(x) unlist(x@forecast$forecasts[[i]])))
	return(forar)
}

setMethod("as.array", signature(x = "ARFIMAmultiforecast"), .arfimamultiforarray)


#----------------------------------------------------------------------------------
# as.data.frame method for distribution object
.arfimadistdf = function(x, row.names = NULL, optional = FALSE, which = "coef", window = 1, ...)
{
	n = x@dist$details$nwindows
	if(window > n) stop("window size greater than actual available", call. = FALSE)
	
	if(which == "rmse"){
		ans = as.data.frame(t(x@dist[[window]]$rmse))
		colnames(ans) = rownames(x@truecoef)
	}
	
	if(which == "stats"){	
		llh = x@dist[[window]]$likelist
		uncmean = x@dist[[window]]$mlongrun
		maxret = x@dist[[window]]$simmaxdata[,1]
		minret = x@dist[[window]]$simmindata[,1]
		meanret = x@dist[[window]]$simmeandata[,1]
		kurtosis = x@dist[[window]]$simmomdata[,1]
		skewness = x@dist[[window]]$simmomdata[,2]
		ans = data.frame(llh = llh, uncmean = uncmean, maxret = maxret, minret = minret, 
				meanret = meanret, kurtosis = kurtosis, skewness  = skewness)
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

setMethod("as.data.frame", signature(x = "ARFIMAdistribution"), .arfimadistdf)

#----------------------------------------------------------------------------------

# as.data.frame method for simulation object
.arfimasimdf = function(x, row.names = NULL, optional = FALSE, which = "series")
{
	# simframe: series, residuals
	#seriesSim=seriesSim, residSim=residSim
	sim = x@simulation
	T = x@model$modeldata$T
	dates = x@model$modeldata$dates[1:T]
	resids = sim$residSim
	m = dim(resids)[2]
	N = dim(resids)[1]
	fwdd = .generatefwd(dates, N = N, dformat = "%Y-%m-%d", periodicity = "days")
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

setMethod("as.data.frame", signature(x = "ARFIMAsim"), .arfimasimdf)
#----------------------------------------------------------------------------------
# as.data.frame method for path simulation object
.arfimapathdf = function(x, row.names = NULL, optional = FALSE, which = "series")
{
	# simframe: sigma, series, residuals
	#sigmaSim=sigmaSim, seriesSim=seriesSim, residSim=residSim
	sim = x@path
	resids = sim$residSim
	m = dim(resids)[2]
	N = dim(resids)[1]
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

setMethod("as.data.frame", signature(x = "ARFIMApath"), .arfimapathdf)
#----------------------------------------------------------------------------------
# as.data.frame method for bootstrap object
.arfimabootdf = function(x, row.names = NULL, optional = FALSE, type = "raw", qtile = c(0.01, 0.099))
{
	n.ahead = x@model$n.ahead
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
	ans
}
#setMethod("as.data.frame", signature(x = "ARFIMAboot"), .arfimabootdf)


#----------------------------------------------------------------------------------
# as.data.frame method for roll object
# valid which = density, fpm, coefs
.arfimarolldf = function(x, row.names = NULL, optional = FALSE, which = "density", n.ahead = 1, refit = 1, aligned = FALSE,
		prepad = FALSE)
{
	n = x@roll$n.ahead
	if(n.ahead>n)
		stop("n.ahead chosen exceeds roll object specification", call. = FALSE)
	if(which == "series"){
		if(refit == "all"){
			nr = x@roll$n.refit
			fn = NULL
			dt = NULL
			for(i in 1:nr){
				fn = c(fn, sapply(x@forecast[[i]]@forecast$forecasts, FUN = function(y) y[n.ahead,1]))
				dt = c(dt, sapply(x@forecast[[i]]@forecast$forecasts, FUN = function(y) rownames(y[n.ahead, 1, drop = FALSE])))
			}
			ans = data.frame(series = fn)
			rownames(ans) = dt
		} else{
			ans = as.data.frame(x@forecast[[refit]], rollframe = "all", aligned = aligned, prepad = prepad)
		}
	}
	else if(which == "coefs"){
		ans = as.data.frame(x@roll$coefs)
		rownames(ans) = paste("refit-", 1:dim(ans)[1], sep = "")
	}
	else if(which == "density"){
		ans =  as.data.frame(x@roll$fdensity[[n.ahead]])
		rownames(ans) = paste("roll-", 1:dim(ans)[1], sep = "")
	}
	else if(which == "coefmat"){
		ans = as.data.frame(x@roll$coefmat[[refit]])
	}
	else if(which == "LLH"){
		ans = as.data.frame(x@roll$LLH)
		rownames(ans) = paste("refit-", 1:dim(ans)[1], sep = "")
		colnames(ans) = "LLH"
	}
	else if(which == "VaR"){
		ans = as.data.frame(x@roll$VaR.out[[n.ahead]])
	} else{
		ans = NA
	}
	return(ans)
}
setMethod("as.data.frame", signature(x = "ARFIMAroll"), .arfimarolldf)

as.ARFIMAforecast = function(object, ...)
{
	setMethod("as.ARFIMAforecast")
}

.roll2arfimaforc = function(object, refit = 1)
{
	n = object@roll$n.refit
	if(refit>n)
		stop("refit chosen exceeds roll object specification", call. = FALSE)
	object@forecast[[refit]]
}

setMethod("as.ARFIMAforecast", signature(object = "ARFIMAroll"), .roll2arfimaforc)

#----------------------------------------------------------------------------------
# residuals method
.arfimafitresids = function(object)
{
	object@fit$residuals
}

setMethod("residuals", signature(object = "ARFIMAfit"), .arfimafitresids)

.arfimafilterresids = function(object)
{
	object@filter$residuals
}

setMethod("residuals", signature(object = "ARFIMAfilter"), .arfimafilterresids)

.arfimamultifitresids = function(object)
{
	sapply(object@fit, FUN = function(x) residuals(x), simplify = TRUE)
}

setMethod("residuals", signature(object = "ARFIMAmultifit"), .arfimamultifitresids)

.arfimamultifilterresids = function(object)
{
	sapply(object@filter, FUN = function(x) residuals(x), simplify = TRUE)
}

setMethod("residuals", signature(object = "ARFIMAmultifilter"), .arfimamultifilterresids)

#----------------------------------------------------------------------------------
# Likelihood method
.arfimafitLikelihood = function(object)
{
	return(object@fit$LLH)
}

setMethod("likelihood", signature(object = "ARFIMAfit"), .arfimafitLikelihood)

.arfimafilterLikelihood = function(object)
{
	return(object@filter$LLH)
}

setMethod("likelihood", signature(object = "ARFIMAfilter"), .arfimafilterLikelihood)


.arfimamultifilterLikelihood = function(object)
{
	sapply(object@filter, FUN = function(x) likelihood(x), simplify = TRUE)
}

setMethod("likelihood", signature(object = "ARFIMAmultifilter"), .arfimamultifilterLikelihood)

.arfimamultifitLikelihood = function(object)
{
	sapply(object@fit, FUN = function(x) likelihood(x), simplify = TRUE)
}

setMethod("likelihood", signature(object = "ARFIMAmultifit"), .arfimamultifitLikelihood)

#----------------------------------------------------------------------------------
# Fitted method
.arfimafitted = function(object)
{
	return(object@fit$fitted.values)
}

setMethod("fitted", signature(object = "ARFIMAfit"), .arfimafitted)

.arfimafiltered = function(object)
{
	res = object@filter$residuals
	T = object@model$modeldata$T
	return(object@model$modeldata$data[1:T] - res)
}

setMethod("fitted", signature(object = "ARFIMAfilter"), .arfimafiltered)

.arfimamultifiltered = function(object)
{
	sapply(object@filter, FUN = function(x) x@model$modeldata$data[1:x@model$modeldata$T] - residuals(x), simplify = TRUE)
}

setMethod("fitted", signature(object = "ARFIMAmultifilter"), .arfimamultifiltered)


.arfimamultifitted = function(object)
{
	sapply(object@fit, FUN = function(x) x@model$modeldata$data[1:x@model$modeldata$T] - residuals(x), simplify = TRUE)
}

setMethod("fitted", signature(object = "ARFIMAmultifit"), .arfimamultifitted)
#----------------------------------------------------------------------------------

# Info Criteria method
.arfimainfocriteria = function(object)
{
	if(is(object, "ARFIMAfit")){
		np = sum(object@fit$ipars[,4])
	} else{
		np = length(coef(object))
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

setMethod("infocriteria", signature(object = "ARFIMAfit"), .arfimainfocriteria)
setMethod("infocriteria", signature(object = "ARFIMAfilter"), .arfimainfocriteria)


#----------------------------------------------------------------------------------
# The mult- methods
#----------------------------------------------------------------------------------
.multispecarfima = function( speclist )
{
	# first create a spec which goes through validation process
	tp = 1
	ans = new("ARFIMAmultispec", 
			spec = speclist,
			type = "equal")
	# then check type
	n = length(speclist)
	for(i in 2:n){
		modelnames1 = rownames( speclist[[i]]@model$pars[speclist[[i]]@model$pars[,3]==1, ] )
		modelnames2 = rownames( speclist[[i-1]]@model$pars[speclist[[i-1]]@model$pars[,3]==1, ] )
		if(length(modelnames1) != length(modelnames2)){
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
	ans = new("ARFIMAmultispec",
			spec = speclist,
			type = type)
	return(ans)
}

setMethod("multifit", signature(multispec = "ARFIMAmultispec"),  definition = .multifitarfima)
setMethod("multifilter", signature(multifitORspec = "ARFIMAmultifit"),  definition = .multifilterarfima1)
setMethod("multifilter", signature(multifitORspec = "ARFIMAmultispec"),  definition = .multifilterarfima2)
setMethod("multiforecast", signature(multifitORspec = "ARFIMAmultifit"),  definition = .multiforecastarfima1)
setMethod("multiforecast", signature(multifitORspec = "ARFIMAmultispec"),  definition = .multiforecastarfima2)

#----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------
# univariate plot method / seperate for fit,sim and forecast
#----------------------------------------------------------------------------------
#setMethod(f = "plot", signature(x = "ARFIMAfit", y = "missing"), .plotarfimafit)

#setMethod(f = "plot", signature(x = "ARFIMAfilter", y = "missing"), .plotarfimafilter)

#setMethod(f = "plot", signature(x = "ARFIMAsim", y = "missing"), .plotarfimasim)

#setMethod(f = "plot", signature(x = "ARFIMAforecast", y = "missing"), .plotarfimaforecast)

#setMethod(f = "plot", signature(x = "ARFIMApath", y = "missing"), .plotarfimapath)

#setMethod(f = "plot", signature(x = "ARFIMAroll", y = "missing"), .plotarfimaroll)

#setMethod(f = "plot", signature(x = "ARFIMAdistribution", y = "missing"), .plotarfimadist)


# Unconditional Mean
.unconditionalmean11 = function(object, method = c("analytical", "simulation"), n.sim = 20000, rseed = NA)
{
	method = method[1]
	N = object@model$modeldata$T
	if(is(object, "ARFIMAfit")) pars = object@fit$ipars[,1] else pars = object@filter$ipars[,1]
	if(method == "analytical"){
		modelinc = object@model$modelinc
		idx = object@model$pidx
		if(modelinc[6]>0){
			mxreg = pars[idx["mxreg",1]:idx["mxreg",2]]
			mexdata = matrix(object@model$modeldata$mexdata[1:N, ], ncol = modelinc[6])
			meanmex = apply(mexdata, 2, "mean")
			umeanmex = sum(mxreg*meanmex)
		} else{
			umeanmex = 0
		}
		if(modelinc[1]>0) mu = pars[idx["mu",1]] else mu=0
		umean = (mu + umeanmex)
		return(umean)
	} else{
		if(is(object, "ARFIMAfit")){
			sim = arfimasim(object, n.sim = n.sim, n.start = 1000, startMethod = "sample", rseed = rseed)
			umean = mean(sim@simulation$seriesSim)
			return(umean)
		} else{
			stop("\nuncmean by simulation not available for ARFIMAfilter class object (used spec instead).")
		}
	}
}

.unconditionalmean21 = function(object, method = c("analytical", "simulation"), n.sim = 20000, rseed = NA)
{
	method = method[1]
	if(is.null(object@model$fixed.pars)) stop("uncmean with ARFIMAspec requires fixed.pars list", call. = FALSE)
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
			mxreg = pars[idx["mxreg",1]:idx["mxreg",2]]
			meanmex = apply(object@model$modeldata$mexdata, 2, "mean")
			umeanmex = sum(mxreg*meanmex)
		} else{
			umeanmex = 0
		}
		if(modelinc[1]>0) mu = pars[idx["mu",1]] else mu=0
		umean = (mu + umeanmex)
		return(umean)
	}  else{
		sim = arfimapath(object, n.sim = n.sim, n.start = 1000, rseed = rseed)
		umean = mean(sim@path$seriesSim)
		return(umean)
	}
}
setMethod("uncmean", signature(object = "ARFIMAfit"),    definition = .unconditionalmean11)
setMethod("uncmean", signature(object = "ARFIMAfilter"), definition = .unconditionalmean11)
setMethod("uncmean", signature(object = "ARFIMAspec"),   definition = .unconditionalmean21)

.fpm11 = function(object, summary = TRUE)
{
	n.ahead = object@forecast$n.ahead
	if(n.ahead == 1){
		n.roll = object@forecast$n.roll
		N = object@forecast$N
		ns = object@forecast$n.start
		if(n.roll == 0 | n.roll<4) stop("\nfpm-->error: Forecast Performance Measures require at least 5 out of sample data points (n.roll>3).")
		forecast = sapply(object@forecast$forecasts, FUN = function(x) x[1,])
		# get only the forecasts for which out.of.sample data is available
		forecast = forecast[1:ns]
		#dates = sapply(object@forecast$forecasts, FUN = function(x) rownames(x[1,]))
		#dates = dates[1:ns]
		actual = object@model$modeldata$data[(N - ns + 1):(N - ns + n.roll+1)]
		DAC = apply(cbind(actual, forecast), 1, FUN = function(x) as.integer(sign(x[1]) == sign(x[2])))
		if(summary){
			ans = data.frame(MSE = mean((forecast - actual)^2), MAE = mean(abs(forecast - actual)), DAC = mean(DAC))
		} else{
			ans = data.frame(SE = (forecast - actual)^2, AE = abs(forecast - actual), DAC = DAC)
			rownames(ans) = object@model$modeldata$dates[(N - ns + 1):(N - ns + n.roll+1)]
		}
	} else{
		if(object@model$n.start<5) stop("\nfpm-->error: Forecast Performance Measures require at least 5 out of sample data points (out.sample>4).")
		if(summary){
			n.roll = object@forecast$n.roll
			actual = object@model$modeldata$data
			dates = as.character( object@model$modeldata$dates )
			ans = matrix(NA, ncol = n.roll+1, nrow = 4)
			colnames(ans) = paste("roll-", seq(0, n.roll))
			rownames(ans) = c("MSE", "MAE", "DAC", "N")
			for(i in 1:(n.roll+1)){
				tmp = as.data.frame(object, aligned = FALSE, rollframe = i-1, type = 1)
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
				ans = matrix(NA, nrow = n.ahead, ncol = 3)
				colnames(ans) = c("SE", "AE", "DAC")
				tmp = as.data.frame(object, aligned = FALSE, rollframe = i-1, type = 1)
				tmp = tmp[,2, drop = FALSE]
				tmp = tmp[which(!is.na(tmp)), , drop = FALSE]
				dt = rownames(tmp)
				actd = actual[match(dt, dates)]
				if(length(actd)>0 && length(actd)>4){
					ans[,1] = (tmp[,1] - actd)^2
					ans[,2] = abs(tmp[,1] - actd)
					ans[,3] = apply(cbind(tmp[,1], actd), 1, FUN = function(x) as.integer(sign(x[1]) == sign(x[2])))
					ans[3,i] = mean(DAC)
				}
				sol[[i]] = ans
			}
			ans = sol
		}
	}
	return( ans )
}


.fpm21 = function(object, summary = TRUE)
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
			forecast = as.data.frame(object,  which = "series", refit = "all", n.ahead = i)
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

setMethod("fpm", signature(object = "ARFIMAforecast"),  definition = .fpm11)
setMethod("fpm", signature(object = "ARFIMAroll"),  definition = .fpm21)


# forecast performance measures
.arfimarollreport = function(object, type = "VaR", n.ahead = 1, VaR.alpha = 0.01, conf.level = 0.95)
{
	switch(type,
			VaR = .rollVaRreport1(object, n.ahead, VaR.alpha, conf.level),
			fpm = .rollfpmreport1(object))
	invisible(object)
}

.rollfpmreport1 = function(object)
{
	cat(paste("\nARFIMA Roll Mean Forecast Performance Measures", sep = ""))
	cat(paste("\n---------------------------------------------", sep = ""))
	cat(paste("\nno.refits : ", object@roll$n.refit, sep = ""))
	cat(paste("\nn.ahead   : ", object@roll$n.ahead, sep = ""))
	cat(paste("\nn.rolls   : ", object@roll$n.refit*object@roll$refit.every, sep = ""))
	cat("\n\n")
	tmp = fpm(object)
	print(signif(tmp, 4))
	cat("\n\n")
}

.rollVaRreport1 = function(object, n.ahead = 1, VaR.alpha = 0.01, conf.level = 0.95)
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

setMethod("report", signature(object = "ARFIMAroll"), .arfimarollreport)