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


rugarch.test1a = function(parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2)){
	tic = Sys.time()
	# simulated parameter distribution
	spec = arfimaspec( mean.model = list(armaOrder = c(2,2), include.mean = TRUE, arfima = FALSE), 
			distribution.model = "norm", fixed.pars = list(ar1=0.6, ar2=0.21, ma1=-0.7, ma2=0.3, mu = 0.02,
					sigma = 0.02))
	dist = arfimadistribution(spec, n.sim = 2000, n.start = 100, 
			m.sim = 100, recursive = TRUE, recursive.length = 10000, recursive.window = 1000,
			parallel = parallel, parallel.control = parallel.control)
	
	save(dist, file = "test1a.rda")
	options(width=150)
	zz <- file("test1a.txt", open="wt")
	sink(zz)
	# slots:
	slotNames(dist)
	# methods:
	# summary
	show(dist)
	# as.data.frame(...., window, which=c("rmse", "stats", "coef", "coefse"))
	# default
	as.data.frame(dist)
	as.data.frame(dist, window = 1, which = "rmse")
	as.data.frame(dist, window = 1, which = "stats")
	as.data.frame(dist, window = 1, which = "coef")
	as.data.frame(dist, window = 1, which = "coefse")
	as.data.frame(dist, window = 8, which = "rmse")
	as.data.frame(dist, window = 8, which = "stats")
	as.data.frame(dist, window = 8, which = "coef")
	as.data.frame(dist, window = 8, which = "coefse")
	sink(type="message")
	sink()
	close(zz)
	
	# create some plots
	nwindows = dist@dist$details$nwindows
	# 2000/3000/4000/5000/6000/7000/8000/9000/10000
	
	# expected reduction factor in RMSE for sqrt(N) consistency
	expexcted.rmsegr = sqrt(2000/seq(3000,10000,by=1000))
	
	# actual RMSE reduction
	actual.rmsegr = matrix(NA, ncol = 8, nrow = 6)
	rownames(actual.rmsegr) = c("mu", "ar1", "ar2", "ma2", "ma2", "sigma")
	# start at 2000 (window 1)
	rmse.start = as.data.frame(dist, window = 1, which = "rmse")
	for(i in 2:nwindows) actual.rmsegr[,i-1] = as.numeric(as.data.frame(dist, window = i, which = "rmse")/rmse.start)
	postscript("test1a.eps", bg = "white", width = 800, height = 800)
	par(mfrow = c(2,3))
	for(i in 1:6){
		plot(seq(3000,10000,by=1000), actual.rmsegr[i,], type = "l", lty = 2, ylab = "RMSE Reduction", xlab = "N (sim)", 
				main = rownames(actual.rmsegr)[i])
		lines(seq(3000,10000,by=1000), expexcted.rmsegr, col = 2)
		legend("topright", legend = c("Actual", "Expected"), col = 1:2, bty = "m", lty = c(2,1))
	}
	dev.off()
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}


rugarch.test1b = function(parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2)){
	# fit/filter
	tic = Sys.time()
	data(sp500ret)	
	fit = vector(mode = "list", length = 9)
	dist = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "jsu")
	for(i in 1:9){
		spec = arfimaspec( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, arfima = FALSE), 
				distribution.model = dist[i])
		fit[[i]] = arfimafit(spec = spec, data = sp500ret, solver = "solnp", fit.control = list(scale = 1))
	}
	cfmatrix = matrix(NA, nrow = 9, ncol = 7)
	colnames(cfmatrix) = c("mu", "ar1", "ma1", "sigma", "skew", "shape", "ghlambda")
	rownames(cfmatrix) = dist
	
	for(i in 1:9){
		cf = coef(fit[[i]])
		cfmatrix[i, match(names(cf), colnames(cfmatrix))] =  cf
	}
	sk = ku = rep(0, 9)
	for(i in 1:9){
		cf = coef(fit[[i]])
		if(fit[[i]]@model$modelinc[16]>0) sk[i] = dskewness(distribution = dist[i], 
					skew = cf["skew"], shape = cf["shape"], lambda = cf["ghlambda"])		
		if(fit[[i]]@model$modelinc[17]>0) ku[i] = dkurtosis(distribution = dist[i], 
					skew = cf["skew"], shape = cf["shape"], lambda = cf["ghlambda"])
	}
	hq = sapply(fit, FUN = function(x) infocriteria(x)[4])
	cfmatrix = cbind(cfmatrix, sk, ku, hq)
	colnames(cfmatrix) = c(colnames(cfmatrix[,1:7]), "skewness", "ex.kurtosis","HQIC")
	
	
	# filter the data to check results
	filt = vector(mode = "list", length = 9)
	for(i in 1:9){
		spec = arfimaspec( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, arfima = FALSE), 
				distribution.model = dist[i])
		setfixed(spec) = as.list(coef(fit[[i]]))
		filt[[i]] = arfimafilter(spec = spec, data = sp500ret)
	}
	
	options(width = 120)
	zz <- file("test1b.txt", open="wt")
	sink(zz)
	print(cfmatrix, digits = 4)
	cat("\nARFIMAfit and ARFIMAfilter residuals check:\n")
	print(head(sapply(filt, FUN = function(x) residuals(x))) == head(sapply(fit, FUN = function(x) residuals(x))))
	cat("\nas.data.frame method:\n")
	print(cbind(head(as.data.frame(filt[[1]])), head(as.data.frame(fit[[1]]))))
	cat("\ncoef method:\n")
	print(cbind(coef(filt[[1]]), coef(fit[[1]])))
	cat("\nfitted method:\n")
	print(cbind(head(fitted(filt[[1]])), head(fitted(fit[[1]]))))
	cat("\ninfocriteria method:\n")
	# For filter, it assumes estimation of parameters else does not make sense!
	print(cbind(infocriteria(filt[[1]]), infocriteria(fit[[1]])))
	cat("\nlikelihood method:\n")
	print(cbind(likelihood(filt[[1]]), likelihood(fit[[1]])))
	cat("\nresiduals method:\n")
	# Note that we the package will always return the full length residuals and fitted
	# object irrespective of the lags (i.e. since this is an ARMA(1,1) i.e. max lag = 1,
	# the first row is zero and should be discarded.
	print(cbind(head(residuals(filt[[1]])), head(residuals(fit[[1]]))))
	cat("\nuncmean method:\n")
	print(cbind(uncmean(filt[[1]]), uncmean(fit[[1]])))
	cat("\nuncmean method (by simulation):\n")
	# For spec and fit
	spec = arfimaspec( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, arfima = FALSE), distribution.model = dist[1])
	setfixed(spec) = as.list(coef(fit[[1]]))
	print(cbind(uncmean(spec, method = "simulation", n.sim = 100000, rseed = 100), 
					uncmean(fit[[1]], method = "simulation", n.sim = 100000, rseed = 100)))
	cat("\nsummary method:\n")
	print(show(filt[[1]]))
	print(show(fit[[1]]))
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test1c = function(parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2)){
	# unconditional forecasting
	tic = Sys.time()
	
	data(sp500ret)	
	fit = vector(mode = "list", length = 9)
	dist = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "jsu")
	for(i in 1:9){
		spec = arfimaspec( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, arfima = FALSE), 
				distribution.model = dist[i])
		fit[[i]] = arfimafit(spec = spec, data = sp500ret, solver = "solnp", fit.control = list(scale = 1))
	}
	cfmatrix = matrix(NA, nrow = 9, ncol = 7)
	colnames(cfmatrix) = c("mu", "ar1", "ma1", "sigma", "skew", "shape", "ghlambda")
	rownames(cfmatrix) = dist
	
	for(i in 1:9){
		cf = coef(fit[[i]])
		cfmatrix[i, match(names(cf), colnames(cfmatrix))] =  cf
	}
	
	umean = rep(0, 9)
	for(i in 1:9){
		umean[i] = uncmean(fit[[i]])
	}
	
	forc = vector(mode = "list", length = 9)
	for(i in 1:9){
		forc[[i]] = arfimaforecast(fit[[i]], n.ahead = 100)
	}
	
	lmean40 = sapply(forc, FUN = function(x) as.numeric(as.data.frame(x)[40,1]))
	cfmatrix1 = cbind(cfmatrix, umean, lmean40)
	colnames(cfmatrix1) = c(colnames(cfmatrix1[,1:7]), "uncmean", "forecast40")
	
	# forecast with spec to check results
	forc2 = vector(mode = "list", length = 9)
	for(i in 1:9){
		spec = arfimaspec( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, arfima = FALSE), 
				distribution.model = dist[i])
		setfixed(spec) = as.list(coef(fit[[i]]))
		forc2[[i]] = arfimaforecast(spec, data = sp500ret, n.ahead = 100)
	}
	lmean240 = sapply(forc2, FUN = function(x) as.numeric(as.data.frame(x)[40,1]))
	cfmatrix2 = cbind(cfmatrix, umean, lmean240)
	colnames(cfmatrix2) = c(colnames(cfmatrix2[,1:7]), "uncmean", "forecast40")
	
	# Test Methods  on object
	
	options(width = 120)
	zz <- file("test1c.txt", open="wt")
	sink(zz)
	cat("\nARFIMAforecast from ARFIMAfit and ARFIMAspec check:")
	cat("\nFit\n")	
	print(cfmatrix1, digits = 4)
	cat("\nSpec\n")	
	print(cfmatrix2, digits = 4)
	slotNames(forc[[1]])
	showMethods(classes="ARFIMAforecast")
	# summary
	print(show(forc[[1]]))
	# Extractor Functions
	# as array (array dimension [3] is 1 since n.roll = 0 i.e. no rolling beyond the first)
	print(as.array(forc[[1]]))
	# as.data.frame
	print(as.data.frame(forc[[1]]))
	# as.list
	print(as.list(forc[[1]]))
	sink(type="message")
	sink()
	close(zz)
	
	nforc = sapply(forc, FUN = function(x) t(unlist(as.data.frame(x, aligned = FALSE))))
	postscript("test1c.eps", width = 12, height = 5)
	dd = c(rownames(tail(sp500ret, 100)), rownames(as.data.frame(forc[[1]])))
	clrs = rainbow(9, alpha = 1, start = 0.4, end = 0.95)
	plot(as.Date(dd), c(tail(sp500ret[,1], 100), nforc[,1]), type = "l", ylim = c(-0.02, 0.02), col = "lightgrey",
			ylab = "", xlab = "", main = "100-ahead Unconditional Forecasts")
	for(i in 1:9){
		tmp = c(tail(sp500ret[,1], 100), rep(NA, 100))
		tmp[101:200] = nforc[1:100,i]
		lines(as.Date(dd), c(rep(NA, 100), tmp[-(1:100)]), col = clrs[i])
	}
	legend("topleft", legend = dist, col = clrs, fill = clrs, bty = "n")
	dev.off()
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test1d = function(parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2)){
	# rolling forecast
	tic = Sys.time()
	
	data(sp500ret)
	fit = vector(mode = "list", length = 9)
	dist = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "jsu")
	for(i in 1:9){
		spec = arfimaspec( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, arfima = FALSE), 
				distribution.model = dist[i])
		fit[[i]] = arfimafit(spec = spec, data = sp500ret, solver = "solnp", 
				out.sample = 1000, fit.control = list(scale = 1))
	}
	cfmatrix = matrix(NA, nrow = 9, ncol = 7)
	colnames(cfmatrix) = c("mu", "ar1", "ma1", "sigma", "skew", "shape", "ghlambda")
	rownames(cfmatrix) = dist
	
	for(i in 1:9){
		cf = coef(fit[[i]])
		cfmatrix[i, match(names(cf), colnames(cfmatrix))] =  cf
	}
	
	
	forc = vector(mode = "list", length = 9)
	for(i in 1:9){
		forc[[i]] = arfimaforecast(fit[[i]], n.ahead = 1, n.roll = 999)
	}
	rollforc = sapply(forc, FUN = function(x) t(unlist(as.data.frame(x, rollframe = "all", aligned = FALSE))))
	
	# forecast performance measures:
	fpmlist = vector(mode = "list", length = 9)
	for(i in 1:9){
		fpmlist[[i]] = fpm(forc[[i]], summary = FALSE)
	}
	
	postscript("test1d.eps", width = 16, height = 5)
	par(mfrow = c(1,2))
	dd = rownames(tail(sp500ret, 1250))
	clrs = rainbow(9, alpha = 1, start = 0.4, end = 0.95)
	plot(as.Date(dd), tail(sp500ret[,1], 1250), type = "l", ylim = c(-0.02, 0.02), col = "lightgrey",
			ylab = "", xlab = "", main = "Rolling 1-ahead Forecasts\nvs Actual")
	for(i in 1:9){
		tmp = tail(sp500ret[,1], 1250)
		tmp[251:1250] = rollforc[1:1000,i]
		lines(as.Date(dd), c(rep(NA, 250), tmp[-(1:250)]), col = clrs[i])
	}
	legend("topleft", legend = dist, col = clrs, fill = clrs, bty = "n")
	
	# plot deviation measures and range
	tmp = vector(mode = "list", length = 9)
	for(i in 1:9){
		tmp[[i]] = fpmlist[[i]][,"AE"]
		names(tmp[[i]]) = dist[i]
	}
	boxplot(tmp, col = clrs, names = dist, range  = 6, notch = TRUE, 
			main = "Rolling 1-ahead Forecasts\nAbsolute Deviation Loss")
	dev.off()
	
	# fpm comparison
	compm = matrix(NA, nrow = 3, ncol = 9)
	compm = sapply(fpmlist, FUN = function(x) c(mean(x[,"SE"]), mean(x[,"AE"]), mean(x[,"DAC"])))
	colnames(compm) = dist
	rownames(compm) = c("MSE", "MAD", "DAC")
	
	zz <- file("test1d.txt", open="wt")
	sink(zz)
	cat("\nRolling Forecast FPM\n")
	print(compm, digits = 4)
	cat("\nMethods Check\n")
	print(as.data.frame(forc[[1]], rollframe = 0))
	print(as.data.frame(forc[[1]], rollframe = 999))
	print(t(as.data.frame(forc[[1]], rollframe = "all", aligned = FALSE)))
	print(fpm(forc[[1]], summary = TRUE))
	print(show(forc[[1]]))
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test1e = function(parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2)){
	# Multi-Methods
	tic = Sys.time()
	
	data(dji30ret)
	Dat = dji30ret[, 1:3, drop = FALSE]
	
	#------------------------------------------------
	# Unequal Spec
	# Fit
	spec1 = arfimaspec(mean.model = list(armaOrder = c(2,1)))
	spec2 = arfimaspec(mean.model = list(armaOrder = c(2,2)))
	spec3 = arfimaspec(mean.model = list(armaOrder = c(1,1)), 
			distribution.model = "sstd")
	speclist = as.list(c(spec1, spec2, spec3))
	mspec = multispec( speclist )
	mfit1 = multifit(multispec = mspec, data = Dat, 
			fit.control = list(stationarity=1), parallel = parallel,
			parallel.control = parallel.control)
	# Filter
	fspec = vector(mode = "list", length = 3)
	fspec[[1]] = spec1
	fspec[[2]] = spec2
	fspec[[3]] = spec3
	for(i in 1:3){
		setfixed(fspec[[i]])<-as.list(coef(mfit1)[[i]])
	}
	mspec1 = multispec( fspec )	
	mfilt1 = multifilter(multifitORspec = mspec1, data = Dat, parallel = parallel,
			parallel.control = parallel.control)
	# Forecast from Fit
	mforc1 = multiforecast(mfit1, n.ahead = 10, parallel = parallel,
			parallel.control = parallel.control)
	# Forecast from Spec
	mforc11 = multiforecast(mspec1, data = Dat, n.ahead = 10, parallel = parallel,
			parallel.control = parallel.control)	
	#------------------------------------------------
	
	#------------------------------------------------
	# Equal Spec
	# Fit
	spec1 = arfimaspec(mean.model = list(armaOrder = c(1,1)))
	mspec = multispec( replicate(3, spec1) )
	mfit2 = multifit(multispec = mspec, data = Dat, parallel = parallel,
			parallel.control = parallel.control)
	# Filter
	fspec = vector(mode = "list", length = 3)
	fspec = replicate(3, spec1)
	for(i in 1:3){
		setfixed(fspec[[i]])<-as.list(coef(mfit2)[,i])
	}
	mspec2 = multispec( fspec )
	mfilt2 = multifilter(multifitORspec = mspec2, data = Dat, parallel = parallel,
			parallel.control = parallel.control)
	# Forecast From Fit
	mforc2 = multiforecast(mfit2, n.ahead = 10)
	# Forecast From Spec
	mforc21 = multiforecast(mspec2, data = Dat, n.ahead = 10, parallel = parallel,
			parallel.control = parallel.control)
	#------------------------------------------------

	#------------------------------------------------
	# Equal Spec/Same Data
	# Fit
	spec1 = arfimaspec(mean.model = list(armaOrder = c(1,1)))
	spec2 = arfimaspec(mean.model = list(armaOrder = c(2,1)))
	spec3 = arfimaspec(mean.model = list(armaOrder = c(3,1)))
	speclist = as.list(c(spec1, spec2, spec3))
	mspec = multispec( speclist )
	mfit3 = multifit(multispec = mspec, data = cbind(Dat[,1], Dat[,1], Dat[,1]))
	# Forecast
	mforc3 = multiforecast(mfit3, n.ahead = 10, parallel = parallel,
			parallel.control = parallel.control)
	#------------------------------------------------
	
	zz <- file("test1e.txt", open="wt")
	sink(zz)
	cat("\nMultifit Evaluation\n")
	cat("\nUnequal Spec\n")
	print(show(mfit1))
	print(likelihood(mfit1))
	print(coef(mfit1))
	print(head(fitted(mfit1)))
	print(head(residuals(mfit1)))
	print(show(mfilt1))
	print(likelihood(mfilt1))
	print(coef(mfilt1))
	print(head(fitted(mfilt1)))
	print(head(residuals(mfilt1)))
	print(show(mforc1))
	print(as.array(mforc1))
	print(as.list(mforc1))
	print(show(mforc11))
	print(as.array(mforc11))
	print(as.list(mforc11))
	cat("\nEqual Spec\n")
	print(show(mfit2))
	print(likelihood(mfit2))
	print(coef(mfit2))
	print(head(fitted(mfit2)))
	print(head(residuals(mfit2)))
	print(show(mfilt2))
	print(likelihood(mfilt2))
	print(coef(mfilt2))
	print(head(fitted(mfilt2)))
	print(head(residuals(mfilt2)))
	print(show(mforc2))
	print(as.array(mforc2))
	print(as.list(mforc2))
	print(show(mforc21))
	print(as.array(mforc21))
	print(as.list(mforc21))
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test1f = function(parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2)){
	# rolling fit/forecast
	tic = Sys.time()
	
	data(sp500ret)
	spec = arfimaspec()
	roll1 = arfimaroll(spec,  data = sp500ret, n.ahead = 1, forecast.length = 500, 
			refit.every = 25, refit.window = "moving", parallel = parallel,
			parallel.control = parallel.control,
			solver = "solnp", fit.control = list(), solver.control = list() ,
			calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05))

	roll2 = arfimaroll(spec,  data = sp500ret, n.ahead = 2, forecast.length = 500, 
			refit.every = 25, refit.window = "moving", parallel = parallel,
			parallel.control = parallel.control,
			solver = "solnp", fit.control = list(), solver.control = list() ,
			calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05))
	
	# as.ARFIMAforecast
	# as.data.frame
	
	zz <- file("test1f.txt", open="wt")
	sink(zz)
	cat("\nForecast Evaluation\n")
	report(roll1, "VaR")
	report(roll1, "fpm")
	report(roll2, "VaR")
	# equivalent to:
	report(roll2, "VaR", n.ahead = 1, VaR.alpha = 0.01, conf.level = 0.95)
	report(roll2, "VaR", n.ahead = 2, VaR.alpha = 0.05, conf.level = 0.95)
	report(roll2, "fpm")
	# Extractor Functions:
	# default:
	print(t(as.data.frame(roll2, which = "series", refit = 1)))
	print(head(as.data.frame(roll2, which = "series", refit = "all", n.ahead = 1)))
	print(head(as.data.frame(roll2, which = "series", refit = "all", n.ahead = 2)))
	print(as.data.frame(roll2, which = "coefs"))
	print(as.data.frame(roll2, which = "coefmat", refit = 1))
	print(as.data.frame(roll2, which = "coefmat", refit = 20))
	print(as.data.frame(roll2, which = "LLH"))
	print(head(as.data.frame(roll2, which = "density", n.ahead = 1)))
	print(head(as.data.frame(roll2, which = "density", n.ahead = 2)))
	print(head(as.data.frame(roll2, which = "VaR", n.ahead = 1)))
	print(head(as.data.frame(roll2, which = "VaR", n.ahead = 2)))
	sink(type="message")
	sink()
	close(zz)
	
	# coefficient plots
	coefx = as.data.frame(roll2, which = "coefs")
	coefv = matrix(NA, ncol = 4, nrow = 20)
	coefse = matrix(NA, ncol = 4, nrow = 20)
	parnames = c("mu", "ar1", "ma1", "sigma")
	for(i in 1:20){
		coefv[i, 1] = coefx[i,1]
		coefv[i, 2] = coefx[i,2]
		coefv[i, 3] = coefx[i,3]
		coefv[i, 4] = coefx[i,4]
		coefxse = as.data.frame(roll2, which = "coefmat", refit = i)
		coefse[i, 1] = coefxse[1,2]
		coefse[i, 2] = coefxse[2,2]
		coefse[i, 3] = coefxse[3,2]
		coefse[i, 4] = coefxse[4,2]
	}
	postscript("test1f.eps", width = 16, height = 5)
	par(mfrow = c(2,2))
	for(i in 1:4){
		plot(1:20, coefv[,i], type = "l", main = paste("coef: ", parnames[i], " with s.e. bands", sep = ""),
				ylab = "", xlab = "refit", ylim = c(min(coefv[,i]-1.96*coefse[, i]), max(coefv[,i]+1.96*coefse[, i])))
		lines(1:20, coefv[,i]+1.96*coefse[, i], lty=2, col = "red")
		lines(1:20, coefv[,i]-1.96*coefse[, i], lty=2, col = "red")
	}
	dev.off()
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test1g = function(parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2)){
	# simulation
	tic = Sys.time()
	
	spec1 = arfimaspec( mean.model = list(armaOrder = c(2,1), include.mean = TRUE, arfima = TRUE), 
			distribution.model = "std", fixed.pars = list(mu = 0.02, ar1 = 0.6, ar2 = 0.01, ma1 = -0.7, arfima = 0,
					shape = 5, sigma = 0.0123))
	spec2 = arfimaspec( mean.model = list(armaOrder = c(2,1), include.mean = TRUE, arfima = FALSE), 
			distribution.model = "std", fixed.pars = list(mu = 0.02, ar1 = 0.6, ar2 = 0.01, ma1 = -0.7,
					shape = 5,sigma = 0.0123))
	sim1 = arfimapath(spec1, n.sim = 1000, m.sim = 1, rseed = 100)
	sim2 = arfimapath(spec2, n.sim = 1000, m.sim = 1, rseed = 100)
	
	zz <- file("test1g-1.txt", open="wt")
	sink(zz)
	cat("\nARFIMA and ARMA simulation tests:\n")
	print(head(as.data.frame(sim1)), digits = 5)
	print(head(as.data.frame(sim2)), digits = 5)
	sink(type="message")
	sink()
	close(zz)
	
	# Now the rugarch simulation of ARFIMA/ARMA with arima.sim of R
	# Note that arima.sim simulates the residuals (i.e no mean):
	#  ARMA(2,2)
	set.seed(33)
	inn = rdist("std", 1000, mu = 0, sigma = 1, lambda = 0, skew = 0, shape = 5)
	
	spec1 = arfimaspec( mean.model = list(armaOrder = c(2,2), include.mean = FALSE, arfima = TRUE), 
			distribution.model = "std", fixed.pars = list(ar1 = 0.6, ar2 = 0.21, arfima = 0, ma1 = -0.7,
					ma2 = 0.3, shape = 5, sigma = 0.0123))
	spec2 = arfimaspec( mean.model = list(armaOrder = c(2,2), include.mean = FALSE, arfima = FALSE), 
			distribution.model = "std", fixed.pars = list(ar1 = 0.6, ar2 = 0.21, ma1 = -0.7,
					ma2 = 0.3, shape = 5,sigma = 0.0123))
	
	sim1 = arfimapath(spec1, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			custom.dist = list(name = "sample",
					distfit = matrix(inn, ncol = 1), type = "z"))
	sim2 = arfimapath(spec2, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			custom.dist = list(name = "sample",
					distfit = matrix(inn, ncol = 1), type = "z"))
	
	# Note that we pass the non-standardized innovations to arima.sim (i.e. multiply by sigma)
	sim3 = arima.sim(model = list(ar = c(0.6, 0.21), ma = c(-0.7, 0.3)), n = 1000,
			n.start = 4, start.innov = c(0,0,0,0), innov = inn*0.0123)
	
	tst = cbind(head(as.data.frame(sim1)), head(as.data.frame(sim2)), head(as.data.frame(sim3)))
	colnames(tst) = c("ARFIMA(d = 0)", "ARMA", "arima.sim")
	
	zz <- file("test1g-2.txt", open="wt")
	sink(zz)
	cat("\nARFIMA, ARMA arima.sim simulation tests:\n")
	print(tst, digits = 6)
	sink(type="message")
	sink()
	close(zz)
	
	#  ARMA(2,1)
	set.seed(33)
	inn = rdist("std", 1000, mu = 0, sigma = 1, lambda = 0, skew = 0, shape = 5)
	
	spec1 = arfimaspec( mean.model = list(armaOrder = c(2,1), include.mean = FALSE, arfima = TRUE), 
			distribution.model = "std", fixed.pars = list(ar1 = 0.6, ar2 = 0.21, arfima = 0, ma1 = -0.7,
					shape = 5, sigma = 0.0123))
	spec2 = arfimaspec( mean.model = list(armaOrder = c(2,1), include.mean = FALSE, arfima = FALSE), 
			distribution.model = "std", fixed.pars = list(ar1 = 0.6, ar2 = 0.21, ma1 = -0.7,
					shape = 5,sigma = 0.0123))
	
	sim1 = arfimapath(spec1, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			custom.dist = list(name = "sample",
					distfit = matrix(inn, ncol = 1), type = "z"))
	sim2 = arfimapath(spec2, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			custom.dist = list(name = "sample",
					distfit = matrix(inn, ncol = 1), type = "z"))
	
	# Note that we pass the non-standardized innovations to arima.sim (i.e. multiply by sigma)
	sim3 = arima.sim(model = list(ar = c(0.6, 0.21), ma = c(-0.7)), n = 1000,
			n.start = 3, start.innov = c(0,0,0), innov = inn*0.0123)
	
	tst = cbind(head(as.data.frame(sim1)), head(as.data.frame(sim2)), head(as.data.frame(sim3)))
	colnames(tst) = c("ARFIMA(d = 0)", "ARMA", "arima.sim")
	
	zz <- file("test1g-3.txt", open="wt")
	sink(zz)
	cat("\nARFIMA, ARMA arima.sim simulation tests:\n")
	print(tst, digits = 6)
	sink(type="message")
	sink()
	close(zz)
	
	
	#  Pure AR
	set.seed(33)
	inn = rdist("std", 1000, mu = 0, sigma = 1, lambda = 0, skew = 0, shape = 5)
	
	spec1 = arfimaspec( mean.model = list(armaOrder = c(2,0), include.mean = FALSE, arfima = TRUE), 
			distribution.model = "std", fixed.pars = list(ar1 = 0.6, ar2 = 0.21, arfima = 0, ma1 = -0.7,
					shape = 5, sigma = 0.0123))
	spec2 = arfimaspec( mean.model = list(armaOrder = c(2,0), include.mean = FALSE, arfima = FALSE), 
			distribution.model = "std", fixed.pars = list(ar1 = 0.6, ar2 = 0.21, ma1 = -0.7,
					shape = 5,sigma = 0.0123))
	
	sim1 = arfimapath(spec1, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			custom.dist = list(name = "sample",
					distfit = matrix(inn, ncol = 1), type = "z"))
	sim2 = arfimapath(spec2, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			custom.dist = list(name = "sample",
					distfit = matrix(inn, ncol = 1), type = "z"))
	
	# Note that we pass the non-standardized innovations to arima.sim (i.e. multiply by sigma)
	sim3 = arima.sim(model = list(ar = c(0.6, 0.21), ma = NULL), n = 1000,
			n.start = 2, start.innov = c(0,0), innov = inn*0.0123)
	
	tst = cbind(head(as.data.frame(sim1)), head(as.data.frame(sim2)), head(as.data.frame(sim3)))
	colnames(tst) = c("ARFIMA(d = 0)", "ARMA", "arima.sim")
	
	zz <- file("test1g-4.txt", open="wt")
	sink(zz)
	cat("\nARFIMA, ARMA arima.sim simulation tests:\n")
	print(tst, digits = 6)
	sink(type="message")
	sink()
	close(zz)
	
	
	#  Pure MA
	set.seed(33)
	inn = rdist("std", 1000, mu = 0, sigma = 1, lambda = 0, skew = 0, shape = 5)
	
	spec1 = arfimaspec( mean.model = list(armaOrder = c(0,2), include.mean = FALSE, arfima = TRUE), 
			distribution.model = "std", fixed.pars = list(ma1 = 0.6, ma2 = -0.21, arfima = 0,
					shape = 5, sigma = 0.0123))
	spec2 = arfimaspec( mean.model = list(armaOrder = c(0,2), include.mean = FALSE, arfima = FALSE), 
			distribution.model = "std", fixed.pars = list(ma1 = 0.6, ma2 = -0.21,
					shape = 5,sigma = 0.0123))
	
	sim1 = arfimapath(spec1, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			custom.dist = list(name = "sample",
					distfit = matrix(inn, ncol = 1), type = "z"))
	sim2 = arfimapath(spec2, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			custom.dist = list(name = "sample",
					distfit = matrix(inn, ncol = 1), type = "z"))
	
	# Note that we pass the non-standardized innovations to arima.sim (i.e. multiply by sigma)
	sim3 = arima.sim(model = list(ar = NULL, ma = c(0.6, -0.21)), n = 1000,
			n.start = 2, start.innov = c(0,0), innov = inn*0.0123)
	
	tst = cbind(head(as.data.frame(sim1)), head(as.data.frame(sim2)), head(as.data.frame(sim3)))
	colnames(tst) = c("ARFIMA(d = 0)", "ARMA", "arima.sim")
	
	zz <- file("test1g-5.txt", open="wt")
	sink(zz)
	cat("\nARFIMA, ARMA arima.sim simulation tests:\n")
	print(tst, digits = 6)
	sink(type="message")
	sink()
	close(zz)
	
	# arfimasim + exogenous regressors + custom innovations
	data(dji30ret)
	Dat = dji30ret[,1, drop = FALSE]
	T = dim(Dat)[1]
	Bench = as.matrix(cbind(apply(dji30ret[,2:10], 1, "mean"), apply(dji30ret[,11:20], 1, "mean")))	

	spec = arfimaspec( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, arfima = FALSE, 
					external.regressors = Bench), distribution.model = "std")
	fit = arfimafit(spec = spec, data = Dat, solver = "solnp", out.sample = 500)
	
	# lag1 Benchmark
	BenchF = Bench[(T-500):(T-500+9), , drop = FALSE]
	
	exsim = vector(mode = "list", length = 10000)
	for(i in 1:10000) exsim[[i]] = as.matrix(BenchF)
	# simulated residuals
	res = residuals(fit)
	ressim = matrix(NA, ncol = 10000, nrow = 10)
	set.seed(10000)
	for(i in 1:10000) ressim[,i] = sample(res, 10, replace = TRUE)
	
	sim = arfimasim(fit, n.sim = 10, m.sim = 10000, startMethod="sample", 
			custom.dist = list(name = "sample", distfit = ressim, type = "res"), mexsimdata = exsim)
	
		
	forc = as.data.frame(arfimaforecast(fit, n.ahead = 10, external.forecasts = list(mregfor = BenchF)))
	simx = as.data.frame(sim)
	actual10 = Dat[(T-500+1):(T-500+10), 1, drop = FALSE]
	
	simm = apply(simx, 1 ,"mean")
	simsd = apply(simx, 1 ,"sd")
	
	zz <- file("test1g-6.txt", open="wt")
	sink(zz)
	print(round(cbind(actual10, forc, simm, simsd),5), digits = 4)
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

# ARFIMA benchmark tests
rugarch.test1h = function(parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2)){
	tic = Sys.time()	
	# ARFIMA(2,d,1)
	truecoef1 = list(mu = 0.005, ar1 = 0.6, ar2 = 0.01, ma1 = -0.7, arfima = 0.3, sigma = 0.0123)
	spec1 = arfimaspec( 
	mean.model = list(armaOrder = c(2,1), include.mean = TRUE, arfima = TRUE), 
	distribution.model = "norm", 
	fixed.pars = truecoef1)
	sim1 = arfimapath(spec1, n.sim = 5000, n.start = 100, m.sim = 1)
	data1 = as.data.frame(sim1)
	#write.csv(data1[,1], file = "D:/temp1.csv")
	spec1 = arfimaspec( 
	mean.model = list(armaOrder = c(2,1), include.mean = TRUE, arfima = TRUE), 
	distribution.model = "norm", 
	start.pars = truecoef1)
	fit1 = arfimafit(spec1, data = data1)
	# Commercial Implementation Program Fit (NLS-with imposed stationarity):
	commcheck1 = c(0.005906822, 0.463559, -0.0337927, -0.562876, 0.310268, 0.0124086)
	chk1 = cbind(coef(fit1), commcheck1, unlist(truecoef1))
	colnames(chk1) = c("rugarch", "check", "true")
	chk1lik = c(likelihood(fit1),  14852.6215)
	
	# ARFIMA(2,d,0)
	truecoef2 = list(mu = 0.005, ar1 = 0.6, ar2 = 0.01, arfima = -0.3, sigma = 0.0123)
	spec2 = arfimaspec( 
	mean.model = list(armaOrder = c(2,0), include.mean = TRUE, arfima = TRUE), 
	distribution.model = "norm", 
	fixed.pars = truecoef2)
	sim2 = arfimapath(spec2, n.sim = 5000, n.start = 100, m.sim = 1)	
	data2 = as.data.frame(sim2)
	#write.csv(data2[,1], file = "D:/temp2.csv")
	spec2 = arfimaspec( 
	mean.model = list(armaOrder = c(2,0), include.mean = TRUE, arfima = TRUE), 
	distribution.model = "norm", 
	start.pars = truecoef2)
	fit2 = arfimafit(spec2, data = data2)
	commcheck2 = c(0.00498073,  0.632362,  -0.0225118, -0.308095,0.0122404)
	chk2 = cbind(coef(fit2), commcheck2 , unlist(truecoef2))
	colnames(chk2) = c("rugarch", "check", "true")
	chk2lik = c(likelihood(fit2), 14920.8979)

	# ARFIMA(0,d,2)
	truecoef3 = list(mu = 0.005, ma1 = 0.3, ma2 = 0.2, arfima = 0.1, sigma = 0.0123)
	spec3 = arfimaspec( 
	mean.model = list(armaOrder = c(0,2), include.mean = TRUE, arfima = TRUE), 
	distribution.model = "norm", 
	fixed.pars = truecoef3)
	sim3 = arfimapath(spec3, n.sim = 5000, n.start = 100, m.sim = 1)
	data3 = as.data.frame(sim3)
	#write.csv(data3[,1], file = "D:/temp3.csv")
	spec3 = arfimaspec( 
	mean.model = list(armaOrder = c(0,2), include.mean = TRUE, arfima = TRUE), 
	distribution.model = "norm", 
	start.pars = truecoef3)
	fit3 = arfimafit(spec3, data = data3)	
	commcheck3 = c(0.00498073,   0.283231, 0.195064, 0.0987104, 0.0122155)
	chk3 = cbind(coef(fit3), commcheck3 , unlist(truecoef3))
	colnames(chk3) = c("rugarch", "check", "true")
	chk3lik = c(likelihood(fit3),  14931.0414)
	
	
	# ARFIMA(2,d,1) simulation:
	truecoef = list(mu = 0.005, ar1 = 0.6, ar2 = 0.01, ma1 = -0.7, arfima = 0.3, sigma = 0.0123)
	spec = arfimaspec( 
	mean.model = list(armaOrder = c(2,1), include.mean = TRUE, arfima = TRUE), 
	distribution.model = "norm", fixed.pars = truecoef)
	sim = arfimapath(spec, n.sim = 5000, n.start = 100, m.sim = 50, rseed = 1:50)
	Data = as.data.frame(sim)
	spec = arfimaspec( 
	mean.model = list(armaOrder = c(2,1), include.mean = TRUE, arfima = TRUE), 
	distribution.model = "norm", 
	start.pars = truecoef)
	coefx = matrix(NA, ncol = 6, nrow = 50)
	for(i in 1:50){
		fit = arfimafit(spec, data = Data[,i])
		if(fit@fit$convergence == 0) coefx[i,] = coef(fit)
	}
	
	zz <- file("test1h.txt", open="wt")
	sink(zz)
	cat("\nARFIMA(2,d,1)\n")
	print(chk1)
	print(chk1lik)
	cat("\nARFIMA(2,d,0)\n")
	print(chk2)
	print(chk2lik)
	cat("\nARFIMA(0,d,2)\n")
	print(chk3)
	print(chk3lik)
	cat("\nARFIMA(2,d,1) mini-simulation/fit\n")
	# small sample/simulation use median:
	print(	cbind(round(apply(coefx, 2, "median"),5), unlist(truecoef)) )	
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}