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
.plotgarchfit = function(x, which="ask",...)
{
	old.par <- par(no.readonly = TRUE)
	choices = c(
			"Series with 2 Conditional SD Superimposed",
			"Series with 2.5% VaR Limits (with unconditional mean)",
			"Conditional SD",
			"ACF of Observations",
			"ACF of Squared Observations",
			"ACF of Absolute Observations",
			"Cross Correlation",
			"Empirical Density of Standardized Residuals",
			"QQ-Plot of Standardized Residuals",
			"ACF of Standardized Residuals",
			"ACF of Squared Standardized Residuals",
			"News-Impact Curve")
	.intergarchfitPlot(x, choices = choices, plotFUN = paste(".plot.garchfit", 1:12, sep = "."), which = which, ...)
	# Return Value:
	par(old.par)
	invisible(x)
}

.intergarchfitPlot = function(x, choices, plotFUN, which, ...)
{
	if (is.numeric(which)) {
		if(which>length(choices)) stop("Not a valid choice. Plots choices are 1-12.\n",call. = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x)
	}
	if(is.character(which))
	{
		if(which!="all" && which!="ask") stop("Not a valid choice.\n",call. = FALSE)
		if (which[1] == "all") {
			#Which = rep(TRUE, times = length(choices))
			par(mfrow=c(3,4))
			for(i in 1:12){
				FUN = match.fun(plotFUN[i])
				FUN(x)
			}
		} else{
			.multgarchfitPlot(x, choices, plotFUN, ...)
		}
	}
	
	invisible(x)
}

.multgarchfitPlot = function(x, choices, ...)
{
	
	pick = 1
	while (pick > 0) {
		pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
		switch (pick,
				.plot.garchfit.1(x),  .plot.garchfit.2(x),  .plot.garchfit.3(x),
				.plot.garchfit.4(x),  .plot.garchfit.5(x),  .plot.garchfit.6(x),
				.plot.garchfit.7(x),  .plot.garchfit.8(x),  .plot.garchfit.9(x),
				.plot.garchfit.10(x), .plot.garchfit.11(x), .plot.garchfit.12(x))
	}
}

# Series with 2 Conditional SD Superimposed
.plot.garchfit.1 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xdates  = x@model$modeldata$dates[insample]
	xseries = x@model$modeldata$data[insample]
	xsigma  = x@fit$sigma
	ci = 2
	plot(xdates, xseries, type = "l", col = "steelblue", ylab = "Returns", xlab="Time", 
			main = "Series with 2 Conditional SD Superimposed", cex.main = 0.8)
	lines(xdates, +ci * xsigma, col = "tomato1")
	lines(xdates, -ci * xsigma, col = "tomato1")
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	abline(h = 0, col = "grey", lty = 3)
	grid()
}

# Series with 2.5% VaR Limits
.plot.garchfit.2 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	xdates  = x@model$modeldata$dates[insample]
	xsigma 	= x@fit$sigma
	distribution = x@model$modeldesc$distribution
	xcmu 	= fitted(x)
	idx = x@model$pidx
	pars  = x@fit$ipars[,1]
	skew  = pars[idx["skew",1]]
	shape = pars[idx["shape",1]]
	if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]
	z1 	= 0.025
	z2 	= 0.975
	cat("\nplease wait...calculating quantiles...\n")
	q025 	= fitted(x) + sigma(x)* qdist(distribution, z1, 0, 1, lambda = ghlambda, skew, shape)
	q975 	= fitted(x) + sigma(x)* qdist(distribution, z2, 0, 1, lambda = ghlambda, skew, shape)
	plot(xdates, xseries, type = "l", col = "steelblue", ylab = "Returns", xlab="Time", main = "Series with with 2.5% VaR Limits", cex.main = 0.8)
	lines(xdates, q025, col = "tomato1")
	lines(xdates, q975, col = "green")
	mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	abline(h = 0, col = "grey", lty = 3)
	grid()
}

# Conditional SD
.plot.garchfit.3 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xdates  = x@model$modeldata$dates[insample]
	xsigma 	= x@fit$sigma
	plot(xdates, xsigma, type = "l", col = "steelblue", ylab = "Volatility", xlab="Time", main = "Conditional SD", cex.main = 0.8)
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	abline(h = 0, col = "grey", lty = 3)
	grid()
}

# ACF of observations
.plot.garchfit.4 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	lag.max = as.integer(10*log10(T))
	acfx 	= acf(xseries, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode = "character", length = lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim = 1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Observations", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# ACF of squared observations
.plot.garchfit.5 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	lag.max = as.integer(10*log10(T))
	acfx 	= acf(xseries^2, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode="character",length=lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim = 1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Squared Observations", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# ACF of absolute observations
.plot.garchfit.6 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	lag.max = as.integer(10*log10(T))
	acfx 	= acf(abs(xseries), lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode="character",length=lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
	barplot(height=as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim = 1.2*ylim, col = clx,
			ylab = "ACF", xlab = "lag", main = "ACF of Absolute Observations", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# Cross Correlations
.plot.garchfit.7 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]	
	lag.max = as.integer(10*log10(T))
	ccfx	= ccf(xseries^2, xseries, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(ccfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(ccfx$acf)))
	clx 	= vector(mode="character",length=(2*lag.max)+1)
	clx[which(as.numeric(ccfx$acf)>=0)]="steelblue"
	clx[which(as.numeric(ccfx$acf)<0)]="orange"
	barplot(height=as.numeric(ccfx$acf), names.arg=as.numeric(ccfx$lag),ylim=1.2*ylim,col=clx,
			ylab = "ACF", xlab="lag", main = "Cross-Correlations of \nSquared vs Actual Observations", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black")
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# Standardized Residuals Empirical Denisty
.plot.garchfit.8 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	zseries = x@fit$z
	distribution = x@model$modeldesc$distribution
	fd 	= .getDensity(distribution)
	idx = x@model$pidx
	pars  = x@fit$ipars[,1]
	skew  = pars[idx["skew",1]]
	shape = pars[idx["shape",1]]
	if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]	
	xmean 	= mean(zseries)
	xmedian = median(zseries)
	xsd 	= sd(zseries)	
	xlim 	= c(min(zseries), max(zseries))
	result 	= hist(x = zseries, col = "grey", border = "white",
			breaks = "Scott", main = "Empirical Density of Standardized Residuals", xlim = xlim, ylim = c(0,0.6),
			probability = TRUE, ylab="Probability", cex.main = 0.8, ...)
	box()
	#s 	= seq(xlim[1], xlim[2], length = 201)
	s = result$breaks
	y	= fd(z = s, hh = 1, lambda = ghlambda, skew = skew, shape = shape)
	lines(s, dnorm(s, 0, 1), lwd = 2, col = "blue")
	lines(s, y, lwd = 2, col = "orange")
	abline(v = xmean, lwd = 2, col = "red")
	abline(v = xmedian, lwd = 2, col = "darkgreen")
	Text = paste("Median: ", round(xmedian, 2), "| Mean: ",signif(xmean, 3))
	mtext(Text, side = 3, adj = 0, col = "darkgrey", cex = 0.7)
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
	}
	abline(h = 0, col = "grey")
	lg.txt = c("normal Density", paste(distribution," (0,1) Fitted Density",sep=""))
	legend("topleft", legend = lg.txt, col = c("blue","orange"), pch = 1, cex = 0.8, bty = "n")
	grid()
}

# QQ-Plot of Standardized Residuals
.plot.garchfit.9 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	zseries = x@fit$z
	distribution = x@model$modeldesc$distribution
	fd 	= .getDensity(distribution)
	idx = x@model$pidx
	pars  = x@fit$ipars[,1]
	skew  = pars[idx["skew",1]]
	shape = pars[idx["shape",1]]
	if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]
	.qqDist(y = zseries, dist = distribution, lambda = ghlambda, skew = skew, shape = shape)
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
	}
}

### Continue Here

# ACF of standardized residuals
.plot.garchfit.10 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	zseries = x@fit$z
	zseries[is.na(zseries)] = 0
	n 		= length(zseries)
	lag.max = as.integer(10*log10(n))
	acfx 	= acf(zseries, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode = "character", length = lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim=1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Standardized Residuals", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# ACF of squared standardized residuals
.plot.garchfit.11 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	zseries = x@fit$z
	zseries[is.na(zseries)] = 0
	n 		= length(zseries)
	lag.max = as.integer(10*log10(n))
	acfx 	= acf(zseries^2, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode="character",length=lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)]="steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)]="orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim=1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Squared Standardized Residuals", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# News Impact Curve
.plot.garchfit.12 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	if(vmodel == "iGARCH"){
		warning(paste("\nplot-->: iGARCH newsimpact not available"))
	} else{
		ni = newsimpact(z = NULL, x)
		ni.y = ni$zy
		ni.x = ni$zx
		xf = ni$xexpr
		yf  = ni$yexpr
		plot( ni.x, ni.y, ylab = yf, xlab = xf, type = "l", lwd = 2, col = "steelblue", main = "News Impact Curve", cex.main = 0.8)
		#mtext(yf, side = 2, adj = 0.5, padj = -2.5, cex = 0.85)
		mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
		if(vmodel == "fGARCH"){
			mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
		}
		grid()
	}
}

#-------------------------------------------------------------------------------
# SECTION GARCH filter plots
#-------------------------------------------------------------------------------

.plotgarchfilter = function(x, which="ask",...)
{
	old.par <- par(no.readonly = TRUE)
	choices = c(
			"Series with 2 Conditional SD Superimposed",
			"Series with 2.5% VaR Limits (with unconditional mean)",
			"Conditional SD",
			"ACF of Observations",
			"ACF of Squared Observations",
			"ACF of Absolute Observations",
			"Cross Correlation",
			"Empirical Density of Standardized Residuals",
			"QQ-Plot of Standardized Residuals",
			"ACF of Standardized Residuals",
			"ACF of Squared Standardized Residuals",
			"News-Impact Curve")
	.intergarchfilterPlot(x, choices = choices, plotFUN = paste(".plot.garchfilter", 1:12, 
					sep = "."), which = which, ...)
	# Return Value:
	par(old.par)
	invisible(x)
}

.intergarchfilterPlot = function(x, choices, plotFUN, which, ...)
{
	if (is.numeric(which)) {
		if(which>length(choices)) 
			stop("Not a valid choice. Plots choices are 1-12.\n",call. = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x)
	}
	if(is.character(which))
	{
		if(which!="all" && which!="ask") stop("Not a valid choice.\n",call. = FALSE)
		if (which[1] == "all") {
			#Which = rep(TRUE, times = length(choices))
			par(mfrow=c(3,4))
			for(i in 1:12){
				FUN = match.fun(plotFUN[i])
				FUN(x)
			}
		} else{
			.multgarchfilterPlot(x, choices, plotFUN, ...)
		}
	}
	
	invisible(x)
}

.multgarchfilterPlot = function(x, choices, ...)
{
	
	pick = 1
	while (pick > 0) {
		pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
		switch (pick,
				.plot.garchfilter.1(x),  .plot.garchfilter.2(x),  .plot.garchfilter.3(x),
				.plot.garchfilter.4(x),  .plot.garchfilter.5(x),  .plot.garchfilter.6(x),
				.plot.garchfilter.7(x),  .plot.garchfilter.8(x),  .plot.garchfilter.9(x),
				.plot.garchfilter.10(x), .plot.garchfilter.11(x), .plot.garchfilter.12(x))
	}
}

# Series with 2 Conditional SD Superimposed
.plot.garchfilter.1 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	xdates  = x@model$modeldata$dates[insample]
	xsigma  = x@filter$sigma
	ci = 2
	plot(xdates, xseries, type = "l", col = "steelblue", ylab = "Returns", xlab="Time", 
			main = "Series with 2 Conditional SD Superimposed", cex.main = 0.8)
	lines(xdates, +ci * xsigma, col = "tomato1")
	lines(xdates, -ci * xsigma, col = "tomato1")
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	abline(h = 0, col = "grey", lty = 3)
	grid()
}

# Series with 2.5% VaR Limits
.plot.garchfilter.2 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	xdates  = x@model$modeldata$dates[insample]
	xsigma 	= x@filter$sigma
	distribution = x@model$modeldesc$distribution
	xcmu 	= fitted(x)
	idx = x@model$pidx
	pars  = x@filter$ipars[,1]
	skew  = pars[idx["skew",1]]
	shape = pars[idx["shape",1]]
	if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]
	z1 	= rep(0.025, length(xsigma))
	z2 	= rep(0.975, length(xsigma))
	spars 	= .scaledist(distribution, xcmu, xsigma, ghlambda, skew, shape)
	cat("\nplease wait...calculating quantiles...\n")
	q025 	= .qdensity(z1, spars[,"mu"], spars[,"sigma"], lambda = ghlambda, spars[,"skew"], spars[,"shape"], distribution = distribution)
	q975 	= .qdensity(z2, spars[,"mu"], spars[,"sigma"], lambda = ghlambda, spars[,"skew"], spars[,"shape"], distribution = distribution)
	plot(xdates, xseries, type = "l", col = "steelblue", ylab = "Returns", xlab="Time", main = "Series with with 2.5% VaR Limits", cex.main = 0.8)
	lines(xdates, q025, col = "tomato1")
	lines(xdates, q975, col = "green")
	mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	abline(h = 0, col = "grey", lty = 3)
	grid()
}

# Conditional SD
.plot.garchfilter.3 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xdates  = x@model$modeldata$dates[insample]
	xsigma 	= x@filter$sigma
	plot(xdates, xsigma, type = "l", col = "steelblue", ylab = "Volatility", xlab="Time", main = "Conditional SD", cex.main = 0.8)
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	abline(h = 0, col = "grey", lty = 3)
	grid()
}

# ACF of observations
.plot.garchfilter.4 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	lag.max = as.integer(10*log10(T))
	acfx 	= acf(xseries[insample], lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode = "character", length = lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim = 1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Observations", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# ACF of squared observations
.plot.garchfilter.5 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	xdates  = x@model$modeldata$dates[insample]
	lag.max = as.integer(10*log10(T))
	acfx 	= acf(xseries^2, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode="character",length=lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim = 1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Squared Observations", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# ACF of absolute observations
.plot.garchfilter.6 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	lag.max = as.integer(10*log10(T))
	acfx 	= acf(abs(xseries), lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode="character",length=lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
	barplot(height=as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim = 1.2*ylim, col = clx,
			ylab = "ACF", xlab = "lag", main = "ACF of Absolute Observations", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# Cross Correlations
.plot.garchfilter.7 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	lag.max = as.integer(10*log10(T))
	ccfx	= ccf(xseries^2, xseries, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(ccfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(ccfx$acf)))
	clx 	= vector(mode="character",length=(2*lag.max)+1)
	clx[which(as.numeric(ccfx$acf)>=0)]="steelblue"
	clx[which(as.numeric(ccfx$acf)<0)]="orange"
	barplot(height=as.numeric(ccfx$acf), names.arg=as.numeric(ccfx$lag),ylim=1.2*ylim,col=clx,
			ylab = "ACF", xlab="lag", main = "Cross-Correlations of \nSquared vs Actual Observations", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black")
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# Standardized Residuals Empirical Denisty
.plot.garchfilter.8 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	zseries = x@filter$z
	distribution = x@model$modeldesc$distribution
	fd 	= .getDensity(distribution)
	idx = x@model$pidx
	pars  = x@filter$ipars[,1]
	skew  = pars[idx["skew",1]]
	shape = pars[idx["shape",1]]
	if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]
	xmean 	= mean(zseries)
	xmedian = median(zseries)
	xsd 	= sd(zseries)	
	xlim 	= c(min(zseries), max(zseries))
	result 	= hist(x = zseries, col = "grey", border = "white",
			breaks = "Scott", main = "Empirical Density of Standardized Residuals", xlim = xlim, ylim = c(0,0.6),
			probability = TRUE, ylab="Probability", cex.main = 0.8, ...)
	box()
	#s 	= seq(xlim[1], xlim[2], length = 201)
	s = result$breaks
	y	= fd(z = s, hh = 1, lambda = ghlambda, skew = skew, shape = shape)
	lines(s, dnorm(s, 0, 1), lwd = 2, col = "blue")
	lines(s, y, lwd = 2, col = "orange")
	abline(v = xmean, lwd = 2, col = "red")
	abline(v = xmedian, lwd = 2, col = "darkgreen")
	Text = paste("Median: ", round(xmedian, 2), "| Mean: ",signif(xmean, 3))
	mtext(Text, side = 3, adj = 0, col = "darkgrey", cex = 0.7)
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
	}
	abline(h = 0, col = "grey")
	lg.txt = c("normal Density", paste(distribution," (0,1) Fitted Density",sep=""))
	legend("topleft", legend = lg.txt, col = c("blue","orange"), pch = 1, cex = 0.8, bty = "n")
	grid()
}

# QQ-Plot of Standardized Residuals
.plot.garchfilter.9 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	zseries = x@filter$z
	distribution = x@model$modeldesc$distribution
	fd 	= .getDensity(distribution)
	idx = x@model$pidx
	pars  = x@filter$ipars[,1]
	skew  = pars[idx["skew",1]]
	shape = pars[idx["shape",1]]
	if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]
	.qqDist(y = zseries, dist = distribution, lambda = ghlambda, skew = skew, shape = shape)
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
	}
}



# ACF of standardized residuals
.plot.garchfilter.10 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	zseries = x@filter$z
	zseries[is.na(zseries)] = 0
	n 		= length(zseries)
	lag.max = as.integer(10*log10(n))
	acfx 	= acf(zseries, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode = "character", length = lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim=1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Standardized Residuals", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# ACF of squared standardized residuals
.plot.garchfilter.11 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	zseries = x@filter$z
	zseries[is.na(zseries)] = 0
	n 		= length(zseries)
	lag.max = as.integer(10*log10(n))
	acfx 	= acf(zseries^2, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode="character",length=lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)]="steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)]="orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim=1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Squared Standardized Residuals", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# News Impact Curve
.plot.garchfilter.12 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	if(vmodel == "iGARCH"){
		warning(paste("\nplot-->: iGARCH newsimpact not available"))
	} else{
		ni	= newsimpact(z = NULL, x)
		ni.y = ni$zy
		ni.x = ni$zx
		xf	= ni$xexpr
		yf  = ni$yexpr
		plot( ni.x, ni.y, ylab = yf, xlab = xf, type = "l", lwd = 2, col = "steelblue", main = "News Impact Curve", cex.main = 0.8)
		#mtext(yf, side = 2, adj = 0.5, padj = -2.5, cex = 0.85)
		mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
		if(vmodel == "fGARCH"){
			mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
		}
		grid()
	}
}

#-------------------------------------------------------------------------------
# SECTION GARCH sim plots
#-------------------------------------------------------------------------------
.plotgarchsim<-function(x, which="ask", m.sim = 1, ...)
{
	
	choices = c(
			"Conditional SD Simulation Path",
			"Return Series Simulation Path",
			"Conditional SD Simulation Density",
			"Return Series Simulation Path Density")
	.intergarchsimPlot(x,choices=choices, plotFUN = paste(".plot.garchsim", 1:4, sep = "."), which = which, 
			m.sim = m.sim, ...)
	# Return Value:
	invisible(x)
}

.intergarchsimPlot<-function(x, choices, plotFUN, which, m.sim = 1, ...)
{
	if (is.numeric(which)) {
		if(which>length(choices)) stop("Not a valid choice. Plots choices are 1-4.\n", call. = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x, m.sim)
	}
	if(is.character(which))
	{
		if(which!="all" && which!="ask") stop("Not a valid choice.\n", call. = FALSE)
		if (which[1] == "all") {
			par(mfrow=c(2,2))
			for(i in 1:4){
				FUN = match.fun(plotFUN[i])
				FUN(x, m.sim = m.sim)
			}
		}
		if (which[1] == "ask") {
			.multgarchsimPlot(x, choices, m.sim = m.sim, ...)
		}
	}
	invisible(x)
}

.multgarchsimPlot<-function(x, choices, m.sim = 1,  ...)
{
	n = dim(x@simulation$sigmaSim)[2]
	if(m.sim>n | m.sim<=0) stop("\ninvalid m.sim input\n", call. = FALSE)
	pick = 1
	while (pick > 0) {
			pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
				switch (pick,
					.plot.garchsim.1(x, m.sim),  .plot.garchsim.2(x, m.sim),  
					.plot.garchsim.3(x, m.sim),  .plot.garchsim.4(x, m.sim))
	}
}

.plot.garchsim.1<-function(x, m.sim,...)
{
	vmodel  = x@model$modeldesc$vmodel
	xseed = x@seed[m.sim]
	simsigma = matrix(x@simulation$sigmaSim[, m.sim], ncol=1)
	N2 = dim(simsigma)[1]
	T = x@model$modeldata$T
	insample = 1:T
	xdates = x@model$modeldata$dates[insample]
	fwddates = .generatefwd(xdates, N = N2, dformat="%Y-%m-%d", periodicity="days")
	sigma = x@model$modeldata$sigma
	Ns = length(sigma)
	nx = ifelse(Ns>500, 500, min(Ns, 200))
	sigma = c(sigma[(Ns-nx):Ns], rep(NA, N2))
	ylim = c(0.95*min(simsigma, sigma, na.rm=TRUE),1.05*max(simsigma, sigma, na.rm=TRUE))
	simsigma = rbind(matrix(NA, ncol = 1, nrow = nx + 1), simsigma)
	plot(c(xdates[(Ns-nx):Ns], fwddates), sigma, type = "l", col = "steelblue", main = "Simulated Conditional Sigma",
			ylab = "Sigma", xlab = "Time/Horizon", ylim = ylim, cex.main = 0.8)
	abline(h = 0, col = "grey", lty = 3)
	abline(v = nx, col = "red", lty = 3)
	lines(c(xdates[(Ns-nx):Ns], fwddates), as.numeric(simsigma), col = 2)
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
		mtext(paste("seed: ", xseed, sep = ""), side = 3, adj = 0, padj=0, col = "gray", cex = 0.5)
	} else{
		mtext(paste("seed: ", xseed, sep = ""), side = 3, adj = 0, padj=0, col = "gray", cex = 0.5)
	}
	lg.txt = c("Actual", "Simulated")
	legend("topleft", legend = lg.txt, col = c("steelblue", "red"), y.intersp = 1.5, pch = 21, cex = 0.7)
	box()
	grid()
}

.plot.garchsim.2<-function(x, m.sim,...)
{
	vmodel  = x@model$modeldesc$vmodel
	xseed = x@seed[m.sim]
	simseries = matrix(x@simulation$seriesSim[,m.sim],ncol=1)
	N2 = dim(simseries)[1]
	simseries = as.numeric(simseries[,1])
	xdates = x@model$modeldata$dates
	Ns = x@model$modeldata$T
	insample = 1:Ns
	fwddates = .generatefwd(xdates[insample], N = N2, dformat = "%Y-%m-%d", periodicity = "days")
	series = x@model$modeldata$data
	nx = ifelse(Ns>500, 500, min(Ns, 200))
	simseries = c(series[(Ns-nx):Ns], simseries)
	series = c(series[(Ns-nx):Ns], rep(NA, N2))
	ylim = c(0.95*min(simseries,na.rm=TRUE), 1.05*max(simseries,na.rm=TRUE))
	plot(c(xdates[(Ns-nx):Ns],fwddates), simseries, type = "l", col = "red", main = "Simulated Series", 
			ylab = "Series", xlab = "Time/Horizon", ylim = ylim, cex.main = 0.8)
	abline(h = 0, col = "grey", lty = 3)
	abline(v = nx+1, col = "red", lty = 3)
	lines(c(xdates[(Ns-nx):Ns], fwddates),series,col="steelblue")
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
		mtext(paste("seed: ", xseed, sep = ""), side = 3, adj = 0, padj=0, col = "gray", cex = 0.5)
	} else{
		mtext(paste("seed: ", xseed, sep = ""), side = 3, adj = 0, padj=0, col = "gray", cex = 0.5)
	}
	lg.txt = c("Actual", "Simulated")
	legend("topleft", legend = lg.txt, col = c("steelblue", "red"), y.intersp = 1.5, pch = 21, cex = 0.7)
	box()
	grid()
}

.plot.garchsim.3<-function(x, m.sim,...)
{
	vmodel  = x@model$modeldesc$vmodel
	xseed = x@seed[m.sim]
	sigma = x@model$modeldata$sigma
	simsigma = matrix(x@simulation$sigmaSim[,m.sim],ncol=1)
	s1 = density(sigma, kernel="epanechnikov")
	s2 = density(as.numeric(simsigma[,1]), kernel="epanechnikov")
	s1x = s1$x
	s1y = s1$y/length(s1$y)
	s2x = s2$x
	s2y = s2$y/length(s2$y)
	ylim = c(0, max(s1y,s2y))
	plot(s1x, s1y, ylim = ylim, main = "Conditional Sigma\n Kernel Density", type = "l", 
			ylab = "Probability", xlab = "Sigma", cex.main = 0.8)
	lines(s2x, s2y, col = "steelblue")
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
		mtext(paste("seed: ", xseed, sep = ""), side = 3, adj = 0, padj=0, col = "gray", cex = 0.5)
	} else{
		mtext(paste("seed: ", xseed, sep = ""), side = 3, adj = 0, padj=0, col = "gray", cex = 0.5)
	}
	lg.txt = c("Actual", "Simulated")
	legend("topleft", legend = lg.txt, col = c("steelblue", "red"), y.intersp = 1.5, pch = 21, cex = 0.7)
	box()
	grid()
}


.plot.garchsim.4<-function(x, m.sim,...)
{
	vmodel  = x@model$modeldesc$vmodel
	xseed = x@seed[m.sim]
	xseries = x@model$modeldata$data
	T = x@model$modeldata$T
	insample = 1:T
	simseries = matrix(x@simulation$seriesSim[,m.sim],ncol=1)
	s2 = density(as.numeric(simseries[,1]), kernel="epanechnikov")
	n = length(s2$y)
	s1 = density(xseries[insample], kernel="epanechnikov",n=n)
	s1x = s1$x
	s1y = s1$y/length(s1$y)
	s2x = s2$x
	s2y = s2$y/length(s2$y)
	ylim = c(0, max(s1y,s2y))
	plot(s1x, s1y, ylim = ylim, main = "Time Series\n Kernel Density", type = "l", 
			ylab = "Probability", xlab = "Returns", cex.main = 0.8)
	lines(s2x, s2y, col = "steelblue")
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
		mtext(paste("seed: ", xseed, sep = ""), side = 3, adj = 0, padj=0, col = "gray", cex = 0.5)
	} else{
		mtext(paste("seed: ", xseed, sep = ""), side = 3, adj = 0, padj=0, col = "gray", cex = 0.5)
	}
	lg.txt = c("Actual", "Simulated")
	legend("topleft", legend = lg.txt, col = c("steelblue", "red"), y.intersp = 1.5, pch = 21, cex = 0.7)
	box()
	grid()
}

#-------------------------------------------------------------------------------
# SECTION GARCH prediction plots
#-------------------------------------------------------------------------------
.plotgarchforecast<-function(x, which="ask", n.roll = 0, ...)
{
	
	choices = c(
			"Time Series Prediction (unconditional)",
			"Time Series Prediction (rolling)",
			"Conditional SD Prediction")
	.intergarchforecastPlot(x,choices=choices, plotFUN = paste(".plot.garchforecast", 1:3, sep = "."), 
			which = which, n.roll = n.roll, ...)
	# Return Value:
	invisible(x)
}

.intergarchforecastPlot<-function(x, choices, plotFUN, which,  n.roll = 0, ...)
{
	if (is.numeric(which)) {
		if(which>length(choices)) stop("Not a valid choice. Plots choices are 1-3.\n", call. = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x, n.roll)
	}
	if(is.character(which))
	{
		if(which!="all" && which!="ask") stop("Not a valid choice.\n", call. = FALSE)
		if (which[1] == "all") {
			#Which = rep(TRUE, times = length(choices))
			par(mfrow=c(1,3))
			for(i in 1:3){
				FUN = match.fun(plotFUN[i])
				FUN(x, n.roll)
			}
		}
		if (which[1] == "ask") {
			.multgarchforecastPlot(x, choices, n.roll, ...)
		}
	}
	
	invisible(x)
}

.multgarchforecastPlot<-function(x, choices, n.roll = 0, ...)
{
	# A function originally written by Diethelm Wuertz
	
	# Description:
	#   Internal plot function
	
	pick = 1
	while (pick > 0) {
		pick = menu(
				### choices = paste("plot:", choices),
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
		# up to 19 plot functions ...
		switch (pick,.plot.garchforecast.1(x, n.roll),  .plot.garchforecast.2(x), 
				.plot.garchforecast.3(x, n.roll))
	}
}

# Time Series Plot
.plot.garchforecast.1 = function(x, n.roll = 0, ...)
{
	vmodel = x@model$modeldesc$vmodel
	# 1. Time Series:
	nr = x@forecast$n.roll
	if(n.roll > nr) stop("plot-->error: n.roll choice is invalid", call. = FALSE)
	n = x@forecast$n.ahead
	N = x@forecast$N - x@forecast$n.start
	fdata = x@forecast$forecast[[n.roll+1]]
	xdates = x@model$modeldata$dates[(N+n.roll-min(N,100)):(N+n.roll)]
	fwddates = x@forecast$fdates[[n.roll+1]]
	forseries = fdata[,"series"]
	series = x@model$modeldata$data[(N+n.roll-min(N,100)):(N+n.roll)]
	forseries = c(series, forseries)
	#forseries = cumsum(c(1,forseries))[-1]
	#series = cumsum(c(1, series))[-1]
	series = c(series, rep(NA, n))
	ylim=c(0.95*min(forseries,na.rm=TRUE), 1.05*max(forseries,na.rm=TRUE))
	plot(c(xdates, fwddates), forseries , type="l", col="tomato2", 
			main = paste("Forecast Series\n (n.roll = ", n.roll,")", sep = "") ,
			ylab="Series",xlab="Time/Horizon", ylim = ylim, cex.main = 0.7)
	abline(h = 0, col = "grey", lty = 3)
	lines(c(xdates, fwddates), series, , col="steelblue2")
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel=="fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
		mtext(paste("Horizon: ",n,sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
	} else{
		mtext(paste("Horizon: ",n,sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
	}
	lg.txt = c("Actual","Forecast")
	legend("topleft", legend = lg.txt, col = c("steelblue2", "tomato2"), y.intersp = 1.5, pch = 21, cex = 0.7)
	box()
	grid()
}

# rolling forecast plot (actual vs forecast)
.plot.garchforecast.2 = function(x, ...)
{
	vmodel = x@model$modeldesc$vmodel
	# 1. Time Series:
	nr = x@forecast$n.roll
	N = x@forecast$N - x@forecast$n.start
	fdata = sapply(x@forecast$forecast, FUN = function(x) x[1,"series"])
	fsigma = sapply(x@forecast$forecast, FUN = function(x) x[1,"sigma"])
	xdata = x@model$modeldata$data[((N+1)-min(N, 100)):(N+nr)]
	xsigma = x@model$modeldata$sigma[((N+1)-min(N, 100)):(N+nr)]
	xdates = x@model$modeldata$dates[((N+1)-min(N, 100)):(N+nr)]
	ns = length(xdata)
	xplus =  xdata + 2*xsigma
	xminus = xdata - 2*xsigma
	fplus = c(rep(NA, (ns - nr)), fdata[-length(fdata)] + 2*fsigma[-length(fsigma)])
	fminus = c(rep(NA, (ns - nr)), fdata[-length(fdata)] - 2*fsigma[-length(fsigma)])
	fdata = c(rep(NA, (ns - nr)), fdata[-length(fdata)])
	
	ylim=c(0.95*min(xminus,na.rm=TRUE), 1.2*max(xplus,na.rm=TRUE))
	plot(xdates, xdata, type="l", col="black", 
			main = paste("Rolling Forecast vs Actual Series\n w/th 2 conditional SD bands", sep = "") ,
			ylab="Series",xlab="Time/Horizon", ylim = ylim, cex.main = 0.7)
	lines(xdates, xplus, col = "lightgrey", lwd = 0.5)
	lines(xdates, xminus, col = "lightgrey", lwd = 0.5)
	abline(h = 0, col = "grey", lty = 3)
	lines(xdates, fdata, col = "tomato2", lwd = 2.5)
	lines(xdates, fplus, col = "brown", lwd = 0.5)
	lines(xdates, fminus, col = "brown", lwd = 0.5)
	
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel=="fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
		mtext(paste("Horizon: ", nr, sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
	} else{
		mtext(paste("Horizon: ", nr, sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
	}
	lg.txt = c("Actual","Forecast")
	legend("topleft", legend = lg.txt, col = c("black", "tomato2"), y.intersp = 1.5, pch = 21, cex = 0.7)
	box()
	grid()
}

.plot.garchforecast.3 = function(x, n.roll = 0, ...)
{
	vmodel = x@model$modeldesc$vmodel
	nr = x@forecast$n.roll
	if(n.roll > nr) stop("plot-->error: n.roll choice is invalid", call. = FALSE)
	n = x@forecast$n.ahead
	N = x@forecast$N - x@forecast$n.start
	fdata = x@forecast$forecast[[n.roll+1]]
	xdates = x@model$modeldata$dates[(N+n.roll-min(N,100)):(N+n.roll)]
	fwddates = x@forecast$fdates[[n.roll+1]]
	forsigma = fdata[,"sigma"]
	sigma = x@model$modeldata$sigma[(N+n.roll-min(N,100)):(N+n.roll)]
	forsigma = c(sigma,forsigma)
	sigma = c(sigma, rep(NA,n))
	ylim=c(0.95*min(forsigma,na.rm=TRUE),1.05*max(forsigma,na.rm=TRUE))
	plot(c(xdates,fwddates), forsigma, type = "l", col = "tomato2", 
			main = paste("Forecast Conditional Sigma\n (n.roll = ", n.roll,")", sep = "") , 
			ylab = "Sigma", xlab = "Time/Horizon", ylim = ylim, cex.main = 0.7)
	abline(h = 0, col = "grey", lty = 3)
	lines(c(xdates,fwddates), sigma, col = "steelblue")
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
		mtext(paste("Horizon: ",n,sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
	} else{
		mtext(paste("Horizon: ",n,sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
	}
	lg.txt<-c("Actual","Forecast")
	legend("topleft", legend=lg.txt,col=c("steelblue2","tomato2"), y.intersp=1.5,pch=21,cex=0.7)
	box()
	grid()
}

#-------------------------------------------------------------------------------
# SECTION GARCH path plots
#-------------------------------------------------------------------------------
.plotgarchpath = function(x, which="ask", m.sim = 1, ...)
{
	choices = c(
			"Conditional SD Simulation Path",
			"Return Series Simulation Path",
			"Conditional SD Simulation Density",
			"Return Series Simulation Path Density")
	.intergarchpathPlot(x, choices = choices, plotFUN = paste(".plot.garchpath", 1:4, sep = "."), 
			which = which, m.sim = m.sim,...)
	# Return Value:
	invisible(x)
}

.intergarchpathPlot = function(x, choices, plotFUN, which, m.sim = 1, ...)
{
	if (is.numeric(which)) {
		if(which>length(choices)) stop("Not a valid choice. Plots choices are 1-4.\n", call. = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x, m.sim = m.sim)
	}
	if(is.character(which))
	{
		if(which!="all" && which!="ask") stop("Not a valid choice.\n", call. = FALSE)
		if (which[1] == "all") {
			par(mfrow=c(2,2))
			for(i in 1:4){
				FUN = match.fun(plotFUN[i])
				FUN(x, m.sim = m.sim)
			}
		}
		if (which[1] == "ask") {
			.multgarchpathPlot(x, choices, m.sim = m.sim, ...)
		}
	}
	invisible(x)
}

.multgarchpathPlot<-function(x, choices, m.sim = 1, ...)
{
	n = dim(x@path$sigmaSim)[2]
	if(m.sim>n | m.sim<=0) stop("\ninvalid m.sim input\n", call. = FALSE)
	pick = 1
	npick = 0
	while (pick > 0) {
		pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
			switch (pick,
					.plot.garchpath.1(x, m.sim), .plot.garchpath.2(x, m.sim),  
					.plot.garchpath.3(x, m.sim), .plot.garchpath.4(x, m.sim))
	}
}

.plot.garchpath.1 = function(x, m.sim = 1, ...)
{
	model = strsplit(class(x),"path")
	xseed = x@seed[m.sim]
	simsigma = matrix(x@path$sigmaSim[,m.sim],ncol=1)
	plot(simsigma, type = "l", col = "steelblue", main = "Simulated Path of Conditional Sigma", 
			ylab="Sigma", xlab="Time/Horizon", cex.main = 0.7)
	mtext(paste("GARCH model :",model), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	box()
	grid()
}

.plot.garchpath.2 = function(x, m.sim = 1, ...)
{
	model = strsplit(class(x),"sim")
	xseed = x@seed[m.sim]
	simseries = matrix(x@path$seriesSim[,m.sim],ncol=1)
	simseries = as.numeric(simseries[,1])
	simseries = simseries
	plot(simseries,type="l",col="red", main = "Simulated Path of Series", ylab = "Series", 
			xlab="Time/Horizon", cex.main = 0.7)
	mtext(paste("GARCH model :",model), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	box()
	grid()
}

.plot.garchpath.3 = function(x, m.sim = 1,...)
{
	model = strsplit(class(x),"sim")
	xseed = x@seed[m.sim]
	simsigma = matrix(x@path$sigmaSim[,m.sim], ncol = 1)
	s2 = density(as.numeric(simsigma[,1]), kernel = "epanechnikov")
	plot(s2, main="Simulated Conditional Sigma\n Kernel Density",type="l",ylab="Probability", xlab="Sigma",
			col="steelblue", cex.main = 0.7)
	mtext(paste("GARCH model :",model), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	box()
	grid()
}


.plot.garchpath.4 = function(x, m.sim = 1,...)
{
	model = strsplit(class(x),"sim")
	xseed = x@seed[m.sim]
	simseries = matrix(x@path$seriesSim[,m.sim], ncol = 1)
	s2 = density(as.numeric(simseries[,1]), kernel = "epanechnikov")
	plot(s2, main = "Simulated Conditional Time Series\n Kernel Density",type="l",ylab="Probability", xlab="Returns",
			col="red", cex.main = 0.7)
	mtext(paste("GARCH model :",model), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	box()
	grid()
}

#-------------------------------------------------------------------------------
# SECTION GARCH roll plots
#-------------------------------------------------------------------------------
.plotgarchroll = function(x, which="ask", n.ahead = 1, VaR.alpha = 0.01, 
		density.support=c(-0.15, 0.15), ...)
{
	
	choices = c(
			"Density Forecast",
			"Sigma Forecast",
			"Series Forecast",
			"VaR Forecast",
			"Fit Coefficients (with s.e. bands)")
	.intergarchrollPlot(x, choices=choices, plotFUN = paste(".plot.garchroll", 1:5, sep = "."), 
			which = which, n.ahead = n.ahead, VaR.alpha = VaR.alpha, 
			density.support = density.support, ...)
	# Return Value:
	invisible(x)
}

.intergarchrollPlot = function(x, choices, plotFUN, which, n.ahead = 1, VaR.alpha = 0.01, 
		density.support=c(-0.15, 0.15), ...)
{
	if (is.numeric(which)) {
		if(which>length(choices)) stop("Not a valid choice. Plots choices are 1-5.\n", call. = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x, n.ahead = n.ahead, VaR.alpha = VaR.alpha, density.support = density.support, ...)
	}
	if(is.character(which))
	{
		if(which!="all" && which!="ask") stop("Not a valid choice.\n", call. = FALSE)
		if (which[1] == "all") {
			par(mfrow=c(2,2))
			.plot.garchroll.1(x, n.ahead = n.ahead, density.support = density.support, ...)
			.plot.garchroll.2(x, n.ahead = n.ahead, ...)
			.plot.garchroll.3(x, n.ahead = n.ahead, ...)
			.plot.garchroll.4(x, n.ahead = n.ahead, VaR.alpha = VaR.alpha, ...)
		}
		if (which[1] == "ask") {
			.multgarchrollPlot(x, choices, n.ahead, VaR.alpha, 
					density.support,...)
		}
	}
	invisible(x)
}

.multgarchrollPlot = function(x, choices, n.ahead = 1, VaR.alpha = 0.01, 
		density.support = c(-0.15, 0.15), ...)
{
	pick = 1
	while (pick > 0) {
		pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
		switch (pick,
				.plot.garchroll.1(x, n.ahead, VaR.alpha, density.support,...),  
				.plot.garchroll.2(x, n.ahead, VaR.alpha, density.support,...),  
				.plot.garchroll.3(x, n.ahead, VaR.alpha, density.support,...),
				.plot.garchroll.4(x, n.ahead, VaR.alpha, density.support,...),
				.plot.garchroll.5(x, n.ahead, VaR.alpha, density.support,...))
	}
}

# rolling sigma forecast comparison plot
.plot.garchroll.1 = function(x,  n.ahead = 1, VaR.alpha = 0.01, 
		density.support = c(-0.15, 0.15), ...)
{
	n = x@roll$n.ahead
	fdensity = x@roll$fdensity
	distribution = x@model$modeldesc$distribution
	forecast.length = x@roll$forecast.length
	# we plot a maximum of 500 forecasts
	esd = floor(seq(1,forecast.length, length.out = min(500, forecast.length)))
	nesd = length(esd)
	colr = topo.colors(nesd, alpha = 1)
	if(is.vector(fdensity)){
		xden = fdensity[[n.ahead]][,-1]
		xdate =  as.Date(as.character(fdensity[[n.ahead]][,"fdate"]), format = "%Y-%m-%d")
		if(is.na(xdate[1])) xdate = as.character(fdensity[[n.ahead]][,"fdate"])
	} else{
		xden = fdensity[,-1]
		xdate =  as.Date(as.character(fdensity[,"fdate"]), format = "%Y-%m-%d")
		if(is.na(xdate[1])) xdate = as.character(fdensity[,"fdate"])
	}
	xseq = seq(density.support[1], density.support[2], length.out=1000)
	yseq = apply(as.data.frame(2:nesd), 1, FUN=function(i)  
				.ddensity(xseq, mu = xden[esd[i],"fmu"], sigma = xden[esd[i],"fsigma"], 
						lambda = xden[esd[i],"fdlambda"], skew = xden[esd[i],"fskew"], shape = xden[esd[i],"fshape"], 
						distribution = distribution))
	plot(xseq, .ddensity(xseq, mu = xden[1,"fmu"], sigma = xden[1,"fsigma"], 
					lambda = xden[1,"fdlambda"], skew = xden[1,"fskew"], shape = xden[1,"fshape"], 
					distribution = distribution), type="l", col = "steelblue", 
			main = paste("n.ahead-", n.ahead," Forecast Density (time varying)",sep=""), ylab="", cex.main = 0.7)
	for(i in 1:(nesd-1)){
		lines(xseq, yseq[,i], col = colr[i+1])
	}
	xesd = floor(seq(1,length(esd), length.out = 5))
	legend("topright", legend = as.character(xdate[esd[xesd]]), fill = colr[xesd], col = colr[xesd], bty = "n")
	invisible(x)
}

# rolling sigma forecast comparison plot
.plot.garchroll.2= function(x,  n.ahead = 1, VaR.alpha = 0.01, 
		density.support = c(-0.15, 0.15), ...)
{
	n = x@roll$n.ahead	
	fdensity = x@roll$f01density
	if(n.ahead > n) stop("\nplot-->error: n.ahead chosen is not valid for object\n", call.=FALSE)
	if(is.vector(fdensity)){
		xsig = fdensity[[n.ahead]][,"f01sigma"]
		xdate = as.Date(as.character(fdensity[[n.ahead]][,"f01date"]), format = "%Y-%m-%d")
		if(is.na(xdate[1])) xdate = as.character(fdensity[[n.ahead]][,"f01date"])
		
	} else{
		xsig = fdensity[,"f01sigma"]
		xdate =  as.Date(as.character(fdensity[,"f01date"]), format = "%Y-%m-%d")
		if(is.na(xdate[1])) xdate =as.character(fdensity[,"f01date"])
		
	}
	plot(xdate, xsig, type="l", col = "steelblue", main = paste("n.ahead-", n.ahead," Sigma Forecast", sep = ""), 
			ylab = "forecast sigma", cex.main = 0.7)
	invisible(x)
}

# rolling sigma forecast comparison plot
.plot.garchroll.3 = function(x,  n.ahead = 1, VaR.alpha = 0.01, 
		density.support = c(-0.15, 0.15), ...)
{
	n = x@roll$n.ahead	
	fdensity = x@roll$f01density
	if(n.ahead > n) stop("\nplot-->error: n.ahead chosen is not valid for object\n", call. = FALSE)
	if(is.vector(fdensity)){
		xmu = fdensity[[n.ahead]][,"f01mu"]
		xdate =  as.Date(as.character(fdensity[[n.ahead]][,"f01date"]), format = "%Y-%m-%d")
		if(is.na(xdate[1])) xdate = as.character(fdensity[[n.ahead]][,"f01date"])
	} else{
		xmu = fdensity[,"f01mu"]
		xdate =  as.Date(as.character(fdensity[,"f01date"]), format = "%Y-%m-%d")
		if(is.na(xdate[1])) xdate =as.character(fdensity[,"f01date"])
	}
	actseries = x@model$modeldata$filterseries
	plot(xdate, xmu, type="l", col = "steelblue", main = paste("n.ahead-",n.ahead," Series Forecast",sep=""), ylab="return", ylim=c(min(min(actseries),min(xmu)), max(max(actseries),max(xmu))), cex.main = 0.7)
	points(xdate, actseries, cex=0.5, col=2)
	legend("topleft", legend=c("Forecast", "Actual"), fill = c("steelblue", "red"), col = c("steelblue", "red"), bty = "n" )
	invisible(x)
}

# rolling VaR backtest plot
.plot.garchroll.4 = function(x,  n.ahead = 1, VaR.alpha = 0.01, 
		density.support = c(-0.15, 0.15), ...)
{
	n = x@roll$n.ahead
	v.a = x@roll$VaR.alpha
	
	if(is.null(x@roll$VaR.out)) stop("\nplot-->error: VaR was not calculated for this object\n", call.=FALSE)
	if(n.ahead > n) stop("\nplot-->error: n.ahead chosen is not valid for object\n", call.=FALSE)
	if(!is.null(v.a) && !any(v.a==VaR.alpha[1])) stop("\nplot-->error: VaR.alpha chosen is invalid for the object\n", call.=FALSE)
	if(is.list(x@roll$VaR.out)){
		dvar = x@roll$VaR.out[[n.ahead]]
		dates = rownames(dvar)
		m = dim(dvar)[2]
		idx = which(colnames(dvar) == paste("alpha(", round(VaR.alpha,2)*100, "%)",sep=""))
		.VaRplot(x@roll$dataname, p = VaR.alpha, actual = dvar[, m], dates = dates, VaR = dvar[, idx])		
	} else{
		dvar = x@roll$VaR.out
		dates = rownames(dvar)
		m = dim(dvar)[2]
		idx = which(colnames(dvar) == paste("alpha(", round(VaR.alpha,2)*100, "%)",sep=""))
		.VaRplot(x@roll$dataname, p = VaR.alpha, actual = dvar[, m], dates = dates, VaR = dvar[, idx])		
	}
	invisible(x)	
}

.plot.garchroll.5 = function(x,  n.ahead = 1, VaR.alpha = 0.01, 
		density.support = c(-0.15, 0.15), ...)
{
	# get the no. of coef of a fit
	coefs = x@roll$coefs
	coefmat = x@roll$coefmat
	N = dim(coefs)[1]
	if(N<2){
		stop("\n 1 refit does not a coef evolution plot make!\n", call.=FALSE)
	}
	n = dim(coefs)[2]
	cnames = colnames(coefs)
	np = .divisortable(n)
	par(mfrow = c(np[1], np[2]))
	for(i in 1:n){
	cse = sapply(coefmat, FUN = function(x) x[i,2], simplify = TRUE)
	plot(coefs[,i], type="l", ylim = c(min(coefs[,i]-cse), max(coefs[,i] + cse)), 
			ylab = "value" , xlab = "", main = "")
	lines(coefs[,i] + cse, col=2)
	lines(coefs[,i] - cse, col=2)
	title(cnames[i], line = 0.4, cex = 0.9)
	}
	title(main = list(paste("uGARCH fit coefficients (across ",N," refits) with s.e. bands",sep=""), 
					cex = 1.5, col = "steelblue", font=2), outer = TRUE, line = -1.4)
	invisible(x)
}



#-------------------------------------------------------------------------------
# SECTION GARCH Distribution (Parameter Uncertainty) plots
#-------------------------------------------------------------------------------
.plotgarchdist = function(x, which="ask", window = 1, ...)
{
	choices = c(
			"Parameter Density Plots",
			"Bivariate Plots",
			"Other Density Plots (Persistence, Likelihood, ...)",
			"Asymptotic Efficiency Plots")
	
	.intergarchdistPlot(x, choices = choices, plotFUN = paste(".plot.garchdist", 1:4, sep = "."), 
			which = which, window = window, ...)
	# Return Value:
	invisible(x)
}

.intergarchdistPlot = function(x, choices, plotFUN, which, window, ...)
{
	if (is.numeric(which)) {
		if(which>length(choices)) 
			stop("Not a valid choice. Plots choices are 1-4.\n", call. = FALSE)
		if(which == 4 && !x@dist$detail$recursive)
			stop("Asymptotic Efficiency Plots only available for option recursive TRUE", .call = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x, window, ...)
	}
	if(is.character(which))
	{
		if(which!="ask") stop("Not a valid choice.\n", call. = FALSE)
		if (which[1] == "ask") {
			.multgarchdistPlot(x, choices, window, ...)
		}
	}
	invisible(x)
}

.multgarchdistPlot = function(x, choices, window, ...)
{
	pick = 1
	while (pick > 0) {
		pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
		switch (pick,
				.plot.garchdist.1(x, window, ...),  
				.plot.garchdist.2(x, window, ...),  
				.plot.garchdist.3(x, ...),
				.plot.garchdist.4(x, ...))
	}
}

# Parameter Density Plots
.plot.garchdist.1 = function(x, window = 1,  ...)
{
	cf = as.data.frame(x, window = window)
	n = dim(cf)[2]
	vmodel = x@model$modeldesc$vmodel
	vsubmodel = x@model$modeldesc$vsubmodel
	truecf = unlist(x@truecoef[,1])
	cnames = names(truecf)
	z = .divisortable(n)
	par(mfrow=c(z[1], z[2]))
	for(i in 1:n){
		str = paste("Parameter ", cnames[i],"\nTrue Value: ", signif(truecf[i], digits = 3))
		plot(density(cf[,i], kernel = "gaussian", na.rm = TRUE), col = "steelblue4", main = str, cex.main = 0.7)
		mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
				cex = 0.5)
		if(vmodel == "fGARCH"){
			mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
					padj = 1.5, col = "gray", cex = 0.5)
		}
	}
}

# Bivariate Plots
.plot.garchdist.2= function(x, window = 1,...)
{
	.dist2dplot(x, window = window, ...)
	invisible()
}

# stat plots
.plot.garchdist.3 = function(x, ...)
{
	.statsplot(x, ...)
}

# rmse backtest plot
.plot.garchdist.4 = function(x,  ...)
{
	if(!x@dist$details$recursive){
		print("no rmse plots for non recursive option")
		return()
	} else{
		.rmseplots(x, ...)
	}
}



#-------------------------------------------------------------------------------
# SECTION GARCH boot plots
#-------------------------------------------------------------------------------
.plotgarchboot = function(x, which="ask", ...)
{
	if(x@model$type!="full"){
		choices = c(
			"Series Standard Error Plots",
			"Sigma Standard Error Plots")
		nc = 2:3
	} else{
		choices = c(
				"Parameter Density Plots",
				"Series Standard Error Plots",
				"Sigma Standard Error Plots")
		nc = 1:3
	}
	.intergarchbootPlot(x, choices = choices, plotFUN = paste(".plot.garchboot", nc, sep = "."), 
			which = which, window = window, ...)
	# Return Value:
	invisible(x)
}

.intergarchbootPlot = function(x, choices, plotFUN, which,  ...)
{
	if (is.numeric(which)) {
		if(which>length(choices)) 
			stop("Not a valid choice. Plots choices are 1-3.\n", call. = FALSE)
		if(which == 1 && x@model$type!="full")
			stop("Parameter Density Plots only available for full bootstrap method", .call = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x, window, ...)
	}
	if(is.character(which))
	{
		if(which!="ask" & which!="all") stop("Not a valid choice.\n", call. = FALSE)
		if (which[1] == "ask") {
			.multgarchbootPlot(x, choices,  ...)
		} else{
				par(mfrow = c(2,1))
				.plot.garchboot.2(x, ...)
				.plot.garchboot.3(x, ...)
		}
	}
	invisible(x)
}

.multgarchbootPlot = function(x, choices,  ...)
{
	pick = 1
	if(x@model$type!="full"){
		while (pick > 0) {
					pick = menu (
							choices = paste(" ", choices),
							title = "\nMake a plot selection (or 0 to exit):")
					switch (pick,
							.plot.garchboot.2(x, ...),  
							.plot.garchboot.3(x, ...))
				}
	} else{
	while (pick > 0) {
		pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
		switch (pick,
				.plot.garchboot.1(x, ...),  
				.plot.garchboot.2(x, ...),  
				.plot.garchboot.3(x, ...))
	}}
}

# Parameter Density Plots (for full method)
.plot.garchboot.1 = function(x,  ...)
{	
	cf = x@bcoef
	vmodel = x@model$modeldesc$vmodel
	vsubmodel = x@model$modeldesc$vsubmodel
	truecf = x@model$truecoef
	n = length(truecf)
	cnames = names(truecf)
	z = .divisortable(n)
	par(mfrow=c(z[1], z[2]))
	for(i in 1:n){
		str = paste("Parameter ", cnames[i],"\nTrue Value: ", signif(truecf[i], digits = 3))
		plot(density(cf[,i], kernel = "gaussian", na.rm = TRUE), col = "steelblue4", main = str,
				cex.main = 0.85)
		mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
				cex = 0.5)
		if(vmodel == "fGARCH"){
			mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
					padj = 1.5, col = "gray", cex = 0.5)
		}
	}
}

# Series Error Plots
.plot.garchboot.2 = function(x, ...)
{
	vmodel = x@model$modeldesc$vmodel
	vsubmodel = x@model$modeldesc$vsubmodel
	n.ahead = x@model$n.ahead
	#fs = rep(NA, n.ahead)
	#rx = rep(NA, n.ahead)
	#if(!is.null(x@model$realized.x)){
	#	n = length(x@model$realized.x)
	#	if(n.ahead>n) rx[1:n] = x@model$realized.x else rx = x@model$realized.x[1:n.ahead]
	#}
	seriesfor = x@forc@forecast$forecast[[1]][,"series"]
	xser = as.data.frame(x, which = "series")
	N = dim(xser)[1]
	serp = as.data.frame(x, which = "series", type = "q", qtile = c(0.05, 0.25, 0.75, 0.95))
	miny = min(serp[1,])
	maxy = max(serp[4,])
	meanser = apply(xser, 2, FUN = function(x) mean(x))
	plot(seriesfor, type = "l", col = "red", ylim = c(miny, maxy), main = "Series Forecast
					with Bootstrap Error Bands\n (q: 5%, 25%, 75%, 95%)", cex.main = 0.7,
			ylab = "returns", xlab = "n.ahead")
	lines(as.numeric(meanser), col = "black")
	#lines(as.numeric(rx), col = "green")
	points(as.numeric(serp[1,]), col = "steelblue1", pch = 19, cex = 0.5)
	points(as.numeric(serp[2,]), col = "steelblue2", pch = 19, cex = 0.5)
	points(as.numeric(serp[3,]), col = "steelblue3", pch = 19, cex = 0.5)
	points(as.numeric(serp[4,]), col = "steelblue4", pch = 19, cex = 0.5)
	n.overlays = round(0.1*n.ahead)
	if(n.overlays == 0) n.overlays = 1
	dens1 = apply(xser, 2, FUN = function(x) density(x, kernel = "gaussian", 
						n = max(512, round(0.01*N))))
	
	.densityoverlay(as.numeric(meanser), dens1, n.overlays = n.overlays)
	mtext(paste("GARCH model :",vmodel), side = 4, adj = 0, padj=0, col = "gray", 
			cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
				padj = 1.5, col = "gray", cex = 0.5)
	}
	legend.txt = c("forecast", "bootstrapped")
	legend("bottomleft", legend = legend.txt, fill = c("red", "black"), col = c("red", "black"),
			bg = "n", bty = "n", cex = 0.6)
}

# Sigma Error Plots
.plot.garchboot.3 = function(x, ...)
{
	vmodel = x@model$modeldesc$vmodel
	vsubmodel = x@model$modeldesc$vsubmodel
	n.ahead = x@model$n.ahead
	fs = rep(NA, n.ahead)
	colr = NULL
	namr = NULL
	if(!is.null(x@model$modeldata$filtered.s)){
		n = length(x@model$modeldata$filtered.s)
		if(n.ahead>n) fs[1:n] = x@model$modeldata$filtered.s else fs = x@model$modeldata$filtered.s[1:n.ahead]
		colr = "green"
		namr = "filtered"
	}
	sigmafor = x@forc@forecast$forecast[[1]][,"sigma"]
	sigdist = vector(mode = "list", length = n.ahead)
	sigp = as.data.frame(x, which = "sigma", type = "q", qtile = c(0.05, 0.25, 0.75, 0.95))
	miny = min(sigp[1,])
	maxy = max(sigp[4,])
	meansig = apply(x@fsigma, 2, FUN = function(x) mean(x))
	plot(sigmafor, type = "l", col = "red", ylim = c(miny, maxy), main = "Sigma Forecast
					with Bootstrap Error Bands\n (q: 5%, 25%, 75%, 95%)", cex.main = 0.7,
			ylab = "sigma", xlab = "n.ahead")
	lines(as.numeric(meansig), col = "black")
	lines(as.numeric(fs), col = colr)
	points(as.numeric(sigp[1,]), col = "steelblue1", pch = 19, cex = 0.5)
	points(as.numeric(sigp[2,]), col = "steelblue2", pch = 19, cex = 0.5)
	points(as.numeric(sigp[3,]), col = "steelblue3", pch = 19, cex = 0.5)
	points(as.numeric(sigp[4,]), col = "steelblue4", pch = 19, cex = 0.5)
	mtext(paste("GARCH model :",vmodel), side = 4, adj = 0, padj=0, col = "gray", 
			cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
				padj = 1.5, col = "gray", cex = 0.5)
	}
	legend.txt = c("forecast", "bootstrapped", namr)
	legend("bottomleft", legend = legend.txt, fill = c("red", "black", colr), col = c("red", "black", colr),
			bg = "n", bty = "n", cex = 0.6)
}