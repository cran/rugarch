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

# BDS Test
# Dechert, Scheinkman and LeBaron (1996)
# BDS i.i.d test of 
# Brock and Potter (1993) and de Lima (1996) show that is
# applying the BDS test to the log of the squared standardized residuals
# log(z^2) the bias in the BDS test is almost corrected because the logarithmic
# transformation transforms the GARCH model into a linear additive model which
# satisfies the nuisance free parameter conditions in de Lima (1996) for the BDS
# test (does not apply to more richly parametrized models e.g. with leverage since
# they cannot be cast into linear additive models).

#.bds.test = function(z){
#	fNonlinear::bdsTest(log(z[is.finite(z)]^2), m = 5)
#}

.information.test = function(LLH, nObs, nPars)
{
	AIC  = (-2*LLH)/nObs + 2 * nPars/nObs
	BIC  = (-2*LLH)/nObs + nPars * log(nObs)/nObs
	SIC  = (-2*LLH)/nObs + log((nObs+2*nPars)/nObs)
	HQIC = (-2*LLH)/nObs + (2*nPars*log(log(nObs)))/nObs
	informationTests = list(AIC = AIC, BIC = BIC, SIC = SIC, HQIC = HQIC)
	return(informationTests)
}

# Q-Statistics on Standardized Residuals
.box.test = function(stdresid, p=1, df = 0)
{
	if(any(!is.finite(stdresid))) stdresid[!is.finite(stdresid)]=0
	# p=1 normal case, p=2 squared std. residuals
	# Q-Statistics on Standardized Residuals
	#H0 : No serial correlation ==> Accept H0 when prob. is High [Q < Chisq(lag)]
	box10 = Box.test(stdresid^p, lag = 10, type = "Ljung-Box", fitdf = df)
	box15 = Box.test(stdresid^p, lag = 15, type = "Ljung-Box", fitdf = df)
	box20 = Box.test(stdresid^p, lag = 20, type = "Ljung-Box", fitdf = df)
	LBSR<-matrix(NA,ncol=2,nrow=3)
	LBSR[1:3,1] = c(box10$statistic[[1]],box15$statistic[[1]],box20$statistic[[1]])
	LBSR[1:3,2] = c(box10$p.value[[1]],box15$p.value[[1]],box20$p.value[[1]])
	rownames(LBSR) = c("Lag10","Lag15","Lag20")
	#LBSR<-as.data.frame(cbind(LBSR,.stars(as.vector(1-LBSR[,2]))))
	colnames(LBSR) = c("statistic","p-value")
	return(LBSR)
}


.archlmtest = function (x, lags, demean = FALSE)
{
	if(any(!is.finite(x))) x[!is.finite(x)] = 0
	x = as.vector(x)
	if(demean) x = scale(x, center = TRUE, scale = FALSE)
	lags = lags + 1
	mat = embed(x^2, lags)
	arch.lm = summary(lm(mat[, 1] ~ mat[, -1]))
	STATISTIC = arch.lm$r.squared * length(resid(arch.lm))
	names(STATISTIC) = "Chi-squared"
	PARAMETER = lags - 1
	names(PARAMETER) = "df"
	PVAL = 1 - pchisq(STATISTIC, df = PARAMETER)
	METHOD = "ARCH LM-test"
	result = list(statistic = STATISTIC, parameter = PARAMETER,
			p.value = PVAL, method = METHOD)
	class(result) = "htest"
	return(result)
}

.nyblomTest = function(object)
{
	#pnames = rownames(object@fit$ipars[object@fit$ipars[,"Estimate"]==1,])
	grad = object@fit$scores
	pnames = colnames(grad)
	# special case when fixed parameters exist but the fixed.se was used in fit.control options
	#if( length(pnames) != dim(grad)[2] && dim(grad)[2] == length(object@fit$ipars[object@fit$ipars[,"Include"]==1,]) ){
	#	pnames = rownames(object@fit$ipars[object@fit$ipars[,"Include"]==1,])
	#}
	if(is(object, "ARFIMAfit")) res = object@fit$residuals/object@model$pars["sigma", 1] else res = object@fit$residuals/object@fit$sigma
	res[!is.finite(res) | is.nan(res) | is.na(res)] = 0
	nn = length(res)
	hes = t(grad)%*%(grad)
	shes = try(solve(hes), silent = TRUE)
	if(inherits(shes,"try-error")){
		shes = try( solve(object@fit$hessian) )
		if(inherits(shes,"try-error")){
			IndividualStat = matrix(rep(NA, length(pnames)), ncol=1)
			IndividualCritical = .nyblomCritical(1)
			JointCritical = .nyblomCritical(length(pnames))
			rownames(IndividualStat) = pnames
			JointStat = NA
		} else{
			zx = matrix(cbind(res, res^2, res^3, res^4), ncol = 4)
			zs = apply(zx, 2, FUN = function(x) scale(x, center = TRUE, scale = FALSE))
			x = as.matrix(apply(grad, 2, FUN = function(x) cumsum(x)))
			xx = t(x)%*%x
			nyblomj = sum(diag(xx%*%shes))/nn
			nyblomt = diag(xx)/(diag(hes)*nn)
			IndividualStat = matrix(nyblomt, ncol=1)
			IndividualCritical = .nyblomCritical(1)
			JointCritical = .nyblomCritical(length(pnames))
			rownames(IndividualStat) = pnames
			JointStat = nyblomj
		}
	} else{
		zx = matrix(cbind(res, res^2, res^3, res^4), ncol = 4)
		zs = apply(zx, 2, FUN = function(x) scale(x, center = TRUE, scale = FALSE))
		x = as.matrix(apply(grad, 2, FUN = function(x) cumsum(x)))
		xx = t(x)%*%x
		nyblomj = sum(diag(xx%*%shes))/nn
		nyblomt = diag(xx)/(diag(hes)*nn)
		IndividualStat = matrix(nyblomt, ncol=1)
		IndividualCritical = .nyblomCritical(1)
		JointCritical = .nyblomCritical(length(pnames))
		rownames(IndividualStat) = pnames
		JointStat = nyblomj
	}
	return(list(IndividualStat=IndividualStat,JointStat=JointStat,
					IndividualCritical=IndividualCritical,JointCritical=JointCritical))
}
.nyblomCritical = function(n){
	# Test for Constancy of Parameters: Hansen Tests (1990 mimeo paper)
	# Null Hypothesis: constant parameters
	cval = matrix(c(0.353, 0.470, 0.748,
					0.610, 0.749, 1.07,
					0.846, 1.01, 1.35,
					1.07,  1.24, 1.60,
					1.28,  1.47, 1.88,
					1.49,  1.68, 2.12,
					1.69,  1.90, 2.35,
					1.89,  2.11, 2.59,
					2.10,  2.32, 2.82,
					2.29,  2.54, 3.05,
					2.49,  2.75, 3.27,
					2.69,  2.96, 3.51,
					2.89,  3.15, 3.69,
					3.08,  3.34, 3.90,
					3.26,  3.54, 4.07,
					3.46,  3.75, 4.30,
					3.64,  3.95, 4.51,
					3.83,  4.14, 4.73,
					4.03,  4.33, 4.92,
					4.22,  4.52, 5.13), byrow=T, ncol=3)
	colnames(cval)=c("10%","5%","1%")
	if(n<=20){
		ans = cval[n,]
		names(ans)=c("10%","5%","1%")
	} else {
		ans="too many parameters"
	}
	return(ans)
}

.signbiasTest = function(object)
{
	if(is(object, "uGARCHfilter")) z = object@filter$z else z = z = object@fit$z
	res = residuals(object)
	z2 = z^2
	n = length(z)
	zminus = as.integer(z<0)
	zplus = 1-zminus
	czminus = zminus*res
	czplus = zplus*res
	cz = cbind(rep(1, n), zminus, czminus, czplus)
	cz = cz[1:(n-1),]
	z2 = matrix(z2[2:n],ncol = 1)
	cm = data.frame(y = z2, const = cz[,1], zminus = cz[,2], czminus = cz[,3], czplus = cz[,4])
	fitA = lm(y~const+zminus+czminus+czplus-1, data = cm)
	resA = residuals(fitA)
	sgmatrix = matrix(ncol = 3, nrow = 4)
	rownames(sgmatrix) = c("Sign Bias","Negative Sign Bias","Positive Sign Bias","Joint Effect")
	colnames(sgmatrix) = c("t-value","prob","sig")
	sgmatrix[1:3,1:2] = abs(summary(fitA)$coef[2:4,3:4])
	jeffect = .linear.hypothesis(fitA, c("zminus = 0", "czminus = 0","czplus  =  0"), test  =  "Chisq")
	sgmatrix[4,1] = jeffect[5]$Chisq[2]
	sgmatrix[4,2] = jeffect[6][2,]
	sgmatrix = as.data.frame(sgmatrix)
	sgmatrix[,3] = .stars(sgmatrix[,2])
	return(sgmatrix)
}

.gofTest = function(object, groups)
{
	modelinc = object@model$modelinc
	idx = object@model$pidx
	if(is(object, "uGARCHfilter")){
		ipars = object@filter$ipars
		z = object@filter$z
	} else{
		ipars = object@fit$ipars
		z = object@fit$z
	}
	# must remove fixed parameters
	dist = object@model$modeldesc$distribution
	cdf = .getDistribution(dist)
	cdfv = cdf(z = sort(z), hh = 0, lambda = ipars[idx["ghlambda",1],1], skew = ipars[idx["skew",1],1], 
			ipars[idx["shape",1],1])
	j = length(groups)
	gofmat = matrix(NA, ncol = 3, nrow = j)
	gofmat[,1] = groups
	for(i in 1:j){
		sq = seq(1/groups[i], 1, by = 1/groups[i])
		ni = tabulate(findInterval(cdfv, c(0, sq), rightmost.closed = TRUE, 
						all.inside = FALSE)) 
		ExpValue = length(cdfv)/groups[i]
		gofmat[i,2] = sum( ( ( ni - ExpValue ) ^2 )/ExpValue )
		gofmat[i,3] = pchisq(q=gofmat[i,2], df = groups[i]-1, lower.tail = FALSE)
	}
	colnames(gofmat) = c("group", "statistic", "p-value(g-1)")
	rownames(gofmat) = 1:j
	return(gofmat)
}

# returns various loss functions (vectors)
lossfn.var = function(realized, forecast)
{
	x = cbind(realized, forecast)
	#exc = which(is.na(x))
	#if(length(exc)>0) x = x[-exc,]
	mse  = (x[,1] - x[,2])^2
	# qlike: quasi likelihood (implied by gaussian likelihood)
	qlike = log(x[,2]) + x[,1]/x[,2]
	# r2log: penalizes forecasts asymmetrically in low and high volatility periods
	r2log = ( log(x[,1]/x[,2]) )^2
	mad = abs(x[,1] - x[,2])
	#hmse = (x[,1]*(1/x[,2]))^2
	lossdf = data.frame(cbind(mse, mad, qlike, r2log))
	return(lossdf)
}

lossfn.ret = function(realized, forecast)
{
	x = cbind(realized, forecast)
	#exc = which(is.na(x), arr.ind = TRUE)
	#if(length(exc)>0) x = x[-exc[,1],]
	mse  = (x[,1] - x[,2])^2
	mae = abs(x[,1] - x[,2])
	#hmse = (x[,1]*(1/x[,2]))^2
	# we count zero as positive for investment purposes
	dac = as.integer(.signpluszero(x[,1])==.signpluszero(x[,2]))
	lossdf = data.frame(cbind(mse, mae, dac))
	return(lossdf)
}

.meanlossfn = function(lossdf)
{
	n = dim(lossdf)[1]
	nm = colnames(lossdf)
	ans = apply(lossdf, 2, FUN = function(x) 1/(n-1) * sum(x[!is.na(x)]))
	names(ans) = nm
	return(ans)
}

.medianlossfn = function(lossdf)
{
	n = dim(lossdf)[1]
	xn = dim(lossdf)[2]
	nm = colnames(lossdf)
	ans = apply(lossdf[,-xn], 2, FUN = function(x) median(x, na.rm = TRUE))
	# obviously we cannot have the median of a 0/1 series.
	ans = c(ans, .meanlossfn(lossdf[,xn]))
	names(ans) = nm
	return(ans)
}

# Diebold Mariano Test (JBES 1995)
#.DM.test = function(x, lagq)
#{
	# x = vector of differences btw the loss function of the benchmark model and the loss fn of the competing model
	# lagq = lag considered to calculate the newey-west var-cov matrix
	#n.ahead = forecast horizon
	# Returns:
	# an anova class object from the linear.hypothesis test of library car
	# x = rnorm(100, 0, 0.01) - rnorm(100,0.0001,0.01)
#	exc = which(is.na(x))
#	if(length(exc)>0) x = x[-exc]
#	dm = lm(x~1)
#	dm.test = .linear.hypothesis(dm,  c("(Intercept) = 0"), test = "Chisq", vcov  =  NeweyWest(dm, lag = lagq, 
#					order.by = NULL, prewhite = TRUE, adjust = TRUE))
#	return(dm.test)
#}

# From his 2001 paper. Also added the Jarque Bera Test as reccomended by Dowd
# since the test does not really account for residuals being from the Normal distribution.

BerkowitzLR = function(data, lags = 1, significance = 0.05)
{
	if( lags < 1 ) stop("\nlags must be 1 or greater!") else lags = as.integer(lags)
	x = as.numeric(data) - mean(data)
	n = length(data)
	xlag = NULL
	for(i in 1:lags) xlag = cbind(xlag, rugarch:::.lagx(x, n.lag = i, pad = 0))
	ans = lm(x~xlag-1)
	uLL = sum(dnorm(residuals(ans)[-c(1:lags)], sd = summary(ans)$sigma, log = TRUE))
	rLL = sum(dnorm(data[-c(1:lags)], log = TRUE))
	LR = 2*(uLL - rLL)
	chid = 1-pchisq(LR, 2 + lags)
	if(chid < significance) res = paste("reject NULL") else res = paste("fail to reject NULL")
	H0 = paste("Normal(0,1) with no autocorrelation")
	m1 = sum(x)/n
	xm = (x - m1)
	m2 = sum(xm^2)/n
	m3 = sum(xm^3)/n
	m4 = sum(xm^4)/n
	k1 = (m3/m2^(3/2))^2
	k2 = (m4/m2^2)
	JB = n * k1/6 + n * (k2 - 3)^2/24
	JBp = 1 - pchisq(JB, df = 2)
	return(list(uLL = uLL, rLL = rLL, LR = LR, LRp = chid, H0 = H0, Test = res, 
					mu = mean(data), sigma = summary(ans)$sigma, rho =  coef(ans)[1:(lags)],
			JB = JB, JBp = JBp))
}


# Tests of Directional Accuracy
DACTest = function(forecast, actual, test = c("PT", "AG"), conf.level = 0.95)
{
  n = length(actual)
  if( length(forecast) != n ) stop("Length of forecast and actual must be the same")
  if( test == "PT"){
    x_t = z_t = y_t = rep(0, n)
    x_t[which(actual>0)] = 1
    y_t[which(forecast>0)] = 1
    p_y = mean(y_t)
    p_x = mean(x_t)
    z_t[which( (forecast*actual)>0 )] = 1
    p_hat = mean(z_t)
    p_star = p_y*p_x + (1 - p_y)*(1-p_x)
    p_hat_var = (p_star*(1-p_star))/n
    p_star_var = ((2*p_y-1)^2*(p_x*(1-p_x)))/n + ((2*p_x-1)^2*(p_y*(1-p_y)))/n + (4*p_x*p_y*(1-p_x)*(1-p_y))/n^2
    s_n = (p_hat - p_star)/sqrt(p_hat_var - p_star_var)
    ans = list(Test = "Pesaran and Timmermann", Stat = s_n, p.value = 1-pnorm(s_n), H0 = "Independently Distributed",
    Decision = if( s_n >  qnorm(conf.level) ) "Reject  H0" else "Fail to Reject H0",
    DirAcc = p_hat)
  } else{
    r_t=sign(forecast)*actual
    A_t = mean(r_t)
    B_t = mean(sign(forecast))*mean(actual)
    p_y = 0.5 * (1 + mean(sign(forecast)))
    
    V_EP = (4/(n^2))*p_y*(1-p_y)*sum((actual-mean(actual))^2)
    EP = (A_t-B_t)/sqrt(V_EP)
    ans = list(Test = "Anatolyev and Gerko", Stat = EP, p.value = 1-pnorm(EP), H0 = "No Predictability",
    Decision = if( EP >  qnorm(conf.level) ) "Reject  H0" else "Fail to Reject H0",
    DirAcc = sum(r_t>0)/n)
  }
  return( ans )
}

lossfn.ret = function(realized, forecast)
{
	x = cbind(realized, forecast)
	#exc = which(is.na(x), arr.ind = TRUE)
	#if(length(exc)>0) x = x[-exc[,1],]
	mse  = (x[,1] - x[,2])^2
	mad = abs(x[,1] - x[,2])
	#hmse = (x[,1]*(1/x[,2]))^2
	# we count zero as positive for investment purposes
	dac = as.integer(.signpluszero(x[,1])==.signpluszero(x[,2]))
	lossdf = data.frame(cbind(mse, mad, dac))
	return(lossdf)
}

.meanlossfn = function(lossdf)
{
	n = dim(lossdf)[1]
	nm = colnames(lossdf)
	ans = apply(lossdf, 2, FUN = function(x) 1/(n-1) * sum(x[!is.na(x)]))
	names(ans) = nm
	return(ans)
}

.medianlossfn = function(lossdf)
{
	n = dim(lossdf)[1]
	xn = dim(lossdf)[2]
	nm = colnames(lossdf)
	ans = apply(lossdf[,-xn], 2, FUN = function(x) median(x, na.rm = TRUE))
	# obviously we cannot have the median of a 0/1 series.
	ans = c(ans, .meanlossfn(lossdf[,xn]))
	names(ans) = nm
	return(ans)
}



.signpluszero = function(x)
{
	z = as.numeric(x)
	ans = rep(0, length(z))
	ans[z>=0] = 1
	ans
}




.VaRplot = function(varname , p, actual, dates, VaR)
{
	
	y.actual = actual
	y.VaR = VaR
	x.dates = dates
	zd = .makedate(as.character(dates))
	if(zd$status == 0) dates = 1:length(dates) else dates = zd$dates
	title = paste("Daily Returns & Value-at-Risk Exceedances\n","(Series: ", varname,", alpha=", p,")",sep="")
	if(zd$status){
		plot(x = as.Date(dates, format = zd$dformat), y = y.actual, type = "n",
				main = title, ylab = "Daily Log Return", xlab = "time", 
				ylim = c(min(y.actual, y.VaR), max(y.actual, y.VaR)))
	} else{
		plot(as.numeric(dates), y.actual, type = "n",
				main = title, ylab = "Daily Log Return", xlab = "time", 
				ylim = c(min(y.actual, y.VaR), max(y.actual, y.VaR)))
	}
	abline(h = 0, col = 2, lty = 2)
	sel  =  which(y.actual>0)
	points(dates[sel], y.actual[sel], pch = 18, col = "green")
	sel  =  which(y.actual<0)
	points(dates[sel], y.actual[sel], pch = 18, col = "orange")
	sel  =  which(y.actual<y.VaR)
	points(dates[sel], y.actual[sel], pch = 18, col = "red")
	if(zd$status){
		lines(x = as.Date(dates, format = "%Y-%m-%d"), y = y.VaR, lwd = 2, col = "black")
	} else{
		lines(as.numeric(dates), y.VaR, lwd = 2, col = "black")
	}
	legend("topleft", max(actual),c("return >= 0","VaR <= return < 0","return < VaR","VaR"),
			col=c("green","orange","red","black"), cex=0.75,
			pch = c(18,18,18,-1), lty=c(0,0,0,1), lwd=c(0,0,0,2), bty = "n")
}



.VaRreport = function(varname, garchmodel, distribution, p, actual, VaR, conf.level = 0.95)
{
	actual<-as.matrix(actual)
	VaR<-as.matrix(VaR)
	result <- LR.cc.test(p=p,actual=actual,VaR=VaR,conf.level=conf.level)
	cat("VaR Backtest Report\n")
	cat("===========================================\n")
	cat(paste("Model:\t\t\t\t","",garchmodel,"-",distribution,"\n",sep=""))
	cat(paste("Backtest Length:\t",result$TN,"\n",sep=""))
	cat(paste("Data:\t\t\t\t","",varname,"\n",sep=""))
	cat("\n==========================================\n")		
	cat(paste("alpha:\t\t\t\t",round(100*p,1),"%\n",sep=""))
	cat(paste("Expected Exceed:\t",round(p*result$TN,1),"\n",sep=""))
	cat(paste("Actual VaR Exceed:\t",result$N,"\n",sep=""))
	cat(paste("Actual %:\t\t\t",round(100*result$N/result$TN,1),"%\n",sep=""))
	cat("\nUnconditional Coverage (Kupiec)")
	cat("\nNull-Hypothesis:\tCorrect Exceedances\n")
	cat(paste("LR.uc Statistic:\t",round(result$stat.uc,3),"\n",sep=""))
	cat(paste("LR.uc Critical:\t\t",round(result$crit.val.uc,3),"\n",sep=""))
	cat(paste("LR.uc p-value:\t\t",round(result$p.value.uc,3),"\n",sep=""))
	cat(paste("Reject Null:\t\t",ifelse(round(result$p.value.uc,3)< (1-conf.level), "YES","NO"),"\n",sep=""))
	cat("\nConditional Coverage (Christoffersen)")
	cat("\nNull-Hypothesis:\tCorrect Exceedances &\n\t\t\t\t\tIndependence of Failures\n")
	cat(paste("LR.cc Statistic:\t",round(result$stat.cc,3),"\n",sep=""))
	cat(paste("LR.cc Critical:\t\t",round(result$crit.val.cc,3),"\n",sep=""))
	cat(paste("LR.cc p-value:\t\t",round(result$p.value.cc,3),"\n",sep=""))
	cat(paste("Reject Null:\t\t",ifelse(result$reject,"YES","NO"),"\n",sep=""))
}


########################################################################
# Available in a number of locations/textbooks.
# This code originally presented in a webcast by Guy Yollin Feb 2006 (from Insightful)
# Functions to perform Hypothesis test
# on VaR models based on # of exceedances
# calc LR.uc statistic

LR.uc = function(p, TN, N)
{
	stat.uc = -2*log((1-p)^(TN-N)*p^N )+2*log((1-N/TN)^(TN-N)*(N/TN)^N)
	
	return(stat.uc)
}

# calc LR.cc statistic

LR.cc = function(p, actual, VaR)
{
	VaR.ind = ifelse(actual < VaR,1,0)
	N = sum(VaR.ind)
	TN = length(VaR.ind)
	T00 = sum(c(0,ifelse(VaR.ind[2:TN]==0 & VaR.ind[1:(TN-1)]==0,1,0)))
	T11 = sum(c(0,ifelse(VaR.ind[2:TN]==1 & VaR.ind[1:(TN-1)]==1,1,0)))
	T01 = sum(c(0,ifelse(VaR.ind[2:TN]==1 & VaR.ind[1:(TN-1)]==0,1,0)))
	T10 = sum(c(0,ifelse(VaR.ind[2:TN]==0 & VaR.ind[1:(TN-1)]==1,1,0)))
	
	T0 = T00+T01
	T1 = T10+T11
	pi0 = T01/T0
	pi1 = T11/T1
	pe = (T01+T11)/(T0+T1)
	stat.ind =  -2*log((1-pe)^(T00+T10)*pe^(T01+T11))+2*log((1-pi0)^T00*pi0^T01*(1-pi1)^T10*pi1^T11)
	stat.uc = LR.uc(p=p,TN=TN,N=N)
	stat.cc = stat.uc + stat.ind
	return(list(stat.cc = stat.cc, stat.uc = stat.uc, N = N, TN = TN))
}

# perform LR.cc confidence test
LR.cc.test = function(p, actual, VaR, conf.level = 0.95)
{
	result = LR.cc(p = p, actual = actual,VaR = VaR)
	crit.val.uc = qchisq(conf.level, df = 1)
	crit.val.cc = qchisq(conf.level, df = 2)
	p.value.cc = 1-pchisq(result$stat.cc, df = 2)
	p.value.uc = 1-pchisq(result$stat.uc, df = 1)
	reject = ifelse(p.value.cc<1-conf.level, TRUE, FALSE)
	return(list(stat.cc = result$stat.cc, stat.uc = result$stat.uc,
					p.value.cc = p.value.cc, p.value.uc = p.value.uc,
					conf.level=conf.level, reject = reject,
					N = result$N, TN = result$TN, crit.val.uc = crit.val.uc,
					crit.val.cc = crit.val.cc))
}
