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

# function to deal with the numerous data formats present
# we want to extract the date from the data
.extractdata<-function(data)
{
	tsclass = class(data)[1]
	valid.choices = c("timeSeries", "zoo", "zooreg", "numeric", "data.frame", "xts", "matrix")
	if(!any(valid.choices == tsclass)) stop("\nrgarch-->error: class of data object not recognized")
	x = switch(tsclass,
			timeSeries = .xseries.timeSeries(data),
			zoo = .xseries.zoo(data),
			zooreg = .xseries.zoo(data),
			xts = .xseries.xts(data),
			numeric = .xseries.numeric(data),
			data.frame = .xseries.dataframe(data),
			matrix = .xseries.matrix(data))
	return(x)
}

.xseries.timeSeries = function(data){
	x = unclass(data)
	if(!is.null(dim(data)[2]) && dim(data)[2]>1) stop("only univariate dataset supported")
	xdata = as.numeric(x)
	rdates = .makedate(as.character(time(data)))
	if(rdates$status){
		xdates = rdates$dates
		dformat = rdates$dformat
	} else{
		xdates = 1:length(x)
		dformat = "numeric"
	}
	return(list(data = xdata, pos = xdates, dformat = dformat))
}

.xseries.zoo = function(data){
	x = unclass(data)
	if(!is.null(dim(data)[2]) && dim(data)[2]>1) stop("only univariate dataset supported")
	xdata = as.numeric(x)
	rdates = .makedate(as.character(index(data)))
	if(rdates$status){
		xdates = rdates$dates
		dformat = rdates$dformat
	} else{
		xdates = 1:length(x)
		dformat = "numeric"
	}
	return(list(data = xdata, pos = xdates, dformat = dformat))
}

.xseries.xts = function(data){
	x = unclass(data)
	if(!is.null(dim(data)[2]) && dim(data)[2]>1) stop("only univariate dataset supported")
	xdata = as.numeric(x)
	rdates = .makedate(as.character(index(data)))
	if(rdates$status){
		xdates = rdates$dates
		dformat = rdates$dformat
	} else{
		xdates = 1:length(x)
		dformat = "numeric"
	}
	return(list(data = xdata, pos = xdates, dformat = dformat))
}

.xseries.numeric = function(data){
	x = unclass(data)
	xdata = x
	if(!is.null(names(data))){
		rdates = .makedate(names(data))
		if(rdates$status){
			xdates = rdates$dates
			dformat = rdates$dformat
		} else{
			xdates = 1:length(x)
			dformat = "numeric"
		}
	} else{
		xdates = 1:length(x)
		dformat = "numeric"
	}
	return(list(data = xdata, pos = xdates, dformat = dformat))
}

.xseries.dataframe = function(data){
	xdata = as.numeric(data[,1])
	if(!is.null(dim(data)[2]) && dim(data)[2]>1) stop("only univariate dataset supported")
	ow <- options("warn")
	options(warn = (-1))
	if(!is.null(rownames(data))){
		if(!is.na(as.numeric(rownames(data)[1]))){
			xdates = as.numeric(rownames(data))
			dformat = "numeric"
		} else{
			rdates = .makedate(rownames(data))
			if(rdates$status){
				xdates = rdates$dates
				dformat = rdates$dformat
			} else{
				xdates = 1:length(xdata)
				dformat = "numeric"
			}
		}
	} else{
		xdates = 1:length(xdata)
		dformat = "numeric"
	}
	options(ow)
	return(list(data = xdata, pos = xdates, dformat = dformat))
}

.xseries.matrix<-function(data){
	xdata = as.numeric(data[,1])
	if(!is.null(dim(data)[2]) && dim(data)[2]>1) stop("only univariate dataset supported")
	if(!is.null(rownames(data))){
		rdates = .makedate(rownames(data))
		if(rdates$status){
			xdates = rdates$dates
			dformat = rdates$dformat
		} else{
			xdates = as.character(1:length(xdata))
			dformat = "numeric"
		}
	} else{
		xdates = as.character(1:length(xdata))
		dformat = "numeric"
	}
	return(list(data = xdata, pos = xdates, dformat = dformat))
}

# make generatefwd search for available dates from original set
.generatefwd = function(Dates, N, dformat, periodicity = "days")
{
	if(is.numeric(Dates[1])){
		n = length(Dates)
		fwd = (Dates[n]+1):(n+N)
	} else{
		n = length(Dates)
		# reformat data
		sdx = format(Dates[n], "%m/%d/%y")
		fwd = seq.dates(from=sdx, by = periodicity, length.=N*8)
		# generate enough to dates to not cause problem with weekend exclusion later
		z1 = which(chron::is.weekend(fwd))
		fwd = as.character(fwd)
		fwd = fwd[-z1]
		fwd = fwd[2:(N+1)]
		fwd = as.character(format(strptime(as.character(fwd), "%m/%d/%y"), dformat))
		fwd = as.Date(fwd, format = dformat)
	}
	return(fwd)
}

ForwardDates = function(Dates, n.ahead, date.format, periodicity="days")
{
	UseMethod("ForwardDates")
}

.ForwardDates = function(Dates, n.ahead, date.format, periodicity="days")
{
	if(is.numeric(Dates[1])){
		n = length(Dates)
		fwd = (Dates[n]+1):(n+n.ahead)
	} else{
		n = length(Dates)
		# reformat data
		sdx = format(as.Date(Dates[n], format = date.format), format = "%m/%d/%y")
		fwd = seq.dates(from = sdx, by = periodicity, length. = 8 * n.ahead)
		# generate enough to dates to not cause problem with weekend exclusion later
		z1 = which(chron::is.weekend(fwd))
		fwd = format(as.Date(as.character(fwd), format = "%m/%d/%y" ), format = date.format)
		if(length(z1)>0) fwd = fwd[-z1]
		fwd = fwd[2:(n.ahead+1)]
		fwd = as.Date(fwd, format = date.format)
	}
	return(fwd)
}
	
setMethod("ForwardDates", signature(Dates = "character"), definition=.ForwardDates)

WeekDayDummy = function(Dates, date.format, weekday = "Monday")
{
	UseMethod("WeekDayDummy")
}

.WeekDayDummy = function(Dates, date.format, weekday = "Monday")
{
	valid.weekdays = c("Monday","Tuesday","Wednesday","Thursday","Friday")
	wd = match(weekday, valid.weekdays)
	if(is.na(wd)) stop("not a valid weekday")
	xdates = .makedate(Dates)
	if(xdates$status == 0) stop("\nForwardDates-->error: date format not recognized\n", call. = FALSE)
	ddates = xdates$dates
	dformat = xdates$dformat
	ans = rep(0, length(Dates))
	zz = which(weekdays(ddates) == weekday)
	ans[zz]=1
	ans
}

setMethod(f = "WeekDayDummy", signature(Dates = "character"), definition = .WeekDayDummy)

# create the forecast date set which depends on out.sample data
.forcdates = function( origdates, n.ahead, N, i, ns , dformat)
{
	# need to return both a timevector and character vector
	if( ns == 0 ){
		fwwdx = .generatefwd(origdates, N = n.ahead, dformat = dformat, 
				periodicity = "days")
	} else{
		if( (n.ahead + i - 1) <= ns ){
			fwwdx = origdates[(N + i):(N + i + n.ahead -1)]
		} else{
			if( i <= ns ){
				nt = max(0, ns - i + 1)
				fwwdx1 = origdates[(N+i):(N+ns)]
				fwwdx2 = .generatefwd(origdates[1:(N+ns)], N = n.ahead-nt, 
								dformat = dformat, periodicity = "days")
				fwwdx =c(fwwdx1, fwwdx2)		
			} else{
				fwwdx = .generatefwd(origdates[1:(N+ns)], N = n.ahead, 
								dformat = dformat, periodicity = "days")
			}
			
		}
	}
	fwdd = fwwdx
	return(fwdd)
}


.makedate = function(x)
{
	# find the divisor: 4 cases "-", "/", ".", and no divisor
	allc = strsplit(x[1], "")
	
	if(any(allc[[1]] == "-")){
		dt = "-"
		ld = length(which(diff(which(allc[[1]]!="-"))==1))+3
		dte = t(apply(as.data.frame(x), 1, FUN=function(z) as.numeric(strsplit(z, dt)[[1]]) ))
	} else if(any(allc[[1]] == "/")){
		dt = "/"
		ld = length(which(diff(which(allc[[1]]!="/"))==1))+3
		dte = t(apply(as.data.frame(x), 1, FUN=function(z) as.numeric(strsplit(z, dt)[[1]]) ))
	} else if(any(allc[[1]] == ".")){
		dt = "."
		dte = t(apply(as.data.frame(x), 1, FUN=function(z) as.numeric(strsplit(z, dt)[[1]]) ))
	} else{
		# this is a little more complicated
		ld = length(allc[[1]])
		if(ld==6){
			dte = t(apply(as.data.frame(x), 1, FUN=function(z) 
								as.numeric(c(substr(z, 1,2), substr(z, 3,4), substr(z, 5,6)))))
		} else if(ld==8){
			# 2 cases either the 4 digit year is at the beginning or else at the end
			dte.1 = as.vector(t(apply(as.data.frame(x), 1, FUN=function(z) 
								as.numeric(c(substr(z, 1,2))))))
			dte.2 = as.vector(t(apply(as.data.frame(x), 1, FUN=function(z) 
										as.numeric(c(substr(z, 5,6))))))
			if(all(dte.1>18)){
				dte = t(apply(as.data.frame(x), 1, FUN=function(z) 
									as.numeric(c(substr(z, 1,4), substr(z, 5,6), substr(z, 7,8)))))
			} else if(all(dte.2>18)){
				dte = t(apply(as.data.frame(x), 1, FUN=function(z) 
									as.numeric(c(substr(z, 1,2), substr(z, 3,4), substr(z, 5,8)))))
			} else{
				return(list(status=0))
			}
		} else{
			return(list(status=0))	
		}
	}
	m = 0
	for(i in 1:3){
		if(all(dte[,i]<=12)) m = i
	}
	if(m==0) return(list(status=0))
	sq = 1:3
	sq = sq[-m]
	y = 0 
	for(i in sq){
		if(any(dte[,i]>31)) y = i
	}
	if(y==0) return(list(status=0))
	d = sq[-y]
	dmatrix = cbind(dte[,d], dte[,m], dte[,y])
	colnames(dmatrix) = c("d","m","y")
	if(ld==6){
		ddates = as.Date(paste(dmatrix[,3], dmatrix[,2], dmatrix[,1], sep = "-"), format="%y-%m-%d")
		dformat = "%y-%m-%d"
	} else{
		ddates = as.Date(paste(dmatrix[,3], dmatrix[,2], dmatrix[,1], sep = "-"), format="%Y-%m-%d")
		dformat = "%Y-%m-%d"
	}
	

	return(list(datesmat = dmatrix, dates = ddates, dformat = dformat, status=1))
}