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

#################################################################################
# GARCH rolling test
#################################################################################


rugarch.test7a = function(parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2))
{
	#cat("\nrugarch-->test7-1: Roll Test (apARCH)\n")
	tic = Sys.time()
	
	data(sp500ret)
	spec = ugarchspec(
			variance.model = list(model = "eGARCH"), distribution.model = "jsu")
	roll = ugarchroll(spec,  data = sp500ret, n.ahead = 1, forecast.length = 500,
			refit.every = 25, refit.window = "recursive", 
			parallel = parallel, parallel.control = parallel.control,
			solver = "solnp", fit.control = list(scale = 1), 
			solver.control = list(tol = 1e-5, delta = 1e-6, trace=0),
			calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05))
	
	postscript("test7a1.eps", width = 12, height = 8)
	plot(roll, which = "all")
	dev.off()
	
	postscript("test7a2.eps", width = 12, height = 8)
	plot(roll, which = 5)
	dev.off()
	
	z1 <- file("test7a1.txt", open="wt")
	sink(z1)
	report(roll, type = "VaR", n.ahead = 1, VaR.alpha = 0.01, conf.level = 0.95)
	sink(type="message")
	sink()
	close(z1)
	
	z2 <- file("test7a2.txt", open="wt")
	sink(z2)
	report(roll, type = "fpm")
	sink(type="message")
	sink()
	close(z2)
	
	z3 <- file("test7a3.txt", open="wt")
	sink(z3)
	print(as.data.frame(roll, which = "coefs"))
	sink(type="message")
	sink()
	close(z3)
	
	z4 <- file("test7a4.txt", open="wt")
	sink(z4)
	print(as.data.frame(roll, which = "density"))
	sink(type="message")
	sink()
	close(z4)
	
	z5 <- file("test7a5.txt", open="wt")
	sink(z5)
	print(as.data.frame(roll, which = "coefmat", refit = 1))
	sink(type="message")
	sink()
	close(z5)
	
	z6 <- file("test7a6.txt", open="wt")
	sink(z6)
	print(as.data.frame(roll, which = "VaR"))
	sink(type="message")
	sink()
	close(z6)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test7b = function(parallel = FALSE, parallel.control = list(pkg = c("multicore", "snowfall"), cores = 2))
{
	#cat("\nrugarch-->test7-2: Roll Test (apARCH)\n")
	tic = Sys.time()
	
	data(sp500ret)
	spec = ugarchspec(
			variance.model = list(model = "apARCH"), distribution = "sstd")
	roll = ugarchroll(spec,  data = sp500ret, n.ahead = 5, forecast.length = 500, 
			refit.every = 25, refit.window = "moving", parallel = parallel, 
			parallel.control = parallel.control, solver = "solnp", 
			fit.control = list(scale = 1), solver.control = list(scale = 1),
			calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05))
	
	postscript("test7b1.eps", width = 12, height = 8)
	plot(roll, n.ahead = 5, which = "all")
	dev.off()
	
	postscript("test7b2.eps", width = 12, height = 8)
	plot(roll, which = 5)
	dev.off()
	
	zz <- file("test7b1.txt", open="wt")
	sink(zz)
	report(roll, type = "VaR", n.ahead = 1, VaR.alpha = 0.01, conf.level = 0.95)
	report(roll, type = "VaR", n.ahead = 5, VaR.alpha = 0.01, conf.level = 0.95)
	sink(type="message")
	sink()
	close(zz)
	
	z1 <- file("test7b2.txt", open="wt")
	sink(z1)
	report(roll, type = "fpm")
	sink(type="message")
	sink()
	close(z1)
	
	z2 <- file("test7b3.txt", open="wt")
	sink(z2)
	print(as.data.frame(roll, which = "coefs"))
	sink(type="message")
	sink()
	close(z2)
	
	z3 <- file("test7b4.txt", open="wt")
	sink(z3)
	print(as.data.frame(roll, which = "density", n.ahead = 1))
	print(as.data.frame(roll, which = "density", n.ahead = 5))
	sink(type="message")
	sink()
	close(z3)
	
	z4 <- file("test7b5.txt", open="wt")
	sink(z4)
	print(as.data.frame(roll, which = "coefmat", refit = 1))
	print(as.data.frame(roll, which = "coefmat", refit = 8))
	sink(type="message")
	sink()
	close(z4)
	
	z5 <- file("test7b6.txt", open="wt")
	sink(z5)
	print(as.data.frame(roll, which = "VaR", n.ahead = 1))
	print(as.data.frame(roll, which = "VaR", n.ahead = 5))
	sink(type="message")
	sink()
	close(z5)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}
