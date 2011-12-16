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
# implements nlminb, lbgs and solnp
# only solnp implements true constraint (stationarity) optimization
.garchsolver = function(solver, pars, fun, Ifn, ILB, IUB, gr, hessian, parscale, 
		control, LB, UB, ux=NULL, ci=NULL, mu=NULL, ...)
{
	gocontrol = control
	control = .getcontrol(solver, control)
	retval = switch(solver,
			nlminb = .nlminbsolver(pars, fun, gr, hessian, parscale, control, LB, UB, ...),
			solnp = .solnpsolver(pars, fun, Ifn, ILB, IUB, control, LB, UB, ...),
			gosolnp = .gosolnpsolver(pars, fun, Ifn, ILB, IUB, gocontrol, LB, UB, ...),
			lbfgs = .lbfgssolver(pars, fun, gr, parscale, control, LB, UB, ...))
	return(retval)
}

.nlminbsolver = function(pars, fun, gr, hessian, parscale, control, LB, UB,...){
	ans = try(nlminb(start = pars, objective = fun, gradient = gr, hessian = hessian,
			..., scale = 1/parscale, control = control, lower = LB, upper = UB), silent = TRUE)
	if(inherits(ans, "try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = ans
		sol$par = rep(NA, length(pars))
		names(sol$par) = names(pars)
	} else{
		sol = ans
	}
	hess = NULL
	return(list(sol = sol,hess = hess))
}

.solnpsolver = function(pars, fun, Ifn, ILB, IUB, control, LB, UB, ...){
	ans = try(solnp(pars, fun = fun, eqfun = NULL, eqB = NULL, ineqfun = Ifn, ineqLB = ILB, 
					ineqUB = IUB, LB = LB, UB = UB, control = control, ...), silent = TRUE)
	if(inherits(ans,"try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = ans
		sol$par = rep(NA, length(pars))
		names(sol$par) = names(pars)
		warning("\nrgarch-->warning: no convergence...\n")
	} else{
		sol = ans
	}
	hess = NULL
	return(list(sol = sol, hess = hess))
}

.gosolnpsolver = function(pars, fun, Ifn, ILB, IUB, gocontrol, LB, UB, ...){
	control = .solnpctrl(gocontrol)
	gocontrol = .gosolnpctrl(gocontrol)
	n.restarts = gocontrol$n.restarts
	parallel = gocontrol$parallel
	parallel.control = gocontrol$parallel.control
	rseed = gocontrol$rseed
	n.sim = gocontrol$n.sim
	op <- options()
	options(warn = 0)
	
	# use the truncated normal distribution
	distr.opt = vector(mode = "list", length = length(pars))
	for(i in 1:length(pars)){
		distr.opt[[i]]$mean = pars[i]
		distr.opt[[i]]$sd = sqrt(pars[i]^2)*2
	}
	# ok parallel will work with snowfall without changing the fun and Ifn to rugarch:::fun and rugarch:::Ifn
	ans = try(gosolnp(pars = pars, fixed = NULL, fun = fun, eqfun = NULL, 
			eqB = NULL, ineqfun = Ifn, ineqLB = ILB, 
			ineqUB = IUB, LB = LB, UB = UB, control = control, distr = rep(2, length(LB)), distr.opt = distr.opt, 
			n.restarts = n.restarts, n.sim = n.sim, parallel = parallel, parallel.control = parallel.control, 
			rseed = rseed, ...),
	silent = TRUE)
	if(inherits(ans,"try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = ans
		sol$par = rep(NA, length(pars))
		names(sol$par) = names(pars)
	} else{
		sol = ans
	}
	hess = NULL
	options(op)
	return(list(sol = sol, hess = hess))
}

.lbfgssolver = function(pars, fun, gr, parscale, control, LB, UB, ...){
	control$parscale = parscale
	ans = try(optim(par = pars, fn = fun, gr = gr, ...,
			method = "L-BFGS-B", lower = LB, upper = UB, control = control, 
			hessian = TRUE),silent=TRUE)
	if(inherits(ans, "try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = ans
		sol$par = rep(NA, length(pars))
		names(sol$par) = names(pars)
	} else{
		sol = ans
	}
	hess = sol$hessian
	return(list(sol = sol, hess = hess))
}

# default control for solvers:
.getcontrol = function(solver, control)
{
	ans = switch(solver,
		nlminb = .nlminbctrl(control),
		solnp = .solnpctrl(control),
		gosolnp = .gosolnpctrl(control),
		lbfgs = .lbfgsctrl(control))
	return(ans)
}

.nlminbctrl = function(control)
{
	if(is.null(control$eval.max)) control$eval.max = 2000
	if(is.null(control$iter.max)) control$iter.max = 1500
	if(is.null(control$abs.tol)) control$abs.tol = 1e-20
	if(is.null(control$rel.tol)) control$rel.tol = 1e-10
	if(is.null(control$x.tol)) control$x.tol = 1.5e-8
	if(is.null(control$step.min)) control$step.min = 2.2e-14
	return(control)
}

.lbfgsctrl = function(control)
{
	if(is.null(control$REPORT)) control$REPORT = 10
	if(is.null(control$lmm)) control$lmm = 15
	if(is.null(control$pgtol)) control$pgtol = 1e-8
	if(is.null(control$factr)) control$factr = 1e-8
	return(control)
}

.solnpctrl = function(control){
	# parameters check is now case independent
	ans = list()
	params = unlist(control)
	if(is.null(params)) {
		ans$rho = 1
		ans$outer.iter = 50
		ans$inner.iter = 1800
		ans$delta = 1.0e-8
		ans$tol = 1.0e-8
		ans$trace = 1
	} else{
		npar = tolower(names(unlist(control)))
		names(params) = npar
		if(any(substr(npar, 1, 3) == "rho")) ans$rho = as.numeric(params["rho"]) else ans$rho = 1
		if(any(substr(npar, 1, 5) == "outer.iter")) ans$outer.iter = as.numeric(params["outer.iter"]) else ans$outer.iter = 50
		if(any(substr(npar, 1, 5) == "inner.iter")) ans$inner.iter = as.numeric(params["inner.iter"]) else ans$inner.iter = 1000
		if(any(substr(npar, 1, 5) == "delta")) ans$delta = as.numeric(params["delta"]) else ans$delta = 1.0e-8
		if(any(substr(npar, 1, 3) == "tol")) ans$tol = as.numeric(params["tol"]) else ans$tol = 1.0e-8
		if(any(substr(npar, 1, 5) == "trace")) ans$trace = as.numeric(params["trace"]) else ans$trace = 1
	}
	return(ans)
}

.gosolnpctrl = function(control){
	# parameters check is now case independent
	ans = list()
	params = unlist(control)
	if(is.null(params)) {
		ans$parallel = FALSE
		ans$parallel.control = list(pkg = "snowfall", cores = 2)
		ans$n.restarts = 1
		ans$rseed
		ans$n.sim = 500
	} else{
		npar = tolower(names(unlist(control)))
		names(params) = npar
		ans$parallel.control = list()
		if(any(substr(npar, 1, 8) == "parallel")) ans$parallel = as.logical(params["parallel"]) else ans$parallel = FALSE
		if(any(substr(npar, 1, 20) == "parallel.control.pkg")) ans$parallel.control$pkg = params["parallel.control.pkg"] else ans$parallel.control$pkg = "snowfall"
		if(any(substr(npar, 1, 22) == "parallel.control.cores")) ans$parallel.control$cores = params["parallel.control.cores"] else ans$parallel.control$cores = "snowfall"
		if(any(substr(npar, 1, 10) == "n.restarts")) ans$n.restarts = as.numeric(params["n.restarts"]) else ans$n.restarts = 1
		if(any(substr(npar, 1, 5) == "rseed")) ans$rseed = as.numeric(params["rseed"]) else ans$rseed = NULL
		if(any(substr(npar, 1, 5) == "n.sim")) ans$n.sim = as.numeric(params["n.sim"]) else ans$n.sim = 500		
	}
	return(ans)
}

#----------------------------------------------------------------------------------
.garchconbounds = function(){
	return(list(LB = eps,UB = 0.999))
}

.sgarchcon = function(pars, data, returnType, garchenv){
	ipars = get("ipars", garchenv)
	estidx = get("estidx", garchenv)
	idx = get("model", garchenv)$pidx
	ipars[estidx, 1] = pars
	distribution = get("model", garchenv)$modeldesc$distribution
	con = .persistsgarch1(pars = ipars[,1], idx, distribution = distribution)
	return(con)
}

.igarchcon = function(pars, data, returnType, garchenv){
	# this is an equality constraint
	ipars = get("ipars", garchenv)
	estidx = get("estidx", garchenv)
	idx = get("model", garchenv)$pidx
	modelinc = get("model", garchenv)$modelinc
	ipars[estidx, 1] = pars
	con = ifelse(modelinc[9]>1, sum(ipars[idx["alpha", 1]:idx["alpha", 2], 1]) + sum(ipars[idx["beta", 1]:(idx["beta", 2]-1), 1]) -
					ipars[idx["beta", 2],1], sum(ipars[idx["alpha", 1]:idx["alpha", 2], 1]) - ipars[idx["beta", 2],1])
	return(con)
}

.aparchcon = function(pars, data, returnType, garchenv){
	ipars = get("ipars", garchenv)
	estidx = get("estidx", garchenv)
	idx = get("model", garchenv)$pidx
	ipars[estidx, 1] = pars
	distribution = get("model", garchenv)$modeldesc$distribution	
	con = .persistaparch1(pars = ipars[,1], idx = idx, distribution = distribution)
	if(is.na(con)) con = 1
	return(con)
}

.fgarchcon = function(pars, data, returnType, garchenv){
	
	ipars = get("ipars", garchenv)
	estidx = get("estidx", garchenv)
	idx = get("model", garchenv)$pidx
	ipars[estidx, 1] = pars
	distribution = get("model", garchenv)$modeldesc$distribution
	vsubmodel = get("model", garchenv)$modeldesc$vsubmodel
	con = .persistfgarch1(pars = ipars[,1], idx = idx, distribution = distribution, submodel = vsubmodel)
	if(is.na(con)) con = 1
	return(con)
}

.gjrgarchcon = function(pars, data, returnType, garchenv){
	ipars = get("ipars", garchenv)
	estidx = get("estidx", garchenv)
	idx = get("model", garchenv)$pidx
	ipars[estidx, 1] = pars
	distribution = get("model", garchenv)$modeldesc$distribution
	con = .persistgjrgarch1(pars = ipars[,1], idx = idx, distribution = distribution)
	if(is.na(con)) con = 1
	return(con)
}

.egarchcon = function(pars, data, returnType, garchenv){
	ipars = get("ipars", garchenv)
	estidx = get("estidx", garchenv)
	idx = get("model", garchenv)$pidx
	ipars[estidx, 1] = pars
	distribution = get("model", garchenv)$modeldesc$distribution	
	con = .persistegarch1(pars = ipars[,1], idx = idx, distribution = distribution)
	if(is.na(con)) con = 1
	return(con)
}