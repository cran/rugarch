/*################################################################################
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
#################################################################################*/
#include "numderiv.h"
using namespace Rcpp;

SEXP hessian2sided(SEXP fun, SEXP x, SEXP H, SEXP deps, SEXP gminus, SEXP gplus)
{
	try {
		Rcpp::NumericVector rx(x);
		Rcpp::NumericMatrix xeps(deps);
		Rcpp::NumericVector gp(gplus);
		Rcpp::NumericVector gm(gminus);
		Rcpp::NumericMatrix xH(H);
		Rcpp::Function f(fun);
		int n = rx.size();
		int nx = xeps.nrow(), mx = xeps.ncol(), i, j;
		Rcpp::NumericVector tmp(n);
		Rcpp::NumericVector xtmp(1);
		Rcpp::NumericVector fx(1);
		fx = f(x);
		Rcpp::NumericMatrix hp(nx, mx);
		Rcpp::NumericMatrix hm(nx, mx);
		for(i=0;i<n;i++){
			for(j=i;j<n;j++){
				tmp = rx+xeps(_,i)+xeps(_,j);
				xtmp = f(tmp);
				hp(i,j) = xtmp[0];
				hp(j,i) = hp(i,j);
				tmp = rx-xeps(_,i)-xeps(_,j);
				xtmp = f(tmp);
				hm(i,j) = xtmp[0];
				hm(j,i) = hm(i,j);
			}
		}
		for(i=0;i<n;i++){
			for(j=i;j<n;j++){
				xH(i,j) = ((hp(i,j)-gp(i)-gp(j)+fx[0]+fx[0]-gm(i)-gm(j)+hm(i,j))/xH(i,j))/2;
				xH(j,i) = xH(i,j);
			}
		}
		return wrap(xH);
	} catch( std::exception &ex ) {
			forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rugarch-->hessian c++ exception (unknown reason)" );
	 }
	 return R_NilValue;
}

SEXP colMaxRcpp(SEXP X_){
	Rcpp::NumericMatrix X(X_);
	int n = X.ncol();
	Rcpp::NumericVector V(n);
	for (int i=0; i<n; i++) {
	Rcpp::NumericVector W = X.column(i);
	V[i] = *std::max_element(W.begin(), W.end());
	}
	return(V);
}
