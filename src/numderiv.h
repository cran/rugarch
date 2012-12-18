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
#ifndef _numderiv_H
#define _numderiv_H
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rcpp.h>
RcppExport SEXP hessian2sided(SEXP fun, SEXP x, SEXP H, SEXP ee, SEXP gm, SEXP gp);
#endif
