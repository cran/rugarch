/*################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2015.
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
#ifndef FRACDIFF_H
#define FRACDIFF_H
void fracdiff(int *n, double *d, double *p, double *x, double *ydiff);
void c_binexpansion(int *n, double *d, double *ans);
void c_figarchcons(double *alpha1, double *delta, double *beta1, double *g, double *psi, int* truncLag);
#endif /* FRACDIFF_H */


