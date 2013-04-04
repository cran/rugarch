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
#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H
double* paramghskt(const double , const double );
double dghsktstd(const double , const double , const double );
double rghyp(const double , const double , const double );
double rsnig(const double , const double );
double rstd(const double );
double rsstd(const double , const double );
double qjsu(const double , const double , const double );
double rjsu(const double , const double );
double rnig(const double , const double , const double , const double );
double rsnorm(const double );
double rged(const double );
double rsged(const double , const double );
double signum(const double );
double dhuge(void);
double depsilon(void);
double kappagh(const double, const double);
double deltakappagh(const double, const double);
double* paramgh(const double, const double, const double);
double dnormstd(const double);
double dsnormstd(const double, const double);
double dgedstd(const double, const double);
double dsgedstd(const double, const double, const double);
double xdt(const double, const double);
double dstdstd(const double, const double);
double dsstdstd(const double, const double, const double);
double dhyp(const double , const double , const double , const double , const double , const int );
double dhypstd(const double ,  const double , const double , const int );
double dgh(const double, const double, const double, const double, const double, const double, const int);
double dghstd(const double, const double, const double, const double, const int);
double dnig(const double, const double, const double, const double, const double, const int);
double dnigstd(const double, const double, const double, const int);
double djsustd(const double , const double , const double );
double garchdistribution(const double, const double, const double , const double , const double, const int);
double rgarchdist(const double , const double , const double , const int );
double pgarchdist(const double , const double , const double , const double , const double , const double , const int );
void xdnormstd(double *x, double *pdf);
void xdsnormstd(double *x, double *xi, double *pdf);
void xdstdstd(double *x, double *nu, double *pdf);
void xdsstdstd(double *x, double *nu, double *xi, double *pdf);
void xdgedstd(double *x, double *nu, double *pdf);
void xdsgedstd(double *x, double *nu, double *xi, double *pdf);
void xdgh(double *x, double *alpha, double *beta, double *delta, double *mu, double *lambda, int *logr, double *pdf);
void xdghstd(double *x, double *zeta, double *rho, double *lambda, int *logr, double *pdf);
void xdnig(double *x, double *alpha, double  *beta, double *delta, double  *mu , int *logr, double *pdf);
void xdnigstd(double *x, double *zeta, double *rho, int *logr, double *pdf);
void xrnig(double *alpha, double *beta, double *delta, double *mu, double *ans, int *n);
void xrghyp(double *mu, double *sigma, double *rho, double *zeta, double *lambda, double *ans, int *n);
void distributionsample(double *shape, double *skew, double *lambda, int *ndis, int *n, double *rvec);
double Heaviside(const double, const double);
double psnorm(const double , const double , const double , const double );
double pged(const double , const double , const double , const double );
double psged(const double , const double , const double , const double , const double );
double pstd(const double , const double , const double , const double );
double psstd(const double , const double , const double , const double , const double );
double pjsu(const double , const double , const double , const double , const double );
void xpsnorm(double *q, double *mu, double *sd, double *xi, double *p, int *n);
void xpged(double *q, double *mu, double *sd, double *nu, double *p, int *n);
void xpsged(double *q, double *mu, double *sd, double *nu, double *xi, double *p, int *n);
void xpstd(double *q, double *mu, double *sd, double *nu, double *p, int *n);
void xpsstd(double *q, double *mu, double *sd, double *nu, double *xi, double *p, int *n);
void xpjsu(double *q, double *mu, double *sd, double *nu, double *tau, double *p, int *n);
#endif /* DISTRIBUTIONS_H */
