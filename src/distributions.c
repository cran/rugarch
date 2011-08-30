/*################################################################################
##
##   R package rgarch by Alexios Ghalanos Copyright (C) 2009, 2010, 2011
##   This file is part of the R package rgarch.
##
##   The R package rgarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rgarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################*/
# include "distributions.h"
# include <R.h>
# include <limits.h>
# include "gig.h"
# include "nig.h"
# include <math.h>
# include <Rmath.h>

double rghyp(const double zeta, const double rho, const double lambda)
{
	double *param;
	param = paramgh(zeta, rho, lambda);
	double alpha = param[0];
	double beta = param[1];
	double delta = param[2];
	double mu = param[2];
	double chi = delta*delta;
	double psi = (alpha * alpha) - (beta * beta);
	double W = rgig(lambda, chi, psi);
	double ans = mu + W * beta + sqrt(W) * rnorm(0, 1);
	free(param);
	return(ans);
}
double rsnig(const double zeta, const double rho)
{
	double gamma, v, y, x0, x1, x2, x3, p1, u, x, ans;
	double *param;
	param = paramgh(zeta, rho, -0.5);
	double alpha = param[0];
	double beta = param[1];
	double delta = param[2];
	double mu = param[3];
	gamma = sqrt(alpha*alpha - beta*beta);
	v = rchisq(1);
	x0 = delta/gamma;
	x1 = sqrt(4 * gamma * delta * v + pow(v,2));
	x2 = x0 + 1/(2 * gamma*gamma) * (v + x1);
	x3 = x0 + + 1/(2 * gamma*gamma) * (v - x1);
	p1 = delta/(delta + gamma * x2);
	u = runif(0,1);
	if(u<p1){
		x = x2;
	} else{
		x = x3;
	}
	y = rnorm(0, 1);
	ans = sqrt(x) * y + mu + beta * x;
	free(param);
	return(ans);
}

double rstd(const double nu)
{
	double ans = 0;
	if(nu > 2.0) 
	{
	double s = sqrt(nu/(nu-2));
	ans = rt(nu) * 1.0 / s;
	}
	return(ans);
}

double rsstd(const double nu, const double xi)
{
	double weight, z, rr, m1, mu, sigma, xx, ans;
	ans = 0.0;
	weight = xi / (xi + 1.0/xi);
	z = runif(-1.0 * weight, 1.0 - weight);
	xx = (z < 0)? 1.0/xi : xi;
	rr = -1.0 * fabs(rstd(nu))/xx * sign(z);
	m1 = 2.0 * sqrt(nu - 2.0) / (nu - 1.0) / beta(0.5, 0.5 * nu);
	mu = m1 * (xi - 1.0/xi);
	sigma =  sqrt((1.0 - (m1 * m1)) * ((xi * xi) + 1.0/(xi * xi)) + 2 * (m1 * m1) - 1.0);
	ans =  (rr - mu ) / sigma;
	return(ans);
}

double qjsu(const double x, const double nu, const double tau)
{
	double rtau, rr, z, w, omega, cc, ans;
	ans = 0.0;
	rtau = 1.0/tau;
	rr = qnorm(x, 0.0, 1.0, 1, 0);
	z = sinh(rtau * (rr + nu));
	w = (rtau<0.0000001) ? 1 : exp(rtau * rtau);
	omega = -1.0 * nu * rtau;
	cc = sqrt(1/(0.5 * (w - 1.0)*(w * cosh(2.0 * omega) + 1)));
	ans = (cc * sqrt(w) * sinh(omega)) + cc * z;
	return(ans);   
}

double rjsu(const double nu, const double tau)
{
	double x, ans;
	ans = 0.0;
	x = runif(0, 1);
	ans = qjsu(x, nu, tau);
	return(ans);
}

double rnig(const double alpha, const double beta, const double delta, const double mu)
{
	double gamma, v, y, x0, x1, x2, x3, p1, u, x, ans;
	gamma = sqrt(alpha*alpha - beta*beta);
	v = rchisq(1);
	x0 = delta/gamma;
	x1 = sqrt(4 * gamma * delta * v + pow(v,2));
	x2 = x0 + 1/(2 * gamma*gamma) * (v + x1);
	x3 = x0 + + 1/(2 * gamma*gamma) * (v - x1);
	p1 = delta/(delta + gamma * x2);
	u = runif(0,1);
	x = (u < p1)? x2 : x3;
	y = rnorm(0, 1);
	ans = sqrt(x) * y + mu + beta * x;
	return(ans);
}

double rsnorm(const double xi)
{
	double weight, z, rr, m1, mu, sigma, xx, ans;
	weight = xi / (xi + 1.0/xi);
	z = runif(-weight, 1.0 - weight);
	xx = (z < 0)? 1.0/xi : xi;
	rr = -1.0 * fabs(rnorm(0, 1))/xx * sign(z);  
	m1 = 2.0/sqrt(2.0 * PI);
	mu = m1 * (xi - 1.0/xi);
	sigma = sqrt((1 - (m1 * m1)) * ( (xi * xi) + 1.0/(xi* xi) ) + 2 * (m1 * m1) - 1.0);
	ans = (rr - mu ) / sigma;
	return(ans);
}

double rged(const double nu)
{
	double lambda, rr, ans;
	ans = 0.0;
	lambda = sqrt ( pow(0.5, 2.0/nu) * gammafn(1.0/nu) / gammafn(3.0/nu) );
	rr = rgamma(1.0/nu, 1.0);
	ans =  lambda * pow(2*rr, 1/nu) * sign(runif(0, 1) - 0.5);
	return(ans);
}

double rsged(const double nu, const double xi)
{
	double weight, lambda, g, z, rr, m1, mu, sigma, xx, ans;
	weight = xi / (xi + 1.0/xi);
	z = runif(-1.0 * weight, 1.0 - weight);
	xx = (z < 0)? 1.0/xi : xi;
	rr = -1.0 * fabs(rged(nu))/xx * sign(z);
	lambda = sqrt ( pow(0.5, 2.0/nu) * gammafn(1.0/nu) / gammafn(3.0/nu) );
	g  = nu / ( lambda * (pow(2, 1.0 +1.0/nu)) * gammafn(1.0/nu) );
	m1 = pow(2, 1.0/nu) * lambda * gammafn(2.0/nu) / gammafn(1.0/nu);
	mu = m1 * (xi - 1.0/xi);
	sigma = sqrt((1 - (m1 * m1)) * ( (xi * xi) + 1.0/(xi* xi) ) + 2 * (m1 * m1) - 1.0);
	ans = (rr - mu ) / sigma;
	return(ans);
}

double signum(const double x)
{
	double res=-(x<0)+(x>0);
	return  res;
}
double dhuge(void)
{
  return HUGE_VAL;
}
double depsilon(void)
{
  double r;
  r = 1.0;
  while( 1.0 < (double)(1.0 + r)){
	r = r / 2.0;
  }
  return (2.0 * r);
}
double kappagh(const double x, const double lambda)
{
	double kappa=0;
	if(lambda == -0.5){
	kappa = 1/x;
	} else{
	kappa = (bessel_k(x, lambda+1, 2) / bessel_k(x, lambda, 2)) / x;
	}
	return kappa;
}
double deltakappagh(const double x, const double lambda)
{
	double deltakappa=0;
	if(lambda == -0.5){
		deltakappa = kappagh(x, lambda+1) - kappagh(x, lambda);
	} else{
		deltakappa = kappagh(x, lambda+1) - kappagh(x, lambda);
	}
	return deltakappa;
}
double* paramgh(const double zeta, const double rho, const double lambda)
{
	double *param = malloc(4*sizeof(double));
	if(param == NULL){
		exit(EXIT_FAILURE);
	}
	double rho2 = 1 - pow(rho,2);
	double alpha = pow(zeta,2) * kappagh(zeta, lambda)/rho2;
	alpha = alpha * ( 1 + pow(rho,2) * pow(zeta,2) * deltakappagh(zeta, lambda)/rho2);
	alpha = sqrt(alpha);
	double beta = alpha * rho;
	double delta = zeta / ( alpha * sqrt(rho2) );
	double mu = -beta * pow(delta,2) * kappagh(zeta, lambda);
	param[0]=(double) alpha;
	param[1]=(double) beta;
	param[2]=(double) delta;
	param[3]=(double) mu;
	return param;
}
double dnormstd(const double x)
{
  double pdf;
  pdf = exp ( -0.5 * x * x ) / sqrt ( 2.0 * PI );
  if(pdf == 0.0) pdf = 0.0 + 2.22507e-128;
  return pdf;
}
double dsnormstd(const double x, const double xi)
{
	double pdf;
	double m1, mu, sigma,z, xxi, g;
	xxi=xi;
	m1 = 2/sqrt(2*PI);
	mu = m1*(xi-1/xi);
	sigma = sqrt((1-pow(m1,2))*(pow(xi,2)+1/pow(xi,2))+2*pow(m1,2)-1);
	z = x*sigma+mu;
	if(z==0){
		xxi=1;
	}
	if(z<0){
		xxi = 1/xi;
	}
	g = 2.0/(xi + 1.0/xi);
	pdf = g * dnormstd(z/xxi)*sigma;
	return pdf;
}
double dgedstd(const double x, const double nu)
{
	double lambda, g, pdf;
	lambda = sqrt(pow(1.0/2.0,2.0/nu)*gammafn( 1.0/nu )/gammafn( 3.0/nu));
    g = nu/(lambda*(pow(2.0,1.0+(1.0/nu)))*gammafn( 1.0/nu));
	pdf = g*exp(-0.5*pow(fabs(x/lambda),nu));
	return pdf;
}
double dsgedstd(const double x, const double nu, const double xi)
{
	double lambda, m1, mu, sigma, z, g, pdf, xxi;
	xxi=xi;
	lambda = sqrt(pow(1.0/2.0,(2/nu))*gammafn( 1.0/nu )/gammafn( 3.0/nu));
	g = nu/(lambda*(pow(2.0,1.0+(1.0/nu)))*gammafn( 1.0/nu));
	m1 = pow(2.0, (1.0/nu))*lambda*gammafn(2.0/nu)/gammafn(1.0/nu);
	mu = m1*(xi-1.0/xi);
	sigma = (1 - pow(m1,2.0))*(pow(xi,2.0)+1.0/(pow(xi,2.0))) + 2.0*(pow(m1,2))-1.0;
	sigma = sqrt(sigma);
	z = x*sigma+mu;
	if(z==0){
		xxi=1;
	}
	if(z<0){
		xxi = 1/xi;
	}
	g = 2.0/(xi + 1.0/xi);
	pdf = g*dgedstd(z/xxi, nu)*sigma;
	return pdf;
}
double xdt(const double x, const double nu)
{
	double a, b, pdf;
	a = gammafn((nu+1.0)/2.0)/sqrt(PI*nu);
    	b = gammafn(nu/2.0)*pow((1.0+(x*x)/nu),((nu+1.0)/2.0));
    	pdf = a/b;
	return pdf;
}
double dstdstd(const double x, const double nu)
{
	double pdf, s;
	if(nu<=2){
		pdf = 999;
	} else{
		s = sqrt(nu/(nu-2.0));
		pdf = s*xdt(x*s,nu);
	}
	return pdf;
}
double dsstdstd(const double x, const double nu, const double xi)
{
	double mu, m1,beta, sigma, z, g,pdf,a,b, xxi;
	xxi=xi;
	a = 1.0/2.0;
	b = nu/2.0;
	beta = (gammafn(a)/gammafn(a+b))*gammafn(b);
	m1 = 2.0*sqrt(nu-2.0)/(nu-1.0)/beta;
	mu = m1*(xi-1.0/xi);
	sigma = sqrt((1.0-pow(m1,2))*(pow(xi,2)+1.0/(pow(xi,2)))+2.0*pow(m1,2)-1.0);
	z = x*sigma + mu;
	if(z==0){
		xxi=1;
	}
	if(z<0){
		xxi = 1/xi;
	}
	g = 2.0/(xi+1.0/xi);
	pdf = g*dstdstd(z/xxi,nu)*sigma;
	return pdf;
}
double dgh(const double x, const double alpha, const double beta, const double delta, const double mu, const double lambda, const int logr)
{
	double pdf=0;
	if(alpha<=0){
		return pdf=0;
	}
	if(delta <= 0){
		return pdf=0;
	}
	if(fabs(beta) >= alpha){
		return pdf=0;
	}
	double arg = delta*sqrt(pow(alpha,2)-pow(beta,2));
	double a = (lambda/2.0)*log(pow(alpha,2)-pow(beta,2)) - (log(sqrt(2*PI)) + (lambda-0.5)*log(alpha)
			+ lambda*log(delta) + log(bessel_k(arg, lambda, 2)) - arg );
	double f = ((lambda-0.5)/2.0)*log(pow(delta,2)+pow((x - mu),2));
	arg = alpha * sqrt(pow(delta,2)+pow((x-mu),2));
	double k = log(bessel_k(arg, lambda-0.5, 2)) - arg;
	double e = beta*(x-mu);
	pdf = a + f + k + e;
	if(logr==1){
		pdf = pdf;
	} else{
		pdf = exp(pdf);
	}

	return pdf;
}
double dghstd(const double x, const double zeta, const double rho, const double lambda, const int logr)
{
	double pdf;
	double *param;
	param = paramgh(zeta, rho, lambda);
	double alpha=param[0];
	double beta=param[1];
	double delta=param[2];
	double mu=param[3];
	pdf = dgh(x, alpha, beta, delta, mu, lambda, logr);
	free(param);
	return pdf;
}

double dnig(const double x, const double alpha, const double  beta, const double delta, const double mu, const int logr)
{
	double pdf=0;
	double lambda=-0.5;
	pdf = dgh(x, alpha, beta, delta, mu, lambda, logr);
	return pdf;
}
double dnigstd(const double x, const double zeta, const double rho, const int logr)
{
	double pdf=0;
	double lambda=-0.5;
	pdf = dghstd(x, zeta, rho, lambda, logr);
	return pdf;
}
double djsustd(const double x, const double shp, const double skw)
{
		double w, z, r, omega, c, pdf=0.0;
		double rtau=1.0/shp;
		if(rtau<0.0000001){
			w = 1.0;
		} else{
			w = exp(rtau*rtau);
		}
		omega=-skw*rtau;
		c=sqrt(1/(0.5*(w-1)*(w*cosh(2*omega)+1)));
		z=(x-(c*sqrt(w)*sinh(omega)))/c;
		r=-skw + asinh(z)/rtau;
		pdf= -log(c)-log(rtau)-0.5*log(z*z+1)-0.5*log(2*PI)-0.5*r*r;
		pdf=exp(pdf);
		return pdf;
}
double garchdistribution(const double zz, const double hh, const double skew, const double shape, const double dlambda, const int ndis)
{
	/*ndist: 1- normal
			 2 - skew-normal
			 3 - student
			 4 - skew-student
			 5 - ged
			 6 - skew-ged
			 7 - Normal Inverse Gaussian (NIG)
			 8 - Generalized Hyperbolic (GHYP)
	distEq: [0]: Skew, [1]: Shape, ... additionally can be used for normal mixture
	(not yet implemented)
	*/
	double pdf=0;
	if(ndis==1)
	{
		pdf=dnormstd(zz)/hh;
	}
	if(ndis==2)
	{
		pdf=dsnormstd(zz, skew)/hh;
	}
	if(ndis==3)
	{
		pdf=dstdstd(zz, shape)/hh;
	}
	if(ndis==4)
	{
		pdf=dsstdstd(zz, shape, skew)/hh;
	}
	if(ndis==5)
	{
		pdf=dgedstd(zz, shape)/hh;
	}
	if(ndis==6)
	{
		pdf=dsgedstd(zz, shape, skew)/hh;
	}
	if(ndis==7)
	{
		pdf=dnigstd(zz, shape, skew, 0)/hh;
	}
	if(ndis==8)
	{
		pdf=dghstd(zz, shape, skew, dlambda, 0)/hh;
	}
	if(ndis==9)
	{
		pdf=djsustd(zz,shape,skew)/hh;
	}
	return pdf;
}

double rgarchdist(const double shape, const double skew, const double lambda, const int ndis)
{
	/*ndist: 1- normal
			 2 - skew-normal
			 3 - student
			 4 - skew-student
			 5 - ged
			 6 - skew-ged
			 7 - Normal Inverse Gaussian (NIG)
			 8 - Generalized Hyperbolic (GHYP)
			 9 - Johnson's SU (JSU)
	distEq: [0]: Skew, [1]: Shape, ... additionally can be used for normal mixture
	(not yet implemented)
	*/
	double ans=0.0;
	if(ndis==1)
	{
		ans=rnorm(0, 1);
	}
	if(ndis==2)
	{
		ans=rsnorm(skew);
	}
	if(ndis==3)
	{
		ans=rstd(shape);
	}
	if(ndis==4)
	{
		ans=rsstd(shape, skew);
	}
	if(ndis==5)
	{
		ans=rged(shape);
	}
	if(ndis==6)
	{
		ans=rsged(shape, skew);
	}
	if(ndis==7)
	{
		ans = rsnig(shape, skew);
	}
	if(ndis==8)
	{
		ans = rghyp(shape, skew, lambda);
	}
	if(ndis==9)
	{
		ans = rjsu(shape, skew);
	}
	return(ans);
}


double pgarchdist(const double q, const double mu, const double sigma, const double shape, const double skew, const double lambda, const int ndis)
{
	/*ndist: 1- normal
			 2 - skew-normal
			 3 - student
			 4 - skew-student
			 5 - ged
			 6 - skew-ged
			 7 - Normal Inverse Gaussian (NIG)
			 8 - Generalized Hyperbolic (GHYP)
			 9 - Johnson's SU (JSU)
	distEq: [0]: Skew, [1]: Shape, ... additionally can be used for normal mixture
	(not yet implemented)
	*/
	double ans=0.0;
	if(ndis==1)
	{
		ans=pnorm(q, mu, sigma, 1, 0);
	}
	if(ndis==2)
	{
		ans=psnorm(q, mu, sigma, skew);
	}
	if(ndis==3)
	{
		ans=pstd(q, mu, sigma, shape);
	}
	if(ndis==4)
	{
		ans=psstd(q, mu, sigma, shape, skew);
	}
	if(ndis==5)
	{
		ans=pged(q, mu, sigma, shape);
	}
	if(ndis==6)
	{
		ans=psged(q, mu, sigma, shape, skew);
	}
	if(ndis==7)
	{
		// Not Yet Implemented in C code
		ans = 0.5;
	}
	if(ndis==8)
	{
		// Not Yet Implemented in C code
		ans = 0.5;
	}
	if(ndis==9)
	{
		ans = pjsu(q, mu, sigma, skew, shape);
	}
	return(ans);
}

void xdnormstd(double *x, double *pdf)
{
	*pdf = dnormstd(*x);
}
void xdsnormstd(double *x, double *xi, double *pdf)
{
	*pdf = dsnormstd(*x,*xi);
}

void xdstdstd(double *x, double *nu, double *pdf)
{
	*pdf = dstdstd(*x,*nu);
}
void xdsstdstd(double *x, double *nu, double *xi, double *pdf)
{
	*pdf = dsstdstd(*x,*nu,*xi);
}
void xdgedstd(double *x, double *nu, double *pdf)
{
	*pdf = dgedstd(*x, *nu);
}
void xdsgedstd(double *x, double *nu, double *xi, double *pdf)
{
	*pdf = dsgedstd(*x, *nu, *xi);
}
void xdgh(double *x, double *alpha, double *beta, double *delta, double *mu, double *lambda, int *logr, double *pdf)
{
	*pdf = dgh(*x, *alpha, *beta, *delta, *mu, *lambda, *logr);
}
void xdghstd(double *x, double *zeta, double *rho, double *lambda, int *logr, double *pdf)
{
	*pdf = dghstd(*x, *zeta, *rho, *lambda, *logr);
}
void xdnig(double *x, double *alpha, double  *beta, double *delta, double  *mu , int *logr, double *pdf)
{
	*pdf = dnig(*x, *alpha, *beta, *delta, *mu, *logr);
}
void xdnigstd(double *x, double *zeta, double *rho, int *logr, double *pdf)
{
	*pdf = dnigstd(*x, *zeta, *rho, *logr);
}

void xrnig(double *alpha, double *beta, double *delta, double *mu, double *ans, int *n)
{
	int i;
	for(i=0;i<*n;i++)
	{
		ans[i] = rnig(*alpha, *beta, *delta, *mu);
	}
}

void xrghyp(double *mu, double *sigma, double *rho, double *zeta, double *lambda, double *ans, int *n)
{
	int i;
	for(i=0;i<*n;i++)
	{
		ans[i] = rghyp(*zeta, *rho, *lambda);
	}
}

void distributionsample(double *shape, double *skew, double *lambda, int *ndis, int *n, double *rvec)
{
    GetRNGstate();
	int i;
	switch(*ndis){
	case 1:
		for(i=0;i<*n;i++){
			rvec[i] = rnorm(0, 1);
		}
		break;
	case 2:
		for(i=0;i<*n;i++){
			rvec[i] = rsnorm(*skew);
		}
		break;
	case 3:
		for(i=0;i<*n;i++){
			rvec[i] = rstd(*shape);
		}
		break;
	case 4:
		for(i=0;i<*n;i++){
			rvec[i] = rsstd(*shape, *skew);
		}
		break;
	case 5:
		for(i=0;i<*n;i++){
			rvec[i] = rged(*shape);
		}
		break;
	case 6:
		for(i=0;i<*n;i++){
			rvec[i] = rsged(*shape, *skew);
		}
		break;
	case 7:
		for(i=0;i<*n;i++){
			rvec[i] = rsnig(*shape, *skew);
		}
		break;
	case 8:
		for(i=0;i<*n;i++){
			rvec[i] = rghyp(*shape, *skew, *lambda);
		}
		break;
	case 9:
		for(i=0;i<*n;i++){
			rvec[i] = rjsu(*shape, *skew);
		}
		break;
	default:
		for(i=0;i<*n;i++){
			rvec[i] = rnorm(0, 1);
		}
		break;
	}
    PutRNGstate();
}

double Heaviside(const double x, const double a){
	return( (signum(x-a) + 1.0)/2.0 );
}

double psnorm(const double q, const double mu, const double sd, const double xi)
{
	double qx = (q-mu)/sd;
	double m1 = 2.0/sqrt(2*PI);
	double mux = m1 * (xi - 1.0/xi);
	double sigma = sqrt((1.0-m1*m1)*(xi*xi+1.0/(xi*xi)) + 2.0*m1*m1 - 1.0);
	double z = qx*sigma + mux;
	double Xi = (z<0)?1.0/xi:xi;
	double g = 2.0/(xi + 1.0/xi);
	double p = Heaviside(z, 0) - signum(z) * g * Xi * pnorm(-fabs(z)/Xi, 0, 1, 1, 0);
	return( p );
}

double pged(const double q, const double mu, const double sd, const double nu)
{
	double qx = (q-mu)/sd;
	double lambda = sqrt ( 1.0/pow(2.0, (2.0/nu)) * gammafn(1.0/nu) / gammafn(3.0/nu) );
	double g  = nu / ( lambda * (pow(2,(1+1.0/nu))) * gammafn(1.0/nu) );
	double h = pow(2.0, (1.0/nu)) * lambda * g * gammafn(1.0/nu) / nu;
	double s = 0.5 * pow( fabs(qx) / lambda , nu );
	double p = 0.5 + signum(qx) * h * pgamma(s, 1.0/nu, 1, 1, 0);
	return( p );
}

double psged(const double q, const double mu, const double sd, const double nu, const double xi)
{
	double qx = (q-mu)/sd;
	double lambda = sqrt ( 1.0/pow(2.0, 2.0/nu) * gammafn(1.0/nu) / gammafn(3.0/nu) );
	double m1 = pow(2.0, 1.0/nu) * lambda * gammafn(2.0/nu) / gammafn(1.0/nu);
	double mux = m1*(xi-1.0/xi);
	double sigma =  sqrt((1.0-m1*m1)*(xi*xi+1/(xi*xi)) + 2.0*m1*m1 - 1);
	double z = qx*sigma + mux;
	double Xi = (z<0)?1.0/xi:xi;
	double g = 2.0/(xi + 1.0/xi);
	double p = Heaviside(z, 0) - signum(z) * g * Xi * pged(-fabs(z)/Xi, 0, 1, nu);
	return( p );
}

double pstd(const double q, const double mu, const double sd, const double nu)
{
	double s = sqrt(nu/(nu-2.0));
	double z = (q-mu)/sd;
	double p = pt(z*s, nu, 1, 0);
	return( p );
}

double psstd(const double q, const double mu, const double sd, const double nu, const double xi)
{
	double qx = (q-mu)/sd;
	double m1 = 2.0 * sqrt(nu-2.0) / (nu-1.0) / beta(0.5, nu/2.0);
	double mux = m1*(xi-1.0/xi);
	double sigma =  sqrt((1-m1*m1)*(xi*xi+1/(xi*xi)) + 2*m1*m1 - 1);
	double z = qx*sigma+mux;
	double Xi = (z<0)?1.0/xi:xi;
	double g = 2.0 / (xi + 1.0/xi);
	double p = Heaviside(z, 0) - signum(z) * g * Xi * pstd(-fabs(z)/Xi, 0, 1, nu);
	return( p );
}

double pjsu(const double q, const double mu, const double sd, const double nu, const double tau)
{
	double rtau = 1.0/tau;
	double w = (rtau<0.0000001)?1.0:exp(rtau*rtau);
	double omega = -1.0*nu*rtau;
	double c = 1/sqrt(0.5*(w-1.0)*(w*cosh(2*omega)+1));
	double z = (q-(mu+c*sd*sqrt(w)*sinh(omega)))/(c*sd);
	double r = -1.0*nu + asinh(z)/rtau;
	double p = pnorm(r,0,1,1,0);
	return( p );
}

void xpsnorm(double *q, double *mu, double *sd, double *xi, double *p, int *n){
	int i;
	for(i=0;i<*n;i++){
		p[i] = psnorm(q[i], *mu, *sd, *xi);
	}
}

void xpged(double *q, double *mu, double *sd, double *nu, double *p, int *n)
{
	int i;
	for(i=0;i<*n;i++){
		p[i] = pged(q[i], *mu, *sd, *nu);
	}
}

void xpsged(double *q, double *mu, double *sd, double *nu, double *xi, double *p, int *n)
{
	int i;
	for(i=0;i<*n;i++){
		p[i] = psged(q[i], *mu, *sd, *nu, *xi);
	}
}

void xpstd(double *q, double *mu, double *sd, double *nu, double *p, int *n)
{
	int i;
	for(i=0;i<*n;i++){
		p[i] = pstd(q[i], *mu, *sd, *nu);
	}
}

void xpsstd(double *q, double *mu, double *sd, double *nu, double *xi, double *p, int *n)
{
	int i;
	for(i=0;i<*n;i++){
		p[i] = psstd(q[i], *mu, *sd, *nu, *xi);
	}
}

void xpjsu(double *q, double *mu, double *sd, double *nu, double *tau, double *p, int *n)
{
	int i;
	for(i=0;i<*n;i++){
		p[i] = pjsu(q[i], *mu, *sd, *nu, *tau);
	}
}
