#include <math.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>
#include "cosmology.h"

/////////////////////////////////////////////////////////
// Routines to perform selected cosmological calculations
// Adopted from "Cosmological Routine Library in Fortran"
// by E. Komatsu
// http://gyudon.as.utexas.edu/~komatsu/CRL/
/////////////////////////////////////////////////////////

// "Common block"
double omega0,baryon0,hubble0,an0,pnorm0;

// function prototypes
double sgm2k(double , void *);
double transeh(double);

float tau0(float zmax, float omegam, float omegab, float h){  

  float thompson = 6.65e-25;
  float YHe = 0.25;
  float hubble0 = 100./3.086e19*h;
  float ne0=omegab*h*h*1.88e-29/1.67e-24;
  float c=3e10;

  float dtau0=thompson*ne0*c/hubble0;

  float tau=0;
  float NHe;
  float epsilon=1e-3;
  int   nz = 1000;
  float dz = zmax / nz;
  for(int i=1;i<nz-1;i++){

    float z=(i-0.5)*dz;

    NHe=1;
    if(z<3) NHe=2;

    float dtau = dtau0*pow((1+z),2)/sqrt(omegam*pow((1+z),3)+1-omegam)*
           (1-YHe+NHe*YHe/4)*dz;
    tau += dtau;

  }

  return tau;

}

float tau(float *za, float *xa, int nz, float omegam, float omegab, float h){  

  float thompson = 6.65e-25;
  float YHe = 0.25;
  float hubble0 = 100./3.086e19*h;
  float ne0=omegab*h*h*1.88e-29/1.67e-24;
  float c=3e10;

  float dtau0=thompson*ne0*c/hubble0;

  float tau=0;
  float NHe;
  float epsilon=1e-3;
  for(int i=1;i<nz-1;i++){

    float z=za[i];
    float x=xa[i];

    float dz=fabs(za[i+1]-za[i]);

    NHe=1;
    if(z<3) NHe=2;

    float dtau = x*dtau0*pow((1+z),2)/sqrt(omegam*pow((1+z),3)+1-omegam)*
           (1-YHe+NHe*YHe/4)*dz;
    tau += dtau;

  }

  return tau;

}

float dtb(float *za, float *dtba, int nz, float omegam, float omegab, float h){  

  float hoc = h * 1e2 / 3e5; // H0/c in units of 1/Mpc
  float Tb0   = 7.3e3*(omegab*h*h/0.02)*sqrt((0.15/omegam/h/h)); // Tb in muK
  float dzdchi0 = pow(omegam,0.5) * hoc;

  float dtb=0;
  for(int i=1;i<nz-1;i++){

    float z=za[i];
    float x=dtba[i];

    float dz=fabs(za[i+1]-za[i]);

    float ddtb = x * Tb0 * hoc * (1+z) * (1+z) * dz;
    dtb += ddtb;

  }

  return dtb;

}

float ksz(float *za, float *dk, int nz, float omegam, float omegab, float h){  

  float thompson = 6.65e-25;
  float YHe = 0.25;
  float hubble0 = 100./3.086e19*h;
  float ne0=omegab*h*h*1.88e-29/1.67e-24;
  float c=3e10;

  float dtau0=thompson*ne0*c/hubble0;

  float tau=0;
  float NHe;
  float ksz=0;
  for(int i=0;i<nz-1;i++){

    float z=za[i];
    float dkszdtau=dk[i];

    float dz=fabs(za[i+1]-za[i]);

    NHe=1;
    if(z<3) NHe=2;

    float dtau = dtau0*pow((1+z),2)/sqrt(omegam*pow((1+z),3)+1-omegam)*
           (1-YHe+NHe*YHe/4)*dz;

    ksz += dkszdtau * dtau;
  }

  return ksz;

}

float deltac(float z, float omegam, float omegal, float w, float *dcoll)
{

  float x,x2,x3;
  float omega;           // Omega_m(z)
  float eta;
  float pi= 3.14159265358979;

  float dc;

  x  = 1+z;
  x2 = x*x;
  x3 = x*x*x;

  omega = omegam*x3/( omegam*x3 + (1-omegam-omegal)*x2 + omegal );
  eta= log( (2/omega-1) + sqrt(pow((2/omega-1),2)-1) );

  if ( (omegam<1) && (omegal==0) ) {    
    dc= 1.5*(3*sinh(eta)*(sinh(eta)-eta)
	       /pow((cosh(eta)-1),2)-2)
               *(1+pow(2*pi/(sinh(eta)-eta),2./3.));
    *dcoll= 4*pi*pi*pow((cosh(eta)-1),3)/pow((sinh(eta)-eta),2);
  }else{
    dc= 3*pow((12*pi),2./3.)/20*(1+0.0123*log10(omega));    
    *dcoll= 18*pi*pi*(1+0.4093*pow((1/omega-1),0.9052));
  }

  dc/=growth(z, omegam, omegal, w);

  return dc;

}

float growth(float z, float omegam, float omegal, float w)
{

  float x,x2,x3,x3w;
  float omega, lambda;    // Omega_m(z, Omega_lambda(z)
  float g, g0, t, t0;

  float D;

  x  = 1+z;
  x2 = x*x;
  x3 = x*x*x;
  x3w= pow(x3,w);

  omega  = omegam*x3/( omegam*x3 + (1-omegam-omegal)*x2 + omegal );
  lambda = omegal*x3*x3w/(omegam*x3+(1-omegam-omegal)*x2+omegal*x3*x3w);

  // Fitting formulae for a w=-1 universe [Carrol, Press & Turner (1992)]

  g  = 2.5*omega /(pow(omega ,4./7.)-lambda+(1+omega /2)*(1+lambda/70));
  g0 = 2.5*omegam/(pow(omegam,4./7.)-omegal+(1+omegam/2)*(1+omegal/70));
		 
  // Fitting formulae for a w/=-1 universe [Ma, Caldwell, Bode & Wang (1999)]
  // << for spatially flat universe only >>

  t = -( 0.255+0.305*w+0.0027/w )*( 1-omega )
      -( 0.366+0.266*w-0.0700/w )*log(omega);
  t0= -( 0.255+0.305*w+0.0027/w )*( 1-omegam )
      -( 0.366+0.266*w-0.0700/w )*log(omegam);
  g  *= pow(-w,t );
  g0 *= pow(-w,t0);

  // Calculate the growth factor

  D= (g/x)/g0;

  return D;

}

void cdmehint(float omega, float baryon, float hubble, float an, 
	      float pnorm, float mass, float *sgm2)
{

  double R, dens;
  double sgm2d, error;

  float rhoc=2.775e11;
  float pi=3.141592654, fourpi=4*pi;

  omega0  = omega;
  baryon0 = baryon;
  hubble0 = hubble;
  an0     = an;
  pnorm0  = pnorm;

  dens = omega*rhoc*hubble*hubble;
  R    = pow( 3*mass/fourpi/dens, 1./3.);

  gsl_integration_workspace * work = gsl_integration_workspace_alloc(1000);

  gsl_function SGM2K;
  SGM2K.function = &sgm2k;
  SGM2K.params = &R;

  gsl_integration_qags(&SGM2K, log(1e-5/R), log(1e5/R), 0, 1e-7, 
		       1000, work, &sgm2d, &error);

  *sgm2=sgm2d;

}

void normalize(float omega, float baryon, float hubble, float an, 
	       float sigma8, float *pnorm)
{

  double R;
  double sgm2d, error;

  omega0  = omega;
  baryon0 = baryon;
  hubble0 = hubble;
  an0     = an;
  pnorm0  = 1;

  R = 8/hubble;

  gsl_integration_workspace * work = gsl_integration_workspace_alloc(1000);

  gsl_function SGM2K;
  SGM2K.function = &sgm2k;
  SGM2K.params = &R;

  gsl_integration_qags(&SGM2K, log(1e-5/R), log(1e5/R), 0, 1e-7, 
		       1000, work, &sgm2d, &error);

  *pnorm=sigma8/sqrt(sgm2d);

}

double sgm2k(double lk, void * params){

  double k  = exp(lk);  
  double R  = *(double *)params;
  double kR =k*R;

  double Pin = pow(k,3+an0);
  double T=transeh(k);
  double W=3*(sin(kR)-kR*cos(kR))/pow(kR,3);

  double sgm2k = Pin*T*T*pow(W,2)*pow(pnorm0,2);

  return sgm2k;

}

double transeh(double k)
{

  double alpha, sound, shape, q;

  alpha = 1
    -0.328*log(431*omega0*hubble0*hubble0)*(baryon0/omega0)  
    +0.38*log(22.3*omega0*hubble0*hubble0)*pow(baryon0/omega0,2);
  sound = 44.5*log(9.83/(omega0*hubble0*hubble0))
    /sqrt(1+10*pow(baryon0*hubble0*hubble0,3./4.));
  shape = omega0*hubble0
    *(alpha+(1-alpha)/(1+pow(0.43*k*sound,4)));
  q = (k/hubble0)*pow(2.728/2.7,2)/shape;
  
  double T= log(2*exp(1)+1.8*q)/(log(2*exp(1)+1.8*q)
			    +(14.2+731/(1+62.5*q))*q*q);
  
  return T;

}

float fcoll(float z, float zInit, float omegam, float omegal, float w, 
	    float sigma_min){

  float arg, dcoll;

  float f;

  arg = deltac(z, omegam, omegal, w, &dcoll)/sqrt(2.)/sigma_min;

  f = erfc(arg);

  return f;

}

