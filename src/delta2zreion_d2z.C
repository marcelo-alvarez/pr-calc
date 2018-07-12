#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "proto_d2z.h"
#include "allvars_d2z.h"

void delta2zreion()
{

  int i,j,k,r,index;

  float pnorm,mass,sgm2;
  float rho0 = 2.775e11*Parameters.Omegam*Parameters.h*Parameters.h;

  // Normalize power spectrum
  normalize(Parameters.Omegam,Parameters.Omegab,Parameters.h,Parameters.ns,
	    Parameters.sigma8,&pnorm);

  // Get erfizeta
  erfizeta = inverse_erfc(1./clParameters.zeta);

  // Get sigma_min
  cdmehint(Parameters.Omegam, Parameters.Omegab, Parameters.h, 
	     Parameters.ns, pnorm, clParameters.M_min, &sgm2);
  sigma_min = sqrt(sgm2);

  // Local copies of parameters
  int N=Parameters.N;
  int Nscales=Parameters.Nscales;

  float Rmin=Parameters.Rmin;
  float BoxSize=Parameters.BoxSize;

  // Range of scales to smooth over
  float Rinitial=BoxSize*Rmin/N;
  float Rfinal=clParameters.Rmax/Parameters.h; // This is where we convert from Mpc/h to Mpc
  float dScale=(log10(Rfinal)-log10(Rinitial))/(Nscales-1);

  // Cell sizes and volumes in k space
  float dk=2.*M_PI/BoxSize;
  float d3k=pow(dk,3);
  float kmax=dk*N;

  // Loop over all scales to find zreion from delta
  if(myid==0){
    printf("\n");
    printf("\n Nscales = %2d    Rmin=%5.3f",Nscales,Rinitial*Parameters.h);
    printf("\n                 Rmax=%5.3f\n",Rfinal*Parameters.h);
    printf("\n      Scale       R [Mpc/h]        sigma(R)      sigma_grid");
    printf("\n ==========================================================\n");
  }
  for(r=0;r<Nscales;r++){

    float R = pow(10,log10(Rinitial)+r*dScale);
    
    // Smooth delta on this scale
    for(i=0;i<Nlocal;i++){
      float kx=(i+Nlocal*myid)*dk;
      if(i+Nlocal*myid>=N/2) kx=kx-kmax;
      float kx2=kx*kx;
      
      for(j=0;j<N;j++){
	float ky=j*dk;
	if(j>=N/2) ky=ky-kmax;
	float ky2=ky*ky;
	
	for(k=0;k<=N/2;k++){
	  float kz=k*dk;
	  if(k>=N/2) kz=kz-kmax;
	  float kz2=kz*kz;
	  
	  int index=i*N*(N+2)+j*(N+2)+2*k;

	  float rk=sqrt(kx2+ky2+kz2);	  
	  float rkR=rk*R;
	  float WofrkR=3.*(sin(rkR)-rkR*cos(rkR))/pow(rkR,3);
	  
	  if(rk!=0){
	    delta_smooth[index]   = delta[index]   * WofrkR;
	    delta_smooth[index+1] = delta[index+1] * WofrkR;
	  } else{
	    delta_smooth[index]   = 0;
	    delta_smooth[index+1] = 0;
	  }
	  
	}
      }
    }
    
    // FFT from k to real space    
    fftwf_execute(iplan);
    NormalizeArray(delta_smooth, N, size_fftw, 3);
    float sigma_grid=GetMeanSigma(delta_smooth, 0, N, Nlocal);

    // Get variance on this scale
    mass = 4*M_PI/3*rho0*R*R*R;
    cdmehint(Parameters.Omegam, Parameters.Omegab, Parameters.h, 
    	     Parameters.ns, pnorm, mass, &sgm2);
    float sigma=sqrt(sgm2);

    if(myid==0) printf("\n %10d %15.3f %15.3f %15.3f",
		       r, R*Parameters.h, sigma, sigma_grid);

    // For each cell find zreion
    if(Parameters.RecordBubbleSizes == 1){

      for(i=0;i<Nlocal;i++){
	for(j=0;j<N;j++){
	  for(k=0;k<N;k++){
	    
	    index = i*N*(N+2) + j*(N+2) + k;
	    
	    float dcur = delta_smooth[index];
	    float zcur = d2z(dcur,sigma);
	    
	    if(zcur>zreion[index]){
	      zreion[index]=zcur;
	      rbubble[index]=R;
	    }
	    
	  }
	}
      }

    }
    else{

      for(i=0;i<Nlocal;i++){
	for(j=0;j<N;j++){
	  for(k=0;k<N;k++){
	    
	    index = i*N*(N+2) + j*(N+2) + k;
	    
	    float dcur = delta_smooth[index];
	    float zcur = d2z(dcur,sigma);
	    
	    if(zcur>zreion[index]){
	      zreion[index]=zcur;
	    }

	  }
	}
      }
      
    }

    // If in lean mode, FFT back to k-space and deconvolve 
    if(clParameters.lean == 1){
      
      // FFT
      
      fftwf_execute(plan);
      
      // Deconvolution
      for(i=0;i<Nlocal;i++){
	float kx=(i+Nlocal*myid)*dk;
	if(i+Nlocal*myid>=N/2) kx=kx-kmax;
	float kx2=kx*kx;
	
	for(j=0;j<N;j++){
	  float ky=j*dk;
	  if(j>=N/2) ky=ky-kmax;
	  float ky2=ky*ky;
	  
	  for(k=0;k<=N/2;k++){
	    float kz=k*dk;
	    if(k>=N/2) kz=kz-kmax;
	    float kz2=kz*kz;
	    
	    int index=i*N*(N+2)+j*(N+2)+2*k;
	    
	    float rk=sqrt(kx2+ky2+kz2);	  
	    float rkR=rk*R;
	    float WofrkR=3.*(sin(rkR)-rkR*cos(rkR))/pow(rkR,3);
	    
	    if(rk!=0){
	      delta[index]   = delta_smooth[index]   / WofrkR;
	      delta[index+1] = delta_smooth[index+1] / WofrkR;
	    } else{
	      delta[index] = 0;
	      delta[index+1] = 0;
	    }
	    
	  }
	}
      }
      
    }
       
    // End loop over smoothing scales
    
  }
  
  if(myid==0){
    printf("\n");
  }

  free(cfftwa_r2c);
  if(clParameters.lean==0) free(cfftwa_c2r);

  return;
  
}

float d2z(float dcur, float sigma)
{

  float deltac=1.686;

  float z=(1+zHigh)*DHigh/deltac*
    (dcur+sqrt(2.*(sigma_min*sigma_min-sigma*sigma))*erfizeta)-1;

  return z;

}

