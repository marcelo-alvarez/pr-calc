#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "proto_ics.h"
#include "allvars_ics.h"

void GenerateNoise()
{

  int i,j,k,r,index;

  // Local copies of parameters
  int N=clParameters.N;

  long seed = (myid + 123) * clParameters.Seed;

  srand48(seed);

  for(int i=0;i<size_fftw;i++){
    delta[i] = gaussdev();
  }

  GetMeanSigma(delta,1,N,Nlocal);

  // FFT of noise to k-space
  fftwf_execute(plan);
  NormalizeArray(delta,N,size_fftw,1.5);
 
  if(myid==0) printf("\n Noise generation done...");

}

void ConvolveNoise()
{

  int i,j,k,r,index;

  // Local copies of parameters
  int N=clParameters.N;
  float BoxSize=clParameters.BoxSize;

  // Cell sizes and volumes in k space
  float dk=2.*M_PI/BoxSize;
  float d3k=pow(dk,3);
  float kmax=dk*N;

  // Convert P(k) from power per d3k to mean amplitude per k-mode
  for(int i=0;i<ps_nk;i++) Wavenumber2PdeltaTable[i]=sqrt(Wavenumber2PdeltaTable[i]*d3k);

  // Convolve noise with the power spectrum
  for(i=0;i<local_nx;i++){
    float kx=(i+local_nx*myid)*dk;
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
	
	if(rk!=0){
	  double ps = Wavenumber2FloatLogSpace(ps_nk, ps_kmin, ps_kmax, rk, Wavenumber2PdeltaTable);
	  delta[index]   *= ps;
	  delta[index+1] *= ps;
	} else{
	  delta[index] = 0;
	  delta[index+1] = 0;
	}
	
      }
    }
  }
  
  // FFT of delta to real space
  fftwf_execute(iplan);

  GetMeanSigma(delta,1,N,Nlocal);

  if(myid==0) printf("\n Convolution done...");
  
  return;
  
}

double gaussdev(){

  // Marsaglia polar method

  double x,y;
  double s = 2;

  while(s>=1){
    x = drand48()*2-1;
    y = drand48()*2-1;
    s = x*x+y*y;
  }

  double val = x*sqrt(-2*log(s)/s);

  return val;

}
