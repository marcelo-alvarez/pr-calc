#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "proto_asm.h"
#include "allvars_asm.h"

void Delta2Velocity()
{

  int i,j,k,r,index;

  // Local copies of parameters
  int N=Parameters.N;
  float BoxSize=Parameters.BoxSize;

  // Cell sizes and volumes in k space
  float dk=2.*M_PI/BoxSize;
  float d3k=pow(dk,3);
  float kmax=dk*N;

  // FFT of delta to k-space
  rfftwnd_mpi(plan, 1, delta, NULL, FFTW_NORMAL_ORDER);
  NormalizeArray(delta, N, size_fftw, 1.5);

  // Convert to velocity
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
	
	float rk2=kx2+ky2+kz2;
	float rk=sqrt(rk2);

	if(i+Nlocal*myid==N/2){
	  xdv_x[index]   = 0;
	  xdv_x[index+1] = 0;
	} else{
	  xdv_x[index]   = -kx/rk2*delta[index+1];
	  xdv_x[index+1] =  kx/rk2*delta[index];
	}

	if(j==N/2){
	  xdv_y[index]   = 0;
	  xdv_y[index+1] = 0;
	} else{
	  xdv_y[index]   = -ky/rk2*delta[index+1];
	  xdv_y[index+1] =  ky/rk2*delta[index];
	}

	if(k==N/2){
	  xdv_z[index]   = 0;
	  xdv_z[index+1] = 0;
	} else{
	  xdv_z[index]   = -kz/rk2*delta[index+1];
	  xdv_z[index+1] =  kz/rk2*delta[index];
	}

	if(rk2==0){
	  xdv_x[index]=0; xdv_x[index+1]=0;
	  xdv_y[index]=0; xdv_y[index+1]=0;
	  xdv_z[index]=0; xdv_z[index+1]=0;
	}

      }
    }
  }

  // Convert density back to real space
  rfftwnd_mpi(iplan, 1, delta, NULL, FFTW_NORMAL_ORDER);
  NormalizeArray(delta, N, size_fftw, 1.5);

  // Convert velocity back to real space
  rfftwnd_mpi(iplan, 1, xdv_x, NULL, FFTW_NORMAL_ORDER);
  rfftwnd_mpi(iplan, 1, xdv_y, NULL, FFTW_NORMAL_ORDER);
  rfftwnd_mpi(iplan, 1, xdv_z, NULL, FFTW_NORMAL_ORDER);
  NormalizeArray(xdv_x, N, size_fftw, 1.5);
  NormalizeArray(xdv_y, N, size_fftw, 1.5);
  NormalizeArray(xdv_z, N, size_fftw, 1.5);

  if(myid==0) printf("\n Velocities calculated...");
  
  return;
  
}

int ijk2index(int i, int j, int k){
   
   int index;
   long N=Parameters.N;

   index = i*N*(N+2)+j*(N+2)+k;

   return index;

}
