#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "proto_lmb.h"
#include "allvars_lmb.h"

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
	  prp_x[index]   = 0;
	  prp_x[index+1] = 0;
	} else{
	  prp_x[index] = -HCurr*aCurr*kx*DCurr/DInit*delta[index+1]/rk2;
	  prp_x[index+1] =  HCurr*aCurr*kx*DCurr/DInit*delta[index]  /rk2;
	}

	if(j==N/2){
	  prp_y[index]   = 0;
	  prp_y[index+1] = 0;
	} else{
	  prp_y[index] = -HCurr*aCurr*ky*DCurr/DInit*delta[index+1]/rk2;
	  prp_y[index+1] =  HCurr*aCurr*ky*DCurr/DInit*delta[index]  /rk2;
	}

	if(k==N/2){
	  prp_z[index]   = 0;
	  prp_z[index+1] = 0;
	} else{
	  prp_z[index] = -HCurr*aCurr*kz*DCurr/DInit*delta[index+1]/rk2;
	  prp_z[index+1] =  HCurr*aCurr*kz*DCurr/DInit*delta[index]  /rk2;
	}

	if(rk2==0){
	  prp_x[index]=0; prp_x[index+1]=0;
	  prp_y[index]=0; prp_y[index+1]=0;
	  prp_z[index]=0; prp_z[index+1]=0;
	}

      }
    }
  }

  // Convert density back to real space
  rfftwnd_mpi(iplan, 1, delta, NULL, FFTW_NORMAL_ORDER);
  NormalizeArray(delta, N, size_fftw, 1.5);

  // Convert velocity back to real space
  rfftwnd_mpi(iplan, 1, prp_x, NULL, FFTW_NORMAL_ORDER);
  rfftwnd_mpi(iplan, 1, prp_y, NULL, FFTW_NORMAL_ORDER);
  rfftwnd_mpi(iplan, 1, prp_z, NULL, FFTW_NORMAL_ORDER);
  NormalizeArray(prp_x, N, size_fftw, 1.5);
  NormalizeArray(prp_y, N, size_fftw, 1.5);
  NormalizeArray(prp_z, N, size_fftw, 1.5);

  if(myid==0){
    long index = i*N*(N+2)+j*(N+2)+k;
    FILE *frho=fopen("delta.bin","wb");
    FILE *fvel=fopen("velocities.bin","wb");
    for(long i=0;i<Nlocal;i++){for(long j=0;j<N;j++){for(long k=0;k<N;k++){
	  long index = i*N*(N+2)+j*(N+2)+k;
	  float val=delta[index]*DCurr/DInit;
	  fwrite(&val,sizeof(float),1,frho);
	}}}
    double vrms=0;
    for(long i=0;i<Nlocal;i++){for(long j=0;j<N;j++){for(long k=0;k<N;k++){
	  long index = i*N*(N+2)+j*(N+2)+k;
	  fwrite(&prp_x[index],sizeof(float),1,fvel);
	  vrms+=prp_x[index]*prp_x[index];
	}}}
    for(long i=0;i<Nlocal;i++){for(long j=0;j<N;j++){for(long k=0;k<N;k++){
	  long index = i*N*(N+2)+j*(N+2)+k;
	  fwrite(&prp_y[index],sizeof(float),1,fvel);
	  vrms+=prp_y[index]*prp_y[index];
	}}}
    for(long i=0;i<Nlocal;i++){for(long j=0;j<N;j++){for(long k=0;k<N;k++){
	  long index = i*N*(N+2)+j*(N+2)+k;
	  fwrite(&prp_z[index],sizeof(float),1,fvel);
	  vrms+=prp_z[index]*prp_z[index];
	}}}
    fclose(frho);
    fclose(fvel);
    printf(" vrms at output is %lf\n",sqrt(vrms/pow(N,3))*2.99e5);
  }

  if(myid==0) printf("\n Velocities calculated...");
  
  return;
  
}

void GenerateMomentum(){

  // Generate 'specific momentum' field from density, zreion, and velocity field
  // Note that xdv_[x,y,z] are currently the proper peculiar velocity in units of
  // the speed of light

  char fname[256];
  float *tmap, *ltmap;

  int i,j,k,r,index;

  float z=1./aCurr-1.;
  long N=Parameters.N;
  float dl=Parameters.BoxSize/N;

  float s = Redshift2Float(z, Redshift2RadiusTable);
  float fov = Parameters.BoxSize / s;

  //  sprintf(fname,"%s_%.2f.tmap",clParameters.BaseOut,z);
  //  FILE *fout=fopen(fname, "wb");  

  tmap = new float[N*N]();
  ltmap = new float[N*N]();

  double vrmsl = 0;
  double vrms;

  for(long i=0;i<Nlocal;i++){    
    for(long j=0;j<N;j++){
      for(long k=0;k<N;k++){
	long index = i*N*(N+2)+j*(N+2)+k;
	long index_zr = i*N*N+j*N+k;

	float d=DCurr/DInit*delta[index];

	float zr=zreion[index_zr];

	float vx=prp_x[index];
	float vy=prp_y[index];
	float vz=prp_z[index];

	float xion=0;
	if(zr>z) xion=1; // ne = xion * ne0 where ne0=nH+nHe (assuming once ionized helium)
	if(clParameters.uniformu==1){
	  xion=Redshift2Float(z, Redshift2HistoryTable);
	  d=0;
	}
	else if(clParameters.uniformd==1){
	  d=0;
	}
	else if(clParameters.uniformx==1){
	  xion=Redshift2Float(z, Redshift2HistoryTable);
	}

	prp_x[index]*=xion*(1+d);
	prp_y[index]*=xion*(1+d);
	prp_z[index]*=xion*(1+d);

	prp_x[index]*=xion;
	prp_y[index]*=xion;
	prp_z[index]*=xion;

	vrmsl += pow(prp_x[index],2)+pow(prp_y[index],2)+pow(prp_z[index],2);

	ltmap[j*N+k]+=A*dl*prp_x[index];

      }
    }
  } 

  

  // sum map over all processors ccumulate

  MPI_Allreduce(ltmap, tmap, N*N, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&vrmsl, &vrms, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(myid == 0) printf(" rms of velocity is %lf\n",sqrt(vrms/N/N/N)*2.99e5);

  /*
  if(myid==0){

    fwrite(&N,sizeof(int),1,fout);
    fwrite(&N,sizeof(int),1,fout);

    fwrite(&fov,sizeof(float),1,fout);
    fwrite(&fov,sizeof(float),1,fout);

    fwrite(tmap, sizeof(float), N*N, fout);

  }

  fclose(fout);
  */

  delete [] tmap;
  delete [] ltmap;

  if(myid==0) printf("\n Specific momentum generated...");

}

void SeparateMomentum(){

  int i,j,k,r,index;  

  // Local copies of parameters
  long        N = Parameters.N;
  float BoxSize = Parameters.BoxSize;
  float       z = 1/aCurr - 1;

  // Cell sizes and volumes in k space
  float dk   = 2.*M_PI/BoxSize;
  float kmax = dk*N;

  // Fourier transform momentum field
  rfftwnd_mpi(plan, 1, prp_x, NULL, FFTW_NORMAL_ORDER);
  rfftwnd_mpi(plan, 1, prp_y, NULL, FFTW_NORMAL_ORDER);
  rfftwnd_mpi(plan, 1, prp_z, NULL, FFTW_NORMAL_ORDER);
  NormalizeArray(prp_x, N, size_fftw, 1.5);
  NormalizeArray(prp_y, N, size_fftw, 1.5);
  NormalizeArray(prp_z, N, size_fftw, 1.5);

  float opz2 = pow(1+z,2);

  // Loop over Fourier space and subtract out parallel part
  
  for(int i=0;i<Nlocal;i++){
    float kx=(i+Nlocal*myid)*dk;
    if(i+Nlocal*myid>=N/2) kx=kx-kmax;
    float kx2=kx*kx;
    
    for(int j=0;j<N;j++){
      float ky=j*dk;
      if(j>=N/2) ky=ky-kmax;
      float ky2=ky*ky;
      
      for(int k=0;k<=N/2;k++){
	float kz=k*dk;
	if(k>=N/2) kz=kz-kmax;
	float kz2=kz*kz;
	
	int index=i*N*(N+2)+j*(N+2)+2*k;
	
	float rk2=kx2+ky2+kz2;
	float rk=sqrt(rk2);

	float kxh = kx/rk;
	float kyh = ky/rk;
	float kzh = kz/rk;
	
	par[index]   = prp_x[index]  *kxh + prp_y[index]  *kyh + prp_z[index]  *kzh;
	par[index+1] = prp_x[index+1]*kxh + prp_y[index+1]*kyh + prp_z[index+1]*kzh;

	prp_x[index]   = prp_x[index]   - kxh * par[index];
	prp_y[index]   = prp_y[index]   - kyh * par[index];
	prp_z[index]   = prp_z[index]   - kzh * par[index];

	prp_x[index+1] = prp_x[index+1] - kxh * par[index+1];
	prp_y[index+1] = prp_y[index+1] - kyh * par[index+1];
	prp_z[index+1] = prp_z[index+1] - kzh * par[index+1];


	par[index]   *= opz2;
	par[index+1] *= opz2;

	if(rk2==0){
	  prp_x[index] = 0; prp_x[index+1] = 0;
	  prp_y[index] = 0; prp_y[index+1] = 0;
	  prp_z[index] = 0; prp_z[index+1] = 0;
	  par[index]   = 0; par[index+1]   = 0;
	}

      }
    }
  }

  if(myid==0) printf("\n Momentum separated...");

}

void DifferentiateMomentum()
{

  // Local copies of parameters
  dz = (Parameters.zmax - Parameters.zmin)/Parameters.nz;

  for(int i=0;i<size_fftw;i++){

    float old_value = paro[i];
    float new_value = par[i];

    float derivative = (new_value - old_value)/dz;

    paro[i] = new_value;
    par[i]  = derivative;

  }

}

int ijk2index(int i, int j, int k){
   
   int index;
   long N=Parameters.N;

   if(i<0) i+=N; if(i>=N) i-=N;
   if(j<0) j+=N; if(j>=N) j-=N;
   if(k<0) k+=N; if(k>=N) k-=N;

   index = i*N*(N+2)+j*(N+2)+k;

   return index;

 }
