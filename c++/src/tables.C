#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "tables.h"
#include <mpi.h>

Float2FloatTable::Float2FloatTable(char* fname){
    
  FILE *fd = fopen(fname,"rb");
      
  fread(&minval, 4, 1, fd);
  fread(&maxval, 4, 1, fd);
  fread(&N,      4, 1, fd);
      
  delta = (maxval - minval) / N;
      
  table = new float[N];
      
  fread(table, 4, N, fd);
      
  fclose(fd);
      
} 
  
float Float2FloatTable::Float2Float(float inval){

  if(inval <= minval) return table[0];
  if(inval >= maxval) return table[N-1];
  
  int bin = int((inval - minval) / delta) ;
  float f = (inval - (minval + bin*delta)) / delta ;
  float outval = (1-f)*table[bin] + f*table[bin+1];
    
  return outval;
    
}

//////////////// Redshift to float //////////////////

float Redshift2Float(float redshift, double *table)
{

  float dztable = ((float)ZTABLE_FINAL - ZTABLE_INITIAL) / NZTABLE;

  int bin = (redshift - ZTABLE_INITIAL) / dztable ;
  float f = (redshift - (ZTABLE_INITIAL + bin*dztable)) / dztable;
  float value = (1-f)*table[bin] + f*table[bin+1];

  return value;

}

//////////////// Radius to float //////////////////

float Radius2Float(float radius, double *table)
{

  float drtable = ((float)RTABLE_FINAL - RTABLE_INITIAL) / NRTABLE;

  int bin = (radius - RTABLE_INITIAL) / drtable ;
  float f = (radius - (RTABLE_INITIAL + bin*drtable)) / drtable;
  float value = (1-f)*table[bin] + f*table[bin+1];

  return value;

}

//////////////// Wavenumber to float //////////////////

float Wavenumber2Float(int nk, float k, double *table, double *wntable){

  for(int i=0;i<nk-1;i++){
    if(k>=wntable[i] && k<=wntable[i+1]){
      float f = (log(k) - log(wntable[i]))/(log(wntable[i+1])-log(wntable[i]));
      float p = (1-f)*table[i] + f*table[i+1];
      return p;
    }
  }

  printf("\n Error: entry not found in wavenumber table for k = %e\n",k);
  for(int i=0;i<nk-1;i++) printf("%d %e %e\n",i,wntable[i],table[i]);

  MPI_Finalize();
  exit(0);

}

//////////////// Wavenumber to float (uniform log spacing in k) //////////////////

float Wavenumber2FloatLogSpace(int nk, float kmin, float kmax, float k, double *table){


  k    = log(k);
  kmin = log(kmin);
  kmax = log(kmax);

  float dk = (kmax-kmin)/nk;


  int bin = (k-kmin) / dk ;
  float f = (k - (kmin + bin*dk)) / dk;
  float value = (1-f)*table[bin] + f*table[bin+1];

  return value;

}

//////////////// Radius table //////////////////

void SetRedshift2RadiusTable(float h, float Omegam, float Omegal, double *table){

  double H0 = h*100.; // in km/sec/Mpc
  double c  = 2.99792458e5;      // in km/sec
  float r0  = c/H0;              // in Mpc 

  float dztable = ((float)ZTABLE_FINAL - ZTABLE_INITIAL) / NZTABLE;
  table[0]=0;

  for(int bin=1;bin<NZTABLE;bin++){
    double z  = ZTABLE_INITIAL + (bin-0.5)*dztable;
    double dr = r0 * dztable / sqrt(Omegam*pow((1+z),3)+Omegal);
    table[bin]=table[bin-1]+dr;
  }

}

//////////////// Redshift table //////////////////

void SetRadius2RedshiftTable(float h, float Omegam, float Omegal, double *table){

  double H0 = h*100.; // in km/sec/Mpc
  double c  = 2.99792458e5;      // in km/sec
  float r0  = c/H0;              // in Mpc 

  float drtable = ((float)RTABLE_FINAL - RTABLE_INITIAL) / NRTABLE;
  table[0]=0;

  for(int bin=1;bin<NRTABLE;bin++){
    double r  = RTABLE_INITIAL + (bin-0.5)*drtable;
    double z  = table[bin-1];
    double dz = drtable * sqrt(Omegam*pow((1+z),3)+Omegal) / r0;
    table[bin]=table[bin-1]+dz;
  }

}

//////////////// History table //////////////////

void SetRedshift2HistoryTable(int N, int Nlocal, float *zreion, double *table){

  float dztable = ((float)ZTABLE_FINAL - ZTABLE_INITIAL) / NZTABLE;

  double *history  = new double[NZTABLE]();

  for(int i=0;i<NZTABLE;i++) history[i]=0;

  for(int i=0;i<Nlocal;i++){
    for(int j=0;j<N;j++){
      for(int k=0;k<N;k++){
	
	int index=i*(N+2)*N+j*(N+2)+k;
	float zcur=zreion[index];
	int bin=(int)((zcur-ZTABLE_INITIAL)/dztable);
	if(bin>=0 && bin<NZTABLE) history[bin]++;

      }
    }
  }

  MPI_Allreduce(history, table, NZTABLE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  

  for(int i=0;i<NZTABLE;i++) table[i]/=pow(N,3);
  for(int i=NZTABLE-2;i>=0;i--) table[i]+=table[i+1];

}

//////////////// Dg table //////////////////

void SetRedshift2DgTable(double *history, double *table){

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  float dztable = ((float)ZTABLE_FINAL - ZTABLE_INITIAL) / NZTABLE;

  for(int i=1;i<NZTABLE-1;i++){
    float x1    = history[i-1];
    float x2    = history[i+1];
    float x     = history[i];
    double dxdz = (x2-x1)/2/dztable;
    float z     = ZTABLE_INITIAL + i*dztable; 
    table[i] = 1.5*pow(1+z,0.5)*x + pow(1+z,1.5)*dxdz;
  }

  table[0]=table[1];
  table[NZTABLE-1]=table[NZTABLE-2];

  return;

}

//////////////// Tau table //////////////////

void SetRedshift2TauTable(float h, float Omegab, float Omegam, double *table, double *Redshift2HistoryTable){

  float dztable = ((float)ZTABLE_FINAL - ZTABLE_INITIAL) / NZTABLE;

  float thompson = 6.65e-25;
  float YHe = 0.25;
  float hubble0 = 100./3.086e19*h;
  float ne0=Omegab*h*h*1.88e-29/1.67e-24;
  float c=3e10;
  
  float dtau0=thompson*ne0*c/hubble0;
  
  float NHe;
  for(int i=1;i<NZTABLE-1;i++){
    
    float z = ZTABLE_INITIAL + (i+0.5)*dztable;
    float x = Redshift2HistoryTable[i];

    NHe=1;
    if(z<3) NHe=2;

    float dtau = x*dtau0*pow((1+z),2)/sqrt(Omegam*pow((1+z),3)+1-Omegam)*
           (1-YHe+NHe*YHe/4)*dztable;

    table[i]+=table[i-1]+dtau;

  }

}

//////////////// Lensing kappa table //////////////////

void SetRedshift2WKappaTable(float h, float Omegam, float Omegal, double *table, double *Redshift2RadiusTable){

  // Get distance to LSS at z=1100
  double H0 = h*100.; // in km/sec/Mpc
  double  c = 2.99792458e5;      // in km/sec
  float  r0 = c/H0;              // in Mpc 
  int    nz = 100000; 
 
  float  dztable = 1100. / nz;
  float  chi0    = 0.;

  for(int bin=1;bin<nz;bin++){
    double z  = (bin-0.5)*dztable;
    double dr = r0 * dztable / sqrt(Omegam*pow((1+z),3)+Omegal);
    chi0 += dr;
  }

  // Make the table

  dztable = ((float)ZTABLE_FINAL - ZTABLE_INITIAL) / NZTABLE;

  for(int i=0;i<NZTABLE-1;i++){

    float z   = ZTABLE_INITIAL + (i+0.5)*dztable;
    float chi = Redshift2RadiusTable[i];
    float wkappa = 3*Omegam*(1+z)*chi*(1-chi/chi0)/2/r0/r0;
    float Hofz = H0*pow(Omegam*pow(1+z,3)+Omegal,0.5);
    table[i] = wkappa;

  }

} 

//////////////// P3D table //////////////////

void SetWavenumber2P3DTable(int nk, double *table, double *wntable, double *ps1, double *ps2, double *ps3, double *ook2){

  // Make the table

  for(int i=0;i<nk;i++){
    
    table[i]=ps1[i]+ps2[i]+ps3[i];
    wntable[i]=sqrt(1/ook2[i]);

  }

}

//////////////// P1D table //////////////////

void SetWavenumber2P1DTable(int nk, double *table, double *wntable, double *ps, double *ook2){

  // Make the table

  for(int i=0;i<nk;i++){
    
    table[i]=ps[i];
    wntable[i]=sqrt(1/ook2[i]);

  }

}

