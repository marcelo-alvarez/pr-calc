#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "proto_lmb.h"
#include "allvars_lmb.h"

void PowerSpectrum(fftw_real *field, double *ps, int transform)
{

  double lsigma=0,sigma=0;

  // Local copies of parameters
  int N=Parameters.N;
  float BoxSize=Parameters.BoxSize;

  // Cell sizes and volumes in k space
  float dk=2.*M_PI/BoxSize;
  float kmax=dk*N;
  float d3k=pow(kmax,3); // Total volume in k-space

  // Binning for the power spectrum
  int nk_ps=Parameters.nk;
  float kmin_ps=Parameters.kmin;
  float kmax_ps=Parameters.kmax;
  float dlnk_ps=(log(kmax_ps)-log(kmin_ps))/nk_ps;

  // Local power spectrum arrays
  int    *linbin = new int[nk_ps]();
  double *lps    = new double[nk_ps]();
  double *look2  = new double[nk_ps]();
  double *lkmean = new double[nk_ps]();

  // Initial FFT of field to k-space if transform = 1
  if(transform==1){
    rfftwnd_mpi(plan, 1, field, NULL, FFTW_NORMAL_ORDER);
    NormalizeArray(field, N, size_fftw, 1.5);
  }

  // Initialize power spectrum arrays
  for(int i=0;i<nk_ps;i++){
    linbin[i]=0;
    lps[i]=0;
    look2[i]=0;
    lkmean[i]=0;

    inbin[i]=0;
    ps[i]=0;
    ook2[i]=0;
    kmean[i]=0;
  }

  // Bin up the power spectrum
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

	int bin = (log(rk)-log(kmin_ps))/dlnk_ps;

	float power = pow(field[index],2)+pow(field[index+1],2); // variance of mode 
	lsigma += 2*power;
	if(bin>=0 && bin<nk_ps){
	  linbin[bin]++;
	  lps[bin]+= power/d3k; // variance of mode per dk^3
	  look2[bin]+=1/rk2;
	  lkmean[bin]+=sqrt(rk2);
	}

      }
    }
  }
    
  // Sum up contributions from all processors
  MPI_Allreduce(linbin, inbin, nk_ps, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(lps, ps, nk_ps, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(look2, ook2, nk_ps, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(lkmean, kmean, nk_ps, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&lsigma, &sigma, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  sigma/=pow(N,3);

  // ---------------------------------------------------------------------
  // IMPORTANT!!! IN COSMOLOGICAL CONVENTION, 
  //   variance per dlnk = Delta(k) = k^3*p_cosmo(k)/(2pi^2).
  // HOWEVER, THE CURRENT DEFINTION OF p(k) is 
  //   p(k) = variance per dk^3 (SEE ABOVE LOOP),
  // WHICH LEADS TO THE EXPRESSION
  //   variance per dlnk = Delta(k) = 4pi*k^3*p(k).
  // SO TO GET THE "RIGHT ONE", p_cosmo(k), WE MUST MULTIPLY BY (2*pi)^3,
  //   p_cosmo(k) = (2*pi)^3 * p(k),
  // TO COMPENSATE.
  // ---------------------------------------------------------------------

  for(int i=0; i<nk_ps; i++) ook2[i] /= inbin[i];
  for(int i=0; i<nk_ps; i++) kmean[i] /= inbin[i];
  for(int i=0; i<nk_ps; i++) ps[i] *= pow((2*M_PI),3)/inbin[i];

  float ps_sigma=0;
  for(int i=1; i<nk_ps; i++){
    if(inbin[i]>0 && inbin[i-1]>0){
      float kmx = sqrt(1/ook2[i]);
      float kmn = sqrt(1/ook2[i-1]);
      float pscur = (ps[i]+ps[i-1])/2;
      float kvol = 4*M_PI/3*(pow(kmx,3)-pow(kmn,3))/pow((2*M_PI),3);
      ps_sigma += kvol*pscur;
    }
  }

  // FFT of field back to real-space
  if(transform==1){
    rfftwnd_mpi(iplan, 1, field, NULL, FFTW_NORMAL_ORDER);
    NormalizeArray(field, N, size_fftw, 1.5);
  }

  delete[] linbin;
  delete[] lps;
  delete[] look2;

  return;
  
}

void AccumulateCl(int DoParallelComponent){

  int   nz   = Parameters.nz;
  float zmin = Parameters.zmin;
  float zmax = Parameters.zmax;
  float dz   = (zmax-zmin)/nz;

  float z = 1/aCurr -1;
  float zp = z - dz/2;

  float dll = (log(lmax)-log(lmin))/nl;
  float s  = Redshift2Float(z,  Redshift2RadiusTable);
  float sp = Redshift2Float(zp, Redshift2RadiusTable);

  float z1 = z - dz/2;
  float z2 = z + dz/2;
  float d1 = Redshift2Float(z1, Redshift2HistoryTable)*pow(1+z1,1.5);
  float d2 = Redshift2Float(z2, Redshift2HistoryTable)*pow(1+z2,1.5);
  float Dg = (d2-d1)/dz;

  float doppler_factor = pow(1+zp,1.5)*pow(sp,2)*Dg*Dg;
  if(myid==0) printf("\n%e %e dopplerfac\n",zp,doppler_factor);

  for(int i=0;i<nl;i++){

    float l  = exp(log(lmin)+(i+0.5)*dll);
    float k  = l/s;
    float kp = l/sp;

    // Perpendicular component
    double dl2cl_prp = pow(1+z,2.5)/pow(s,2)*
      Wavenumber2Float(Parameters.nk, k,  Wavenumber2PprpTable, Wavenumbers4PprpTable);
      
    // Parallel component
    double dl2cl_par = pow(1+z,1.5)*
      Wavenumber2Float(Parameters.nk, kp, Wavenumber2PparTable, Wavenumbers4PparTable);

    // 1st order Parallel (i.e. velocity only)
    double dl2cl_dop = pow(1+z,1.5)*pow(s,2)*
      Wavenumber2Float(Parameters.nk, kp, Wavenumber2PdeltaTable, Wavenumbers4PdeltaTable)*Dg*Dg;
    //      pow(Redshift2Float(z, Redshift2DgTable),2);

    if(DoParallelComponent==0) dl2cl_par = 0;

    l2cl_prp[i] += dl2cl_prp * dz;
    l2cl_par[i] += dl2cl_par * dz;
    l2cl_dop[i] += dl2cl_dop * dz;

  }

  if(myid==0) printf("\n l^2Cl/(2pi) accumulated...");

}

