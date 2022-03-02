#include "allvars_lmb.h"
#include "proto_lmb.h"
#include <math.h>

/////////////////////////////////////////////////////////////////////////
// THIS PROGRAM MAKES AN OPTICAL DEPTH MAP FROM A REIONIZATION REDSHIFT
// FIELD, ASSUMING THAT THE MATTER DENSITY IS UNIFORM.
//
//                                                AUTHOR: MARCELO ALVAREZ  
//                                             LAST EDIT:        01.25.13
/////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  
  // Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // Parse command line
  CommandLine(argc, &argv[0]);

  // Read parameter file
  ReadParameterFile();
  int N=Parameters.N;
  Parameters.BoxSize/=Parameters.h;
  nl=Parameters.nl;

  // Growth factor at initial time
  DInit=growth(Parameters.zInit,Parameters.Omegam,Parameters.Omegal,
	       Parameters.w);

  // Growth factor at z>>1
  float zhigh=100;
  float Dhigh=growth(zhigh,Parameters.Omegam,Parameters.Omegal,
	       Parameters.w);
  float Df = (1+zhigh)*Dhigh/DInit;

  // Make FFTW plans
  plan  = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,N, N, N,
				  FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
  iplan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,N, N, N,
				  FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
  
  rfftwnd_mpi_local_sizes(plan, &local_nx, &local_x_start,
			  &local_ny_after_transpose,
			  &local_y_start_after_transpose,
			  &total_local_size);

  // Allocate arrays
  AllocateArrays();

  // Read input files
  char ReionFile[256];
  sprintf(ReionFile,"%s.zreion",clParameters.BaseIn);

  ReadReionFile(ReionFile);
  ReadDeltaFile(Parameters.DeltaFile);

  GetMeanSigma(delta,1,N,Nlocal);
  GetMeanSigma(zreion,1,N,Nlocal);

  // FFT delta
  rfftwnd_mpi(plan, 1, delta, NULL, FFTW_NORMAL_ORDER);
  NormalizeArray(delta, N, size_fftw, 1.5);
  
  // Obtain power spectrum for density contrast at input time
  PowerSpectrum(delta,ps_delta,0);

  // Set history table
  Redshift2HistoryTable = new double[NZTABLE];
  SetRedshift2HistoryTable(Parameters.N, Nlocal, zreion, Redshift2HistoryTable);  

  // Loop over ionized fractions
  dz=-1e-5;
  zCurr = 20;
  float x = Redshift2Float((float)zCurr, Redshift2HistoryTable);
  float dx    = 0.01;
  float xNext = x + dx;
  float xmax  = 0.99;
  while(x < xmax){    
    zCurr += dz;
    x = Redshift2Float((float)zCurr, Redshift2HistoryTable);
    if(x < xNext) continue;    
    x = xNext;
    xNext += dx;
    xCurr = x;
    if(myid==0) printf("\n\n PS for z=%f x=%f --------- \n",zCurr,x);

    // Get ionization field at this redshift
    zreion2xion();
    
    // FFT xion
    rfftwnd_mpi(plan, 1, xion, NULL, FFTW_NORMAL_ORDER);
    NormalizeArray(xion, N, size_fftw, 1.5);

    // Get ionization field power spectrum
    PowerSpectrum(xion,ps_xx,0);

    // Get cross power spectrum
    CrossPower(xion,delta,ps_xdelta,0);

    // Write power spectrum
    WritePS();

    // Accumulate Cls
    AccumulateCl(1);

  }

  // WriteCl
  WriteCl();

  // Finalize and return

  if(myid==0) printf("\n\n");
  fclose(stdout);
  MPI_Finalize();    
  exit(0);

}

