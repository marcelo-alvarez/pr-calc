#include "allvars_ics.h"
#include "proto_ics.h"
#include <math.h>

/////////////////////////////////////////////////////////////////////////
// THIS PROGRAM GENERATES UNIGRID INITIAL GAUSSIAN RANDOM FIELD
//
//                                                AUTHOR: MARCELO ALVAREZ  
//                                             LAST EDIT:        04.02.13
/////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  
  // Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // Parse command line  
  CommandLine(argc, &argv[0]);

  // Read parameter file
  ReadParameters();
  int N=clParameters.N;
  clParameters.BoxSize/=Parameters.h;

  // FFTW Initialization and plan creation
  fftwf_mpi_init();
  total_local_size = fftwf_mpi_local_size_3d(N,N,N/2+1,MPI_COMM_WORLD,
					    &local_nx,&local_x_start);
  // Allocate arrays
  AllocateArrays();

  // Make FFTW plans
  plan  = fftwf_mpi_plan_dft_r2c_3d(N, N, N, delta, cdelta, MPI_COMM_WORLD,
				   FFTW_ESTIMATE);
  iplan = fftwf_mpi_plan_dft_c2r_3d(N, N, N, cdelta, delta, MPI_COMM_WORLD,
				   FFTW_ESTIMATE);

  // Read power spectrum
  ReadPowerSpectrum();

  // Generate noise
  GenerateNoise();
  
  // Convolve noise with power spectrum
  ConvolveNoise();

  // Write output file
  WriteDelta();

  // Finalize and return
  MPI_Finalize();  if(myid==0) printf("\n\n"); 
  
  fclose(stdout);
  return 0;

}

