#include "allvars_lmb.h"
#include "proto_lmb.h"
#include <math.h>

void AllocateArrays()
{

  // Set local slab size

  long N = Parameters.N;
  int nk = Parameters.nk;

  if(Parameters.N % nproc != 0){
    if(myid==0) printf(
	        "\n # of slices %d does not divide into # of procs %d, exiting\n\n",
		Parameters.N,nproc);
    MPI_Finalize();
    return;
  }
  else{
    Nlocal = N / nproc;
    size = Nlocal*N*N;
    size_fftw = Nlocal*N*(N+2);
  }

  zreion  = new float[size];
  delta   = new float[size_fftw];
  xion    = new float[size_fftw]();

  ps_delta  = new double[nk]();
  ps_xdelta = new double[nk]();
  ps_xx     = new double[nk]();

  inbin     = new int[nk]();

  ook2    = new double[nk]();
  kmean   = new double[nk]();

  if(myid==0) printf("\n Arrays allocated...");

}
