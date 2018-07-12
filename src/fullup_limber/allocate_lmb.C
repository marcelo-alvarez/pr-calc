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

  prp_x   = new float[size_fftw]();
  prp_y   = new float[size_fftw]();
  prp_z   = new float[size_fftw]();
  par     = new float[size_fftw]();
  paro    = new float[size_fftw]();

  ps_prpx = new double[nk]();
  ps_prpy = new double[nk]();
  ps_prpz = new double[nk]();
  ps_par  = new double[nk]();

  ps_delta = new double[nk]();

  inbin   = new int[nk]();

  ook2    = new double[nk]();
  kmean   = new double[nk]();

  l2cl_prp = new double[nl]();
  l2cl_par = new double[nl]();
  l2cl_dop = new double[nl]();

  if(myid==0) printf("\n Arrays allocated...");

}
