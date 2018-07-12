#include "allvars_map.h"
#include "proto_map.h"
#include <math.h>

void AllocateArrays()
{

  // Set local slab size

  long N = Parameters.N;
  mapsize = Parameters.NPixels*Parameters.NPixels;

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

  zreion  = new float[size_fftw]();
  delta   = new float[size_fftw]();

  xdv_x   = new float[size_fftw]();
  xdv_y   = new float[size_fftw]();
  xdv_z   = new float[size_fftw]();

  taumap     = new float[mapsize]();
  kszmap     = new float[mapsize]();
  dtbmap     = new float[mapsize]();

  if(myid==0) printf("\n Arrays allocated...");

}

