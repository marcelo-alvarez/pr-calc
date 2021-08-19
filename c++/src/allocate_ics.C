#include "allvars_ics.h"
#include "proto_ics.h"
#include <math.h>

void AllocateArrays()
{

  // Set local slab size

  long N = clParameters.N;

  if(N % nproc != 0){
    if(myid==0) printf(
	        "\n # of slices %ld does not divide into # of procs %d, exiting\n\n",
		N,nproc);
    MPI_Finalize();
    exit(0);
  }

  Nlocal = N / nproc;
  if(local_nx != Nlocal){
    if(myid==0) printf(
		"\n Nlocal = %ld not equal to local_nx = %ld for process %d",
		Nlocal,local_nx,myid);
    MPI_Finalize();
    exit(0);
  }

  size_fftw = 2*total_local_size;

  cdelta = fftwf_alloc_complex(total_local_size);
  delta = (float *)cdelta;

  if(myid==0) printf("\n Arrays allocated...");

  for(int i=0;i<size_fftw;i++) delta[i]=0;
  return;

}
