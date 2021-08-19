#include "allvars_d2z.h"
#include "proto_d2z.h"
#include <math.h>

void AllocateArrays()
{

  // Set local slab size

  long N = Parameters.N;

  if(Parameters.N % nproc != 0){
    if(myid==0) printf(
	        "\n # of slices %d does not divide into # of procs %d, exiting\n\n",
		Parameters.N,nproc);
  }

  Nlocal = N / nproc;
  if(local_nx != Nlocal){
    if(myid==0) printf(
		"\n Nlocal = %ld not equal to local_nx = %ld for process %d",
		Nlocal,local_nx,myid);
    MPI_Finalize();
    exit(0);
  }

  cfftwa_r2c = fftwf_alloc_complex(total_local_size);
  fftwa_r2c  = (float *)cfftwa_r2c;
  cfftwa_c2r = cfftwa_r2c;
  
  if(clParameters.lean == 0) cfftwa_c2r = fftwf_alloc_complex(total_local_size);    

  fftwa_c2r = (float *)cfftwa_c2r;

  delta_smooth = fftwa_c2r;

  zreion = new float[size_fftw]();
  if(Parameters.RecordBubbleSizes == 1) rbubble = new float[size_fftw];

  if(myid==0) printf("\n Arrays allocated...");

  return;

}
