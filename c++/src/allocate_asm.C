#include "allvars_asm.h"
#include "proto_asm.h"
#include <math.h>

void AllocateArrays()
{

  // Set local slab size

  long N = Parameters.N;
  long NSide = Parameters.NSide;
  mapsize = nside2npix(NSide);

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

  int nmesh = 5;
  int nmaps = 3;

  float overhead_usage = (float)ovrt*(float)total_local_size;
  float mesh_usage = (float)total_local_size*nmesh*4.;
  float buffer_usage = (float)Nlocal*(float)N*N*4.*0;
  float map_usage = mapsize*nmaps*4.*2.; // factor of two is for local map buffers

  float estimated_usage = overhead_usage + mesh_usage + buffer_usage + map_usage;
  estimated_usage /= pow(1024.,3); // GB per process
  overhead_usage  /= pow(1024.,3);
  mesh_usage      /= pow(1024.,3);
  buffer_usage    /= pow(1024.,3);
  map_usage       /= pow(1024.,3);

  if(myid==0) printf( "\n  overhead usage is %f GB",overhead_usage);
  if(myid==0) printf( "\n      mesh usage is %f GB",mesh_usage);
  if(myid==0) printf( "\n    buffer usage is %f GB",buffer_usage);
  if(myid==0) printf( "\n       map usage is %f GB",map_usage);
  if(myid==0) printf( "\n estimated usage is %f GB per process",estimated_usage);

  zreion  = new float[size_fftw]();
  delta   = new float[size_fftw]();

  xdv_x   = new float[size_fftw]();
  xdv_y   = new float[size_fftw]();
  xdv_z   = new float[size_fftw]();

  taumap     = (float *)malloc(mapsize*sizeof(float));
  dtbmap     = (float *)malloc(mapsize*sizeof(float));
  kszmap     = (float *)malloc(mapsize*sizeof(float));
  kapmap     = (float *)malloc(mapsize*sizeof(float));

  if(myid==0) printf("\n Arrays allocated...");

}

