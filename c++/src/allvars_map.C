#include <stdio.h>
#include "allvars_map.h"

// Table Stuff
#include "globaltablevars.C"

// IO
FILE *logfile;

// MPI Variables
int myid, nproc;

// FFTW Variables
rfftwnd_mpi_plan plan, iplan; 
int local_nx, local_x_start, 
           local_ny_after_transpose, local_y_start_after_transpose, 
           total_local_size;

// Slab sizes
long Nlocal;        // Local slab dimension 
long int size;    // Local slab size 
long int size_fftw;
int mapsize;

// Arrays
fftw_real *delta, *xdv_x, *xdv_y, *xdv_z;
float *zreion, *taumap, *kszmap, *dtbmap;

// Parameters from parameter file
struct Parameter Parameters;

// Command line parameters
struct clParameter clParameters;

// Derived parameters
float DInit;
