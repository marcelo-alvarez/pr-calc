#include <stdio.h>
#include "allvars_ics.h"

#include "globaltablevars.C"

// IO
FILE *logfile;

// MPI Variables
int myid, nproc;

// FFTW Variables
fftwf_plan plan, iplan; 
ptrdiff_t local_nx, local_x_start,total_local_size;

// Slab sizes
long Nlocal;           // Local slab dimension 
long int size;         // Local slab size 
long int size_fftw;    // Local slab size for FFTW arrays

// Arrays
float *delta;
fftwf_complex *cdelta;

// Filenames
char DeltaFile[256];

// Parameters
struct Parameter Parameters;

// Command line parameters
struct clParameter clParameters;

// Global variables
float ps_kmin, ps_kmax;
int ps_nk;

