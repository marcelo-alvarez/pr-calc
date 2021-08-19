#include <stdio.h>
#include "allvars_d2z.h"

// IO
FILE *logfile;

// MPI Variables
int myid, nproc;

// FFTW Variables
fftwf_plan plan, iplan; 
ptrdiff_t local_nx, local_x_start,total_local_size;

// Slab sizes
long Nlocal;        // Local slab dimension 
long int size;      // Local slab size 
long int size_fftw; // Local slab size for FFTW arrays

// Arrays
fftwf_complex *cfftwa_r2c, *cfftwa_c2r;
float         *fftwa_r2c,  *fftwa_c2r;

float *zreion, *rbubble;
float *delta, *delta_smooth;
float *chistory, *history, *historyl, *mfp;
float *chistory_z, *history_z, *mfp_x, *dfdmfpl, *dfdmfp;

// Parameters from parameter file
struct Parameter Parameters;

// Command line parameters
struct clParameter clParameters;

// Tables
float dzhistory, dzchistory, dxmfp;

// Derived parameters
float DInit, sigma_min, erfizeta, DHigh, zHigh;
