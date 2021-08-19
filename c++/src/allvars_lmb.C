#include <stdio.h>
#include "allvars_lmb.h"

// Table stuff
#include "globaltablevars.C"

// IO
FILE *logfile;

// MPI Variables
int myid, nproc;

// FFTW Variables
rfftwnd_mpi_plan plan, iplan; 
int local_nx, local_x_start, local_ny_after_transpose, 
    local_y_start_after_transpose, total_local_size;

// Slab sizes
long Nlocal;        // Local slab dimension 
long int size;      // Local slab size 
long int size_fftw; // Local slab size for FFTW arrays

// Arrays
fftw_real *zreion, *delta, *prp_x, *prp_y, *prp_z, *par, *paro;
fftw_real *xion;
double    *ps_prpx, *ps_prpy, *ps_prpz, *ps_par, *ps_delta, *ook2, *kmean;
double    *ps_xdelta, *ps_xx;
double    *l2cl_par, *l2cl_prp, *l2cl_dop;
int       *inbin;

// Parameters from parameter file
struct Parameter Parameters;

// Command line parameters
struct clParameter clParameters;

// Derived parameters
float DInit, DCurr, D0, dz, clprefac_prp, clprefac_par, clprefac_dop, clprefac_bnd, A;
double aCurr, HCurr, zCurr, xCurr;
float lmin,lmax;
int nl;
