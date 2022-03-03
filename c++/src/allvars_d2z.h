#define NHISTORY 100000
#define HISTORY_INITIAL 50.
#define HISTORY_FINAL   0.

#define NCHISTORY 10000
#define CHISTORY_INITIAL 50.
#define CHISTORY_FINAL   0.

#define NX_MFP 1000
#define X_MFP_INITIAL 0.01
#define X_MFP_FINAL   0.999

#define NMFP 200
#define MFP_INITIAL_CELL 1
#define MFP_FINAL_BOX    1.0

#define NSKEWERS 20000

#ifndef ALLVARS_H
#define ALLVARS_H

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <fftw3-mpi.h>

#include "cosmology.h"
#include "memorytracking.h"
#include "arrayoperations.h"
#include "tables.h"

// IO
extern FILE *logfile;

// MPI Variables
extern int myid, nproc;

// FFTW Variables
extern fftwf_plan plan, iplan; 
extern ptrdiff_t local_nx, local_x_start,total_local_size;

// Slab sizes
extern long Nlocal;           // Local slab dimension 
extern long int size;         // Local slab size 
extern long int size_fftw;    // Local slab size for FFTW arrays

// Arrays
extern fftwf_complex *cfftwa_r2c, *cfftwa_c2r;
extern float         *fftwa_r2c,  *fftwa_c2r;

extern float *zreion, *rbubble;
extern float *delta, *delta_smooth;
extern float *chistory, *history, *historyl, *mfp;
extern float *chistory_z, *history_z, *mfp_x, *dfdmfpl, *dfdmfp;

// Parameters
extern struct Parameter{
  float Omegam, Omegab, Omegal, h, ns, sigma8, w;
  float BoxSize, zInit, Rmin;
  int N, Nscales;
  char DeltaFile[256];
  int RecordBubbleSizes;
} Parameters;

extern struct clParameter{
  float Rmax, zeta, M_min;
  int verbose, everbose,lean;
  char Paramfile[256], Basename[256];
} clParameters;

// Tables
extern float dzhistory, dzchistory, dxmfp;

// Derived parameters
extern float DInit, sigma_min, erfizeta, zHigh, DHigh;

#endif
