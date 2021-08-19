#ifndef ALLVARS_H
#define ALLVARS_H

#include <mpi.h>
/*
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

#ifdef NOTYPEPREFIX_FFTW
#include        <rfftw_mpi.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h>	
#else
#include     <srfftw_mpi.h>
#endif
#endif

#include "tables.h"
#include "globaltablevars.h"
#include "chealpix.h"

// IO
extern FILE *logfile;

// MPI Variables
extern int myid, nproc;

// FFTW Variables
extern rfftwnd_mpi_plan plan, iplan; 
extern int local_nx, local_x_start, 
           local_ny_after_transpose, local_y_start_after_transpose, 
           total_local_size;

// Slab sizes
extern long Nlocal;             // Local slab dimension 
extern long int size;           // Local slab size 
extern long int size_fftw;
extern long mapsize;

// Arrays
extern fftw_real *delta, *xdv_x, *xdv_y, *xdv_z;
extern float *zreion, *taumap, *kszmap, *dtbmap;

// Parameters
extern struct Parameter{
  float Omegam, Omegab, Omegal, h, ns, Sigma8, w;
  float BoxSize, zInit, fov, theta, phi;
  float InitialRedshift, FinalRedshift;  
  int N, NSide, NPixels, NRedshifts;
  char DeltaFile[256];
} Parameters;

extern struct clParameter{
  int verbose, uniform;
  char Paramfile[256], BaseIn[256], BaseOut[256];
} clParameters;

// Redshift table
extern float dztable;

// Derived parameters
extern float DInit, ovrt, oram;
*/
#endif
