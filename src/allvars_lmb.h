#ifndef ALLVARS_H
#define ALLVARS_H

#include <mpi.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

#ifdef NOTYPEPREFIX_FFTW
#include        <rfftw_mpi.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h>	/* double precision FFTW */
#else
#include     <srfftw_mpi.h>
#endif
#endif

// Table Stuff
#include "tables.h"
#include "globaltablevars.h"

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

// Arrays
extern fftw_real *zreion,   *delta,   *prp_x,   *prp_y,  *prp_z,    *par, *paro;
extern fftw_real *xion;
extern double    *ps_prpx,  *ps_prpy, *ps_prpz, *ps_par, *ps_delta, *ook2, *kmean;
extern double    *ps_xdelta, *ps_xx;
extern double    *l2cl_par, *l2cl_prp, *l2cl_dop;
extern int       *inbin;

// Parameters
extern struct Parameter{
  float Omegam, Omegab, Omegal, h, ns, Sigma8, w;
  float BoxSize, zInit;
  int N, nk, nz, nl;
  float kmin,kmax,zmin,zmax;
  char DeltaFile[256];
} Parameters;

extern struct clParameter{
  int verbose, uniformu, uniformx, uniformd;
  char Paramfile[256], BaseIn[256], BaseOut[256];
} clParameters;

// Derived parameters
extern float DInit, DCurr, D0, dz, clprefac_prp, clprefac_par, clprefac_dop, clprefac_bnd, A;
extern double aCurr, HCurr, zCurr, xCurr;
extern float lmin,lmax;
extern int nl;

#endif
