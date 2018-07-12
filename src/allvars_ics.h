#ifndef ALLVARS_H
#define ALLVARS_H

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <fftw3-mpi.h>

#include "tables.h"
#include "globaltablevars.h"
#include "memorytracking.h"
#include "arrayoperations.h"

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
extern float *delta;
extern fftwf_complex *cdelta; 

// Filenames
extern char DeltaFile[256];

// Parameters
extern struct Parameter{
  float h;
} Parameters;

// Command line parameters
extern struct clParameter{
  float BoxSize;
  int verbose, Seed, N;
  char ParameterFile[256], PowerSpectrumFile[256], Base[256];
} clParameters;

// Global variables
extern float ps_kmin, ps_kmax;
extern int ps_nk;

#endif
