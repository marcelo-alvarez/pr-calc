#ifndef PROTO_H
#define PROTO_H

#include <stdio.h>
#include "cosmology.h"
#include "tables.h"

#ifndef ALLVARS_H
#include "allvars_ics.h"
#endif

// allocate.C
void   AllocateArrays(); 

// commandline.C
void usage();
void CommandLine(int, char **);

// grf.C
void GenerateNoise();
void ConvolveNoise();
double gaussdev();

// io.C
void ReadParameters();
void ReadPowerSpectrum();
void WriteDelta();
int  parallel_write(char *, long, long, void *);
int  parallel_read(char *, long, long, void *);

#endif

