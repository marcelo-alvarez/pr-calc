#ifndef PROTO_H
#define PROTO_H

#include <stdio.h>
#include "cosmology.h"
#include "geometry.h"
#include "arrayoperations.h"
#include "memorytracking.h"
#include "parallel_io.h"

#ifndef ALLVARS_H
#include "allvars_asm.h"
#endif

// arrayoperations.C
void AllocateArrays(); 

// commandline.C
void usage();
void CommandLine(int, char **);

// delta2velocity.C
void Delta2Velocity();
void GenerateMomentum();
void GenerateCurlMomentum();
void GenerateCurlMomentumRealSpace();
int  ijk2index(int, int, int);

// io.C
void ReadReionFile(char *);
void ReadDeltaFile(char *);
void ReadParameterFile();
void WriteMaps();

// makemaps.C
void MakeMapsRayTrace();
void MakeMapsBinCells();

#endif

