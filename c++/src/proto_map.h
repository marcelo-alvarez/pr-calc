#ifndef PROTO_H
#define PROTO_H

#include <stdio.h>
#include "cosmology.h"
#include "arrayoperations.h"

#ifndef ALLVARS_H
#include "allvars_map.h"
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
int  parallel_write(char *, long, long, void *);
int  parallel_read(char *, long, long, void *);
void ReadReionFile(char *);
void ReadDeltaFile(char *);
void ReadParameterFile();
void WriteMaps();

// makemaps.C
void MakeMaps();

#endif

