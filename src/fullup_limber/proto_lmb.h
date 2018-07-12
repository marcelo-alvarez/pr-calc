#ifndef PROTO_H
#define PROTO_H

#include <stdio.h>
#include "cosmology.h"
#include "arrayoperations.h"
#include "tables.h"

#ifndef ALLVARS_H
#include "allvars_lmb.h"
#endif

// allocate.C
void    AllocateArrays();

// CommandLine.C
void usage();
void CommandLine(int, char **);

// Delta2Velocity.C
void Delta2Velocity();
void GenerateMomentum();
void SeparateMomentum();
void DifferentiateMomentum();
int  ijk2index(int, int, int);

// io.C
int  parallel_write(char *, long, long, void *);
int  parallel_read(char *, long, long, void *);
void ReadParameterFile();
void ReadReionFile(char *);
void ReadDeltaFile(char *);
void WritePS();
void WriteCl();

// MakeMap.C
void MakeMap();

// Tables.C
//void  SetRedshift2RadiusTable(float, float, float, double *);
//void  SetRedshift2HistoryTable(int, int, float *, double *);
//void  SetRedshift2DgTable(double *, double *);
//void  SetRedshift2TauTable(float, float, float, double *, double *);

//void  SetWavenumber2P1DTable(int, double *, double *, double *, double *);
//void  SetWavenumber2P3DTable(int, double *, double *, double *, double *, double *, double *);

//float Redshift2Float(float, double *);
//float Wavenumber2Float(int, float, double *, double *);

// PowerSpectrum.C
void PowerSpectrum(fftw_real *, double *, int);
void AccumulateCl(int);

#endif

