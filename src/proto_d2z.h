#ifndef PROTO_H
#define PROTO_H

#include <stdio.h>

#ifndef ALLVARS_H
#include "allvars_d2z.h"
#endif

// allocate.C
void   AllocateArrays(); 

// bubblehistory.C
float GetHistory();
void  GetCHistory();
void  GetBubbleMFP();
float BubbleMFP(float );
float zreioncorrect();
void  CorrectHistory();

// commandline.C
void usage();
void CommandLine(int, char **);
float d2z(float, float);

// delta2zreion.C
void delta2zreion();

// io.C
int  parallel_write(char *, long, long, void *);
int  parallel_read(char *, long, long, void *);
void ReadDeltaFile(char *);
void ReadParameterFile();
void WriteZreionFile();
void WriteRbubbleFile();
void WriteHistory();
void WriteCHistory();
void WriteMFP();
void WriteMFPPDF();
void OutputMap(float *, int);

// math.C
float inverse_erfc(float);

// memorymanagement.C
int parseLine(char *);
int getValue();

#endif

