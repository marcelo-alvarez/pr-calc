#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "proto_ics.h"
#include "allvars_ics.h"
#include <iostream>
#include <fstream>
using namespace std;

#define FLOAT 1
#define DOUBLE 2
#define STRING 3
#define INT 4
#define MAXTAGS 300

void ReadParameters()
{

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  char fname[256];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];

  sprintf(fname,"%s",clParameters.ParameterFile);

  nt=0;
  
  strcpy(tag[nt], "h");
  addr[nt] = &Parameters.h;
  id[nt++] = FLOAT;

  if((fd = fopen(fname, "r")))
    {
      sprintf(buf, "%s%s", fname, "-usedvalues");
      if(!(fdout = fopen(buf, "w")))
	{
	  printf("error opening file '%s' \n", buf);
	}
      else
	{
	  while(!feof(fd))
	    {
	      *buf = 0;
	      fgets(buf, 200, fd);
	      if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		continue;
	      
	      if(buf1[0] == '%' || buf1[0] == '#')
		continue;
	      
	      for(i = 0, j = -1; i < nt; i++)
		if(strcmp(buf1, tag[i]) == 0)
		  {
		    j = i;
		    tag[i][0] = 0;
		    break;
		  }
	      
	      if(j >= 0)
		{
		  switch (id[j])
		    {
		    case FLOAT:
		      *((float *) addr[j]) = atof(buf2);
		      fprintf(fdout, "%-35s%f\n", buf1, *((float *) addr[j]));
		      break;
		    case DOUBLE:
		      *((float *) addr[j]) = atof(buf2);
		      fprintf(fdout, "%-35s%f\n", buf1, *((float *) addr[j]));
		      break;
		    case STRING:
		      strcpy((char *)addr[j], buf2);
		      fprintf(fdout, "%-35s%s\n", buf1, buf2);
		      break;
		    case INT:
		      *((int *) addr[j]) = atoi(buf2);
		      fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
		      break;
		    }
		}
	      else
		{
		  fprintf(stdout, 
			  "Error in %s: Tag '%s' not allowed or multiple defined.\n",
			  fname, buf1);
		}
	    }
	  fclose(fd);
	  fclose(fdout);
	  
	}
    }
  else
    {
      printf("\nParameter file %s not found.\n\n", fname);
    }
  
   if(myid==0) printf("\n Parameter file read...");

#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS
  
}

void ReadPowerSpectrum()
{
    
  ifstream fd;
  fd.open(clParameters.PowerSpectrumFile);
  if(!fd.is_open()){
    if(myid==0) printf("\n failed to open file %s, exiting...\n",clParameters.PowerSpectrumFile);
    MPI_Finalize();
    exit(0);
  }

  float dum;

  ps_nk=0;

  fd >> ps_kmin >> dum;
  while(!fd.eof()){
    fd >> ps_kmax >> dum;
    ps_nk++;
  }

  fd.close();

  Wavenumber2PdeltaTable = new double[ps_nk];

  fd.open(clParameters.PowerSpectrumFile);
  for(int i=0;i<ps_nk;i++) fd >> dum >> Wavenumber2PdeltaTable[i];
  for(int i=0;i<ps_nk;i++) Wavenumber2PdeltaTable[i] /= pow(2*M_PI,3); // Convert to spectral density
  fd.close();

  if(myid==0) printf("\n Input power spectrum read...");

}

void WriteDelta()
{

  char fname[256];

  sprintf(fname,"%s",DeltaFile);

  int N=clParameters.N;

  float *slab = new float[N*N];

  long offset = (long)myid*(long)N*(long)N*(long)Nlocal*(long)sizeof(float);
  long size = N*N*sizeof(float);

  for(int i=0;i<Nlocal;i++){
    
    for(int j=0;j<N;j++){
      for(int k=0;k<N;k++){
	long index = i*N*(N+2)+j*(N+2)+k;
	slab[j*N+k] = delta[index];
      }
    }

    // Write this slab to the output file
    
    long int offset_local = i*N*N*sizeof(float) + offset;
    parallel_write(fname, size, offset_local, slab);

  }
  
  if(myid==0) printf("\n Density contrast delta written...");

}

