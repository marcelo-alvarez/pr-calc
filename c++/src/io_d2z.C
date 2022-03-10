#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "proto_d2z.h"
#include "allvars_d2z.h"

#define FLOAT 1
#define DOUBLE 2
#define STRING 3
#define INT 4
#define MAXTAGS 300

void ReadParameterFile()
{

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  char fname[256];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];

  Parameters.RecordBubbleSizes = 0;

  sprintf(fname,"%s",clParameters.Paramfile);

  nt=0;
  
  strcpy(tag[nt], "Omegam");
  addr[nt] = &Parameters.Omegam;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Omegab");
  addr[nt] = &Parameters.Omegab;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Omegal");
  addr[nt] = &Parameters.Omegal;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "h");
  addr[nt] = &Parameters.h;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "ns");
  addr[nt] = &Parameters.ns;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "sigma8");
  addr[nt] = &Parameters.sigma8;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "w");
  addr[nt] = &Parameters.w;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "BoxSize");
  addr[nt] = &Parameters.BoxSize;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "zInit");
  addr[nt] = &Parameters.zInit;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Rmin");
  addr[nt] = &Parameters.Rmin;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "N");
  addr[nt] = &Parameters.N;
  id[nt++] = INT;

  strcpy(tag[nt], "Nscales");
  addr[nt] = &Parameters.Nscales;
  id[nt++] = INT;

  strcpy(tag[nt], "DeltaFile");
  addr[nt] = &Parameters.DeltaFile;
  id[nt++] = STRING;

  strcpy(tag[nt], "RecordBubbleSizes");
  addr[nt] = &Parameters.RecordBubbleSizes;
  id[nt++] = INT;

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

void ReadDeltaFile(char *fname)
{
    
  int N=Parameters.N;

  float *slab = new float[N*N];

  double mean=0, sigma=0, lmean=0, lsigma=0;
  
  long offset = (long)myid*(long)N*(long)N*(long)Nlocal*(long)sizeof(float);
  long size = N*N*sizeof(float);

  FILE *fd = fopen(fname,"rb");

  for(long i=0;i<Nlocal;i++){
    
    // Read this slab from the input file
    
    long int offset_local = i*N*N*sizeof(float) + offset;
    parallel_read(fname, size, offset_local, slab);

    for(long j=0;j<N;j++){
      for(long k=0;k<N;k++){
	long index = i*N*(N+2)+j*(N+2)+k;
	delta[index] = slab[j*N+k]; // Normalize here
	lmean+=delta[index];
	lsigma+=delta[index]*delta[index];
      }
    }
  }
  
  MPI_Allreduce(&lmean,  &mean,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&lsigma, &sigma, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  mean/= pow((float)N,3);

  sigma/= pow((float)N,3);
  sigma = sqrt(sigma);

  if(myid==0) printf("\n Input data read, mean = %le sigma = %le",mean,sigma);

  delete[] slab;

}

void WriteZreionFile()
{

  char fname[256];

  if(myid==0) printf("\n Writing reionization redshifts...");

  sprintf(fname,"%s.zreion",clParameters.Basename);

  int N=Parameters.N;

  float *slab = new float[N*N];

  long offset = (long)myid*(long)N*(long)N*(long)Nlocal*(long)sizeof(float);
  long size = N*N*sizeof(float);

  for(int i=0;i<Nlocal;i++){
    
    for(int j=0;j<N;j++){
      for(int k=0;k<N;k++){
	long index = i*N*(N+2)+j*(N+2)+k;
	slab[j*N+k] = zreion[index];
      }
    }

    // Write this slab to the output file
    
    long int offset_local = i*N*N*sizeof(float) + offset;
    parallel_write(fname, size, offset_local, slab);

  }
  
  if(myid==0) printf("\n Reionization redshifts written...");

  delete[] slab;

}


void WriteRbubbleFile()
{

  char fname[256];

  if(myid==0) printf("\n Writing bubble sizes...");

  sprintf(fname,"%s.rbubble",clParameters.Basename);

  int N=Parameters.N;

  float *slab = new float[N*N];

  long offset = (long)myid*(long)N*(long)N*(long)Nlocal*(long)sizeof(float);
  long size = N*N*sizeof(float);

  for(int i=0;i<Nlocal;i++){
    
    for(int j=0;j<N;j++){
      for(int k=0;k<N;k++){
	long index = i*N*(N+2)+j*(N+2)+k;
	slab[j*N+k] = rbubble[index];
      }
    }

    // Write this slab to the output file
    
    long int offset_local = i*N*N*sizeof(float) + offset;
    parallel_write(fname, size, offset_local, slab);

  }
  
  if(myid==0) printf("\n Bubble sizes written...");

}

void WriteHistory()
{

  FILE *fd;
  char fname[256];

  if(myid==0){
    printf("\n Writing ionization history...");

    sprintf(fname,"%s.history",clParameters.Basename);

    if(!(fd = fopen(fname, "w"))){
      if(myid==0) printf("\n File %s could not be opened\n\n", fname);
      MPI_Finalize();
      return;
    }

    Float2FloatTable ExternalFcoll("fcoll_table.tab");
    float zInit  = Parameters.zInit;
    float omegam = Parameters.Omegam;
    float omegal = Parameters.Omegal;
    float w      = Parameters.w;
    float smin   = sigma_min;
    for(int i=0;i<NHISTORY;i++) fprintf(fd,"%e %e %e %e\n",history_z[i],history[i],ExternalFcoll.Float2Float(history_z[i]),fcoll(history_z[i], zInit, omegam, omegal, w, smin));

    fclose(fd);

  }


}


void WriteCHistory()
{

  FILE *fd;
  char fname[256];

  if(myid==0){
    printf("\n Writing correct ionization history...");

    sprintf(fname,"%s.chistory",clParameters.Basename);
    
    if(!(fd = fopen(fname, "w"))){
      if(myid==0) printf("\n File %s could not be opened\n\n", fname);
      MPI_Finalize();
      return;
    }
    
    for(int i=0;i<NCHISTORY;i++) fprintf(fd,"%f %f\n",chistory_z[i],chistory[i]);

    fclose(fd);

  }


}

void WriteMFP()
{

  FILE *fd;
  char fname[256];

  if(myid==0){
    printf("\n Writing bubble mfp...");

    sprintf(fname,"%s.mfp",clParameters.Basename);
    
    if(!(fd = fopen(fname, "w"))){
      if(myid==0) printf("\n File %s could not be opened\n\n", fname);
      MPI_Finalize();
      return;
    }
    
    for(int i=0;i<NX_MFP;i++) fprintf(fd,"%f %f\n",mfp_x[i],mfp[i]);
    
  }

}

void WriteMFPPDF()
{

  FILE *fd;
  char fname[256];

  if(myid==0){
    printf("\n Writing bubble mfp...");

    sprintf(fname,"%s.mfppdf",clParameters.Basename);
    
    if(!(fd = fopen(fname, "wb"))){
      if(myid==0) printf("\n File %s could not be opened\n\n", fname);
      MPI_Finalize();
      return;
    }
    
    int         N  = Parameters.N;
    float BoxSize  = Parameters.BoxSize;
    float CellSize = BoxSize / N;

    int   *header1 = new int[2]();
    float *header2 = new float[4]();
    header1[0] = NX_MFP;
    header1[1] = NMFP;
    header2[0] = X_MFP_INITIAL;
    header2[1] = X_MFP_FINAL;
    header2[2] = log10(MFP_INITIAL_CELL * CellSize);
    header2[3] = log10(MFP_FINAL_BOX    *  BoxSize);
    
    fwrite(header1,sizeof(int),2,fd);
    fwrite(header2,sizeof(float),4,fd);
    fwrite(dfdmfp,sizeof(float),NX_MFP*NMFP,fd);

  }

}

void OutputMap(float *array, int snap)
{

  FILE *fd;
  char str[100];

  sprintf(str,"%03d.map",snap);

  fd=fopen(str,"wb");

  int N=Parameters.N;
  int i=Nlocal/2;

  for(int j=0;j<N;j++){
    for(int k=0;k<N;k++){
      int index = i*N*(N+2)+j*(N+2) + k;
      fwrite(&array[index],1,sizeof(float),fd);
    }
  }

  fclose(fd);

}

