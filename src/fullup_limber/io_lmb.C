#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "proto_lmb.h"
#include "allvars_lmb.h"

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

  strcpy(tag[nt], "Sigma8");
  addr[nt] = &Parameters.Sigma8;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "w");
  addr[nt] = &Parameters.w;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "BoxSize");
  addr[nt] = &Parameters.BoxSize;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "N");
  addr[nt] = &Parameters.N;
  id[nt++] = INT;

  strcpy(tag[nt], "zInit");
  addr[nt] = &Parameters.zInit;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "DeltaFile");
  addr[nt] = &Parameters.DeltaFile;
  id[nt++] = STRING;

  strcpy(tag[nt], "kmin");
  addr[nt] = &Parameters.kmin;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "kmax");
  addr[nt] = &Parameters.kmax;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "nk");
  addr[nt] = &Parameters.nk;
  id[nt++] = INT;

  strcpy(tag[nt], "zmin");
  addr[nt] = &Parameters.zmin;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "zmax");
  addr[nt] = &Parameters.zmax;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "nz");
  addr[nt] = &Parameters.nz;
  id[nt++] = INT;

  strcpy(tag[nt], "nl");
  addr[nt] = &Parameters.nl;
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
      MPI_Finalize();
      exit(0);

    }
  
   if(myid==0) printf("\n Parameter file read...");

#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS
  
}

void ReadReionFile(char *fname)
{
    
  int N=Parameters.N;

  float *slab = new float[N*N];

  long offset = (long)myid*(long)N*(long)N*(long)Nlocal*(long)sizeof(float);
  long bsize = N*N*sizeof(float);

  FILE *fd = fopen(fname,"rb");

  for(long i=0;i<Nlocal;i++){
    
    // Read this slab from the input file

    long int offset_local = i*N*N*sizeof(float) + offset;
    parallel_read(fname, bsize, offset_local, slab);

    for(long j=0;j<N;j++){
      for(long k=0;k<N;k++){
	long index = i*N*N+j*N+k;
	zreion[index] = slab[j*N+k]; 
      }
    }

  }
  
  delete[] slab;

  if(myid==0) printf("\n ZReion data read...");

}

void ReadDeltaFile(char *fname)
{
    
  int N=Parameters.N;

  float *slab = new float[N*N];

  long offset = (long)myid*(long)N*(long)N*(long)Nlocal*(long)sizeof(float);
  long bsize = N*N*sizeof(float);

  FILE *fd = fopen(fname,"rb");
  if(fd==NULL){
    printf("could not open file %s\n",fname);
    exit(0);
    MPI_Finalize();
  }

  for(long i=0;i<Nlocal;i++){
    
    // Read this slab from the input file

    long int offset_local = i*N*N*sizeof(float) + offset;
    parallel_read(fname, bsize, offset_local, slab);

    for(long j=0;j<N;j++){
      for(long k=0;k<N;k++){
	long index = i*N*(N+2)+j*(N+2)+k;
	delta[index] = slab[j*N+k]; 
      }
    }

  }
  
  delete[] slab;

  if(myid==0) printf("\n Density data read...");

}

void WritePS()
{

  char fname[256];

  FILE *fout;

  float z=1/aCurr-1;
  float dzbox = HCurr * Parameters.BoxSize;

  if(myid==0) printf("\n dzbox = %e z+dz = %e",dzbox,z+dzbox);

  float s = Redshift2Float(z, Redshift2RadiusTable);

  sprintf(fname,"%s_%.2f.ps",clParameters.BaseOut,z);
  fout=fopen(fname, "w");

  fprintf(fout,"\n#Column 1: k [1/Mpc]");
  fprintf(fout,"\n#Column 2: p(k) [Mpc]^3");
  fprintf(fout,"\n#Column 3: p_prp(k) [Mpc]");
  fprintf(fout,"\n#Column 4: p_par(k) [Mpc]");
  fprintf(fout,"\n#Curl momentum field is defined as the part of [(1+delta)*xion*v]"); 
  fprintf(fout,"\n# whose Fourier components are perpendicular to k where");
  fprintf(fout,"\n# v=velocity in units of c; see Park et al. (2013) eq. 4 for definition");
  fprintf(fout,"\n");

  int nk=Parameters.nk;
  float kmin=Parameters.kmin;
  float kmax=Parameters.kmax;
  for(int i=0;i<nk;i++){

    float k=sqrt(1/ook2[i]);
    k=kmean[i];

    float H=100*Parameters.h*sqrt(Parameters.Omegam*pow((1+z),3)+1-Parameters.Omegam)/2.99e5;

    float ps_dlt = ps_delta[i];
    float ps_prp = ps_prpx[i]+ps_prpy[i]+ps_prpz[i];

    float dps_dlt = pow(k,3)*ps_dlt/2/M_PI/M_PI;
    float dps_prp = pow(k,3)*ps_prp/2/M_PI/M_PI*pow(2.99e5,2); // converted to (km/s)^2 
    float dps_par = pow(k,3)*ps_par[i]/2/M_PI/M_PI;
    float dps_prp_los = sqrt(dps_prp)*k/H;

    if(inbin[i]>0) fprintf(fout,"%e %e %e %e\n",k,ps_dlt,ps_prp,ps_par[i]);

  }

  if(myid==0) printf("\n Power spectrum written...");

  fclose(fout);

}

void WriteCl()
{

  char fname[256];

  FILE *fout;

  sprintf(fname,"%s_lmb.cl",clParameters.BaseOut);
  fout=fopen(fname, "w");

  /*
  fprintf(fout,"\n#Column 1: ell");
  fprintf(fout,"\n#Column 2: l^2Cl/(2pi) [muK^2]");
  fprintf(fout,"\n");
  */

  float dll = (log(lmax)-log(lmin))/nl;

  float zmin = Parameters.zmin;
  float s0 = Redshift2Float(zmin, Redshift2RadiusTable);

  for(int i=0;i<nl;i++){
    
    float l = exp(log(lmin)+(i+0.5)*dll);
    float k0 = l/s0;

    // Do boundary term

    double l2cl_bnd = clprefac_bnd / l *
      Wavenumber2Float(Parameters.nk, k0, Wavenumber2PdeltaTable, Wavenumbers4PdeltaTable)*s0*pow(1+zmin,3);

    l2cl_prp[i] *= clprefac_prp * pow(l,2);
    l2cl_par[i] *= clprefac_par ;
    l2cl_dop[i] *= clprefac_dop / pow(l,2);

    float l2cl = l2cl_prp[i]+l2cl_par[i];    

    fprintf(fout,"%e %le %le %le %le %le\n",l,l2cl,l2cl_prp[i],l2cl_par[i],l2cl_dop[i],l2cl_bnd);

  }

  if(myid==0) printf("\n l^2Cl/(2pi) written...");

  fclose(fout);

}
