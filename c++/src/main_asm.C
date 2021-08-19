#include "allvars_asm.h"
#include "proto_asm.h"
#include <math.h>

/////////////////////////////////////////////////////////////////////////
// THIS PROGRAM MAKES A MAP FROM INPUT REIONIZATION REDSHIFT AND LINEAR
// DENSITY FIELDS
//
//                                                AUTHOR: MARCELO ALVAREZ  
//                                             LAST EDIT:        10.10.16
/////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  
  // Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // Parse command line  
  CommandLine(argc, &argv[0]);

  // Read parameter file
  ReadParameterFile();
  int N=Parameters.N;
  Parameters.BoxSize/=Parameters.h;

  // Growth factor at initial time
  DInit=growth(Parameters.zInit,Parameters.Omegam,Parameters.Omegal,
	       Parameters.w);

  // Make FFTW plans
  plan  = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,N, N, N,
				  FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
  iplan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,N, N, N,
				  FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
  
  rfftwnd_mpi_local_sizes(plan, &local_nx, &local_x_start,
			  &local_ny_after_transpose,
			  &local_y_start_after_transpose,
			  &total_local_size);

  ovrt=0;
  oram=0;
  GetOverhead(total_local_size,&ovrt,&oram);
  ReportMemory("after plan creation",total_local_size,ovrt,oram);

  // Allocate arrays
  AllocateArrays();

  ReportMemory("after array allocation",size_fftw,ovrt,oram);

  // Read input files
  char ReionFile[256];
  sprintf(ReionFile,"%s.zreion",clParameters.BaseIn);

  ReadReionFile(ReionFile);
  ReadDeltaFile(Parameters.DeltaFile);

  // Set history table
  Redshift2HistoryTable = new double[NZTABLE];
  SetRedshift2HistoryTable(Parameters.N, Nlocal, zreion, Redshift2HistoryTable);

  // Use linear theory to go from density to velocity
  Delta2Velocity();
  
  // Make the optical depth map
  MakeMapsRayTrace();
  
  // Write map
  WriteMaps();

  // Finalize and return
  MPI_Finalize();  if(myid==0) printf("\n\n"); 
  
  fclose(stdout);
  exit(0);

}

