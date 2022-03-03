#include "allvars_d2z.h"
#include "proto_d2z.h"
#include <math.h>

/////////////////////////////////////////////////////////////////////////
// THIS PROGRAM CALCULATES THE REIONIZATION REDSHIFT FROM AN INPUT MATTER
// DENSITY FIELD USING A FILTERING TECHNIQUE. THE INPUT IS MATTER DENSITY
// ON A MESH, THE OUTPUT IS THE REIONIZATION REDSHIFT ON THE SAME MESH
//
//
//                                                AUTHOR: MARCELO ALVAREZ  
//                                             LAST EDIT:        07.23.12
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

  // Create Float2FloatTable object for external fcoll 
  Float2FloatTable ExternalFcoll("fcoll_table.tab");
  
  // FFTW Initialization and local sizes
  fftwf_mpi_init();
  total_local_size = fftwf_mpi_local_size_3d(N,N,N/2+1,MPI_COMM_WORLD,
					    &local_nx,&local_x_start);
  size_fftw = 2*total_local_size;

  float ovrt=0, oram=0;
  GetOverhead((float)size_fftw,&ovrt,&oram);
  ReportMemory("after plan creation",(float)size_fftw,ovrt,oram);

  // Allocate arrays
  AllocateArrays();

  // Make FFTW plans
  plan  = fftwf_mpi_plan_dft_r2c_3d(N, N, N, fftwa_r2c, cfftwa_r2c, MPI_COMM_WORLD,
				   FFTW_ESTIMATE);
  iplan = fftwf_mpi_plan_dft_c2r_3d(N, N, N, cfftwa_c2r, fftwa_c2r, MPI_COMM_WORLD,
				   FFTW_ESTIMATE);

  ReportMemory("after array allocation",(float)size_fftw,ovrt,oram);

  // Read input file
  delta = fftwa_r2c;
  ReadDeltaFile(Parameters.DeltaFile);

  // Extrapolate to z=0
  DInit=growth(Parameters.zInit,Parameters.Omegam,Parameters.Omegal,
	       Parameters.w);
  for(int i=0;i<size_fftw;i++) delta[i]/=DInit;

  // Scaling of growth for matter domination at high-z=20
  zHigh=20.0;
  DHigh=growth(zHigh,Parameters.Omegam,Parameters.Omegal,
	       Parameters.w);
  
  // Initial FFT of delta to k-space
  fftwf_execute(plan); 
  delta = (float *)cfftwa_r2c;

  // Find uncorrected zreion field by looping over all scales
  delta2zreion();

  ReportMemory("after delta2zreion",size_fftw,ovrt,oram);
  
  // Bin the ionization history and write
  GetHistory();

  // Get bubble mfp
  GetBubbleMFP();

  MPI_Barrier(MPI_COMM_WORLD);

  // Get corrected ionization history
  GetCHistory();

  // Correct the reionization history
  CorrectHistory();

  // Write Zreion file
  WriteZreionFile();

  ReportMemory("after WriteZreionFile",size_fftw,ovrt,oram);

  // Write Rbubble file
  // if(Parameters.RecordBubbleSizes == 1) WriteRbubbleFile();

  float tau = GetHistory();
  if(myid==0) printf("\n\n tau = %f",tau);

  // Write ionization history
  WriteHistory();

  // Finalize and return
  MPI_Finalize();  if(myid==0) printf("\n\n"); 
  
  fclose(stdout);
  exit(0);

}

