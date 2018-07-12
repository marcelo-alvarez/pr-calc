#include <math.h>
#include <mpi.h>

double GetMeanSigma(float *array, int report, int N, int Nlocal)
{
  
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  double lmean = 0, lsigma = 0;
  double mean =0, sigma=0;

  for(long i=0;i<Nlocal;i++){
    for(long j=0;j<N;j++){
      for(long k=0;k<N;k++){
	long index = i*N*(N+2)+j*(N+2)+k;
	lmean+=array[index];
	lsigma+=array[index]*array[index];
      }
    }
  }

  MPI_Allreduce(&lmean,  &mean,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&lsigma, &sigma, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  mean/= pow((float)N,3);

  sigma/= pow((float)N,3);
  sigma = sqrt(sigma);

  if(report==1 && myid==0) printf("\n mean = %e sigma = %e",mean,sigma);

  return sigma;

}

double GetSigmaComplex(float *array, int report, int N, int Nlocal)
{
  
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  double lsigma = 0;
  double sigma = 0;

  long size_fftw = (long)N*(long)(N+2)*(long)Nlocal;

  for(long i=0;i<size_fftw;i++){
    lsigma+=array[i]*array[i];
  }

  MPI_Allreduce(&lsigma, &sigma, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  sigma/= pow((float)N,3);
  sigma = sqrt(sigma*2);

  if(myid==0 && report==1) printf("\n sigma = %e",sigma);

  return sigma;

}

void NormalizeArray(float *array, int N, long size_fftw, float power)
{

  for(int i=0;i<size_fftw;i++) array[i]/=pow((float)N,power);

}

void CopyArray(float *array1, float *array2, long size_fftw)
{
  for(int i=0;i<size_fftw;i++) array2[i]=array1[i];
}
