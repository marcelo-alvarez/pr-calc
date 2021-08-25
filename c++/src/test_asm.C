#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv){

  int myid,nproc;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if(myid==0) system("[ -d junk ] || mkdir junk");

}


