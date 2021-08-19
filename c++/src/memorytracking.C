#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "memorytracking.h"
#ifdef DARWIN
#include <mach/mach.h>
#endif

void GetOverhead(float N, float *vrtpc, float *rampc){

  int myid;
  int nproc;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  float vrt;
  float ram;

  float vrt_local = getVRT(); // Virtual memory used by local process
  float ram_local = getRAM(); // RAM used by local process

  MPI_Allreduce(&vrt_local,&vrt,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&ram_local,&ram,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

  *vrtpc = vrt/N*1024./nproc;
  *rampc = ram/N*1024./nproc;

  vrt/=nproc*pow(1024.,2); // Mean virtual memory used per node in GB
  ram/=nproc*pow(1024.,2); // RAM used per node in GB

  if(myid==0){
    fprintf(stderr,"\n Overhead:");
    fprintf(stderr,"\n     virtual: %f GB per process (%f bytes per cell)",vrt,*vrtpc);
    fprintf(stderr,"\n    physical: %f GB per process (%f bytes per cell)",ram,*rampc);
    fprintf(stderr,"\n");
  }

  return;

}

void ReportMemory(const char *tag, float N, float ovrt, float oram)
{

  int myid;
  int nproc;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  float vrt;
  float ram;

  float vrt_local = getVRT(); // Virtual memory used by local process
  float ram_local = getRAM(); // RAM used by local process

  MPI_Allreduce(&vrt_local,&vrt,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&ram_local,&ram,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

  float vrtpc = vrt/N*1024./nproc;
  float rampc = ram/N*1024./nproc;

  vrt/=nproc*pow(1024.,2); // Mean virtual memory used per process in GB
  ram/=nproc*pow(1024.,2); // RAM used per node in GB

  if(myid==0){
    fprintf(stderr,"\n %s:",tag);
    fprintf(stderr,"\n     virtual: %f GB per process (%f bytes per cell)",vrt,vrtpc-ovrt);
    fprintf(stderr,"\n    physical: %f GB per process (%f bytes per cell)",ram,rampc-oram);
    fprintf(stderr,"\n");
  }

}

int parseLine(char* line){
  int i = strlen(line);
  while (*line < '0' || *line > '9') line++;
  line[i-3] = '\0';
  i = atoi(line);
  return i;
}
    
long getVRT(){ //Note: this value is in KB!

#ifndef DARWIN
  FILE* file = fopen("/proc/self/status", "r");
  int result = -1;
  char line[128];
    

  while (fgets(line, 128, file) != NULL){
    if (strncmp(line, "VmSize:", 7) == 0){
      result = parseLine(line);
      break;
    }
  }
  fclose(file);
  return result;
#else
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
  
  if (KERN_SUCCESS != task_info(mach_task_self(),
				TASK_BASIC_INFO, (task_info_t)&t_info, 
				&t_info_count)) return -1;
  return t_info.virtual_size/1024;
#endif  
}

long getRAM(){ //Note: this value is in KB!

#ifndef DARWIN
  FILE* file = fopen("/proc/self/status", "r");
  int result = -1;
  char line[128];

  while (fgets(line, 128, file) != NULL){
    if (strncmp(line, "VmRSS:", 6) == 0){
      result = parseLine(line);
      break;
    }
  }
  fclose(file);
  return result;
#else
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
  
  if (KERN_SUCCESS != task_info(mach_task_self(),
				TASK_BASIC_INFO, (task_info_t)&t_info, 
				&t_info_count)) return -1;
  return t_info.resident_size/1024;
#endif  
}
