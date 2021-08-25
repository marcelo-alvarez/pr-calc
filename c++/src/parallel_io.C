#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#define OFFSET_MAX 2199023255552 // 4*8192^3

int parallel_write(char *filename, long size, long offset, void *buffer) {  
 
  int fd; //file descriptor
  int stat; // I/O status
  if (size + offset > OFFSET_MAX) {
    fprintf(stderr,"You are trying to access a file location\n");
    fprintf(stderr,"which is bigger than %ld\n",(long)OFFSET_MAX);
    fprintf(stderr,"Verify your code and/or change OFFSET_MAX\n");
    return(1);
  }
  fd = open(filename,O_RDWR|O_CREAT,S_IRUSR|S_IWUSR);
  stat=pwrite(fd,buffer,size,offset);
  close(fd);
  if (stat != size) return(1);
  else return(0);
}

int parallel_read(char *filename, long size, long offset, void *buffer) {

  int fd; //file descriptor
  int stat; // I/O status
  if (size + offset > OFFSET_MAX) {
    fprintf(stderr,"You are trying to access a file location\n");
    fprintf(stderr,"which is bigger than %ld\n",(long)OFFSET_MAX);
    fprintf(stderr,"Verify your code and/or change OFFSET_MAX\n");
    return(1);
  }
  fd = open(filename,O_RDONLY);
  stat = pread(fd,buffer,size,offset);
  close(fd);
  if (stat != size) return(1);
  else return(0);
}

