#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include "myio.h"

extern"C" void myfwrite(char* filename, float* buffer, long size, long offset){

  int fd   = open(filename,O_RDWR|O_CREAT,S_IRUSR|S_IWUSR);
  int stat = pwrite(fd,buffer,size,offset);
  close(fd);

  return;

}
