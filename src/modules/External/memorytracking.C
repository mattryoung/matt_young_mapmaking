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
#ifdef BGQ
#include <spi/include/kernel/memory.h>
#endif

extern"C" void GetOverhead(){

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

  vrt/=nproc*pow(1024.,2); // Mean virtual memory used per node in GB
  ram/=nproc*pow(1024.,2); // RAM used per node in GB

  if(myid==0){
    fprintf(stderr,"   Overhead:");
    fprintf(stderr,"\n       virtual: %f GB per process",vrt);
    fprintf(stderr,"\n      physical: %f GB per process",ram);
    fprintf(stderr,"\n");
  }

  return;

}

extern"C" void ReportMemory(const char *tag)
{

  int myid;
  int nproc;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

#ifdef BGQ
  if(myid==0){
  uint64_t shared, persist, heapavail, stackavail, stack, heap, guard, mmap;

  Kernel_GetMemorySize(KERNEL_MEMSIZE_GUARD, &guard);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_SHARED, &shared);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_PERSIST, &persist);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &heapavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &stackavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &stack);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_MMAP, &mmap);
  printf("\n   %s:\n",tag);
  printf("   Allocated heap: %.2f MB, avail. heap: %.2f MB\n", (double)heap/(1024*1024),(double)heapavail/(1024*1024));
  printf("   Allocated stack: %.2f MB, avail. stack: %.2f MB\n", (double)stack/(1024*1024), (double)stackavail/(1024*1024));
  printf("   Memory: shared: %.2f MB, persist: %.2f MB, guard: %.2f MB, mmap: %.2f MB\n", (double)shared/(1024*1024), (double)persist/(1024*1024), (double)guard/(1024*1024), (double)mmap/(1024*1024));
  }

  return;
#endif

  float vrt;
  float ram;

  float vrt_local = getVRT(); // Virtual memory used by local process
  float ram_local = getRAM(); // RAM used by local process

  MPI_Allreduce(&vrt_local,&vrt,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&ram_local,&ram,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

  vrt/=nproc*pow(1024.,2); // Mean virtual memory used per process in GB
  ram/=nproc*pow(1024.,2); // RAM used per node in GB

  if(myid==0){
    fprintf(stderr,"\n   %s:",tag);
    fprintf(stderr,"\n       virtual: %f GB per process ",vrt);
    fprintf(stderr,"\n      physical: %f GB per process ",ram);
    fprintf(stderr,"\n\n");
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
