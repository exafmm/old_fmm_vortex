#include "mpi.h"
#include "../misc/constants.h"

extern float *sortd;

extern void memoryuse();
extern void memoryfree();
extern void mpialltoallvf(float*, int*, int*, float*, int*, int*, int);

void mpisortvar(int* nsdsp, int* nscnt, int* nrdsp, int*nrcnt, float* var) {
  int n,i;

  sortd = new float [npmax];
  mem = npmax*4;
  memoryuse();

  n = 0;
  for( i=0; i<nprocs; i++ ) n += nrcnt[i];
  mpialltoallvf(var,nscnt,nsdsp,sortd,nrcnt,nrdsp,npmax);
  for( i=0; i<n; i++ ) var[i] = sortd[i];

  delete[] sortd;
  mem = npmax*4;
  memoryfree();
}

void mpiunsortvar(int* nsdsp, int* nscnt, int* nrdsp, int*nrcnt, float* var) {
  int n,i;

  sortd = new float [npmax];
  mem = npmax*4;
  memoryuse();

  n = 0;
  for( i=0; i<nprocs; i++ ) n += nscnt[i];
  mpialltoallvf(var,nrcnt,nrdsp,sortd,nscnt,nsdsp,npmax);
  for( i=0; i<n; i++ ) var[i] = sortd[i];

  delete[] sortd;
  mem = npmax*4;
  memoryfree();
}
