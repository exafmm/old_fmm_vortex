#include "../misc/constants.h"

extern float *sortd;

extern void memoryuse();
extern void memoryfree();

void sortvar(int n0, int n1, float* var, int* nvar) {
  int i;
  sortd = new float [npmax];
  mem = npmax*4;
  memoryuse();
  for( i=n0; i<n1; i++ ) {
    sortd[i] = var[nvar[i]];
  }
  for( i=n0; i<n1; i++ ) {
    var[i] = sortd[i];
  }
  delete[] sortd;
  mem = npmax*4;
  memoryfree();
}

void unsortvar(int n0, int n1, float* var, int* nvar) {
  int i;
  sortd = new float [npmax];
  mem = npmax*4;
  memoryuse();
  for( i=n0; i<n1; i++ ) {
    sortd[nvar[i]] = var[i];
  }
  for( i=n0; i<n1; i++ ) {
    var[i] = sortd[i];
  }
  delete[] sortd;
  mem = npmax*4;
  memoryfree();
}
