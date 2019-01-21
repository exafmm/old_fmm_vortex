#include "mpi.h"
#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern int *isort,*isdsp,*iscnt,*irdsp,*ircnt;

extern void memoryuse();
extern void mpirange(int, int, int&, int&);

void mpiprei(int ni, int& mi) {
  int ista,iend,i;

  isort = new int [npmax];
  isdsp = new int [nprocs];
  iscnt = new int [nprocs];
  irdsp = new int [nprocs];
  ircnt = new int [nprocs];
  mem = npmax*4+nprocs*4*4;
  memoryuse();

  mpirange(0,ni-1,ista,iend);
  for( i=ista; i<=iend; i++ ) {
    mi = i-ista;
    xi[mi] = xi[i];
    yi[mi] = yi[i];
    zi[mi] = zi[i];
    gxi[mi] = gxi[i];
    gyi[mi] = gyi[i];
    gzi[mi] = gzi[i];
    vi[mi] = vi[i];
  }
  mi++;

}
