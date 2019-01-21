#include "mpi.h"
#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi,*sortd;
extern int *isort,*isdsp,*iscnt,*irdsp,*ircnt;

extern void memoryuse();
extern void memoryfree();

void mpiposti(int mi, int& ni) {
  int ista,iend,i;

  sortd = new float [npmax];
  mem = npmax*4;
  memoryuse();

  MPI_Allgather(&mi,1,MPI_INT,ircnt,1,MPI_INT,MPI_COMM_WORLD);
  ni = 0;
  for( i=0; i<nprocs; i++ ) {
    irdsp[i] = ni;
    ni += ircnt[i];
  }
  MPI_Allgatherv(xi,mi,MPI_FLOAT,sortd,ircnt,irdsp,MPI_FLOAT,MPI_COMM_WORLD);
  for( i=0; i<ni; i++ ) xi[i] = sortd[i];
  MPI_Allgatherv(yi,mi,MPI_FLOAT,sortd,ircnt,irdsp,MPI_FLOAT,MPI_COMM_WORLD);
  for( i=0; i<ni; i++ ) yi[i] = sortd[i];
  MPI_Allgatherv(zi,mi,MPI_FLOAT,sortd,ircnt,irdsp,MPI_FLOAT,MPI_COMM_WORLD);
  for( i=0; i<ni; i++ ) zi[i] = sortd[i];
  MPI_Allgatherv(gxi,mi,MPI_FLOAT,sortd,ircnt,irdsp,MPI_FLOAT,MPI_COMM_WORLD);
  for( i=0; i<ni; i++ ) gxi[i] = sortd[i];
  MPI_Allgatherv(gyi,mi,MPI_FLOAT,sortd,ircnt,irdsp,MPI_FLOAT,MPI_COMM_WORLD);
  for( i=0; i<ni; i++ ) gyi[i] = sortd[i];
  MPI_Allgatherv(gzi,mi,MPI_FLOAT,sortd,ircnt,irdsp,MPI_FLOAT,MPI_COMM_WORLD);
  for( i=0; i<ni; i++ ) gzi[i] = sortd[i];
  MPI_Allgatherv(vi,mi,MPI_FLOAT,sortd,ircnt,irdsp,MPI_FLOAT,MPI_COMM_WORLD);
  for( i=0; i<ni; i++ ) vi[i] = sortd[i];

  delete[] sortd;
  delete[] isort;
  delete[] isdsp;
  delete[] iscnt;
  delete[] irdsp;
  delete[] ircnt;
  mem = npmax*2*4+nprocs*4*4;
  memoryfree();
}
