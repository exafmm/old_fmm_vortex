#include "mpi.h"
#include "../misc/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj,*sortd;
extern int *jsort,*jsdsp,*jscnt,*jrdsp,*jrcnt;

extern void memoryuse();
extern void memoryfree();

void mpipostj(int mj, int& nj) {
  int jsta,jend,i;

  sortd = new float [npmax];
  mem = npmax*4;
  memoryuse();

  MPI_Allgather(&mj,1,MPI_INT,jrcnt,1,MPI_INT,MPI_COMM_WORLD);
  nj = 0;
  for( i=0; i<nprocs; i++ ) {
    jrdsp[i] = nj;
    nj += jrcnt[i];
  }
  MPI_Allgatherv(xj,mj,MPI_FLOAT,sortd,jrcnt,jrdsp,MPI_FLOAT,MPI_COMM_WORLD);
  for( i=0; i<nj; i++ ) xj[i] = sortd[i];
  MPI_Allgatherv(yj,mj,MPI_FLOAT,sortd,jrcnt,jrdsp,MPI_FLOAT,MPI_COMM_WORLD);
  for( i=0; i<nj; i++ ) yj[i] = sortd[i];
  MPI_Allgatherv(zj,mj,MPI_FLOAT,sortd,jrcnt,jrdsp,MPI_FLOAT,MPI_COMM_WORLD);
  for( i=0; i<nj; i++ ) zj[i] = sortd[i];
  MPI_Allgatherv(gxj,mj,MPI_FLOAT,sortd,jrcnt,jrdsp,MPI_FLOAT,MPI_COMM_WORLD);
  for( i=0; i<nj; i++ ) gxj[i] = sortd[i];
  MPI_Allgatherv(gyj,mj,MPI_FLOAT,sortd,jrcnt,jrdsp,MPI_FLOAT,MPI_COMM_WORLD);
  for( i=0; i<nj; i++ ) gyj[i] = sortd[i];
  MPI_Allgatherv(gzj,mj,MPI_FLOAT,sortd,jrcnt,jrdsp,MPI_FLOAT,MPI_COMM_WORLD);
  for( i=0; i<nj; i++ ) gzj[i] = sortd[i];
  MPI_Allgatherv(vj,mj,MPI_FLOAT,sortd,jrcnt,jrdsp,MPI_FLOAT,MPI_COMM_WORLD);
  for( i=0; i<nj; i++ ) vj[i] = sortd[i];
  MPI_Allgatherv(sj,mj,MPI_FLOAT,sortd,jrcnt,jrdsp,MPI_FLOAT,MPI_COMM_WORLD);
  for( i=0; i<nj; i++ ) sj[i] = sortd[i];

  delete[] sortd;
  delete[] jsort;
  delete[] jsdsp;
  delete[] jscnt;
  delete[] jrdsp;
  delete[] jrcnt;
  mem = npmax*2*4+nprocs*4*4;
  memoryfree();
}
