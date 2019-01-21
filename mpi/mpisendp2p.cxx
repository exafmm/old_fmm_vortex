#include "mpi.h"
#include "../misc/constants.h"

extern int **ndj,**nsij,*kscnt,*lsdsp,*lscnt,*lrdsp,*lrcnt;
extern float *fsend,*frecv;

extern void mpialltoallvf(float*, int*, int*, float*, int*, int*, int);

void mpisendp2p(int mjp, int np, float* var, double& tic, double* tfmm) {
  int ic,ii,i,j,k;
  double toc;

  ic = 0;
  for( ii=0; ii<nprocs; ii++ ) {
    for( i=0; i<kscnt[ii]; i++ ) {
      j = nsij[i][ii];
      for( k=ndj[0][j]; k<=ndj[1][j]; k++ ) {
        fsend[ic] = var[k];
        ic++;
      }
    }
  }
  toc = tic;
  tic = get_time();
  tfmm[8] += tic-toc;
  MPI_Barrier(MPI_COMM_WORLD);
  mpialltoallvf(fsend,lscnt,lsdsp,frecv,lrcnt,lrdsp,npmax);
  toc = tic;
  tic = get_time();
  tfmm[0] += tic-toc;
  ic = mjp;
  for( i=0; i<np; i++ ) {
    var[ic] = frecv[i];
    ic++;
  }
}
