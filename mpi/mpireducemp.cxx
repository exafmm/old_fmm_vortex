#include "mpi.h"
#include "../misc/constants.h"

extern int *nlbj;
extern double (*gx)[nspm],(*gy)[nspm],(*gz)[nspm];
extern int *nsend,*nrecv,*nbuf;
extern float *fsend,*frecv,*fbuf;

extern void memoryuse();
extern void memoryfree();
extern void mpisendi(int*, int, int, int);
extern void mpisendrecvi(int, int, int, int);
extern void mpirecvi(int*, int, int, int);
extern void mpisendf(float*, int, int, int);
extern void mpisendrecvf(int, int, int, int);
extern void mpirecvf(float*, int, int, int);

void mpireducem(int nmp, int lbj, int lev, double& tic, double* tfmm) {
  int ic,lbjd,jj,j,k,num,isend,lrecv,lsend,i;
  double toc;

  if( pow(8,lev) < nprocs ) {
    nbuf = new int [npmax];
    fbuf = new float [npmax];
    mem = npmax*2*4;
    memoryuse();

// reduce gx,gy,gz that are on different nodes but belong to the same box
    ic = 0;
    lbjd = lbj;
    for( jj=0; jj<lbj; jj++ ) {
      j = jj+nlbj[lev-1];
      for( k=0; k<nmp; k++ ) {
        fsend[ic] = gx[j][k];
        fsend[ic+lbj*nmp] = gy[j][k];
        fsend[ic+2*lbj*nmp] = gz[j][k];
        ic++;
      }
    }
    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;

    num = int(nprocs/pow(8,lev));
    for( isend=1; isend<num; isend++ ) {
      lrecv = (myrank/num)*num;
      lsend = (myrank/num)*num+isend;
      nsend[0] = lbj;
      mpisendi(nsend,0,0,lsend);
      mpisendrecvi(1,lsend,lrecv,1);
      mpirecvi(nrecv,0,0,lrecv);
      if( myrank == lrecv ) lbj = nrecv[0];
      if( lbj != 0 ) {
        mpisendf(fsend,0,3*lbj*nmp-1,lsend);
        mpisendrecvf(3*lbj*nmp,lsend,lrecv,1);
        mpirecvf(frecv,0,3*lbj*nmp-1,lrecv);
        if( myrank == lrecv ) {
          for( i=0; i<3*lbj*nmp; i++ ) fsend[i] += frecv[i];
        }
      }
    }
    lbj = lbjd;
    toc = tic;
    tic = get_time();
    tfmm[4] += tic-toc;

    ic = 0;
    for( jj=0; jj<lbj; jj++ ) {
      j = jj+nlbj[lev-1];
      for( k=0; k<nmp; k++ ) {
        gx[k][j] = fsend[ic];
        gy[k][j] = fsend[ic+lbj*nmp];
        gz[k][j] = fsend[ic+2*lbj*nmp];
        ic++;
      }
    }

    delete[] nbuf;
    delete[] fbuf;
    mem = npmax*2*4;
    memoryfree();

    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;
  }
}
