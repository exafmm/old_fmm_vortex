#include "mpi.h"
#include "../misc/constants.h"

extern int *nlbj;
extern std::complex<double> (*bx)[mpsym],(*by)[mpsym],(*bz)[mpsym];
extern int *nsend,*nrecv,*nbuf;
extern std::complex<float> *csend,*crecv,*cbuf;

extern void memoryuse();
extern void memoryfree();
extern void mpisendi(int*, int, int, int);
extern void mpisendrecvi(int, int, int, int);
extern void mpirecvi(int*, int, int, int);
extern void mpisendc(std::complex<float>*, int, int, int);
extern void mpisendrecvc(int, int, int, int);
extern void mpirecvc(std::complex<float>*, int, int, int);

void mpireducem(int nmp, int lbj, int lev, double& tic, double* tfmm) {
  int ic,lbjd,jj,j,k,num,isend,lrecv,lsend,i;
  double toc;

  if( pow(8,lev) < nprocs ) {
    nbuf = new int [npmax];
    cbuf = new std::complex<float> [npmax];
    mem = npmax*4+npmax*8;
    memoryuse();

// reduce bx,by,bz that are on different nodes but belong to the same box
    ic = 0;
    lbjd = lbj;
    for( jj=0; jj<lbj; jj++ ) {
      j = jj+nlbj[lev-1];
      for( k=0; k<nmp; k++ ) {
        csend[ic] = bx[j][k];
        csend[ic+lbj*nmp] = by[j][k];
        csend[ic+2*lbj*nmp] = bz[j][k];
        ic++;
      }
    }
    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;

    num = nprocs/int(pow(8,lev));
    for( isend=1; isend<num; isend++ ) {
      lrecv = (myrank/num)*num;
      lsend = (myrank/num)*num+isend;
      nsend[0] = lbj;
      mpisendi(nsend,0,0,lsend);
      mpisendrecvi(1,lsend,lrecv,1);
      mpirecvi(nrecv,0,0,lrecv);
      if( myrank == lrecv ) lbj = nrecv[0];
      if( lbj != 0 ) {
        mpisendc(csend,0,3*lbj*nmp-1,lsend);
        mpisendrecvc(3*lbj*nmp,lsend,lrecv,1);
        mpirecvc(crecv,0,3*lbj*nmp-1,lrecv);
        if( myrank == lrecv ) {
          for( i=0; i<3*lbj*nmp; i++ ) csend[i] += crecv[i];
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
        bx[j][k] = csend[ic];
        by[j][k] = csend[ic+lbj*nmp];
        bz[j][k] = csend[ic+2*lbj*nmp];
        ic++;
      }
    }

    delete[] nbuf;
    delete[] cbuf;
    mem = npmax*4+npmax*8;
    memoryfree();

    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;
  }
}
