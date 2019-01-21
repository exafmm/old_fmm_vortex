#include "mpi.h"
#include "../misc/constants.h"

extern int *nej,*nfj,*nlbj,*nfo,*nij,**neij,**nsij,**npx,**npxd,*ncnt,*nfrom;
extern double (*gx)[nspm],(*gy)[nspm],(*gz)[nspm];
extern int *kscnt,*lsdsp,*lscnt,*lrdsp,*lrcnt;
extern float *fsend,*frecv;

extern void memoryuse();
extern void memoryfree();
extern void mpialltoallvf(float*, int*, int*, float*, int*, int*, int);

void mpisendm2l(int nmp, int nr, int lbjp, int& lbjrp, int lev, double& tic, double* tfmm) {
  int nnmax,nmmax,ic,jj,j,k,ii,i,ieij,npar,ij;
  double toc;

  npxd = new int* [3];
  for( i=0; i<3; i++ ) npxd[i] = new int [nbnp];
  ncnt = new int [nbnes];
  mem = nbnp*3*4+nbnes*4;
  memoryuse();

  nnmax = nbnp*nprocs;
  nmmax = nmp*nnmax;

  if( lev == 1 ) {
// reduce gx-gz
    ic = 0;
    for( jj=0; jj<8; jj++ ) {
      j = jj+nlbj[lev-1];
      for( k=0; k<nmp; k++ ) {
        fsend[ic] = gx[j][k];
        fsend[ic+8*nmp] = gy[j][k];
        fsend[ic+16*nmp] = gz[j][k];
        ic++;
      }
    }
    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;
    if( nprocs == 1 ) {
      for( k=0; k<24*nmp; k++ ) frecv[k] = fsend[k];
    } else {
      MPI_Allreduce(fsend,frecv,24*nmp,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    }
    toc = tic;
    tic = get_time();
    tfmm[4] += tic-toc;
    ic = 0;
    for( jj=0; jj<8; jj++ ) {
      j = jj+nlbj[lev-1];
      for( k=0; k<nmp; k++ ) {
        gx[j][k] = frecv[ic];
        gy[j][k] = frecv[ic+8*nmp];
        gz[j][k] = frecv[ic+16*nmp];
        ic++;
      }
    }
  } else {
// gx
    ic = 0;
    for( ii=0; ii<nprocs; ii++ ) {
      for( i=0; i<kscnt[ii]; i++ ) {
        j = nej[nfj[nsij[i][ii]]]+nlbj[lev-1];
        for( k=0; k<nmp; k++ ) {
          fsend[ic] = gx[j][k];
          ic++;
        }
      }
    }
    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;
    mpialltoallvf(fsend,lscnt,lsdsp,frecv,lrcnt,lrdsp,nmmax);
    toc = tic;
    tic = get_time();
    tfmm[4] += tic-toc;
    for( i=0; i<nbnes; i++ ) ncnt[i] = -1;
    for( i=0; i<nbnp; i++ ) nij[i] = 0;
    for( jj=0; jj<lbjp; jj++ ) {
      j = nej[nfj[jj]]+nlbj[lev-1];
      ncnt[j] = myrank;
      neij[nij[j]][j] = npx[2][jj]*100+npx[1][jj]*10+npx[0][jj];
      nij[j]++;
    }
    jj = lbjp;
    lbjrp = lbjp;
    ic = 0;
    for( i=0; i<nr; i++ ) {
      j = nej[nfj[jj]]+nlbj[lev-1];
      ieij = npx[2][jj]*100+npx[1][jj]*10+npx[0][jj];
      if( ncnt[j] != nfrom[i] ) {
        for( k=0; k<nmp; k++ ) {
          gx[j][k] += frecv[ic];
          ic++;
        }
        ncnt[j] = nfrom[i];
      } else {
        ic += nmp;
      }
      npar = 0;
      for( ij=0; ij<nij[j]; ij++ ) {
        if( neij[ij][j] == ieij ) npar = 1;
      }
      if( npar == 0 ) {
        nfo[lbjrp] = nfj[jj];
        for( k=0; k<3; k++ ) npxd[k][lbjrp] = npx[k][jj];
        neij[nij[j]][j] = ieij;
        nij[j]++;
        lbjrp++;
      }
      jj++;
    }
// gy
    ic = 0;
    for( ii=0; ii<nprocs; ii++ ) {
      for( i=0; i<kscnt[ii]; i++ ) {
        j = nej[nfj[nsij[i][ii]]]+nlbj[lev-1];
        for( k=0; k<nmp; k++ ) {
          fsend[ic] = gy[j][k];
          ic++;
        }
      }
    }
    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;
    mpialltoallvf(fsend,lscnt,lsdsp,frecv,lrcnt,lrdsp,nmmax);
    toc = tic;
    tic = get_time();
    tfmm[4] += tic-toc;
    for( i=0; i<nbnes; i++ ) ncnt[i] = -1;
    for( i=0; i<nbnp; i++ ) nij[i] = 0;
    for( jj=0; jj<lbjp; jj++ ) {
      j = nej[nfj[jj]]+nlbj[lev-1];
      ncnt[j] = myrank;
      neij[nij[j]][j] = npx[2][jj]*100+npx[1][jj]*10+npx[0][jj];
      nij[j]++;
    }
    jj = lbjp;
    lbjrp = lbjp;
    ic = 0;
    for( i=0; i<nr; i++ ) {
      j = nej[nfj[jj]]+nlbj[lev-1];
      ieij = npx[2][jj]*100+npx[1][jj]*10+npx[0][jj];
      if( ncnt[j] != nfrom[i] ) {
        for( k=0; k<nmp; k++ ) {
          gy[j][k] += frecv[ic];
          ic++;
        }
        ncnt[j] = nfrom[i];
      } else {
        ic += nmp;
      }
      npar = 0;
      for( ij=0; ij<nij[j]; ij++ ) {
        if( neij[ij][j] == ieij ) npar = 1;
      }
      if( npar == 0 ) {
        nfo[lbjrp] = nfj[jj];
        for( k=0; k<3; k++ ) npxd[k][lbjrp] = npx[k][jj];
        neij[nij[j]][j] = ieij;
        nij[j]++;
        lbjrp++;
      }
      jj++;
    }
// gz
    ic = 0;
    for( ii=0; ii<nprocs; ii++ ) {
      for( i=0; i<kscnt[ii]; i++ ) {
        j = nej[nfj[nsij[i][ii]]]+nlbj[lev-1];
        for( k=0; k<nmp; k++ ) {
          fsend[ic] = gz[j][k];
          ic++;
        }
      }
    }
    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;
    mpialltoallvf(fsend,lscnt,lsdsp,frecv,lrcnt,lrdsp,nmmax);
    toc = tic;
    tic = get_time();
    tfmm[4] += tic-toc;
    for( i=0; i<nbnes; i++ ) ncnt[i] = -1;
    for( i=0; i<nbnp; i++ ) nij[i] = 0;
    for( jj=0; jj<lbjp; jj++ ) {
      j = nej[nfj[jj]]+nlbj[lev-1];
      ncnt[j] = myrank;
      neij[nij[j]][j] = npx[2][jj]*100+npx[1][jj]*10+npx[0][jj];
      nij[j]++;
    }
    jj = lbjp;
    lbjrp = lbjp;
    ic = 0;
    for( i=0; i<nr; i++ ) {
      j = nej[nfj[jj]]+nlbj[lev-1];
      ieij = npx[2][jj]*100+npx[1][jj]*10+npx[0][jj];
      if( ncnt[j] != nfrom[i] ) {
        for( k=0; k<nmp; k++ ) {
          gz[j][k] += frecv[ic];
          ic++;
        }
        ncnt[j] = nfrom[i];
      } else {
        ic += nmp;
      }
      npar = 0;
      for( ij=0; ij<nij[j]; ij++ ) {
        if( neij[ij][j] == ieij ) npar = 1;
      }
      if( npar == 0 ) {
        nfo[lbjrp] = nfj[jj];
        for( k=0; k<3; k++ ) npxd[k][lbjrp] = npx[k][jj];
        neij[nij[j]][j] = ieij;
        nij[j]++;
        lbjrp++;
      }
      jj++;
    }
    for( jj=lbjp; jj<lbjrp; jj++ ) {
      nfj[jj] = nfo[jj];
      for( k=0; k<3; k++ ) npx[k][jj] = npxd[k][jj];
    }
  }

  for( i=0; i<3; i++ ) delete[] npxd[i];
  delete[] npxd;
  delete[] ncnt;
  mem = nbnp*3*4+nbnes*4;
  memoryfree();
}
