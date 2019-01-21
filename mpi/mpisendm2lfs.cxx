#include "mpi.h"
#include "../misc/constants.h"

extern int *nej,*nfj,*nlbj,*nel,*nlbl,*nem,*nlbm,*nfo,*nij,**neij,**nsij,**npx,**npxd,*ncnt,*nfrom;
extern std::complex<double> (*bx)[mpsym],(*by)[mpsym],(*bz)[mpsym];
extern int *kscnt,*lsdsp,*lscnt,*lrdsp,*lrcnt;
extern std::complex<float> *csend,*crecv;

extern void memoryuse();
extern void memoryfree();
extern void mpialltoallvc(std::complex<float>*, int*, int*, std::complex<float>*, int*, int*, int);

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
// reduce bx-bz
    ic = 0;
    for( jj=0; jj<8; jj++ ) {
      if( nr == 0 ) {
        j = jj+nlbj[lev-1];
      } else if( nr > 0 ) {
        j = jj+nlbl[lev-1];
      } else {
        j = jj+nlbm[lev-1];
      }
      for( k=0; k<nmp; k++ ) {
        csend[ic] = bx[j][k];
        csend[ic+8*nmp] = by[j][k];
        csend[ic+16*nmp] = bz[j][k];
        ic++;
      }
    }
    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;
    if( nprocs == 1 ) {
      for( k=0; k<24*nmp; k++ ) crecv[k] = csend[k];
    } else {
      MPI_Allreduce(csend,crecv,24*nmp,MPI_COMPLEX,MPI_SUM,MPI_COMM_WORLD);
    }
    toc = tic;
    tic = get_time();
    tfmm[4] += tic-toc;
    ic = 0;
    for( jj=0; jj<8; jj++ ) {
      if( nr == 0 ) {
        j = jj+nlbj[lev-1];
      } else if( nr > 0 ) {
        j = jj+nlbl[lev-1];
      } else {
        j = jj+nlbm[lev-1];
      }
      for( k=0; k<nmp; k++ ) {
        bx[j][k] = crecv[ic];
        by[j][k] = crecv[ic+8*nmp];
        bz[j][k] = crecv[ic+16*nmp];
        ic++;
      }
    }
  } else {
// bx
    ic = 0;
    for( ii=0; ii<nprocs; ii++ ) {
      for( i=0; i<kscnt[ii]; i++ ) {
        lbjrp = nsij[i][ii];
        if( npx[2][lbjrp] == 0 ) {
          j = nej[nfj[lbjrp]]+nlbj[lev-1];
        } else if( npx[2][lbjrp] > 0 ) {
          j = nel[nfj[lbjrp]]+nlbl[lev-1];
        } else {
          j = nem[nfj[lbjrp]]+nlbm[lev-1];
        }
        for( k=0; k<nmp; k++ ) {
          csend[ic] = bx[j][k];
          ic++;
        }
      }
    }
    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;
    mpialltoallvc(csend,lscnt,lsdsp,crecv,lrcnt,lrdsp,nmmax);
    toc = tic;
    tic = get_time();
    tfmm[4] += tic-toc;
    for( i=0; i<nbnes; i++ ) ncnt[i] = -1;
    for( i=0; i<nbnp; i++ ) nij[i] = 0;
    for( jj=0; jj<lbjp; jj++ ) {
      if( npx[2][jj] == 0 ) {
        j = nej[nfj[jj]]+nlbj[lev-1];
      } else if( npx[2][jj] > 0 ) {
        j = nel[nfj[jj]]+nlbl[lev-1];
      } else {
        j = nem[nfj[jj]]+nlbm[lev-1];
      }
      ncnt[j] = myrank;
      neij[nij[j]][j] = npx[2][jj]*100+npx[1][jj]*10+npx[0][jj];
      nij[j]++;
    }
    jj = lbjp;
    lbjrp = lbjp;
    ic = 0;
    for( i=0; i<nr; i++ ) {
      if( npx[2][jj] == 0 ) {
        j = nej[nfj[jj]]+nlbj[lev-1];
      } else if( npx[2][jj] > 0 ) {
        j = nel[nfj[jj]]+nlbl[lev-1];
      } else {
        j = nem[nfj[jj]]+nlbm[lev-1];
      }
      ieij = npx[2][jj]*100+npx[1][jj]*10+npx[0][jj];
      if( ncnt[j] != nfrom[i] ) {
        for( k=0; k<nmp; k++ ) {
          bx[j][k] += crecv[ic];
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
// by
    ic = 0;
    for( ii=0; ii<nprocs; ii++ ) {
      for( i=0; i<kscnt[ii]; i++ ) {
        lbjrp = nsij[i][ii];
        if( npx[2][lbjrp] == 0 ) {
          j = nej[nfj[lbjrp]]+nlbj[lev-1];
        } else if( npx[2][lbjrp] > 0 ) {
          j = nel[nfj[lbjrp]]+nlbl[lev-1];
        } else {
          j = nem[nfj[lbjrp]]+nlbm[lev-1];
        }
        for( k=0; k<nmp; k++ ) {
          csend[ic] = by[j][k];
          ic++;
        }
      }
    }
    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;
    mpialltoallvc(csend,lscnt,lsdsp,crecv,lrcnt,lrdsp,nmmax);
    toc = tic;
    tic = get_time();
    tfmm[4] += tic-toc;
    for( i=0; i<nbnes; i++ ) ncnt[i] = -1;
    for( i=0; i<nbnp; i++ ) nij[i] = 0;
    for( jj=0; jj<lbjp; jj++ ) {
      if( npx[2][jj] == 0 ) {
        j = nej[nfj[jj]]+nlbj[lev-1];
      } else if( npx[2][jj] > 0 ) {
        j = nel[nfj[jj]]+nlbl[lev-1];
      } else {
        j = nem[nfj[jj]]+nlbm[lev-1];
      }
      ncnt[j] = myrank;
      neij[nij[j]][j] = npx[2][jj]*100+npx[1][jj]*10+npx[0][jj];
      nij[j]++;
    }
    jj = lbjp;
    lbjrp = lbjp;
    ic = 0;
    for( i=0; i<nr; i++ ) {
      if( npx[2][jj] == 0 ) {
        j = nej[nfj[jj]]+nlbj[lev-1];
      } else if( npx[2][jj] > 0 ) {
        j = nel[nfj[jj]]+nlbl[lev-1];
      } else {
        j = nem[nfj[jj]]+nlbm[lev-1];
      }
      ieij = npx[2][jj]*100+npx[1][jj]*10+npx[0][jj];
      if( ncnt[j] != nfrom[i] ) {
        for( k=0; k<nmp; k++ ) {
          by[j][k] += crecv[ic];
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
// bz
    ic = 0;
    for( ii=0; ii<nprocs; ii++ ) {
      for( i=0; i<kscnt[ii]; i++ ) {
        lbjrp = nsij[i][ii];
        if( npx[2][lbjrp] == 0 ) {
          j = nej[nfj[lbjrp]]+nlbj[lev-1];
        } else if( npx[2][lbjrp] > 0 ) {
          j = nel[nfj[lbjrp]]+nlbl[lev-1];
        } else {
          j = nem[nfj[lbjrp]]+nlbm[lev-1];
        }
        for( k=0; k<nmp; k++ ) {
          csend[ic] = bz[j][k];
          ic++;
        }
      }
    }
    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;
    mpialltoallvc(csend,lscnt,lsdsp,crecv,lrcnt,lrdsp,nmmax);
    toc = tic;
    tic = get_time();
    tfmm[4] += tic-toc;
    for( i=0; i<nbnes; i++ ) ncnt[i] = -1;
    for( i=0; i<nbnp; i++ ) nij[i] = 0;
    for( jj=0; jj<lbjp; jj++ ) {
      if( npx[2][jj] == 0 ) {
        j = nej[nfj[jj]]+nlbj[lev-1];
      } else if( npx[2][jj] > 0 ) {
        j = nel[nfj[jj]]+nlbl[lev-1];
      } else {
        j = nem[nfj[jj]]+nlbm[lev-1];
      }
      ncnt[j] = myrank;
      neij[nij[j]][j] = npx[2][jj]*100+npx[1][jj]*10+npx[0][jj];
      nij[j]++;
    }
    jj = lbjp;
    lbjrp = lbjp;
    ic = 0;
    for( i=0; i<nr; i++ ) {
      if( npx[2][jj] == 0 ) {
        j = nej[nfj[jj]]+nlbj[lev-1];
      } else if( npx[2][jj] > 0 ) {
        j = nel[nfj[jj]]+nlbl[lev-1];
      } else {
        j = nem[nfj[jj]]+nlbm[lev-1];
      }
      ieij = npx[2][jj]*100+npx[1][jj]*10+npx[0][jj];
      if( ncnt[j] != nfrom[i] ) {
        for( k=0; k<nmp; k++ ) {
          bz[j][k] += crecv[ic];
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
