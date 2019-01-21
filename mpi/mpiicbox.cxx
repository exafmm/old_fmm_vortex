#include "mpi.h"
#include "../misc/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern int *nei,*nfi,**ndj,*nej,*nfj,*nlbj,**nsij,**npx;
extern int *ksdsp,*kscnt,*krdsp,*krcnt,*lsdsp,*lscnt,*lrdsp,*lrcnt;
extern int *nsend,*nrecv;
extern float *fsend,*frecv;
extern std::complex<float> *csend,*crecv;

extern void memoryuse();
extern void memoryfree();
extern void mpireducem(int, int, int, double&, double*);
extern void jcbox(int, int&, int, int, int);
extern void ijbox(int, int, int, int, int);
extern void jcnt(int&, int, int, int, int);
extern void mpialltoallvi(int*, int*, int*, int*, int*, int*, int);
extern void mpisendp2p(int, int, float*, double&, double*);
extern void mpisendm2l(int, int, int, int&, int, double&, double*);
extern void boxdatai(int, int, int, int&, double&);

void icbox(int n0, int n1, int& mjp, int nmp, int lbi, int lbj, int& lbjr, int& lbjrp, int lev, int ipb, int npb, double* tfmm) {
  int nnmax,n3max,nmmax,ii,i,lbjp,k,j,nr,ic,jj,np;
  double tic,toc,rb;

  tic = get_time();
  nnmax = nbnp*nprocs;
  n3max = 3*nnmax;
  nmmax = nmp*nnmax;

  ksdsp = new int [nprocs];
  kscnt = new int [nprocs];
  krdsp = new int [nprocs];
  krcnt = new int [nprocs];
  lsdsp = new int [nprocs];
  lscnt = new int [nprocs];
  lrdsp = new int [nprocs];
  lrcnt = new int [nprocs];
  nsend = new int [n3max];
  nrecv = new int [n3max];
  fsend = new float [npmax];
  frecv = new float [npmax];
  csend = new std::complex<float> [nmmax];
  crecv = new std::complex<float> [nmmax];
  mem = nprocs*8*4+n3max*2*4+npmax*2*4+nmmax*2*8;
  memoryuse();

  if( lev == 1 ) {
    mpisendm2l(nmp,nr,lbj,lbjrp,lev,tic,tfmm);
    tic = get_time();
  } else {
    mpireducem(nmp,lbj,lev,tic,tfmm);
    tic = get_time();

    lbi = pow(8,lev);
    for( ii=0; ii<lbi; ii++ ) {
      nei[ii] = ii;
      nfi[ii] = ii;
    }
    for( i=0; i<nbnp; i++ ) {
      for( j=0; j<3; j++ ) {
        npx[j][i] = 0;
      }
    }
    jcbox(lbj,lbjp,lev,ipb,npb);
    ijbox(lbi,lbjp,lev,ipb,npb);
    jcnt(nr,lbi,lbjp,lev,ipb);
// nfj
    ic = 0;
    for( ii=0; ii<nprocs; ii++ ) {
      for( i=0; i<kscnt[ii]; i++ ) {
        jj = nsij[i][ii];
        nsend[ic] = nfj[jj];
        ic++;
      }
    }
    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;
    mpialltoallvi(nsend,kscnt,ksdsp,nrecv,krcnt,krdsp,nnmax);
    toc = tic;
    tic = get_time();
    tfmm[4] += tic-toc;
    lbjrp = lbjp;
    for( i=0; i<nr; i++ ) {
      nfj[lbjrp] = nrecv[i];
      lbjrp++;
    }
// nej
    lbjr = lbj;
    lbjrp = lbjp;
    for( ii=0; ii<nr; ii++ ) {
      if( nej[nfj[lbjrp]] == -1 ) {
        nej[nfj[lbjrp]] = lbjr;
        lbjr++;
      }
      lbjrp++;
    }
// lscnt for npx
    for( i=0; i<nprocs; i++ ) {
      lsdsp[i] = 3*ksdsp[i];
      lscnt[i] = 3*kscnt[i];
    }
    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;
    MPI_Alltoall(lscnt,1,MPI_INT,lrcnt,1,MPI_INT,MPI_COMM_WORLD);
    toc = tic;
    tic = get_time();
    tfmm[4] += tic-toc;
    ic = 0;
    for( i=0; i<nprocs; i++ ) {
      lrdsp[i] = ic;
      ic += lrcnt[i];
    }
// npx
    ic = 0;
    for( ii=0; ii<nprocs; ii++ ) {
      for( i=0; i<kscnt[ii]; i++ ) {
        jj = nsij[i][ii];
        for( k=0; k<3; k++ ) {
          nsend[ic] = npx[k][jj];
          ic++;
        }
      }
    }
    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;
    mpialltoallvi(nsend,lscnt,lsdsp,nrecv,lrcnt,lrdsp,nnmax);
    toc = tic;
    tic = get_time();
    tfmm[4] += tic-toc;
    lbjrp = lbjp;
    ic = 0;
    for( i=0; i<nr; i++ ) {
      for( j=0; j<3; j++ ) {
        npx[j][lbjrp] = nrecv[ic];
        ic++;
      }
      lbjrp++;
    }

    if( ipb == -2 ) {
// iscnt for xj-sj
      for( i=0; i<nprocs; i++ ) lscnt[i] = 0;
      for( ii=0; ii<nprocs; ii++ ) {
        for( i=0; i<kscnt[ii]; i++ ) {
          jj = nsij[i][ii];
          lscnt[ii] += ndj[1][jj]-ndj[0][jj]+1;
        }
      }
      ic = 0;
      for( i=0; i<nprocs; i++ ) {
        lsdsp[i] = ic;
        ic += lscnt[i];
      }
      toc = tic;
      tic = get_time();
      tfmm[8] += tic-toc;
      MPI_Alltoall(lscnt,1,MPI_INT,lrcnt,1,MPI_INT,MPI_COMM_WORLD);
      toc = tic;
      tic = get_time();
      tfmm[0] += tic-toc;
      np = 0;
      for( i=0; i<nprocs; i++ ) {
        lrdsp[i] = np;
        np += lrcnt[i];
      }
// xj-vj
      mpisendp2p(mjp,np,xj,tic,tfmm);
      mpisendp2p(mjp,np,yj,tic,tfmm);
      mpisendp2p(mjp,np,zj,tic,tfmm);
      mpisendp2p(mjp,np,gxj,tic,tfmm);
      mpisendp2p(mjp,np,gyj,tic,tfmm);
      mpisendp2p(mjp,np,gzj,tic,tfmm);
      mpisendp2p(mjp,np,vj,tic,tfmm);
      mpisendp2p(mjp,np,sj,tic,tfmm);
      tic = get_time();
// ndj
      ic = 0;
      for( ii=0; ii<nprocs; ii++ ) {
        for( i=0; i<kscnt[ii]; i++ ) {
          jj = nsij[i][ii];
          nsend[ic] = ndj[1][jj]-ndj[0][jj]+1;
          ic++;
        }
      }
      toc = tic;
      tic = get_time();
      tfmm[8] += tic-toc;
      mpialltoallvi(nsend,kscnt,ksdsp,nrecv,krcnt,krdsp,nnmax);
      toc = tic;
      tic = get_time();
      tfmm[0] += tic-toc;
      ic = mjp;
      lbjrp = lbjp;
      for( i=0; i<nr; i++ ) {
        ndj[0][lbjrp] = ic;
        ic += nrecv[i];
        ndj[1][lbjrp] = ic-1;
        lbjrp++;
      }
      mjp = ic;
    } else {
// lscnt for bx-bz (gx-gz)
      for( i=0; i<nprocs; i++ ) {
        lsdsp[i] = nmp*ksdsp[i];
        lscnt[i] = nmp*kscnt[i];
      }
      toc = tic;
      tic = get_time();
      tfmm[8] += tic-toc;
      MPI_Alltoall(lscnt,1,MPI_INT,lrcnt,1,MPI_INT,MPI_COMM_WORLD);
      toc = tic;
      tic = get_time();
      tfmm[4] += tic-toc;
      ic = 0;
      for( i=0; i<nprocs; i++ ) {
        lrdsp[i] = ic;
        ic += lrcnt[i];
      }
      mpisendm2l(nmp,nr,lbjp,lbjrp,lev,tic,tfmm);
      tic = get_time();

    }

    boxdatai(n0,n1,lev,lbi,rb);
    ijbox(lbi,lbjrp,lev,ipb,npb);
  }

  toc = tic;
  tic = get_time();
  tfmm[8] += tic-toc;

  delete[] ksdsp;
  delete[] kscnt;
  delete[] krdsp;
  delete[] krcnt;
  delete[] lsdsp;
  delete[] lscnt;
  delete[] lrdsp;
  delete[] lrcnt;
  delete[] nsend;
  delete[] nrecv;
  delete[] fsend;
  delete[] frecv;
  delete[] csend;
  delete[] crecv;
  mem = nprocs*8*4+nnmax*2*4+npmax*2*4+nmmax*2*8;
  memoryfree();
}
