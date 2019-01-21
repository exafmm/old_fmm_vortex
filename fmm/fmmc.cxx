#include "../misc/parameters.h"
#include "../misc/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern int *nbi,**ndi,*nei,*nfi;
extern int *nbj,**ndj,*nej,*nfj,*nlbj,*neo,*nfo;
extern int **nxs,*nfn,*nij,*njb,*nnp,**neij,**nsij,**npx,*nfrom,*nc,*nd;
extern float *gxd,*gyd,*gzd,*vd;
extern int *nek,*ixadj,*nxadj,*nadjncy,*nvwgt,*nadjwgt,*nvtxdist;,*npart;
extern int *npartd,*irank,*na,*nb;

extern void memoryuse();
extern void memoryfree();
extern void nlevel(int);
extern void boxallocate(int, int, int, int);
extern void alloc();
extern void dealloc();
extern void trans(int&, int);
extern void boxpart(int, int, int);
extern void sorti(int&);
extern void sortj(int&);
extern void boxdatai(int, int, int, int&, double&);
extern void boxdataj(int, int, int, int&, double&);
extern void jcbox(int, int&, int, int, int);
extern void ijbox(int, int, int, int, int);
extern void icbox(int, int, int&, int, int, int, int&, int&, int, int, int, double*);
extern void potp2p(int, int);
extern void bsp2p(int, int);
extern void stp2p(int, int);
extern void psep2p(int, int, double);
extern void fgtp2p(int, int);
extern void remp2p(int, int);
extern void wallp2p(int, int, double);
extern void p2m(int, int, int, double, int);
extern void m2m1(int, int, int, int, int, int*, double);
extern void m2m2(int, int, int, int, int, double);
extern void m2l(int, int, int, int, int, int, double);
extern void boxc(int, int, int*);
extern void l2l1(int, int, int, double);
extern void l2l2(int, int, int, double);
extern void bsl2p(int, int, int, double);
extern void stl2p(int, int, int, double);
extern void potl2p(int, int, int, double);
extern void unsorti(int&);
extern void unsortj(int&);

void fmm(int mi, int mj, int neqd, double* tfmm){
  int neq,lev,lmin,nebmp,nbnep,i,nmp,lbi,nbyte,ii,lbj,j,lbjp,ic,jj,mjp,lbjr,lbjo,ib,ipb,jb;
  double tic,toc,rb,vist;
  tic = get_time();

  neq = neqd;
  nlevel(mj);
  boxallocate(0,mi,0,mj);
  lev = lmax;
  lmin = 1;
  nbmax = pow(8,lmax);
  nbnes = nbnet+2*npb;
  nebmp = std::max(nebm,nprocs);
  nbnep = std::max(nbne,nprocs);

  nbi = new int [npmax];
  ndi = new int* [2];
  for( i=0; i<2; i++ ) ndi[i] = new int [nbne];
  nei = new int [nbmax];
  nfi = new int [nbmax];
  nbj = new int [npmax];
  ndj = new int* [2];
  for( i=0; i<2; i++ ) ndj[i] = new int [nbnp];
  nej = new int [nbmax];
  nfj = new int [nbnp];
  nlbj = new int [lmax];
  neo = new int [nbmax];
  nfo = new int [nbnp];
  nxs = new int* [3];
  for( i=0; i<3; i++ ) nxs[i] = new int [npmax];
  nfn = new int [npmax];
  nij = new int [nbnp];
  njb = new int [nbnp];
  nnp = new int [nbnp];
  npx = new int* [3];
  for( i=0; i<3; i++ ) npx[i] = new int [nbnp];
  nfrom = new int [nbnp];
  nc  = new int [3];
  nd  = new int [3];
  gxd = new float [npmax];
  gyd = new float [npmax];
  gzd = new float [npmax];
  vd = new float [npmax];
  nek = new int [nbmax];
  ixadj = new int [nbne];
  nxadj = new int [nbne+1];
  nadjncy = new int [npmax];
  nvwgt = new int [nbne];
  nadjwgt = new int [npmax];
  nvtxdist = new int [nprocs+1];
  npart = new int [nbne+1];
  npartd = new int [nbnep];
  irank = new int [nbnep];
  na = new int [npmax];
  nb = new int [npmax];
  mem = npmax*14*4+nbmax*4*4+nbne*8*4+nbnp*11*4+nbnep*2*4+lmax*4+3*2*4;
  memoryuse();

  neij = new int* [nebmp];
  for( i=0; i<nebmp; i++ ) neij[i] = new int [nbnp];
  nsij = new int* [nbnp];
  for( i=0; i<nbnp; i++ ) nsij[i] = new int [nebmp];
  mem = nebmp*nbnp*2*4;
  memoryuse();

  alloc();
  trans(nmp,mp);
  nlbj[lev-1] = 0;

  boxpart(mi,mj,lev);

  sorti(mi);
  sortj(mj);

  boxdatai(0,mi,lev,lbi,rb);
  if( neq == 10 && myrank == 0 ) {
    std::fstream fid;
    fid.open("density.dat",std::ios::out|std::ios::binary|std::ios::ate);
    nbyte = sizeof(int);
    fid.write((char *)(&nbyte),sizeof(int));
    fid.write((char *)(&lbi),sizeof(int));
    fid.write((char *)(&nbyte),sizeof(int));
    for( ii=0; ii<lbi; ii++ ) {
      nbyte = sizeof(int);
      i = ndi[1][ii]-ndi[0][ii]+1;
      fid.write((char *)(&nbyte),sizeof(int));
      fid.write((char *)(&i),sizeof(int));
      fid.write((char *)(&nbyte),sizeof(int));
    }
    fid.close();
  }

// Step 1. P2P-summation @ lev lmax

  boxdataj(0,mj,lev,lbj,rb);

  for( i=0; i<nbnp; i++ ) {
    for( j=0; j<3; j++ ) {
      npx[j][i] = 0;
    }
  }
  jcbox(lbj,lbjp,lev,-2,npb);
  ijbox(lbi,lbjp,lev,-2,npb);

  ic = mj;
  for( jj=lbj; jj<lbjp; jj++ ) {
    ndj[0][jj] = ic;
    for( j=ndj[0][nnp[jj]]; j<=ndj[1][nnp[jj]]; j++ ) {
      xj[ic] = xj[j]+npx[0][jj]*rd;
      yj[ic] = yj[j]+npx[1][jj]*rd;
      zj[ic] = zj[j]+npx[2][jj]*rd;
      gxj[ic] = gxj[j];
      gyj[ic] = gyj[j];
      gzj[ic] = gzj[j];
      vj[ic] = vj[j];
      sj[ic] = sj[j];
      ic++;
    }
    ndj[1][jj] = ic-1;
  }
  mjp = ic;

  toc = tic;
  tic = get_time();
  tfmm[8] += tic-toc;

// Step 1.1 Send and recv particles on other nodes

  icbox(0,mi,mjp,nmp,lbi,lbj,lbjr,lbjp,lev,-2,npb,tfmm);

  toc = tic;
  tic = get_time();

// Step 1.2 P2P-summation for other nodes

  for( i=0; i<mi; i++ ) {
    gxd[i] = 0;
    gyd[i] = 0;
    gzd[i] = 0;
    vd[i] = 0;
  }
  if( neq == 5 ) {
    for( i=0; i<mjp; i++ ) {
      gyd[i] = gxj[i];
    }
  } else if( neq == 6 ) {
    for( i=0; i<mjp; i++ ) {
      gyd[i] = gyj[i];
    }
  } else if( neq == 7 ) {
    for( i=0; i<mjp; i++ ) {
      gyd[i] = gzj[i];
    }
  }
  if( lbi*lbjp != 0 ) {
    if( neq == -1 ) {
      potp2p(lbi,lbjp);
    } else if( neq == 0 ) {
      bsp2p(lbi,lbjp);
    } else if( 1 <= neq && neq <= 3 ) {
      stp2p(lbi,lbjp);
    } else if( neq == 4 ) {
      psep2p(lbi,lbjp,vis);
    } else if( 5 <= neq && neq <= 7 ) {
      fgtp2p(lbi,lbjp);
    } else if( neq == 8 ) {
      remp2p(lbi,lbjp);
    } else if( neq == 9 ) {
      vist = 1/sqrt(2*vis*dt);
      wallp2p(lbi,lbjp,vist);
    }
  }

  toc = tic;
  tic = get_time();
  tfmm[1] += tic-toc;

  lev = lmax;
  if( neq <= 3 ) {

    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;

// Step 2. P2M-expansion @ lev lmax

    boxdataj(0,mj,lev,lbj,rb);
    p2m(nmp,mp,lbj,rb,0);

    toc = tic;
    tic = get_time();
    tfmm[2] += tic-toc;

// Step 3. M2M-translation @ lev lmax-1 to 2

    if(lmax > lmin) {
      for( lev=lmax-1; lev>=lmin; lev-- ) {
        for( jj=0; jj<lbj; jj++ ) {
          nfo[jj] = nfj[jj];
        }
        lbjo = lbj;
        nlbj[lev-1] = nlbj[lev]+lbj;
        boxdataj(0,mj,lev,lbj,rb);
        m2m1(nmp,mp,lev,lbj,lbjo,nlbj,rb);
      }
      lev = lmin;
    }
    for( ii=0; ii<8; ii++ ) {
      nei[ii] = ii;
      nfi[ii] = ii;
      nfj[ii] = ii;
    }
    lbi = 8;
    lbj = 8;

    toc = tic;
    tic = get_time();
    tfmm[3] += tic-toc;

    icbox(0,mi,mj,nmp,lbi,lbj,lbjr,lbjp,lev,0,npb,tfmm);

    toc = tic;
    tic = get_time();

// Step 4. M2M-translation @ lev 1 to sublev 1

    if( npb > 0 ) {
      ib = 1+nbnet;
      m2m2(nmp,mp,ib,nlbj[0],2,rb);
      ib = 2+nbnet;
      m2m2(nmp,mp,ib,nlbj[0],3,rb);
      rb *= 2;
    }

// Step 5. M2M-translation @ sublev 2 to npb

    for( ipb=2; ipb<=npb; ipb++ ) {
      ib = 2*ipb-1+nbnet;
      jb = 2*ipb-3+nbnet;
      m2m2(nmp,mp,ib,jb,6,rb);
      ib = 2*ipb+nbnet;
      jb = 2*ipb-2+nbnet;
      m2m2(nmp,mp,ib,jb,7,rb);
      rb *=2;
    }

    toc = tic;
    tic = get_time();
    tfmm[3] += tic-toc;

// Step 6. M2L-translation @ sublev npb

    jcbox(lbj,lbjp,lev,npb,npb);
    ijbox(lbi,lbjp,lev,npb,npb);
    if( npb != 0 ) {
      for( jj=0; jj<lbjp; jj++ ) {
        boxc(nfj[jj],3,nc);
        njb[jj] = 2*npb-2+nc[2]+nbnet;
      }
    }

    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;

    m2l(nmp,mp,lbi,lbjp,lev,lbi,rb);

    toc = tic;
    tic = get_time();
    tfmm[5] += tic-toc;

// Step 7. L2L-translation @ sublev npb-1 to 1

    for( ipb=npb-1; ipb>=1; ipb-- ) {
      rb *= 0.5;
      l2l2(nmp,mp,lbi,rb);

      toc = tic;
      tic = get_time();
      tfmm[6] += tic-toc;

// Step 8. M2L-translation @ sublev npb-1 to 1

      jcbox(lbj,lbjp,lev,ipb,npb);
      ijbox(lbi,lbjp,lev,ipb,npb);
      for( jj=0; jj<lbjp; jj++ ) {
        boxc(nfj[jj],3,nc);
        njb[jj] = 2*ipb-2+nc[2]+nbnet;
      }

      toc = tic;
      tic = get_time();
      tfmm[8] += tic-toc;

      m2l(nmp,mp,lbi,lbjp,lev,0,rb);

      toc = tic;
      tic = get_time();
      tfmm[5] += tic-toc;

    }
    ipb = 0;

// Step 9. L2L-translation @ sublev 1 to lev 1

    if( npb > 0 ) {
      rb *= 0.5;
      l2l2(nmp,mp,lbi,rb);

      toc = tic;
      tic = get_time();
      tfmm[6] += tic-toc;

// Step 10. M2L-translation @ sublev 1 to lev 1

      jcbox(lbj,lbjp,lev,ipb,npb);
      ijbox(lbi,lbjp,lev,ipb,npb);

      toc = tic;
      tic = get_time();
      tfmm[8] += tic-toc;

      m2l(nmp,mp,lbi,lbjp,lev,0,rb);

      toc = tic;
      tic = get_time();
      tfmm[5] += tic-toc;

    }
    ipb = -1;

// Step 11. L2L-translation @ lev 2 to lmax

    if( lmax > lmin ) {
      for( lev=lmin+1; lev<=lmax; lev++ ) {
        for( ii=0; ii<int(pow(8,lev-1)); ii++ ) neo[ii] = nei[ii];
        boxdatai(0,mi,lev,lbi,rb);

        l2l1(nmp,mp,lbi,rb);

        boxdataj(0,mj,lev,lbj,rb);

        toc = tic;
        tic = get_time();
        tfmm[6] += tic-toc;

// Step 12. M2L-translation @ lev 2 to lmax

        icbox(0,mi,mj,nmp,lbi,lbj,lbjr,lbjp,lev,ipb,npb,tfmm);

        toc = tic;
        tic = get_time();

        m2l(nmp,mp,lbi,lbjp,lev,0,rb);

        toc = tic;
        tic = get_time();
        tfmm[5] += tic-toc;

      }
      lev = lmax;
    }

// Step 13. L2P-expansion @ lev lmax

    if( neq == 0 ) {
      bsl2p(nmp,mp,lbi,rb);
    } else if( 1 <= neq && neq <= 3) {
      stl2p(nmp,mp,lbi,rb);
    } else {
      potl2p(nmp,mp,lbi,rb);
    }

    toc = tic;
    tic = get_time();
    tfmm[7] += tic-toc;

  }

  unsorti(mi);
  unsortj(mj);

  delete[] nbi;
  for( i=0; i<2; i++ ) delete[] ndi[i];
  delete[] ndi;
  delete[] nei;
  delete[] nfi;
  delete[] nbj;
  for( i=0; i<2; i++ ) delete[] ndj[i];
  delete[] ndj;
  delete[] nej;
  delete[] nfj;
  delete[] nlbj;
  delete[] neo;
  delete[] nfo;
  for( i=0; i<3; i++ ) delete[] nxs[i];
  delete[] nxs;
  delete[] nfn;
  delete[] nij;
  delete[] njb;
  delete[] nnp;
  for( i=0; i<nebmp; i++ ) delete[] neij[i];
  delete[] neij;
  for( i=0; i<nbnp; i++ ) delete[] nsij[i];
  delete[] nsij;
  for( i=0; i<3; i++ ) delete[] npx[i];
  delete[] npx;
  delete[] nfrom;
  delete[] nc;
  delete[] nd;
  delete[] gxd;
  delete[] gyd;
  delete[] gzd;
  delete[] vd;
  delete[] nek;
  delete[] ixadj;
  delete[] nxadj;
  delete[] nadjncy;
  delete[] nvwgt;
  delete[] nadjwgt;
  delete[] nvtxdist;
  delete[] npart;
  delete[] npartd;
  delete[] irank;
  delete[] na;
  delete[] nb;
  mem = npmax*14*4+nbmax*4*4+nbne*8*4+nbnp*11*4+nbnep*2*4+lmax*4+3*2*4+nebmp*nbnp*2*4;
  memoryfree();

  dealloc();

  toc = tic;
  tic = get_time();
  tfmm[8] += tic-toc;

}
