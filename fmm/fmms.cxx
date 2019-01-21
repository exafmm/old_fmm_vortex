#include "../misc/parameters.h"
#include "../misc/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern int *nbi,**ndi,*nei,*nfi,*nbj,**ndj,*nej,*nfj,*nlbj;
extern int **ndl,*nel,*nfl,*nlbl,**ndm,*nem,*nfm,*nlbm,*neo,*nfo;
extern int **nxs,*nfn,*nij,*njb,*nnp,**neij,**nsij,**npx,*nfrom,*nc,*nd;
extern float *xd,*gxd,*gyd,*gzd,*vd;
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
extern void boxdatak(int, int, int);
extern void jsbox(int, int, int, int&, int, int, int);
extern void ijbox(int, int, int, int, int);
extern void isbox(int, int, int, int&, int, int, int, int, int, int&, int&, int, int, int, int, double*);
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
  int neq,lev,lmin,nebmp,nbnep,i,nmp,lbi,nbyte,ii,ipb,lbj,j,lbl,lbm,lbjp,ic,jj,mjp,lbjr,lbjo,ib,npbo,npbs,jpb,jz;
  double tic,toc,rb,vist;
  tic = get_time();

  neq = neqd;
  nlevel(mj);
  boxallocate(0,mi,0,mj);
  lev = lmax;
  lmin = 1;
  nbmax = int(pow(8,lmax));
  nbnes = 3*nbnet+int(pow(2,npb+2))+4;
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
  ndl = new int* [2];
  for( i=0; i<2; i++ ) ndl[i] = new int [nbnp];
  nel = new int [nbmax];
  nfl = new int [nbnp];
  nlbl = new int [lmax];
  ndm = new int* [2];
  for( i=0; i<2; i++ ) ndm[i] = new int [nbnp];
  nem = new int [nbmax];
  nfm = new int [nbnp];
  nlbm = new int [lmax];
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
  xd = new float [npmax];
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
  mem = npmax*15*4+nbmax*6*4+nbne*8*4+nbnp*17*4+nbnep*2*4+lmax*3*4+3*2*4;
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
  nlbl[lev-1] = nbnet;
  nlbm[lev-1] = 2*nbnet;

  boxpart(mi,mj,lev);

  sorti(mi);
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

  for( i=0; i<mj; i++ ) xd[i] = xj[i];

// Step 1. P2M-expansion and M2M-translation of lower shear periodic boxes

  if( neq<=3 && npb!=0 ) {

    for( ipb=1; ipb<=int(pow(2,npb-1)); ipb++ ) {

// Step 1.1 Initialize geometry for lower box ipb (bx)

      for( i=0; i<mj; i++ ) {
        xj[i] = xd[i]-stx*(pow(2,npb)-2*ipb+2);
        while( xj[i] < -pi ) {
          xj[i] += 2*pi;
        }
      }
      lev = lmax;

      toc = tic;
      tic = get_time();
      tfmm[8] += tic-toc;

// Step 1.2 P2M-expansion @ lev lmax for lower box ipb (bx)

      sortj(mj);
      boxdataj(0,mj,lev,lbj,rb);
      p2m(nmp,mp,lbj,rb,0);

      toc = tic;
      tic = get_time();
      tfmm[2] += tic-toc;

// Step 1.3 M2M-translation @ lev lmax to 1 for lower box ipb (bx)

      if(lmax > lmin) {
        for( lev=lmax-1; lev>=lmin; lev-- ) {
          for( jj=0; jj<lbj; jj++ ) nfo[jj] = nfj[jj];
          lbjo = lbj;
          nlbj[lev-1] = nlbj[lev]+lbj;
          boxdataj(0,mj,lev,lbj,rb);
          m2m1(nmp,mp,lev,lbj,lbjo,nlbj,rb);
        }
        lev = lmin;
      }
      unsortj(mj);

      toc = tic;
      tic = get_time();
      tfmm[3] += tic-toc;

      isbox(0,mi,mj,mj,nmp,lbi,lbj,lbl,lbm,lbjr,lbjp,lev,ipb,npb,0,tfmm);

      toc = tic;
      tic = get_time();

// Step 1.4 M2M-translation @ lev 1 to sublev 1 for lower box ipb (bx)

      if( ipb != 1 ) {
        ib = 2*ipb+6+3*nbnet;
        m2m2(nmp,mp,ib,nlbj[0],2,rb);
        m2m2(nmp,mp,ib,nlbm[0],1,rb);
      }
      if( ipb == int(pow(2,npb-1)) ) {
        m2m2(nmp,mp,ib,nlbj[0],9,rb);
      }

      toc = tic;
      tic = get_time();
      tfmm[3] += tic-toc;

// Step 1.5 Initialize geometry for lower box ipb (bxm)

      for( i=0; i<mj; i++ ) {
        xj[i] = xd[i]-stx*(pow(2,npb)-2*ipb+1);
        while( xj[i] < -pi ) xj[i] += 2*pi;
      }
      lev = lmax;

      toc = tic;
      tic = get_time();
      tfmm[8] += tic-toc;

// Step 1.6 P2M-expansion @ lev lmax for lower box ipb (bxm)

      sortj(mj);
      boxdataj(0,mj,lev,lbj,rb);
      lbm = lbj;
      for( jj=0; jj<lbm; jj++ ) nfm[jj] = nfj[jj];
      p2m(nmp,mp,lbj,rb,2*nbnet);

      toc = tic;
      tic = get_time();
      tfmm[2] += tic-toc;

// Step 1.7 M2M-translation @ lev lmax to 1 for lower box ipb (bxm)

      if(lmax > lmin) {
        for( lev=lmax-1; lev>=lmin; lev-- ) {
          for( jj=0; jj<lbj; jj++ ) nfo[jj] = nfj[jj];
          lbjo = lbj;
          nlbm[lev-1] = nlbm[lev]+lbj;
          boxdataj(0,mj,lev,lbj,rb);
          m2m1(nmp,mp,lev,lbj,lbjo,nlbm,rb);
        }
        lev = lmin;
      }
      unsortj(mj);

      toc = tic;
      tic = get_time();
      tfmm[3] += tic-toc;

      isbox(0,mi,mj,mj,nmp,lbi,lbj,lbl,lbm,lbjr,lbjp,lev,ipb,npb,-1,tfmm);

      toc = tic;
      tic = get_time();

// Step 1.8 M2M-translation @ lev 1 to sublev 1 for lower box ipb (bxm)

      ib = 2*ipb+7+3*nbnet;
      m2m2(nmp,mp,ib,nlbm[0],2,rb);
      m2m2(nmp,mp,ib,nlbj[0],1,rb);

      toc = tic;
      tic = get_time();
      tfmm[3] += tic-toc;

    }

// Step 2. P2M-expansion and M2M-translation of upper shear periodic boxes

    for( ipb=int(pow(2,npb)); ipb>=int(pow(2,npb-1))+1; ipb-- ) {

// Step 2.1 Initialize geometry for upper box ipb (bx)

      for( i=0; i<mj; i++ ) {
        xj[i] = xd[i]+stx*(2*ipb-pow(2,npb));
        while( xj[i] > pi ) {
          xj[i] -= 2*pi;
        }
      }
      lev = lmax;

      toc = tic;
      tic = get_time();
      tfmm[8] += tic-toc;

// Step 2.2 P2M-expansion @ lev lmax for upper box ipb (bx)

      sortj(mj);
      boxdataj(0,mj,lev,lbj,rb);
      p2m(nmp,mp,lbj,rb,0);

      toc = tic;
      tic = get_time();
      tfmm[2] += tic-toc;

// Step 2.3 M2M-translation @ lev lmax to 1 for upper box ipb (bx)

      if(lmax > lmin) {
        for( lev=lmax-1; lev>=lmin; lev-- ) {
          for( jj=0; jj<lbj; jj++ ) nfo[jj] = nfj[jj];
          lbjo = lbj;
          nlbj[lev-1] = nlbj[lev]+lbj;
          boxdataj(0,mj,lev,lbj,rb);
          m2m1(nmp,mp,lev,lbj,lbjo,nlbj,rb);
        }
        lev = lmin;
      }
      unsortj(mj);

      toc = tic;
      tic = get_time();
      tfmm[3] += tic-toc;

      isbox(0,mi,mj,mj,nmp,lbi,lbj,lbl,lbm,lbjr,lbjp,lev,ipb,npb,0,tfmm);

      toc = tic;
      tic = get_time();

// Step 2.4 M2M-translation @ lev 1 to sublev 1 for upper box ipb (bx)

      if( ipb != int(pow(2,npb)) ) {
        ib = 2*ipb+9+3*nbnet;
        m2m2(nmp,mp,ib,nlbj[0],3,rb);
        m2m2(nmp,mp,ib,nlbl[0],0,rb);
      }
      if( ipb == int(pow(2,npb-1))+1 ) {
        m2m2(nmp,mp,ib,nlbj[0],8,rb);
      }

      toc = tic;
      tic = get_time();
      tfmm[3] += tic-toc;

// Step 2.5 Initialize geometry for upper box ipb (bxl)

      for( i=0; i<mj; i++ ) {
        xj[i] = xd[i]+stx*(2*ipb-pow(2,npb)-1);
        while( xj[i] > pi ) xj[i] -= 2*pi;
      }
      lev = lmax;

      toc = tic;
      tic = get_time();
      tfmm[8] += tic-toc;

// Step 2.6 P2M-expansion @ lev lmax for upper box ipb (bxl)

      sortj(mj);
      boxdataj(0,mj,lev,lbj,rb);
      lbl = lbj;
      for( jj=0; jj<lbl; jj++ ) nfl[jj] = nfj[jj];
      p2m(nmp,mp,lbj,rb,nbnet);

      toc = tic;
      tic = get_time();
      tfmm[2] += tic-toc;

// Step 2.7 M2M-translation @ lev lmax to 1 for upper box ipb (bxl)

      if(lmax > lmin) {
        for( lev=lmax-1; lev>=lmin; lev-- ) {
          for( jj=0; jj<lbj; jj++ ) nfo[jj] = nfj[jj];
          lbjo = lbj;
          nlbl[lev-1] = nlbl[lev]+lbj;
          boxdataj(0,mj,lev,lbj,rb);
          m2m1(nmp,mp,lev,lbj,lbjo,nlbl,rb);
        }
        lev = lmin;
      }
      unsortj(mj);

      toc = tic;
      tic = get_time();
      tfmm[3] += tic-toc;

      isbox(0,mi,mj,mj,nmp,lbi,lbj,lbl,lbm,lbjr,lbjp,lev,ipb,npb,1,tfmm);

      toc = tic;
      tic = get_time();

// Step 2.8 M2M-translation @ lev 1 to sublev 1 for upper box ipb (bxl)

      ib = 2*ipb+8+3*nbnet;
      m2m2(nmp,mp,ib,nlbl[0],3,rb);
      m2m2(nmp,mp,ib,nlbj[0],0,rb);

      toc = tic;
      tic = get_time();
      tfmm[3] += tic-toc;

    }

  } else {

// Step 1. P2M-expansion and M2M-translation of lower shear periodic boxes

// Step 1.1 Initialize geometry for lower box (bxm)

    for( i=0; i<mj; i++ ) {
      xj[i] = xd[i]-stx;
      while( xj[i] < -pi ) xj[i] +=2*pi;
    }
    lev = lmax;

    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;

// Step 1.2 P2M-expansion @ lev lmax for lower box (bxm)

    sortj(mj);
    boxdataj(0,mj,lev,lbj,rb);
    lbm = lbj;
    for( jj=0; jj<lbm; jj++ ) nfm[jj] = nfj[jj];

    if( neq <= 3 ) {
      p2m(nmp,mp,lbj,rb,2*nbnet);

      toc = tic;
      tic = get_time();
      tfmm[2] += tic-toc;
 
// Step 1.3 M2M-translation @ lev lmax to 1 for lower box (bxm)

      if(lmax > lmin) {
        for( lev=lmax-1; lev>=lmin; lev-- ) {
          for( jj=0; jj<lbj; jj++ ) nfo[jj] = nfj[jj];
          lbjo = lbj;
          nlbm[lev-1] = nlbm[lev]+lbj;
          boxdataj(0,mj,lev,lbj,rb);
          m2m1(nmp,mp,lev,lbj,lbjo,nlbm,rb);
        }
        lev = lmin;
      }
    }
    unsortj(mj);

    toc = tic;
    tic = get_time();
    tfmm[3] += tic-toc;

    if( neq <= 3 ) {
      isbox(0,mi,mj,mj,nmp,lbi,lbj,lbl,lbm,lbjr,lbjp,lev,ipb,npb,-1,tfmm);
    }

    toc = tic;
    tic = get_time();

// Step 2. P2M-expansion and M2M-translation of upper shear periodic boxes

// Step 2.1 Initialize geometry for upper box ipb (bxl)

    for( i=0; i<mj; i++ ) {
      xj[i] = xd[i]+stx;
      while( xj[i] > pi ) xj[i] -= 2*pi;
    }
    lev = lmax;

    toc = tic;
    tic = get_time();
    tfmm[8] += tic-toc;

// Step 2.2 P2M-expansion @ lev lmax for upper box ipb (bxl)

    sortj(mj);
    boxdataj(0,mj,lev,lbj,rb);
    lbl = lbj;
    for( jj=0; jj<lbl; jj++ ) nfl[jj] = nfj[jj];

    if( neq <= 3 ) {
      p2m(nmp,mp,lbj,rb,nbnet);

      toc = tic;
      tic = get_time();
      tfmm[2] += tic-toc;

// Step 2.3 M2M-translation @ lev lmax to 1 for upper box ipb (bxl)

      if(lmax > lmin) {
        for( lev=lmax-1; lev>=lmin; lev-- ) {
          for( jj=0; jj<lbj; jj++ ) nfo[jj] = nfj[jj];
          lbjo = lbj;
          nlbl[lev-1] = nlbl[lev]+lbj;
          boxdataj(0,mj,lev,lbj,rb);
          m2m1(nmp,mp,lev,lbj,lbjo,nlbl,rb);
        }
        lev = lmin;
      }
    }
    unsortj(mj);

    toc = tic;
    tic = get_time();
    tfmm[3] += tic-toc;

    if( neq <= 3 ) {
      isbox(0,mi,mj,mj,nmp,lbi,lbj,lbl,lbm,lbjr,lbjp,lev,ipb,npb,1,tfmm);
    }

    toc = tic;
    tic = get_time();

  }

// Step 3. P2P-summation @ lev lmax

// Step 3.1 Shift sheared neighbor particles

  for( i=0; i<mj; i++ ) xj[i] = xd[i];
  lev = lmax;
  ic = 2*mj;
  sortj(mj);
  boxdataj(0,mj,lev,lbj,rb);
  jsbox(lbj,lbl,lbm,lbjp,lev,-2,npb);
  ijbox(lbi,lbjp,lev,-2,npb);
  for( jj=lbj; jj<lbjp; jj++ ) {
    if( npx[2][jj] == 0 ) {
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
  }
  unsortj(mj);
  for( i=0; i<mj; i++ ) {
    xj[i] = xd[i]+stx;
    while( xj[i] > pi ) xj[i] -= 2*pi;
  }
  sortj(mj);
  boxdatak(0,mj,lev);
  for( jj=lbj; jj<lbjp; jj++ ) {
    if( npx[2][jj] == 1 ) {
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
  }
  unsortj(mj);
  for( i=0; i<mj; i++ ) {
    xj[i] = xd[i]-stx;
    while( xj[i] < -pi ) xj[i] += 2*pi;
  }
  sortj(mj);
  boxdatak(0,mj,lev);
  for( jj=lbj; jj<lbjp; jj++ ) {
    if( npx[2][jj] == -1 ) {
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
  }
  unsortj(mj);
  for( i=0; i<mj; i++ ) xj[i] = xd[i];
  sortj(mj);
  boxdatak(0,mj,lev);
  mjp = ic;

  toc = tic;
  tic = get_time();
  tfmm[8] += tic-toc;

// Step 3.1 Send and recv particles on other nodes

  isbox(0,mi,mj,mjp,nmp,lbi,lbj,lbl,lbm,lbjr,lbjp,lev,-2,npb,0,tfmm);

  toc = tic;
  tic = get_time();

// Step 3.2 P2P-summation

  boxdataj(0,mj,lev,lbj,rb);
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
  unsortj(mj);

  toc = tic;
  tic = get_time();
  tfmm[1] += tic-toc;

  lev = lmax;
  if( neq <= 3 ) {

// Step 4. P2M-expansion and M2M-translation of center box

// Step 4.1 P2M-expansion @ lev lmax for center box

    sortj(mj);
    boxdataj(0,mj,lev,lbj,rb);
    p2m(nmp,mp,lbj,rb,0);

    toc = tic;
    tic = get_time();
    tfmm[2] += tic-toc;

// Step 4.2 M2M-translation @ lev lmax to 1 for center box

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
    unsortj(mj);
    boxdatai(0,mi,lev,lbi,rb);
    for( ii=0; ii<8; ii++ ) {
      nfi[ii] = ii;
      nfl[ii] = ii;
      nfj[ii] = ii;
      nfm[ii] = ii;
    }
    lbi = 8;
    lbl = 8;
    lbj = 8;
    lbm = 8;

    toc = tic;
    tic = get_time();
    tfmm[3] += tic-toc;

    isbox(0,mi,mj,mj,nmp,lbi,lbj,lbl,lbm,lbjr,lbjp,lev,0,npb,0,tfmm);

    toc = tic;
    tic = get_time();

    if( npb > 0 ) {

// Step 4.3 M2M-translation @ lev 1 to sublev 1 for lower and center box

      ib = int(pow(2,npb))+8+3*nbnet;
      m2m2(nmp,mp,ib,nlbj[0],2,rb);
      m2m2(nmp,mp,ib,nlbm[0],1,rb);

// Step 4.4 M2M-translation @ lev 1 to sublev 1 for upper and center box

      ib = int(pow(2,npb))+9+3*nbnet;
      m2m2(nmp,mp,ib,nlbj[0],3,rb);
      m2m2(nmp,mp,ib,nlbl[0],0,rb);
      rb *= 2;
      npbo = 8+3*nbnet;
    }

// Step 5. M2M-translation @ sublev 2 to npb

    for( ipb=2; ipb<=npb; ipb++ ) {
      npbs = int(pow(2,npb+2)-pow(2,npb-ipb+3))+8+3*nbnet;
      for( jpb=1; jpb<=int(pow(2,npb-ipb+2)); jpb++ ) {
        ib = jpb+npbs;
        m2m2(nmp,mp,ib,2*jpb+npbo,5,rb);
      }
      npbo = npbs;
      rb *=2;
    }

    toc = tic;
    tic = get_time();
    tfmm[3] += tic-toc;

// Step 6. M2L-translation @ sublev npb

    jsbox(lbj,lbl,lbm,lbjp,lev,npb,npb);
    ijbox(lbi,lbjp,lev,npb,npb);
    if( npb == 0 ) {
      for( jj=0; jj<lbjp; jj++ ) {
        boxc(nfj[jj],3,nc);
        jz = nc[2]+npx[2][jj]*int(pow(2,lev));
        if( jz == -1 ) {
          njb[jj] = nfj[jj]+nlbm[0];
        } else if( jz == 2 ) {
          njb[jj] = nfj[jj]+nlbl[0];
        } else {
          njb[jj] = nfj[jj]+nlbj[0];
        }
      }
    } else {
      for( jj=0; jj<lbjp; jj++ ) {
        boxc(nfj[jj],3,nc);
        njb[jj] = nc[2]+npx[2][jj]*2+1+npbo;
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
      npbs = int(pow(2,npb+2)-3*pow(2,npb-ipb+1))+4+3*nbnet;
      rb *= 0.5;
      l2l2(nmp,mp,lbi,rb);

      toc = tic;
      tic = get_time();
      tfmm[6] += tic-toc;

// Step 8. M2L-translation @ sublev npb-1 to 1

      jsbox(lbj,lbl,lbm,lbjp,lev,ipb,npb);
      ijbox(lbi,lbjp,lev,ipb,npb);
      for( jj=0; jj<lbjp; jj++ ) {
        boxc(nfj[jj],3,nc);
        njb[jj] = nc[2]+npx[2][jj]*2+3+npbs;
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
    boxdatai(0,mi,lev,lbi,rb);

// Step 9. L2L-translation @ sublev 1 to lev 1

    if( npb > 0 ) {
      l2l2(nmp,mp,lbi,rb);
      rb *= 0.5;

      toc = tic;
      tic = get_time();
      tfmm[6] += tic-toc;

      for( i=0; i<mj; i++ ) {
        xj[i] = xd[i]+stx;
        while( xj[i] > pi ) xj[i] -= 2*pi;
      }
      sortj(mj);
      boxdataj(0,mj,lev,lbl,rb);
      unsortj(mj);
      for( jj=0; jj<lbl; jj++ ) nfl[jj] = nfj[jj];
      for( i=0; i<mj; i++ ) {
        xj[i] = xd[i]-stx;
        while( xj[i] < -pi ) xj[i] += 2*pi;
      }
      sortj(mj);
      boxdataj(0,mj,lev,lbm,rb);
      unsortj(mj);
      for( jj=0; jj<lbm; jj++ ) nfm[jj] = nfj[jj];
      for( i=0; i<mj; i++ ) xj[i] = xd[i];
      sortj(mj);
      boxdataj(0,mj,lev,lbj,rb);
      unsortj(mj);

// Step 10. M2L-translation @ sublev 1 to lev 1

      jsbox(lbj,lbl,lbm,lbjp,lev,ipb,npb);
      ijbox(lbi,lbjp,lev,ipb,npb);

      for( jj=0; jj<lbjp; jj++ ) {
        boxc(nfj[jj],3,nc);
        jz = nc[2]+npx[2][jj]*2+4;
        if( jz == 1 || jz == 8 ) {
          njb[jj] = nfj[jj]+3*nbnet;
        } else if ( jz < 4 ) {
          njb[jj] = nfj[jj]+nlbm[lev-1];
        } else if ( jz < 6 ) {
          njb[jj] = nfj[jj]+nlbj[lev-1];
        } else {
          njb[jj] = nfj[jj]+nlbl[lev-1];
        }
      }

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

        toc = tic;
        tic = get_time();
        tfmm[6] += tic-toc;

        for( i=0; i<mj; i++ ) {
          xj[i] = xd[i]+stx;
          while( xj[i] > pi ) xj[i] -= 2*pi;
        }
        sortj(mj);
        boxdataj(0,mj,lev,lbl,rb);
        unsortj(mj);
        for( jj=0; jj<int(pow(8,lev)); jj++ ) nel[jj] = nej[jj];
        for( jj=0; jj<lbl; jj++ ) nfl[jj] = nfj[jj];
        for( i=0; i<mj; i++ ) {
          xj[i] = xd[i]-stx;
          while( xj[i] < -pi ) xj[i] += 2*pi;
        }
        sortj(mj);
        boxdataj(0,mj,lev,lbm,rb);
        unsortj(mj);
        for( jj=0; jj<int(pow(8,lev)); jj++ ) nem[jj] = nej[jj];
        for( jj=0; jj<lbm; jj++ ) nfm[jj] = nfj[jj];
        for( i=0; i<mj; i++ ) xj[i] = xd[i];
        sortj(mj);
        boxdataj(0,mj,lev,lbj,rb);

        toc = tic;
        tic = get_time();
        tfmm[8] += tic-toc;

// Step 12. M2L-translation @ lev 2 to lmax

        isbox(0,mi,mj,mj,nmp,lbi,lbj,lbl,lbm,lbjr,lbjp,lev,ipb,npb,0,tfmm);

        unsortj(mj);

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
  for( i=0; i<2; i++ ) delete[] ndl[i];
  delete[] ndl;
  delete[] nel;
  delete[] nfl;
  delete[] nlbl;
  for( i=0; i<2; i++ ) delete[] ndm[i];
  delete[] ndm;
  delete[] nem;
  delete[] nfm;
  delete[] nlbm;
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
  delete[] xd;
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
  mem = npmax*15*4+nbmax*6*4+nbne*8*4+nbnp*17*4+nbnep*2*4+lmax*3*4+3*2*4+nebmp*nbnp*2*4;
  memoryfree();

  dealloc();

  toc = tic;
  tic = get_time();
  tfmm[8] += tic-toc;
}
