#include "../misc/parameters.h"
#include "../misc/constants.h"

extern float *xi,*yi,*zi,*xj,*yj,*zj;
extern int *nbi,**ndi,*nei,*nfi;
extern int *nbj,**ndj,*nej,*nfj;
extern int **ndl,*nfl,**ndm,*nfm,*nfo;
extern int **nxs,*nfn,*nc,*nd,**npx,*nnp;
extern float *xd;
extern int *na,*nb;

extern void memoryuse();
extern void memoryfree();
extern void boxdatai(int, int, int, int&, double&);
extern void boxdataj(int, int, int, int&, double&);
extern void jpbox(int, int&, int, int, int);
extern void jsbox(int, int, int, int&, int, int, int);
extern void jcbox(int, int&, int, int, int);

void boxallocate(int n0, int n1, int n2, int n3) {
  int i,lbi,lbj,lev,lbjp,lbl,lbm;
  double rb;
  nbmax = int(pow(pow(2,lmax)+4,3));

  nxs = new int* [3];
  for( i=0; i<3; i++ ) nxs[i] = new int [npmax];
  nfn = new int [npmax];
  na  = new int [npmax];
  nb  = new int [npmax];
  nbi = new int [npmax];
  ndi = new int* [2];
  for( i=0; i<2; i++ ) ndi[i] = new int [nbmax];
  nei = new int [nbmax];
  nfi = new int [nbmax];

  mem = npmax*7*4+nbmax*4*4;
  memoryuse();

  boxdatai(n0,n1,lmax,lbi,rb);
  nbne = lbi;

  for( i=0; i<2; i++ ) delete[] ndi[i];
  delete[] nbi;
  delete[] ndi;
  delete[] nei;
  delete[] nfi;
  mem = npmax*4+nbmax*4*4;
  memoryfree();
  nbj = new int [npmax];
  ndj = new int* [2];
  for( i=0; i<2; i++ ) ndj[i] = new int [nbmax];
  nej = new int [nbmax];
  nfj = new int [nbmax];
  mem = npmax*4+nbmax*4*4;
  memoryuse();

  nbnet = 0;
  for( lev=1; lev<=lmax; lev++ ) {
    boxdataj(n2,n3,lev,lbj,rb);
    nbnet += lbj;
  }

  for( i=0; i<2; i++ ) delete[] ndj[i];
  delete[] nbj;
  delete[] ndj;
  delete[] nej;
  mem = npmax*4+nbmax*3*4;
  memoryfree();
  nfo = new int [nbmax];
  nc  = new int [3];
  nd  = new int [3];
  npx = new int* [3];
  for( i=0; i<3; i++ ) npx[i] = new int [nbmax];
  nnp = new int [nbmax];
  mem = nbmax*5*4+3*2*4;
  memoryuse();

  if( ngeo == 1 ) {
    jpbox(lbj,lbjp,lmax,-1,1);
  } else if( ngeo == 2 ) {
    nbj = new int [npmax];
    ndj = new int* [2];
    for( i=0; i<2; i++ ) ndj[i] = new int [nbmax];
    nej = new int [nbmax];
    ndl = new int* [2];
    for( i=0; i<2; i++ ) ndl[i] = new int [nbmax];
    nfl = new int [nbmax];
    ndm = new int* [2];
    for( i=0; i<2; i++ ) ndm[i] = new int [nbmax];
    nfm = new int [nbmax];
    xd  = new float [npmax];
    mem = npmax*2*4+nbmax*9*4;
    memoryuse();
    for( i=n2; i<=n3; i++ ) {
      xd[i] = xj[i];
      xj[i] = xd[i]+stx;
      while( xj[i] > pi ) {
        xj[i] -= 2*pi;
      }
    }
    boxdataj(n2,n3,lmax,lbl,rb);
    for( i=0; i<lbl; i++ ) {
      ndl[0][i] = ndj[0][i];
      ndl[1][i] = ndj[1][i];
      nfl[i] = nfj[i];
    }
    for( i=n2; i<=n3; i++ ) {
      xj[i] = xd[i]-stx;
      while( xj[i] < -pi ) {
        xj[i] += 2*pi;
      }
    }
    boxdataj(n2,n3,lmax,lbm,rb);
    for( i=0; i<lbm; i++ ) {
      ndm[0][i] = ndj[0][i];
      ndm[1][i] = ndj[1][i];
      nfm[i] = nfj[i];
    }
    for( i=n2; i<=n3; i++ ) {
      xj[i] = xd[i];
    }
    boxdataj(n2,n3,lmax,lbj,rb);
    jsbox(lbj,lbl,lbm,lbjp,lmax,-1,1);
    for( i=n2; i<=n3; i++ ) {
      xj[i] = xd[i];
    }
    for( i=0; i<2; i++ ) {
      delete[] ndj[i];
      delete[] ndl[i];
      delete[] ndm[i];
    }
    delete[] nbj;
    delete[] ndj;
    delete[] nej;
    delete[] ndl;
    delete[] nfl;
    delete[] ndm;
    delete[] nfm;
    delete[] xd;
    mem = npmax*2*4+nbmax*9*4;
    memoryfree();
  } else if( ngeo == 3 ) {
    jcbox(lbj,lbjp,lmax,-1,1);
  } else {
    lbjp = lbj;
  }

  if( nbne < 8 ) nbne = 8;
  nbnp = std::max(lbjp,8*8*8);

// temporary

  nbmax = int(pow(8,lmax));
  nbne = int(pow(8,lmax));
  nbnp = int(pow(pow(2,lmax)+4,3));
  nbnet = 0;
  for( lev=1; lev<=lmax; lev++ ) nbnet += int(pow(8,lev));

  for( i=0; i<3; i++ ) {
    delete[] nxs[i];
    delete[] npx[i];
  }
  delete[] nfj;
  delete[] nfo;
  delete[] nxs;
  delete[] nfn;
  delete[] nc;
  delete[] nd;
  delete[] npx;
  delete[] nnp;
  delete[] na;
  delete[] nb;
  mem = npmax*6*4+nbmax*6*4+3*2*4;
  memoryfree();
}
