#include "../misc/constants.h"

extern int *nej,*nfo,*nc;
extern std::complex<double> (*bx)[mpsym],(*by)[mpsym],(*bz)[mpsym];

extern void boxc(int, int, int*);
extern void boxn1(int*, int&, int);
extern void fm2m(std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, double, int, int);

void m2m1(int nmp, int mp, int lev, int lbj, int lbjo, int* nlbj, double rb) {
  int ii,ib,j,jj,nfjp,nfjc,jb,je;
  double rh;
  std::complex<double> bbx[mpsym],bby[mpsym],bbz[mpsym],bbxd[mpsym],bbyd[mpsym],bbzd[mpsym];

  if( lev == 1 ) lbj = 8;
  for( ii=0; ii<lbj; ii++ ) {
    ib = ii+nlbj[lev-1];
    for( j=0; j<nmp; j++ ) {
      bx[ib][j] = 0;
      by[ib][j] = 0;
      bz[ib][j] = 0;
    }
  }
  for( jj=0; jj<lbjo; jj++ ) {
    nfjp = nfo[jj]/8;
    nfjc = nfo[jj]%8;
    if( lev == 1 ) {
      ib = nfjp+nlbj[lev-1];
    } else {
      ib = nej[nfjp]+nlbj[lev-1];
    }
    jb = jj+nlbj[lev];
    boxc(nfjc,3,nc);
    nc[0] = 4-nc[0]*2;
    nc[1] = 4-nc[1]*2;
    nc[2] = 4-nc[2]*2;
    boxn1(nc,je,3);
    rh = rb*sqrt(3.0)/4;
    for( j=0; j<nmp; j++ ) {
      bbxd[j] = bx[jb][j];
      bbyd[j] = by[jb][j];
      bbzd[j] = bz[jb][j];
    }
    fm2m(bbx,bby,bbz,bbxd,bbyd,bbzd,rh,mp,je);
    for( j=0; j<nmp; j++ ) {
      bx[ib][j] += bbx[j];
      by[ib][j] += bby[j];
      bz[ib][j] += bbz[j];
    }
  }
}
