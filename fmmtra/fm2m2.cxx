#include "../misc/constants.h"

extern int *nc;
extern std::complex<double> (*bx)[mpsym],(*by)[mpsym],(*bz)[mpsym];

extern void boxc(int, int, int*);
extern void boxn1(int*, int&, int);
extern void fm2m(std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, double, int, int);

void m2m2(int nmp, int mp, int ip, int jp, int nc3, double rb) {
  int mc3,i,j,ib,jb,je;
  double rh;
  std::complex<double> bbx[mpsym],bby[mpsym],bbz[mpsym],bbxd[mpsym],bbyd[mpsym],bbzd[mpsym];

  mc3 = nc3;
  if( 2 <= nc3 && nc3 <= 7 ) {
    for( i=0; i<nmp; i++ ) {
      bx[ip-1][i] = 0;
      by[ip-1][i] = 0;
      bz[ip-1][i] = 0;
    }
    mc3 = nc3-2;
  } else if ( nc3 == 9 ) {
    for( i=0; i<8; i++ ) {
      for( j=0; j<nmp; j++ ) {
        bx[i+3*nbnet][j] = 0;
        by[i+3*nbnet][j] = 0;
        bz[i+3*nbnet][j] = 0;
      }
    }
  }
  for( i=0; i<8; i++ ) {
    boxc(i,3,nc);
// M2M @ lev 1 to sublev 1
    if( nc[2] == mc3 ) {
      ib = ip-1;
      jb = i+jp;
      nc[0] = nc[0]*2+2;
      nc[1] = nc[1]*2+2;
      nc[2] = nc[2]*2+2;
      boxn1(nc,je,3);
      rh = rb*sqrt(3.0)/2;
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
// M2M @ sublev 2 to npb
    } else if ( 4 <= nc3 && nc3 <= 5 ) {
      ib = ip-1;
      if( nc[2] == nc3-5 ) {
        jb = jp-1;
      } else {
        jb = jp-2;
      }
      nc[0] = nc[0]*2+2;
      nc[1] = nc[1]*2+2;
      nc[2] = nc[2]*2+2;
      boxn1(nc,je,3);
      rh = rb*sqrt(3.0)/2;
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
// M2M @ sublev 2 to npb (channel)
    } else if ( 6 <= nc3 && nc3 <= 7 ) {
      if( nc[2] == nc3-6 ) {
        ib = ip-1;
        jb = jp-1;
        nc[0] = nc[0]*2+2;
        nc[1] = nc[1]*2+2;
        nc[2] = nc[2]*2+2;
        boxn1(nc,je,3);
        rh = rb*sqrt(3.0)/2;
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
// save bx for M2L @ lev 1
    } else if ( 8 <= nc3 && nc3 <= 9 ) {
      if( nc[2] == nc3-8 ) {
        ib = i+3*nbnet;
        jb = i+jp;
        for( j=0; j<nmp; j++ ) {
          bx[ib][j] += bx[jb][j];
          by[ib][j] += by[jb][j];
          bz[ib][j] += bz[jb][j];
        }
      }
    }
  }
}
