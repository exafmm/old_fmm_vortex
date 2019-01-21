#include "../misc/constants.h"

extern int *nfi,*nc;
extern std::complex<double> (*ax)[mpsym],(*axo)[mpsym],(*ay)[mpsym],(*ayo)[mpsym],(*az)[mpsym],(*azo)[mpsym];

extern void boxc(int, int, int*);
extern void boxn1(int*, int&, int);
extern void fl2l(std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, double, int, int);

void l2l2(int nmp, int mp, int lbi, double rb) {
  int ii,i,je;
  double rh;
  std::complex<double> aax[mpsym],aay[mpsym],aaz[mpsym],aaxd[mpsym],aayd[mpsym],aazd[mpsym];

  for( ii=0; ii<lbi; ii++ ) {
    for( i=0; i<nmp; i++ ) {
      axo[ii][i] = ax[nfi[ii]][i];
      ayo[ii][i] = ay[nfi[ii]][i];
      azo[ii][i] = az[nfi[ii]][i];
    }
  }
  for( ii=0; ii<lbi; ii++ ) {
    boxc(nfi[ii],3,nc);
    nc[0] = 4-nc[0]*2;
    nc[1] = 4-nc[1]*2;
    nc[2] = 4-nc[2]*2;
    boxn1(nc,je,3);
    rh = rb*sqrt(3.0)/2;
    for( i=0; i<nmp; i++ ) {
      aaxd[i] = axo[ii][i];
      aayd[i] = ayo[ii][i];
      aazd[i] = azo[ii][i];
    }
    fl2l(aax,aay,aaz,aaxd,aayd,aazd,rh,mp,je);
    for( i=0; i<nmp; i++ ) {
      ax[ii][i] = aax[i];
      ay[ii][i] = aay[i];
      az[ii][i] = aaz[i];
    }
  }
}
