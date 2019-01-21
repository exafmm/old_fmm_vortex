#include "../misc/constants.h"

extern int *nfi,*neo,*nc;
extern std::complex<double> (*ax)[mpsym],(*axo)[mpsym],(*ay)[mpsym],(*ayo)[mpsym],(*az)[mpsym],(*azo)[mpsym];

extern void boxc(int, int, int*);
extern void boxn1(int*, int&, int);
extern void fl2l(std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, double, int, int);

void l2l1(int nmp, int mp, int lbi, double rb) {
  int lbio,ii,i,nfip,nfic,ib,je;
  double rh;
  std::complex<double> aax[mpsym],aay[mpsym],aaz[mpsym],aaxd[mpsym],aayd[mpsym],aazd[mpsym];

  lbio = lbi;
  if( lbio < 8 ) lbio = 8;
  for( ii=0; ii<lbio; ii++ ) {
    for( i=0; i<nmp; i++ ) {
      axo[ii][i] = ax[ii][i];
      ayo[ii][i] = ay[ii][i];
      azo[ii][i] = az[ii][i];
    }
  }
  for( ii=0; ii<lbi; ii++ ) {
    nfip = nfi[ii]/8;
    nfic = nfi[ii]%8;
    boxc(nfic,3,nc);
    nc[0] = nc[0]*2+2;
    nc[1] = nc[1]*2+2;
    nc[2] = nc[2]*2+2;
    boxn1(nc,je,3);
    rh = rb*sqrt(3.0)/2;
    ib = neo[nfip];
    for( i=0; i<nmp; i++ ) {
      aaxd[i] = axo[ib][i];
      aayd[i] = ayo[ib][i];
      aazd[i] = azo[ib][i];
    }
    fl2l(aax,aay,aaz,aaxd,aayd,aazd,rh,mp,je);
    for( i=0; i<nmp; i++ ) {
      ax[ii][i] = aax[i];
      ay[ii][i] = aay[i];
      az[ii][i] = aaz[i];
    }
  }
}
