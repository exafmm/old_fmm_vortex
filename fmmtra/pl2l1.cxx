#include "../misc/constants.h"

extern int *nfi,*neo;
extern double (*px)[nspm],(*pxo)[nspm],(*py)[nspm],(*pyo)[nspm],(*pz)[nspm],(*pzo)[nspm];

extern void pl2l(double*, double*, double*, double*, double*, double*, int, int);

void l2l1(int nmp, int mp, int lbi, double rb) {
  int lbio,ii,i,nfip,nfic,ib;
  double ppx[nspm],ppy[nspm],ppz[nspm],ppxd[nspm],ppyd[nspm],ppzd[nspm];

  lbio = lbi;
  if( lbio < 8 ) lbio = 8;
  for( ii=0; ii<lbio; ii++ ) {
    for( i=0; i<nmp; i++ ) {
      pxo[ii][i] = px[ii][i];
      pyo[ii][i] = py[ii][i];
      pzo[ii][i] = pz[ii][i];
    }
  }
  for( ii=0; ii<lbi; ii++ ) {
    nfip = nfi[ii]/8;
    nfic = nfi[ii]%8;
    ib = neo[nfip];
    for( i=0; i<nmp; i++ ) {
      ppxd[i] = pxo[ib][i];
      ppyd[i] = pyo[ib][i];
      ppzd[i] = pzo[ib][i];
    }
    pl2l(ppx,ppy,ppz,ppxd,ppyd,ppzd,nmp,nfic);
    for( i=0; i<nmp; i++ ) {
      px[ii][i] = ppx[i];
      py[ii][i] = ppy[i];
      pz[ii][i] = ppz[i];
    }
  }
}
