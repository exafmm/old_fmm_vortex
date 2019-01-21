#include "../misc/constants.h"

extern int *nfi;
extern double (*px)[nspm],(*pxo)[nspm],(*py)[nspm],(*pyo)[nspm],(*pz)[nspm],(*pzo)[nspm];

extern void pl2l(double*, double*, double*, double*, double*, double*, int, int);

void l2l2(int nmp, int mp, int lbi, double rb) {
  int ii,i,nfic;
  double ppx[nspm],ppy[nspm],ppz[nspm],ppxd[nspm],ppyd[nspm],ppzd[nspm];

  for( ii=0; ii<lbi; ii++ ) {
    for( i=0; i<nmp; i++ ) {
      pxo[ii][i] = px[nfi[ii]][i];
      pyo[ii][i] = py[nfi[ii]][i];
      pzo[ii][i] = pz[nfi[ii]][i];
    }
  }
  for( ii=0; ii<lbi; ii++ ) {
    nfic = 7-nfi[ii];
    for( i=0; i<nmp; i++ ) {
      ppxd[i] = pxo[ii][i];
      ppyd[i] = pyo[ii][i];
      ppzd[i] = pzo[ii][i];
    }
    pl2l(ppx,ppy,ppz,ppxd,ppyd,ppzd,nmp,nfic);
    for( i=0; i<nmp; i++ ) {
      px[ii][i] = ppx[i];
      py[ii][i] = ppy[i];
      pz[ii][i] = ppz[i];
    }
  }
}
