#include "../misc/constants.h"

extern int *nej,*nfo;
extern double (*gx)[nspm],(*gy)[nspm],(*gz)[nspm];

extern void pm2m(double*, double*, double*, double*, double*, double*, int, int);

void m2m1(int nmp, int mp, int lev, int lbj, int lbjo, int* nlbj, double rb) {
  int ii,ib,j,jj,nfjp,nfjc,jb;
  double ggx[nspm],ggy[nspm],ggz[nspm],ggxd[nspm],ggyd[nspm],ggzd[nspm];

  if( lev == 1 ) lbj = 8;
  for( ii=0; ii<lbj; ii++ ) {
    ib = ii+nlbj[lev-1];
    for( j=0; j<nmp; j++ ) {
      gx[ib][j] = 0;
      gy[ib][j] = 0;
      gz[ib][j] = 0;
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
    for( j=0; j<nmp; j++ ) {
      ggxd[j] = gx[jb][j];
      ggyd[j] = gy[jb][j];
      ggzd[j] = gz[jb][j];
    }
    pm2m(ggx,ggy,ggz,ggxd,ggyd,ggzd,nmp,nfjc);
    for( j=0; j<nmp; j++ ) {
      gx[ib][j] += ggx[j];
      gy[ib][j] += ggy[j];
      gz[ib][j] += ggz[j];
    }
  }
}
