#include "../misc/constants.h"

extern float *xi,*yi,*zi;
extern int **ndi,*nei,*nfi,**nxs,*nfn,*na,*nb;

extern void boxn(int, int, int);
extern void sort(int);

void boxdatai(int n0, int n1, int lev, int& lbi, double& rb) {
  int kl,n,i,j,nbc,ii;

  kl = int(pow(2,lev));
  rb = rd/kl;
  n = 0;
  for( i=n0; i<n1; i++ ) {
    nxs[0][n] = int((xi[i]-xmin)/rb);
    nxs[1][n] = int((yi[i]-ymin)/rb);
    nxs[2][n] = int((zi[i]-zmin)/rb);
    for( j=0; j<3; j++ ) {
      if( nxs[j][n] >= pow(2,lev) ) nxs[j][n] = nxs[j][n]-1;
    }
    n++;
  }

  boxn(n,3,lev);
  for( i=0; i<n; i++ ) {
    na[i] = nfn[i];
    nb[i] = i+n0;
  }
  sort(n);

  lbi = 0;
  nbc = -1;
  for( i=0; i<nbmax; i++ ) nei[i] = -1;
  for( i=0; i<n; i++ ) {
    if( na[i] != nbc ) {
      nei[na[i]] = lbi;
      nfi[lbi] = na[i];
      ndi[0][lbi] = i+n0;
      if( lbi > 0 ) ndi[1][lbi-1] = i+n0-1;
      nbc = na[i]; 
      lbi++;
    }
  }
  if( lbi > 0 ) ndi[1][lbi-1] = n1-1;

  if( lev == 1 ) {
    for( ii=0; ii<lbi; ii++ ) {
      na[ii] = ndi[0][ii];
      nb[ii] = ndi[1][ii];
    }
    for( ii=0; ii<8; ii++ ) {
      if( nei[ii] == -1 ) {
        ndi[0][ii] = 1;
        ndi[1][ii] = 0;
      } else {
        ndi[0][ii] = na[nei[ii]];
        ndi[1][ii] = nb[nei[ii]];
      }
      nei[ii] = ii;
      nfi[ii] = ii;
    }
    lbi = 8;
  }
}
