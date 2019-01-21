#include "../misc/constants.h"

extern float *xj,*yj,*zj;
extern int **ndj,**nxs,*nfn,*na,*nb;

extern void boxn(int, int, int);
extern void sort(int);

void boxdatak(int n2, int n3, int lev) {
  int kl,n,i,j,lbj,nbc;
  double rb;

  kl = int(pow(2,lev));
  rb = rd/kl;
  n = 0;
  for( i=n2; i<n3; i++ ) {
    nxs[0][n] = int((xj[i]-xmin)/rb);
    nxs[1][n] = int((yj[i]-ymin)/rb);
    nxs[2][n] = int((zj[i]-zmin)/rb);
    for( j=0; j<3; j++ ) {
      if( nxs[j][n] >= pow(2,lev) ) nxs[j][n]--;
    }
    n++;
  }

  boxn(n,3,lev);
  for( i=0; i<n; i++ ) {
    na[i] = nfn[i];
    nb[i] = i+n2;
  }
  sort(n);

  lbj = 0;
  nbc = -1;
  for( i=0; i<n; i++ ) {
    if( na[i] != nbc ) {
      ndj[0][lbj] = i+n2;
      if( lbj > 0 ) ndj[1][lbj-1] = i+n2-1;
      nbc = na[i]; 
      lbj++;
    }
  }
  if( lbj > 0 ) ndj[1][lbj-1] = n3-1;
}
