#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi,*gxd,*gyd,*gzd,*vd;
extern int *nbi,**nxs,*nfn,*na,*nb;

extern void boxn(int, int, int);
extern void sort(int);
extern void sortvar(int, int, float*, int*);
extern void unsortvar(int, int, float*, int*);

void sorti(int& mi) {
  int kl,i,j;
  double rb;

  kl = int(pow(2,lmax));
  rb = rd/kl;
  for( i=0; i<mi; i++ ) {
    nxs[0][i] = int((xi[i]-xmin)/rb);
    nxs[1][i] = int((yi[i]-ymin)/rb);
    nxs[2][i] = int((zi[i]-zmin)/rb);
    for( j=0; j<3; j++ ) {
      if( nxs[j][i] >= int(pow(2,lmax)) ) nxs[j][i] = nxs[j][i]-1;
    }
  }

  boxn(mi,3,lmax);
  for( i=0; i<mi; i++ ) {
    na[i] = nfn[i];
    nb[i] = i;
  }
  sort(mi);
  for( i=0; i<mi; i++ ) {
    nbi[i] = nb[i];
  }

  sortvar(0,mi,xi,nbi);
  sortvar(0,mi,yi,nbi);
  sortvar(0,mi,zi,nbi);
  sortvar(0,mi,vi,nbi);
  if( neq < 5 || 7 < neq ) {
    sortvar(0,mi,gxi,nbi);
    sortvar(0,mi,gyi,nbi);
    sortvar(0,mi,gzi,nbi);
  }

}

void unsorti(int& mi) {
  int i;

  if( neq == 5 ) {
    for( i=0; i<npmax; i++ ) gxi[i] = 0;
    for( i=0; i<mi; i++ ) gxi[nbi[i]] = gxd[i];
  } else if( neq == 6 ) {
    for( i=0; i<npmax; i++ ) gyi[i] = 0;
    for( i=0; i<mi; i++ ) gyi[nbi[i]] = gxd[i];
  } else if( neq == 7 ) {
    for( i=0; i<npmax; i++ ) gzi[i] = 0;
    for( i=0; i<mi; i++ ) gzi[nbi[i]] = gxd[i];
  } else {
    for( i=0; i<npmax; i++ ) {
      gxi[i] = 0;
      gyi[i] = 0;
      gzi[i] = 0;
      vi[i] = 0;
    }
    for( i=0; i<mi; i++ ) {
      gxi[nbi[i]] = gxd[i];
      gyi[nbi[i]] = gyd[i];
      gzi[nbi[i]] = gzd[i];
      vi[nbi[i]] = vd[i];
    }
  }

  unsortvar(0,mi,xi,nbi);
  unsortvar(0,mi,yi,nbi);
  unsortvar(0,mi,zi,nbi);
}
