#include "../misc/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern int *nbj,**nxs,*nfn,*na,*nb;

extern void boxn(int, int, int);
extern void sort(int);
extern void sortvar(int, int, float*, int*);
extern void unsortvar(int, int, float*, int*);

void sortj(int& mj) {
  int kl,i,j;
  double rb;

  kl = int(pow(2,lmax));
  rb = rd/kl;
  for( i=0; i<mj; i++ ) {
    nxs[0][i] = int((xj[i]-xmin)/rb);
    nxs[1][i] = int((yj[i]-ymin)/rb);
    nxs[2][i] = int((zj[i]-zmin)/rb);
    for( j=0; j<3; j++ ) {
      if( nxs[j][i] >= int(pow(2,lmax)) ) nxs[j][i] = nxs[j][i]-1;
    }
  }

  boxn(mj,3,lmax);
  for( i=0; i<mj; i++ ) {
    na[i] = nfn[i];
    nb[i] = i;
  }
  sort(mj);
  for( i=0; i<mj; i++ ) {
    nbj[i] = nb[i];
  }

  sortvar(0,mj,xj,nbj);
  sortvar(0,mj,yj,nbj);
  sortvar(0,mj,zj,nbj);
  sortvar(0,mj,gxj,nbj);
  sortvar(0,mj,gyj,nbj);
  sortvar(0,mj,gzj,nbj);
  sortvar(0,mj,vj,nbj);
  sortvar(0,mj,sj,nbj);
}

void unsortj(int& mj) {
  unsortvar(0,mj,xj,nbj);
  unsortvar(0,mj,yj,nbj);
  unsortvar(0,mj,zj,nbj);
  unsortvar(0,mj,gxj,nbj);
  unsortvar(0,mj,gyj,nbj);
  unsortvar(0,mj,gzj,nbj);
  unsortvar(0,mj,vj,nbj);
  unsortvar(0,mj,sj,nbj);
}
