#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;

void dir(int n0, int n1, int n2, int n3) {
  int i,j;
  double uij,vij,wij,dxij,dyij,dzij,rij,sij,cutoff;
  for( i=0; i<npmax; i++ ) {
    gxi[i] = 0;
    gyi[i] = 0;
    gzi[i] = 0;
  }
  for( i=n0; i<=n1; i++ ) {
    uij = 0;
    vij = 0;
    wij = 0;
    for( j=n2; j<=n3; j++ ){
      #include "rij.cxx"
      sij = 2*sj[j]*sj[j];
      cutoff = (rij+2.5*sij)/pow(rij+sij,2.5);
      uij += 0.25/pi*(dyij*gzj[j]-dzij*gyj[j])*cutoff;
      vij += 0.25/pi*(dzij*gxj[j]-dxij*gzj[j])*cutoff;
      wij += 0.25/pi*(dxij*gyj[j]-dyij*gxj[j])*cutoff;
    }
    gxi[i] = uij;
    gyi[i] = vij;
    gzi[i] = wij;
  }
}
