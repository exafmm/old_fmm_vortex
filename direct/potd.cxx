#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern float *xj,*yj,*zj,*gxj,*vj;

void dird(int n0, int n1, int n2, int n3) {
  int i,j;
  double pij,uij,vij,wij,dxij,dyij,dzij,rij,rsij;
  for( i=0; i<npmax; i++ ) {
    gxi[i] = 0;
    gyi[i] = 0;
    gzi[i] = 0;
    vi[i] = 0;
  }
  for( i=n0; i<=n1; i++ ) {
    pij = 0;
    uij = 0;
    vij = 0;
    wij = 0;
    for( j=n2; j<=n3; j++ ){
      #include "rij.cxx"
      rsij = 0.25/pi*gxj[j]/rij/sqrt(rij);
      pij += 0.25/pi/sqrt(rij)*gxj[j];
      uij -= dxij*rsij;
      vij -= dyij*rsij;
      wij -= dzij*rsij;
    }
    gxi[i] = uij;
    gyi[i] = vij;
    gzi[i] = wij;
    vi[i] = pij;
  }
}
