#include "../misc/constants.h"

extern float *xi,*yi,*zi,*vi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj;

void dir(int n0, int n1, int n2, int n3) {
  int i,j;
  double pij,dxij,dyij,dzij,rij,rsij;
  for( i=0; i<npmax; i++ ) {
    vi[i] = 0;
  }
  for( i=n0; i<=n1; i++ ) {
    pij = 0;
    for( j=n2; j<=n3; j++ ){
      #include "rij.cxx"
      rsij = 0.25/pi*vj[j]/rij/sqrt(rij);
      pij += (dxij*gxj[j]+dyij*gyj[j]+dzij*gzj[j])*rsij;
    }
    vi[i] = pij;
  }
}
