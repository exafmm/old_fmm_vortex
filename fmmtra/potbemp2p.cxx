#include "../misc/constants.h"

extern float *xi,*yi,*zi,*xj,*yj,*zj,*gxj,*gyj,*gzj,*vj;
extern int **ndi,**ndj,*nej,*nij,**neij;
extern float *vd;

void potbemp2p(int lbi, int lbj) {
  int ii,ij,jj,i,j;
  double pij,dxij,dyij,dzij,rij,rsij;

  for( ii=0; ii<lbi; ii++ ) {
    for( ij=0; ij<nij[ii]; ij++ ) {
      jj = neij[ij][ii];
      for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
        pij = 0;
        for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
          #include "../direct/rij.cxx"
          rsij = 0.25/pi*vj[j]/sqrt(rij)/rij;
          pij += (dxij*gxj[j]+dyij*gyj[j]+dzij*gzj[j])*rsij;
        }
        vd[i] += pij;
      } 
    }
  }
}
