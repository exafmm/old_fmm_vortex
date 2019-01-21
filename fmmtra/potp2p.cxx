#include "../misc/constants.h"

extern float *xi,*yi,*zi,*xj,*yj,*zj,*gxj,*vj;
extern int **ndi,**ndj,*nej,*nij,**neij;
extern float *gxd,*gyd,*gzd,*vd;

void potp2p(int lbi, int lbj) {
  int ii,ij,jj,i,j;
  double pij,uij,vij,wij,dxij,dyij,dzij,rij,rsij;

  for( ii=0; ii<lbi; ii++ ) {
    for( ij=0; ij<nij[ii]; ij++ ) {
      jj = neij[ij][ii];
      for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
        pij = 0;
        uij = 0;
        vij = 0;
        wij = 0;
        for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
          #include "../direct/rij.cxx"
          rsij = 0.25/pi*gxj[j]/pow(rij,1.5);
          pij += 0.25/pi/sqrt(rij)*gxj[j];
          uij -= dxij*rsij;
          vij -= dyij*rsij;
          wij -= dzij*rsij;
        }
        vd[i] += pij;
        gxd[i] += uij;
        gyd[i] += vij;
        gzd[i] += wij;
      } 
    }
  }
}
