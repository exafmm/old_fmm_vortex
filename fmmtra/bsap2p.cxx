#include "../misc/constants.h"

extern float *xi,*yi,*zi,*vi,*xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern int **ndi,**ndj,*nej,*nij,**neij;
extern float *gxd,*gyd,*gzd;

void bsp2p(int lbi, int lbj) {
  int ii,ij,jj,i,j;
  double uij,vij,wij,dxij,dyij,dzij,rij,sij,cutoff;

  for( ii=0; ii<lbi; ii++ ) {
    for( ij=0; ij<nij[ii]; ij++ ) {
      jj = neij[ij][ii];
      for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
        uij = 0;
        vij = 0;
        wij = 0;
        for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
          #include "../direct/rij.cxx"
          sij = sj[j]*sj[j];
          cutoff = (rij+2.5*sij)/pow(rij+sij,2.5);
          uij += 0.25/pi*(dyij*gzj[j]-dzij*gyj[j])*cutoff;
          vij += 0.25/pi*(dzij*gxj[j]-dxij*gzj[j])*cutoff;
          wij += 0.25/pi*(dxij*gyj[j]-dyij*gxj[j])*cutoff;
        }
        gxd[i] += uij;
        gyd[i] += vij;
        gzd[i] += wij;
      } 
    }
  }
}
