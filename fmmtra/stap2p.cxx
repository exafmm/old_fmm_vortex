#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern int **ndi,**ndj,*nej,*nij,**neij;
extern float *gxd,*gyd,*gzd;

void stp2p(int lbi, int lbj) {
  int ii,ij,jj,i,j;
  double gxij,gyij,gzij,dxij,dyij,dzij,rij,sij,cutoff;

  for( ii=0; ii<lbi; ii++ ) {
    for( ij=0; ij<nij[ii]; ij++ ) {
      jj = neij[ij][ii];
      for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
        gxij = 0;
        gyij = 0;
        gzij = 0;
        for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
          #include "../direct/rij.cxx"
          sij = sj[j]*sj[j];
          cutoff = (rij+2.5*sij)/pow(rij+sij,2.5);
          gxij += 0.25/pi*(gyi[i]*gzj[j]-gzi[i]*gyj[j])*cutoff;
          gyij += 0.25/pi*(gzi[i]*gxj[j]-gxi[i]*gzj[j])*cutoff;
          gzij += 0.25/pi*(gxi[i]*gyj[j]-gyi[i]*gxj[j])*cutoff;
          cutoff = (rij+3.5*sij)/pow(rij+sij,3.5)*(gxi[i]*dxij+gyi[i]*dyij+gzi[i]*dzij);
          gxij += 0.75/pi*(gyj[j]*dzij-gzj[j]*dyij)*cutoff;
          gyij += 0.75/pi*(gzj[j]*dxij-gxj[j]*dzij)*cutoff;
          gzij += 0.75/pi*(gxj[j]*dyij-gyj[j]*dxij)*cutoff;
        }
        gxd[i] += gxij;
        gyd[i] += gyij;
        gzd[i] += gzij;
      } 
    }
  }
}
