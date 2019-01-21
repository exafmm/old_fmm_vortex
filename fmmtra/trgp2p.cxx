#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern int **ndi,**ndj,*nej,*nij,**neij;
extern float *gxd,*gyd,*gzd;

void stp2p(int lbi, int lbj) {
  int ii,ij,jj,i,j;
  double gxij,gyij,gzij,dxij,dyij,dzij,rij,sij,rsij,cutoff;

  for( ii=0; ii<lbi; ii++ ) {
    for( ij=0; ij<nij[ii]; ij++ ) {
      jj = neij[ij][ii];
      for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
        gxij = 0;
        gyij = 0;
        gzij = 0;
        for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
          #include "../direct/rij.cxx"
          sij = 2*sj[j]*sj[j];
          rsij = rij/sij;
          cutoff = 1.0/rij/sqrt(rij)*(erf(sqrt(rsij))-sqrt(4.0/pi*rsij)*exp(-rsij));
          gxij -= 0.25/pi*(gyi[i]*gzj[j]-gzi[i]*gyj[j])*cutoff;
          gyij -= 0.25/pi*(gzi[i]*gxj[j]-gxi[i]*gzj[j])*cutoff;
          gzij -= 0.25/pi*(gxi[i]*gyj[j]-gyi[i]*gxj[j])*cutoff;
          cutoff = (3*erf(sqrt(rsij))-(2*rsij+3)*sqrt(4.0/pi*rsij)*exp(-rsij))/3.0/pow(rij,2.5)*(gxi[i]*(gyj[j]*dzij-gzj[j]*dyij)+gyi[i]*(gzj[j]*dxij-gxj[j]*dzij)+gzi[i]*(gxj[j]*dyij-gyj[j]*dxij));
          gxij += 0.75/pi*dxij*cutoff;
          gyij += 0.75/pi*dyij*cutoff;
          gzij += 0.75/pi*dzij*cutoff;
        }
        gxd[i] += gxij;
        gyd[i] += gyij;
        gzd[i] += gzij;
      } 
    }
  }
}
