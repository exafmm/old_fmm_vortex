#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern int **ndi,**ndj,*nej,*nij,**neij;
extern float *gxd,*gyd,*gzd;

void stp2p(int lbi, int lbj) {
  int ii,ij,jj,i,j;
  double gxij,gyij,gzij,dxij,dyij,dzij,rij,sij,rsij,cutoff,st,tr;

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
          cutoff = (3*erf(sqrt(rsij))-(2*rsij+3)*sqrt(4.0/pi*rsij)*exp(-rsij))/3.0/pow(rij,2.5);
          st = gxi[i]*dxij+gyi[i]*dyij+gzi[i]*dzij;
          tr = gxi[i]*(gyj[j]*dzij-gzj[j]*dyij)+gyi[i]*(gzj[j]*dxij-gxj[j]*dzij)+gzi[i]*(gxj[j]*dyij-gyj[j]*dxij);
          gxij += 0.375/pi*((gyj[j]*dzij-gzj[j]*dyij)*st+dxij*tr)*cutoff;
          gyij += 0.375/pi*((gzj[j]*dxij-gxj[j]*dzij)*st+dyij*tr)*cutoff;
          gzij += 0.375/pi*((gxj[j]*dyij-gyj[j]*dxij)*st+dzij*tr)*cutoff;
        }
        gxd[i] += gxij;
        gyd[i] += gyij;
        gzd[i] += gzij;
      } 
    }
  }
}
