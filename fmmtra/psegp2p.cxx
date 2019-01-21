#include "../misc/parameters.h"
#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi,*xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern int **ndi,**ndj,*nej,*nij,**neij;
extern float *gxd,*gyd,*gzd;

void psep2p(int lbi, int lbj, double vis) {
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
          sij = 2*sj[j]*sj[j];
          if( rij<100*sij ) {
            cutoff = 4.0*vis/sij/pow(pi*sij,1.5)*exp(-rij/sij);
            gxij +=(gxj[j]*vi[i]-gxi[i]*vj[j])*cutoff;
            gyij +=(gyj[j]*vi[i]-gyi[i]*vj[j])*cutoff;
            gzij +=(gzj[j]*vi[i]-gzi[i]*vj[j])*cutoff;
          }
        }
        gxd[i] += gxij;
        gyd[i] += gyij;
        gzd[i] += gzij;
      } 
    }
  }
}
