#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi,*xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern int **ndi,**ndj,*nej,*nij,**neij;
extern float *gxd,*gyd,*gzd;

extern double ierfc(double);

void wallp2p(int lbi, int lbj, double vist) {
  int ii,ij,jj,i,j;
  double gxij,gyij,gzij,dxij,dyij,dzij,rij,sij,cutoff;

  for( ii=0; ii<lbi; ii++ ) {
    for( ij=0; ij<nij[ii]; ij++ ) {
      jj = neij[ij][ii];
      for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
        dz = vi[i]/dx/dy;
        gxij = 0;
        gyij = 0;
        gzij = 0;
        for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
          dxij = xi[i]-xj[j];
          dyij = yi[i]-yj[j];
          dzij = std::abs(zi[i]-zj[j]);
          rij = dxij*dxij+dyij*dyij+dzij*dzij+eps;
          sij = sj[j]*sj[j];
          if( rij < 50*sij ) {
            cutoff = (erfc((dzij-dz/2)*vist)-erfc((dxij+dz/2)*vist))*
                     (0.5/vist*(ierfc((dxij+dx)*vist)+ierfc((dxij-dx)*vist)-2*ierfc(dxij*vist)))*
                     (0.5/vist*(ierfc((dyij+dy)*vist)+ierfc((dyij-dy)*vist)-2*ierfc(dyij*vist)));
            gxij += gxj[j]*cutoff;
            gyij += gyj[j]*cutoff;
            gzij += gzj[j]*cutoff;
          }
        }
        gxd[i] += gxij;
        gyd[i] += gyij;
        gzd[i] += gzij;
      } 
    }
  }
}
