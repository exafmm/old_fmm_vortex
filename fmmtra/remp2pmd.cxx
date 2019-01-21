#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern int **ndi,*nei,**ndj,*nij,**neij;
extern float *gxd,*gyd,*gzd;

void remp2p(int lbi, int lbj) {
  int jj,ij,ii,i,j;
  double gxij,gyij,gzij,xij,yij,zij,rij,w;

  for( jj=0; jj<lbj; jj++ ) {
    for( ij=0; ij<nij[jj]; ij++ ) {
      ii = neij[ij][jj];
      for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
        gxij = 0;
        gyij = 0;
        gzij = 0;
        for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
          xij = std::abs(xi[i]-xj[j])/dx;
          yij = std::abs(yi[i]-yj[j])/dy;
          zij = std::abs(zi[i]-zj[j])/dz;
          rij = xij*xij+yij*yij+zij*zij+eps;
          if( rij < 12 ) {
            if( xij > 2 ) {
              w = 0;
            } else if( xij >= 1 ) {
              w = 0.5*(2-xij)*(2-xij)*(1-xij);
            } else {
              w = 1-2.5*xij*xij+1.5*xij*xij*xij;
            }
            if( yij > 2 ) {
              w = 0;
            } else if( yij >= 1 ) {
              w = 0.5*(2-yij)*(2-yij)*(1-yij)*w;
            } else {
              w = (1-2.5*yij*yij+1.5*yij*yij*yij)*w;
            }
            if( zij > 2 ) {
              w = 0;
            } else if( zij >= 1 ) {
              w = 0.5*(2-zij)*(2-zij)*(1-zij)*w;
            } else {
              w = (1-2.5*zij*zij+1.5*zij*zij*zij)*w;
            }
            gxij += gxj[j]*w;
            gyij += gyj[j]*w;
            gzij += gzj[j]*w;
          }
        }
        gxd[i] += gxij;
        gyd[i] += gyij;
        gzd[i] += gzij;
      } 
    }
  }
}
