#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern float *gxd,*gyd,*gzd;

extern void memoryuse();
extern void memoryfree();

void dir(int n0, int n1, int n2, int n3) {
  int i,j;
  double gxdd,gydd,gzdd,xij,yij,zij,rij,w;

  gxd = new float [npmax];
  gyd = new float [npmax];
  gzd = new float [npmax];
  mem = npmax*3*4;
  memoryuse();

  for( i=n0; i<=n1; i++ ) {
    dz = vi[i]/dx/dy;
    gxdd = 0;
    gydd = 0;
    gzdd = 0;
    for( j=n2; j<=n3; j++ ){
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
        gxdd += gxj[j]*w;
        gydd += gyj[j]*w;
        gzdd += gzj[j]*w;
      }
    }
    gxd[i] = gxdd;
    gyd[i] = gydd;
    gzd[i] = gzdd;
  }
  for( i=0; i<npmax; i++ ) {
    gxi[i] = 0;
    gyi[i] = 0;
    gzi[i] = 0;
  }
  for( i=n0; i<=n1; i++ ) {
    gxi[i] = gxd[i];
    gyi[i] = gyd[i];
    gzi[i] = gzd[i];
  }

  delete[] gxd;
  delete[] gyd;
  delete[] gzd;
  mem = npmax*3*4;
  memoryfree();
}
