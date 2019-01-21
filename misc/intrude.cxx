#include "../misc/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj,*xo,*yo,*zo;
extern double **ptr;
extern int **ntr,*nne;

void intrude(int& np, int& npo, int nge, int nin) {
  int nln,i,n,nz,j,ny,k;
  int l1[4];
  double x1,y1,z1,x2,y2,z2,x3,y3,z3,a,b,c,d;

  nln = 3;
  for( i=0; i<nln-1; i++ ) {
    l1[i] = i+1;
  }
  l1[nln-1] = 0;
  n = 0;
  for( i=0; i<np; i++ ) {
    nz = 0;
    for( j=0; j<nwe; j++ ) {
      ny = 0;
      for( k=0; k<nln; k++ ) {
        x1 = ptr[0][ntr[k][j]];
        y1 = ptr[1][ntr[k][j]];
        x2 = ptr[0][ntr[l1[k]][j]];
        y2 = ptr[1][ntr[l1[k]][j]];
        if( ((x1 <= xj[i] && xj[i] <= x2) || (x2 <= xj[i] && xj[i] <= x1)) &&
            yj[i] < (y1-y2)*(xj[i]-x2)/(x1-x2)+y2 ) {
          ny++;
        }
      }
      if( ny%2 == 1 ) {
        x1 = ptr[0][ntr[0][j]];
        y1 = ptr[1][ntr[0][j]];
        z1 = ptr[2][ntr[0][j]];
        x2 = ptr[0][ntr[1][j]];
        y2 = ptr[1][ntr[1][j]];
        z2 = ptr[2][ntr[1][j]];
        x3 = ptr[0][ntr[2][j]];
        y3 = ptr[1][ntr[2][j]];
        z3 = ptr[2][ntr[2][j]];
        a = (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1);
        b = (z2-z1)*(x3-x1)-(x2-x1)*(z3-z1);
        c = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);
        d = a*x1+b*y1+c*z1;
        if( std::abs(c) < eps ) {
          if( zj[i] < z1 ) nz++;
        } else if( zj[i] < (d-b*yj[i]-a*xj[i])/c ) {
          nz++;
        }
      }
    }
    if( nz%2 == 1 ) {
      if (nin == 1) {
        xj[i] = xo[i];
        yj[i] = yo[i];
        zj[i] = zo[i];
      }
    } else if( nin == 2 ) {
      xj[n] = xj[i];
      yj[n] = yj[i];
      zj[n] = zj[i];
      gxj[n] = gxj[i];
      gyj[n] = gyj[i];
      gzj[n] = gzj[i];
      vj[n] = vj[i];
      sj[n] = sj[i];
      n++;
    }
  }
  npo = np;
  if( nin == 2 ) np = n;
}
