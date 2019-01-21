#include "../misc/constants.h"

extern float *xi,*yi,*zi;
extern int **ndi,*nfi,*nc;
extern std::complex<double> (*ax)[mpsym],*bnm,*bth;
extern float *vd;

extern void boxc(int, int, int*);
extern void multipoled(double, double, double, int);

void potbeml2p(int nmp, int mp, int lbi, double rb) {
  int ii,i,n,nm,nms,m;
  double xic,yic,zic,xiic,yiic,ziic,r,th,ph,gx;
  std::complex<double> aax[mpsym],rs,cnm(0.0,1.0),eim;

  for( ii=0; ii<lbi; ii++ ) {
    boxc(nfi[ii],3,nc);
    xic = xmin+(nc[0]+0.5)*rb;
    yic = ymin+(nc[1]+0.5)*rb;
    zic = zmin+(nc[2]+0.5)*rb;
    for( i=0; i<nmp; i++ ) aax[i] = ax[ii][i];
    for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
      xiic = xi[i]-xic;
      yiic = yi[i]-yic;
      ziic = zi[i]-zic;
      r = sqrt(xiic*xiic+yiic*yiic+ziic*ziic)+eps;
      th = acos(ziic/r);
      if( std::abs(xiic)+std::abs(yiic) < eps ) {
        ph = 0;
      } else if( std::abs(xiic) < eps ) {
        ph = yiic/std::abs(yiic)*pi*0.5;
      } else if( xiic > 0 ) {
        ph = atan(yiic/xiic);
      } else {
        ph = atan(yiic/xiic)+pi;
      }
      gx = 0;
      multipoled(r,th,ph,mp);
      for( n=0; n<mp; n++ ) {
        nm = n*n+n;
        nms = n*(n+1)/2;
        rs = pow(r,n)*bnm[nm];
        gx += real(rs*aax[nms]);
        for( m=1; m<=n; m++ ) {
          nm = n*n+n+m;
          nms = n*(n+1)/2+m;
          eim = exp(m*ph*cnm);
          rs = pow(r,n)*bnm[nm]*eim;
          gx += 2*real(rs*aax[nms]);
        }
      }
      vd[i] += 0.25/pi*(gx);
    }
  }
}
