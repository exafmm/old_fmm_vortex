#include "../misc/constants.h"

extern float *xi,*yi,*zi;
extern int **ndi,*nfi,*nc;
extern std::complex<double> (*ax)[mpsym],*bnm,*bth;
extern float *gxd,*gyd,*gzd,*vd;

extern void boxc(int, int, int*);
extern void multipoled(double, double, double, int);

void potl2p(int nmp, int mp, int lbi, double rb) {
  int ii,i,n,nm,nms,m;
  double xic,yic,zic,xiic,yiic,ziic,r,th,ph,gx,gxr,gxth,gxph;
  double gxx,gxy,gxz;
  std::complex<double> aax[mpsym],rs,rr,rth,rph,cnm(0.0,1.0),eim;

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
      gxr = 0;
      gxth = 0;
      gxph = 0;
      multipoled(r,th,ph,mp);
      for( n=0; n<mp; n++ ) {
        nm = n*n+n;
        nms = n*(n+1)/2;
        rs = pow(r,n)*bnm[nm];
        rr = n*pow(r,n-1)*bnm[nm];
        rth = pow(r,n)*bth[nm];
        gx += real(rs*aax[nms]);
        gxr += real(rr*aax[nms]);
        gxth += real(rth*aax[nms]);
        for( m=1; m<=n; m++ ) {
          nm = n*n+n+m;
          nms = n*(n+1)/2+m;
          eim = exp(m*ph*cnm);
          rs = pow(r,n)*bnm[nm]*eim;
          rr = n*pow(r,n-1)*bnm[nm]*eim;
          rth = pow(r,n)*bth[nm]*eim;
          rph = m*pow(r,n)*bnm[nm]*eim*cnm;
          gx += 2*real(rs*aax[nms]);
          gxr += 2*real(rr*aax[nms]);
          gxth += 2*real(rth*aax[nms]);
          gxph += 2*real(rph*aax[nms]);
        }
      }
      gxx = sin(th)*cos(ph)*gxr+cos(th)*cos(ph)/r*gxth-sin(ph)/r/sin(th)*gxph;
      gxy = sin(th)*sin(ph)*gxr+cos(th)*sin(ph)/r*gxth+cos(ph)/r/sin(th)*gxph;
      gxz = cos(th)*gxr-sin(th)/r*gxth;
      vd[i] += 0.25/pi*gx;
      gxd[i] += 0.25/pi*gxx;
      gyd[i] += 0.25/pi*gxy;
      gzd[i] += 0.25/pi*gxz;
    }
  }
}
