#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi;
extern int **ndi,*nfi,*nc;
extern std::complex<double> (*ax)[mpsym],(*ay)[mpsym],(*az)[mpsym],*bnm,*bth;
extern float *gxd,*gyd,*gzd;

extern void boxc(int, int, int*);
extern void multipoled(double, double, double, int);

void stl2p(int nmp, int mp, int lbi, double rb) {
  int ii,i,n,nm,nms,m;
  double xic,yic,zic,xiic,yiic,ziic,r,th,ph,gxr,gxth,gxph,gyr,gyth,gyph,gzr,gzth,gzph;
  double gxx,gxy,gxz,gyx,gyy,gyz,gzx,gzy,gzz;
  std::complex<double> aax[mpsym],aay[mpsym],aaz[mpsym],rr,rth,rph,cnm(0.0,1.0),eim;

  for( ii=0; ii<lbi; ii++ ) {
    boxc(nfi[ii],3,nc);
    xic = xmin+(nc[0]+0.5)*rb;
    yic = ymin+(nc[1]+0.5)*rb;
    zic = zmin+(nc[2]+0.5)*rb;
    for( i=0; i<nmp; i++ ) {
      aax[i] = ax[ii][i];
      aay[i] = ay[ii][i];
      aaz[i] = az[ii][i];
    }
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
      gxr = 0;
      gxth = 0;
      gxph = 0;
      gyr = 0;
      gyth = 0;
      gyph = 0;
      gzr = 0;
      gzth = 0;
      gzph = 0;
      multipoled(r,th,ph,mp);
      for( n=0; n<mp; n++ ) {
        nm = n*n+n;
        nms = n*(n+1)/2;
        rr = n*pow(r,n-1)*bnm[nm];
        rth = pow(r,n)*bth[nm];
        gxr += real(rr*aax[nms]);
        gxth += real(rth*aax[nms]);
        gyr += real(rr*aay[nms]);
        gyth += real(rth*aay[nms]);
        gzr += real(rr*aaz[nms]);
        gzth += real(rth*aaz[nms]);
        for( m=1; m<=n; m++ ) {
          nm = n*n+n+m;
          nms = n*(n+1)/2+m;
          eim = exp(m*ph*cnm);
          rr = n*pow(r,n-1)*bnm[nm]*eim;
          rth = pow(r,n)*bth[nm]*eim;
          rph = m*pow(r,n)*bnm[nm]*eim*cnm;
          gxr += 2*real(rr*aax[nms]);
          gxth += 2*real(rth*aax[nms]);
          gxph += 2*real(rph*aax[nms]);
          gyr += 2*real(rr*aay[nms]);
          gyth += 2*real(rth*aay[nms]);
          gyph += 2*real(rph*aay[nms]);
          gzr += 2*real(rr*aaz[nms]);
          gzth += 2*real(rth*aaz[nms]);
          gzph += 2*real(rph*aaz[nms]);
        }
      }
      gxx = sin(th)*cos(ph)*gxr+cos(th)*cos(ph)/r*gxth-sin(ph)/r/sin(th)*gxph;
      gyx = sin(th)*cos(ph)*gyr+cos(th)*cos(ph)/r*gyth-sin(ph)/r/sin(th)*gyph;
      gzx = sin(th)*cos(ph)*gzr+cos(th)*cos(ph)/r*gzth-sin(ph)/r/sin(th)*gzph;
      gxy = sin(th)*sin(ph)*gxr+cos(th)*sin(ph)/r*gxth+cos(ph)/r/sin(th)*gxph;
      gyy = sin(th)*sin(ph)*gyr+cos(th)*sin(ph)/r*gyth+cos(ph)/r/sin(th)*gyph;
      gzy = sin(th)*sin(ph)*gzr+cos(th)*sin(ph)/r*gzth+cos(ph)/r/sin(th)*gzph;
      gxz = cos(th)*gxr-sin(th)/r*gxth;
      gyz = cos(th)*gyr-sin(th)/r*gyth;
      gzz = cos(th)*gzr-sin(th)/r*gzth;
      gxd[i] -= 0.25/pi*(gxi[i]*gxx+gyi[i]*gxy+gzi[i]*gxz);
      gyd[i] -= 0.25/pi*(gxi[i]*gyx+gyi[i]*gyy+gzi[i]*gyz);
      gzd[i] -= 0.25/pi*(gxi[i]*gzx+gyi[i]*gzy+gzi[i]*gzz);
    }
  }
}
