#include "../misc/constants.h"

extern float *xi,*yi,*zi;
extern int **ndi,*nfj,*nc,**npx,**neij,*nij,*njb;
extern std::complex<double> (*bx)[mpsym],(*by)[mpsym],(*bz)[mpsym],*bnm,*bth;
extern float *gxd,*gyd,*gzd;

extern void boxc(int, int, int*);
extern void multipoled(double, double, double, int);

void bsm2p(int nmp, int mp, int lbi, int lbj, int lev, int ipb, double rb) {
  int jj,jb,j,jx,jy,jz,ij,ii,i,n,nm,nms,m;
  double xjc,yjc,zjc,xijc,yijc,zijc,r,th,ph,gxr,gxth,gxph,gyr,gyth,gyph,gzr,gzth,gzph;
  double gxy,gxz,gyx,gyz,gzx,gzy;
  std::complex<double> bbx[mpsym],bby[mpsym],bbz[mpsym],rr,rth,rph,cnm(0.0,1.0),eim;

  for( jj=0; jj<lbj; jj++ ) {
    jb = njb[jj];
    for( j=0; j<nmp; j++ ) {
      bbx[j] = bx[jb][j];
      bby[j] = by[jb][j];
      bbz[j] = bz[jb][j];
    }
    boxc(nfj[jj],3,nc);
    jx = nc[0]+npx[0][jj]*int(pow(2,lev));
    jy = nc[1]+npx[1][jj]*int(pow(2,lev));
    jz = nc[2]+npx[2][jj]*int(pow(2,lev));
    xjc = xmin+(jx-0.5+pow(0.5,ipb))*rb;
    yjc = ymin+(jy-0.5+pow(0.5,ipb))*rb;
    zjc = zmin+(jz-0.5+pow(0.5,ipb))*rb;
    for( ij=0; ij<nij[jj]; ij++ ) {
      ii = neij[ij][jj];
      for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
        xijc = xi[i]-xjc;
        yijc = yi[i]-yjc;
        zijc = zi[i]-zjc;
        r = sqrt(xijc*xijc+yijc*yijc+zijc*zijc)+eps;
        th = acos(zijc/r);
        if( std::abs(xijc)+std::abs(yijc) < eps ) {
          ph = 0;
        } else if( std::abs(xijc) < eps ) {
          ph = yijc/std::abs(yijc)*pi*0.5;
        } else if( xijc > 0 ) {
          ph = atan(yijc/xijc);
        } else {
          ph = atan(yijc/xijc)+pi;
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
          rr = -(n+1)*pow(r,-n-2)*bnm[nm];
          rth = pow(r,-n-1)*bth[nm];
          gxr += real(rr*bbx[nms]);
          gxth += real(rth*bbx[nms]);
          gyr += real(rr*bby[nms]);
          gyth += real(rth*bby[nms]);
          gzr += real(rr*bbz[nms]);
          gzth += real(rth*bbz[nms]);
          for( m=1; m<=n; m++ ) {
            nm = n*n+n+m;
            nms = n*(n+1)/2+m;
            eim = exp(m*ph*cnm);
            rr = -(n+1)*pow(r,-n-2)*bnm[nm]*eim;
            rth = pow(r,-n-1)*bth[nm]*eim;
            rph = m*pow(r,-n-1)*bnm[nm]*eim*cnm;
            gxr += 2*real(rr*bbx[nms]);
            gxth += 2*real(rth*bbx[nms]);
            gxph += 2*real(rph*bbx[nms]);
            gyr += 2*real(rr*bby[nms]);
            gyth += 2*real(rth*bby[nms]);
            gyph += 2*real(rph*bby[nms]);
            gzr += 2*real(rr*bbz[nms]);
            gzth += 2*real(rth*bbz[nms]);
            gzph += 2*real(rph*bbz[nms]);
          }
        }
        gxy = sin(th)*sin(ph)*gxr+cos(th)*sin(ph)/r*gxth+cos(ph)/r/sin(th)*gxph;
        gxz = cos(th)*gxr-sin(th)/r*gxth;
        gyx = sin(th)*cos(ph)*gyr+cos(th)*cos(ph)/r*gyth-sin(ph)/r/sin(th)*gyph;
        gyz = cos(th)*gyr-sin(th)/r*gyth;
        gzx = sin(th)*cos(ph)*gzr+cos(th)*cos(ph)/r*gzth-sin(ph)/r/sin(th)*gzph;
        gzy = sin(th)*sin(ph)*gzr+cos(th)*sin(ph)/r*gzth+cos(ph)/r/sin(th)*gzph;
        gxd[i] += 0.25/pi*(gyz-gzy);
        gyd[i] += 0.25/pi*(gzx-gxz);
        gzd[i] += 0.25/pi*(gxy-gyx);
      }
    }
  }
}
