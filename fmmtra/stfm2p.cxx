#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi;
extern int **ndi,*nfj,*nc,**npx,**neij,*nij,*njb;
extern std::complex<double> (*bx)[mpsym],(*by)[mpsym],(*bz)[mpsym],*bnm,*bth;
extern float *gxd,*gyd,*gzd;

extern void boxc(int, int, int*);
extern void multipoled(double, double, double, int);

void stm2p(int nmp, int mp, int lbi, int lbj, int lev, int ipb, double rb) {
  int ii,i,ij,jj,jb,j,jx,jy,jz,n,nm,nms,m;
  double xjc,yjc,zjc,xijc,yijc,zijc,r,th,ph,rn,gxr,gxth,gxph,gyr,gyth,gyph,gzr,gzth,gzph;
  double gxx,gxy,gxz,gyx,gyy,gyz,gzx,gzy,gzz;
  std::complex<double> bbx[mpsym],bby[mpsym],bbz[mpsym],rr,rth,rph,cnm(0.0,1.0),eim;

  for( ii=0; ii<lbi; ii++ ) {
    for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
      for( ij=0; ij<nij[ii]; ij++ ) {
        jj = neij[ij][ii];
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
        rn = 1/r;
        for( n=0; n<mp; n++ ) {
          rn /= r;
          nm = n*n+n;
          nms = n*(n+1)/2;
          rr = -(n+1)*rn*bnm[nm];
          rth = rn*r*bth[nm];
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
            rr = -(n+1)*rn*bnm[nm]*eim;
            rth = rn*r*bth[nm]*eim;
            rph = m*rn*r*bnm[nm]*eim*cnm;
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
}
