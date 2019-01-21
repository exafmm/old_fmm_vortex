#include "../misc/constants.h"

extern float *xi,*yi,*zi;
extern int **ndi,*nfj,*nc,**npx,**neij,*nij,*njb;
extern std::complex<double> (*bx)[mpsym],*bnm,*bth;
extern float *gxd,*gyd,*gzd,*vd;

extern void boxc(int, int, int*);
extern void multipoled(double, double, double, int);

void potm2p(int nmp, int mp, int lbi, int lbj, int lev, int ipb, double rb) {
  int ii,i,ij,jj,jb,j,jx,jy,jz,n,nm,nms,m;
  double xjc,yjc,zjc,xijc,yijc,zijc,r,th,ph,rn,gx,gxr,gxth,gxph,gxx,gxy,gxz;
  std::complex<double> bbx[mpsym],rs,rr,rth,rph,cnm(0.0,1.0),eim;

  for( ii=0; ii<lbi; ii++ ) {
    for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
      for( ij=0; ij<nij[ii]; ij++ ) {
        jj = neij[ij][ii];
        jb = njb[jj];
        for( j=0; j<nmp; j++ ) bbx[j] = bx[jb][j];
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
        gx = 0;
        gxr = 0;
        gxth = 0;
        gxph = 0;
        multipoled(r,th,ph,mp);
        rn = 1/r;
        for( n=0; n<mp; n++ ) {
          rn /= r;
          nm = n*n+n;
          nms = n*(n+1)/2;
          rs = rn*r*bnm[nm];
          rr = -(n+1)*rn*bnm[nm];
          rth = rn*r*bth[nm];
          gx += real(rs*bbx[nms]);
          gxr += real(rr*bbx[nms]);
          gxth += real(rth*bbx[nms]);
          for( m=1; m<=n; m++ ) {
            nm = n*n+n+m;
            nms = n*(n+1)/2+m;
            eim = exp(m*ph*cnm);
            rs = rn*r*bnm[nm]*eim;
            rr = -(n+1)*rn*bnm[nm]*eim;
            rth = rn*r*bth[nm]*eim;
            rph = m*rn*r*bnm[nm]*eim*cnm;
            gx += 2*real(rs*bbx[nms]);
            gxr += 2*real(rr*bbx[nms]);
            gxth += 2*real(rth*bbx[nms]);
            gxph += 2*real(rph*bbx[nms]);
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
}
