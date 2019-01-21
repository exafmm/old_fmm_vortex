#include "../misc/constants.h"

extern int *nfi,*nfj,*nlbj,*nc,**npx,**neij,*nij,*njb;
extern double (*px)[nspm],(*gx)[nspm],(*py)[nspm],(*gy)[nspm],(*pz)[nspm],(*gz)[nspm],*xsp,*ysp,*zsp;

extern void boxc(int, int, int*);

void m2l(int nmp, int mp, int lbi, int lbj, int lev, int ini, double rb) {
  int i,j,ii,ix,iy,iz,ij,jj,jb,jx,jy,jz,k;
  double rsp,xic,yic,zic,xijc,yijc,zijc,ppxd,ppyd,ppzd,r;
  double ggx[nspm],ggy[nspm],ggz[nspm];

  rsp = rb*sqrt(3.0)*0.5;
  if( ini != 0 ) {
    for( i=0; i<ini; i++ ) {
      for( j=0; j<nmp; j++ ) {
        px[i][j] = 0;
        py[i][j] = 0;
        pz[i][j] = 0;
      }
    }
  }
  for( ii=0; ii<lbi; ii++ ) {
    boxc(nfi[ii],3,nc);
    ix = nc[0];
    iy = nc[1];
    iz = nc[2];
    for( ij=0; ij<nij[ii]; ij++ ) {
      jj = neij[ij][ii];
      jb = njb[jj];
      for( j=0; j<nmp; j++ ) {
        ggx[j] = gx[jb][j];
        ggy[j] = gy[jb][j];
        ggz[j] = gz[jb][j];
      }
      boxc(nfj[jj],3,nc);
      jx = nc[0]+npx[0][jj]*int(pow(2,lev));
      jy = nc[1]+npx[1][jj]*int(pow(2,lev));
      jz = nc[2]+npx[2][jj]*int(pow(2,lev));
      xic = (ix-jx)*rb;
      yic = (iy-jy)*rb;
      zic = (iz-jz)*rb;
      for( j=0; j<nmp; j++ ) {
        ppxd = 0;
        ppyd = 0;
        ppzd = 0;
        for( k=0; k<nmp; k++ ) {
          xijc = xic+xsp[j]*rsp-xsp[k]*rsp;
          yijc = yic+ysp[j]*rsp-ysp[k]*rsp;
          zijc = zic+zsp[j]*rsp-zsp[k]*rsp;
          r = sqrt(xijc*xijc+yijc*yijc+zijc*zijc)+eps;
          ppxd += 0.25/pi*ggx[k]/r;
          ppyd += 0.25/pi*ggy[k]/r;
          ppzd += 0.25/pi*ggz[k]/r;
        }
        px[ii][j] += ppxd;
        py[ii][j] += ppyd;
        pz[ii][j] += ppzd;
      }
    }
  }
  for( jj=0; jj<lbj; jj++ ) {
    jb = njb[jj];
    for( j=0; j<nmp; j++ ) {
      gx[jb][j] = 0;
      gy[jb][j] = 0;
      gz[jb][j] = 0;
    }
  }
}
