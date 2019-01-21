#include "../misc/constants.h"

extern float *xi,*yi,*zi;
extern int **ndi,*nfj,*nc,**npx,**neij,*nij,*njb;
extern double (*gx)[nspm],(*gy)[nspm],(*gz)[nspm],*xsp,*ysp,*zsp;
extern float *gxd,*gyd,*gzd;

extern void boxc(int, int, int*);

void bsm2p(int nmp, int mp, int lbi, int lbj, int lev, int ipb, double rb) {
  int ii,ij,jj,jb,j,jx,jy,jz,i;
  double rsp,xjc,yjc,zjc,xijc,yijc,zijc,r;
  double ggx[nspm],ggy[nspm],ggz[nspm];

  rsp = rb*sqrt(3.0)*0.5;
  for( ii=0; ii<lbi; ii++ ) {
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
      xjc = xmin+(jx-0.5+pow(0.5,ipb))*rb;
      yjc = ymin+(jy-0.5+pow(0.5,ipb))*rb;
      zjc = zmin+(jz-0.5+pow(0.5,ipb))*rb;
      for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
        for( j=0; j<nmp; j++ ) {
          xijc = xi[i]-xjc-xsp[j]*rsp;
          yijc = yi[i]-yjc-ysp[j]*rsp;
          zijc = zi[i]-zjc-zsp[j]*rsp;
          r = sqrt(xijc*xijc+yijc*yijc+zijc*zijc)+eps;
          gxd[i] += 0.25/pi*(yijc*ggz[j]-zijc*ggy[j])/pow(r,3);
          gyd[i] += 0.25/pi*(zijc*ggx[j]-xijc*ggz[j])/pow(r,3);
          gzd[i] += 0.25/pi*(xijc*ggy[j]-yijc*ggx[j])/pow(r,3);
        }
      }
    }
  }
}
