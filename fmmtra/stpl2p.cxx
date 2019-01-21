#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi;
extern int **ndi,*nfi,*nc;
extern double (*px)[nspm],(*py)[nspm],(*pz)[nspm],*xsp,*ysp,*zsp,**p2g;
extern float *gxd,*gyd,*gzd;

extern void boxc(int, int, int*);

void stl2p(int nmp, int mp, int lbi, double rb) {
  int ii,i,j;
  double rsp,xic,yic,zic,ggxd,ggyd,ggzd,xiic,yiic,ziic,r,cutoff;
  double ggx[nspm],ggy[nspm],ggz[nspm];

  rsp = rb*sqrt(3.0)*0.5;
  for( ii=0; ii<lbi; ii++ ) {
    boxc(nfi[ii],3,nc);
    xic = xmin+(nc[0]+0.5)*rb;
    yic = ymin+(nc[1]+0.5)*rb;
    zic = zmin+(nc[2]+0.5)*rb;
    for( i=0; i<nmp; i++ ) {
      ggxd = 0;
      ggyd = 0;
      ggzd = 0;
      for( j=0; j<nmp; j++ ) {
        ggxd += p2g[j][i]*px[ii][j]*rsp;
        ggyd += p2g[j][i]*py[ii][j]*rsp;
        ggzd += p2g[j][i]*pz[ii][j]*rsp;
      }
      ggx[i] = ggxd;
      ggy[i] = ggyd;
      ggz[i] = ggzd;
    }
    for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
      for( j=0; j<nmp; j++ ) {
        xiic = xi[i]-xic-xsp[j]*rsp*2;
        yiic = yi[i]-yic-ysp[j]*rsp*2;
        ziic = zi[i]-zic-zsp[j]*rsp*2;
        r = sqrt(xiic*xiic+yiic*yiic+ziic*ziic)+eps;
        gxd[i] += 0.25/pi*(gyi[i]*ggz[j]-gzi[i]*ggy[j])/pow(r,3)+0.75/pi*(gxi[i]*xiic+gyi[i]*yiic+gzi[i]*ziic)*(ggy[j]*ziic-ggz[j]*yiic)/pow(r,5);
        gyd[i] += 0.25/pi*(gzi[i]*ggx[j]-gxi[i]*ggz[j])/pow(r,3)+0.75/pi*(gxi[i]*xiic+gyi[i]*yiic+gzi[i]*ziic)*(ggz[j]*xiic-ggx[j]*ziic)/pow(r,5);
        gzd[i] += 0.25/pi*(gxi[i]*ggy[j]-gyi[i]*ggx[j])/pow(r,3)+0.75/pi*(gxi[i]*xiic+gyi[i]*yiic+gzi[i]*ziic)*(ggx[j]*yiic-ggy[j]*xiic)/pow(r,5);
      }
    }
  }
}
