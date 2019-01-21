#include "../misc/constants.h"

extern float *xi,*yi,*zi;
extern int **ndi,*nfi,*nc;
extern double (*px)[nspm],*xsp,*ysp,*zsp,**p2g;
extern float *gxd,*gyd,*gzd,*vd;

extern void boxc(int, int, int*);

void potl2p(int nmp, int mp, int lbi, double rb) {
  int ii,i,j;
  double rsp,xic,yic,zic,ggxd,ggyd,ggzd,xiic,yiic,ziic,r,rsij;
  double ggx[nspm];

  rsp = rb*sqrt(3.0)*0.5;
  for( ii=0; ii<lbi; ii++ ) {
    boxc(nfi[ii],3,nc);
    xic = xmin+(nc[0]+0.5)*rb;
    yic = ymin+(nc[1]+0.5)*rb;
    zic = zmin+(nc[2]+0.5)*rb;
    for( i=0; i<nmp; i++ ) {
      ggxd = 0;
      for( j=0; j<nmp; j++ ) {
        ggxd += p2g[j][i]*px[ii][j]*rsp;
      }
      ggx[i] = ggxd;
    }
    for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
      for( j=0; j<nmp; j++ ) {
        xiic = xi[i]-xic-xsp[j]*rsp*2;
        yiic = yi[i]-yic-ysp[j]*rsp*2;
        ziic = zi[i]-zic-zsp[j]*rsp*2;
        r = sqrt(xiic*xiic+yiic*yiic+ziic*ziic)+eps;
        rsij = 0.25/pi*ggx[j]/pow(r,1.5);
        vd[i] += 0.25/pi/sqrt(r)*ggx[j];
        gxd[i] -= xiic*rsij;
        gyd[i] -= yiic*rsij;
        gzd[i] -= ziic*rsij;
      }
    }
  }
}
