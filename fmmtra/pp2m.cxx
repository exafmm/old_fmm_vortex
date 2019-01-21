#include "../misc/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj;
extern int **ndj,*nfj,*nc;
extern double (*gx)[nspm],(*gy)[nspm],(*gz)[nspm],*xsp,*ysp,*zsp,**p2g;

extern void boxc(int, int, int*);

void p2m(int nmp, int mp, int lbj, double rb, int nlm) {
  int jj,j,k;
  double rsp,xjc,yjc,zjc,xjjc,yjjc,zjjc,rh,ggxd,ggyd,ggzd;
  double ppx[nspm],ppy[nspm],ppz[nspm];

  rsp = rb*sqrt(3.0)*0.5;
  for( jj=0; jj<lbj; jj++ ) {
    boxc(nfj[jj],3,nc);
    xjc = xmin+(nc[0]+0.5)*rb;
    yjc = ymin+(nc[1]+0.5)*rb;
    zjc = zmin+(nc[2]+0.5)*rb;
    for( j=0; j<nmp; j++ ) {
      ppx[j] = 0;
      ppy[j] = 0;
      ppz[j] = 0;
    }
    for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
      for( k=0; k<nmp; k++ ) {
        xjjc = xjc+2*xsp[k]*rsp-xj[j];
        yjjc = yjc+2*ysp[k]*rsp-yj[j];
        zjjc = zjc+2*zsp[k]*rsp-zj[j];
        rh = sqrt(xjjc*xjjc+yjjc*yjjc+zjjc*zjjc)+eps;
        ppx[k] += 0.25/pi*gxj[j]/rh;
        ppy[k] += 0.25/pi*gyj[j]/rh;
        ppz[k] += 0.25/pi*gzj[j]/rh;
      }
    }
    for( j=0; j<nmp; j++ ) {
      ggxd = 0;
      ggyd = 0;
      ggzd = 0;
      for( k=0; k<nmp; k++ ) {
        ggxd += p2g[k][j]*ppx[k]*rsp;
        ggyd += p2g[k][j]*ppy[k]*rsp;
        ggzd += p2g[k][j]*ppz[k]*rsp;
      }
      gx[jj+nlm][j] = ggxd;
      gy[jj+nlm][j] = ggyd;
      gz[jj+nlm][j] = ggzd;
    }
  }
}
