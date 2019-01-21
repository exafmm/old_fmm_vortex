#include "../misc/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern double (*gx)[nspm],(*gy)[nspm],(*gz)[nspm],*xsp,*ysp,*zsp,**p2g;
extern int **ndj,*nfj,*nc;
extern int *istagp,*iendgp,*nvecd;
extern float *xig,*yig,*zig,*gxig,*gyig,*gzig,*vig;
extern float *xjg,*yjg,*zjg,*gxjg,*gyjg,*gzjg,*vjg,*sjg;
extern float *gxdg,*gydg,*gzdg,*vdg;

extern void boxc(int, int, int*);
extern "C" void p2pgpu_(int*, double*, double*, double*, double*, double*, double*, double*,
        float*, float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*);

void p2m(int nmp, int mp, int lbj, double rb, int nlm) {
  int idev,mblok,nj,ncall,jj,icall,iblok,jc,jbase,j,jsize,isize,is,i,iiblok;
  double rsp,op,xjc,yjc,zjc,visd,epsd,dxd,dyd,dzd,dxyzd,gxdd,gydd,gzdd;

  idev = myrank%ngpu;
  mblok = 3;

  rsp = rb*sqrt(3.0)*0.5;

  nj = 0;
  ncall = 0;
  istagp[0] = 0;
  for( jj=0; jj<lbj; jj++ ) {
    nj += ndj[1][jj]-ndj[0][jj]+1;
    if( nj > njmax ) {
      iendgp[ncall] = jj-1;
      ncall++;
      istagp[ncall] = jj;
      nj = ndj[1][jj]-ndj[0][jj]+1;
    }
  }
  iendgp[ncall] = lbj-1;
  ncall++;

  for( icall=0; icall<ncall; icall++ ) {
    iblok = 0;
    jc = 0;
    op = 0;
    for( jj=istagp[icall]; jj<=iendgp[icall]; jj++ ) {
      jbase = jc;
      for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
        xjg[jc] = xj[j];
        yjg[jc] = yj[j];
        zjg[jc] = zj[j];
        gxjg[jc] = gxj[j];
        gyjg[jc] = gyj[j];
        gzjg[jc] = gzj[j];
        vjg[jc] = vj[j];
        sjg[jc] = sj[j];
        jc++;
      }
      jsize = jc-jbase;
      boxc(nfj[jj],3,nc);
      xjc = xmin+(nc[0]+0.5)*rb;
      yjc = ymin+(nc[1]+0.5)*rb;
      zjc = zmin+(nc[2]+0.5)*rb;
      isize = nmp;
      for( is=0; is<isize; is+=nblok0 ) {
        for( i=0; i<std::min(isize-is,nblok0); i++ ) {
          xig[iblok*nblok0+i] = xjc+xsp[is+i]*rsp*2;
          yig[iblok*nblok0+i] = yjc+ysp[is+i]*rsp*2;
          zig[iblok*nblok0+i] = zjc+zsp[is+i]*rsp*2;
          gxig[iblok*nblok0+i] = 0;
          gyig[iblok*nblok0+i] = 0;
          gzig[iblok*nblok0+i] = 0;
          vig[iblok*nblok0+i] = 0;
        }
        for( i=isize-is; i<nblok0; i++ ) {
          xig[iblok*nblok0+i] = 0;
          yig[iblok*nblok0+i] = 0;
          zig[iblok*nblok0+i] = 0;
          gxig[iblok*nblok0+i] = 0;
          gyig[iblok*nblok0+i] = 0;
          gzig[iblok*nblok0+i] = 0;
          vig[iblok*nblok0+i] = 0;
        }
        nvecd[iblok*mblok+10] = 1;
        nvecd[iblok*mblok+11] = jbase;
        nvecd[iblok*mblok+12] = jsize;
        op += nblok0*jsize;
        iblok++;
      }
    }
    nvecd[0] = idev;
    nvecd[1] = iblok;
    nvecd[2] = mblok;
    nvecd[3] = jc;
    nvecd[4] = 24;
    nvecd[5] = myrank;
    visd = 0;
    epsd = eps;
    dxd = dx;
    dyd = dy;
    dzd = dz;
    dxyzd = dxyz;
    p2pgpu_(nvecd,&op,&visd,&epsd,&dxd,&dyd,&dzd,&dxyzd,xig,yig,zig,gxig,gyig,gzig,vig,
            xjg,yjg,zjg,gxjg,gyjg,gzjg,vjg,sjg,gxdg,gydg,gzdg,vdg);
    iiblok = 0;
    for( jj=istagp[icall]; jj<=iendgp[icall]; jj++ ) {
      for( j=0; j<nmp; j++ ) {
        gxdd = 0;
        gydd = 0;
        gzdd = 0;
        isize = nmp;
        iblok = iiblok;
        for( is=0; is<isize; is+=nblok0 ) {
          for( i=0; i<std::min(isize-is,nblok0); i++ ) {
            gxdd += p2g[j][is+i]*gxdg[iblok*nblok0+i]*rsp;
            gydd += p2g[j][is+i]*gydg[iblok*nblok0+i]*rsp;
            gzdd += p2g[j][is+i]*gzdg[iblok*nblok0+i]*rsp;
          }
          iblok++;
        }
        gx[jj+nlm][j] = gxdd;
        gy[jj+nlm][j] = gydd;
        gz[jj+nlm][j] = gzdd;
      }
      iiblok = iblok;
    }
  }
}
