#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern int **ndi,*nfi,*nc;
extern double (*px)[nspm],(*py)[nspm],(*pz)[nspm],*xsp,*ysp,*zsp,**p2g;
extern float *gxd,*gyd,*gzd;
extern int *istagp,*iendgp,*nvecd;
extern float *xig,*yig,*zig,*gxig,*gyig,*gzig,*vig;
extern float *xjg,*yjg,*zjg,*gxjg,*gyjg,*gzjg,*vjg,*sjg;
extern float *gxdg,*gydg,*gzdg,*vdg;

extern void boxc(int, int, int*);
extern "C" void p2pgpu_(int*, double*, double*, double*, double*, double*, double*, double*,
        float*, float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*);

void stl2p(int nmp, int mp, int lbi, double rb) {
  int idev,mblok,ni,nj,ncall,ii,icall,iblok,jc,jbase,j,i,jsize,ibase,isize,is;
  double rsp,op,xic,yic,zic,ggxd,ggyd,ggzd,visd,epsd,dxd,dyd,dzd,dxyzd;

  idev = myrank%ngpu;
  mblok = 3;

  rsp = rb*sqrt(3.0)*0.5;

  ni = 0;
  nj = 0;
  ncall = 0;
  istagp[0] = 0;
  for( ii=0; ii<lbi; ii++ ) {
    ni += ((ndi[1][ii]-ndi[0][ii]+nblok0)/nblok0+1)*nblok0;
    nj += nmp;
    if( ni > nimax || nj > njmax ) {
      iendgp[ncall] = ii-1;
      ncall++;
      istagp[ncall] = ii;
      ni = ((ndi[1][ii]-ndi[0][ii]+nblok0)/nblok0+1)*nblok0;
      nj = nmp;
    }
  }
  iendgp[ncall] = lbi-1;
  ncall++;

  for( icall=0; icall<ncall; icall++ ) {
    iblok = 0;
    jc = 0;
    op = 0;
    for( ii=istagp[icall]; ii<=iendgp[icall]; ii++ ) {
      boxc(nfi[ii],3,nc);
      xic = xmin+(nc[0]+0.5)*rb;
      yic = ymin+(nc[1]+0.5)*rb;
      zic = zmin+(nc[2]+0.5)*rb;
      jbase = jc;
      for( j=0; j<nmp; j++ ) {
        ggxd = 0;
        ggyd = 0;
        ggzd = 0;
        for( i=0; i<nmp; i++ ) {
          ggxd += p2g[j][i]*px[ii][i]*rsp;
          ggyd += p2g[j][i]*py[ii][i]*rsp;
          ggzd += p2g[j][i]*pz[ii][i]*rsp;
        }
        xjg[jc] = xic+xsp[j]*rsp*2;
        yjg[jc] = yic+ysp[j]*rsp*2;
        zjg[jc] = zic+zsp[j]*rsp*2;
        gxjg[jc] = ggxd;
        gyjg[jc] = ggyd;
        gzjg[jc] = ggzd;
        vjg[jc] = 0;
        sjg[jc] = eps;
        jc++;
      }
      jsize = jc-jbase;
      ibase = ndi[0][ii];
      isize = ndi[1][ii]-ibase+1;
      for( is=0; is<isize; is+=nblok0 ) {
        for( i=0; i<std::min(isize-is,nblok0); i++ ) {
          xig[iblok*nblok0+i] = xi[ibase+is+i];
          yig[iblok*nblok0+i] = yi[ibase+is+i];
          zig[iblok*nblok0+i] = zi[ibase+is+i];
          gxig[iblok*nblok0+i] = gxi[ibase+is+i];
          gyig[iblok*nblok0+i] = gyi[ibase+is+i];
          gzig[iblok*nblok0+i] = gzi[ibase+is+i];
          vig[iblok*nblok0+i] = vi[ibase+is+i];
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
    nvecd[4] = 21;
    nvecd[5] = myrank;
    visd = 0;
    epsd = eps;
    dxd = dx;
    dyd = dy;
    dzd = dz;
    dxyzd = dxyz;
    p2pgpu_(nvecd,&op,&visd,&epsd,&dxd,&dyd,&dzd,&dxyzd,xig,yig,zig,gxig,gyig,gzig,vig,
            xjg,yjg,zjg,gxjg,gyjg,gzjg,vjg,sjg,gxdg,gydg,gzdg,vdg);
    iblok = 0;
    for( ii=istagp[icall]; ii<=iendgp[icall]; ii++ ) {
      ibase = ndi[0][ii];
      isize = ndi[1][ii]-ibase+1;
      for( is=0; is<isize; is+=nblok0 ) {
        for( i=0; i<std::min(isize-is,nblok0); i++ ) {
          gxd[ibase+is+i] += gxdg[iblok*nblok0+i];
          gyd[ibase+is+i] += gydg[iblok*nblok0+i];
          gzd[ibase+is+i] += gzdg[iblok*nblok0+i];
        }
        iblok++;
      }
    }
  }
}
