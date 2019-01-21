#include "../misc/constants.h"

extern int *nfi,*nfj,*nlbj,*nc,**npx,**neij,*nij,*njb;
extern double (*px)[nspm],(*gx)[nspm],(*py)[nspm],(*gy)[nspm],(*pz)[nspm],(*gz)[nspm],*xsp,*ysp,*zsp;
extern int *jbase,*jsize,*istagp,*iendgp,*nicall,*nvecd,*njj;
extern float *xig,*yig,*zig,*gxig,*gyig,*gzig,*vig;
extern float *xjg,*yjg,*zjg,*gxjg,*gyjg,*gzjg,*vjg,*sjg;
extern float *gxdg,*gydg,*gzdg,*vdg;

extern void boxc(int, int, int*);
extern "C" void p2pgpu_(int*, double*, double*, double*, double*, double*, double*, double*,
        float*, float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*);

void m2l(int nmp, int mp, int lbi, int lbj, int lev, int ini, double rb) {
  int idev,mblok,i,j,ni,nj,ncall,jj,ii,ij,icall,iblok,jc,jjd,jb,isize,is,jjdd;
  double rsp,op,xjc,yjc,zjc,xic,yic,zic,visd,epsd,dxd,dyd,dzd,dxyzd;

  idev = myrank%ngpu;
  mblok = 2*nebm+1;

  rsp = rb*sqrt(3.0)*0.5;
  if( ini != 0 ) {
    for( i=0; i<ini; i++ ) {
      for( j=0; j<nmp; j++ ) {
        px[j][i] = 0;
        py[j][i] = 0;
        pz[j][i] = 0;
      }
    }
  }

  ni = 0;
  nj = 0;
  ncall = 0;
  istagp[0] = 0;
  for( jj=0; jj<nbnp; jj++ ) njj[jj] = 0;
  for( ii=0; ii<lbi; ii++ ) {
    if( nij[ii] != 0 ) {
      ni += nmp;
      for( ij=0; ij<nij[ii]; ij++ ) {
        jj = neij[ij][ii];
        if( njj[jj] == 0 ) {
          nj += nmp;
          njj[jj] = 1;
        }
      }
      if ( ni > nimax || nj > njmax ) {
        iendgp[ncall] = ii-1;
        ncall++;
        istagp[ncall] = ii;
        ni = nmp;
        nj = 0;
        for( jj=0; jj<nbnp; jj++ ) njj[jj] = 0;
        for( ij=0; ij<nij[ii]; ij++ ) {
          jj = neij[ij][ii];
          nj += nmp;
          njj[jj] = 1;
        }
      }
    }
  }
  iendgp[ncall] = lbi-1;
  ncall++;

  for( icall=0; icall<ncall; icall++ ) {
    iblok = 0;
    jc = 0;
    jjd = 0;
    op = 0;
    for( jj=0; jj<nbnp; jj++ ) njj[jj] = 0;
    for( ii=istagp[icall]; ii<=iendgp[icall]; ii++ ) {
      for( ij=0; ij<nij[ii]; ij++ ) {
        jj = neij[ij][ii];
        jb = njb[jj];
        if( njj[jj] == 0 ) {
          boxc(nfj[jj],3,nc);
          xjc = (nc[0]+npx[0][jj]*pow(2,lev)+0.5)*rb;
          yjc = (nc[1]+npx[1][jj]*pow(2,lev)+0.5)*rb;
          zjc = (nc[2]+npx[2][jj]*pow(2,lev)+0.5)*rb;
          jbase[jjd] = jc;
          for( j=0; j<nmp; j++ ) {
            xjg[jc] = xjc+xsp[j]*rsp;
            yjg[jc] = yjc+ysp[j]*rsp;
            zjg[jc] = zjc+zsp[j]*rsp;
            gxjg[jc] = gx[jb][j];
            gyjg[jc] = gy[jb][j];
            gzjg[jc] = gz[jb][j];
            vjg[jc] = 0;
            sjg[jc] = 1;
            jc++;
          }
          jsize[jjd] = jc-jbase[jjd];
          njj[jj] = jjd;
          jjd++;
        }
      }
      boxc(nfi[ii],3,nc);
      xic = (nc[0]+0.5)*rb;
      yic = (nc[1]+0.5)*rb;
      zic = (nc[2]+0.5)*rb;
      isize = nmp;
      for( is=0; is<isize; is+=nblok1 ) {
        for( i=0; i<std::min(isize-is,nblok0); i++ ) {
          xig[iblok*nblok0+i] = xic+xsp[is+i]*rsp;
          yig[iblok*nblok0+i] = yic+ysp[is+i]*rsp;
          zig[iblok*nblok0+i] = zic+zsp[is+i]*rsp;
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
        nvecd[iblok*mblok+10] = nij[ii];
        for( ij=0; ij<nij[ii]; ij++ ) {
          jj = neij[ij][ii];
          jjdd = njj[jj];
          nvecd[iblok*mblok+2*ij+11] = jbase[jjdd];
          nvecd[iblok*mblok+2*ij+12] = jsize[jjdd];
          op += (double) nblok0*jsize[jjdd];
        }
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
    iblok = 0;
    for( ii=istagp[icall]; ii<=iendgp[icall]; ii++ ) {
      isize = nmp;
      for( is=0; is<isize; is+=nblok1 ) {
        for( i=0; i<std::min(isize-is,nblok1); i++ ) {
          px[ii][is+i] += gxdg[iblok*nblok0+i];
          py[ii][is+i] += gydg[iblok*nblok0+i];
          pz[ii][is+i] += gzdg[iblok*nblok0+i];
        }
        iblok++;
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
