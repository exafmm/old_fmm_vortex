#include "../misc/constants.h"

extern float *xi,*yi,*zi;
extern int **ndi,*nfi,*nc;
extern std::complex<double> (*ax)[mpsym],(*ay)[mpsym],(*az)[mpsym];
extern float *fac,*gxd,*gyd,*gzd;
extern int *istagp,*iendgp,*nvecd;
extern float *xig,*yig,*zig,*gxig,*gyig,*gzig,*vig;
extern float *arex,*arey,*arez,*aimx,*aimy,*aimz;

extern "C" void l2pgpu_(int*, double*, double*, double*, double*, double*, double*,
        float*, float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*, float*, float*, float*);

void bsl2p(int nmp, int mp, int lbi, double rb) {
  int idev,mblok,ni,nj,ncall,ii,icall,iblok,jc,jbase,j,jsize,ibase,isize,is,i;
  double op,rbd,epsd,xmind,ymind,zmind;

  idev = myrank%ngpu;
  mblok = 3;

  ni = 0;
  nj = 0;
  ncall = 0;
  istagp[0] = 0;
  for( ii=0; ii<lbi; ii++ ) {
    ni += ((ndi[1][ii]-ndi[0][ii]+nblok1)/nblok1+1)*nblok1;
    nj += nmp;
    if( ni > nimax || nj > njmax ) {
      iendgp[ncall] = ii-1;
      ncall++;
      istagp[ncall] = ii;
      ni = ((ndi[1][ii]-ndi[0][ii]+nblok1)/nblok1+1)*nblok1;
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
      jbase = jc;
      for( j=0; j<nmp; j++ ) {
        arex[jc] = std::real(ax[ii][j]);
        arey[jc] = std::real(ay[ii][j]);
        arez[jc] = std::real(az[ii][j]);
        aimx[jc] = std::imag(ax[ii][j]);
        aimy[jc] = std::imag(ay[ii][j]);
        aimz[jc] = std::imag(az[ii][j]);
        jc++;
      }
      jsize = jc-jbase;
      ibase = ndi[0][ii];
      isize = ndi[1][ii]-ibase+1;
      for( is=0; is<isize; is+=nblok1 ) {
        for( i=0; i<std::min(isize-is,nblok1); i++ ) {
          xig[iblok*nblok1+i] = xi[ibase+is+i];
          yig[iblok*nblok1+i] = yi[ibase+is+i];
          zig[iblok*nblok1+i] = zi[ibase+is+i];
          gxig[iblok*nblok1+i] = 0;
          gyig[iblok*nblok1+i] = 0;
          gzig[iblok*nblok1+i] = 0;
        }
        for( i=isize-is; i<nblok1; i++ ) {
          xig[iblok*nblok1+i] = 0;
          yig[iblok*nblok1+i] = 0;
          zig[iblok*nblok1+i] = 0;
          gxig[iblok*nblok1+i] = 0;
          gyig[iblok*nblok1+i] = 0;
          gzig[iblok*nblok1+i] = 0;
        }
        nvecd[iblok*mblok+10] = nfi[ii];
        nvecd[iblok*mblok+11] = jbase;
        nvecd[iblok*mblok+12] = jsize;
        op += nblok1*jsize;
        iblok++;
      }
    }
    nvecd[0] = idev;
    nvecd[1] = iblok;
    nvecd[2] = mblok;
    nvecd[3] = jc;
    nvecd[4] = 0;
    nvecd[5] = myrank;
    nvecd[6] = mp;
    rbd = rb;
    epsd = eps;
    xmind = xmin;
    ymind = ymin;
    zmind = zmin;
    l2pgpu_(nvecd,&op,&rbd,&epsd,&xmind,&ymind,&zmind,xig,yig,zig,gxig,gyig,gzig,vig,
            arex,arey,arez,aimx,aimy,aimz,fac);
    iblok = 0;
    for( ii=istagp[icall]; ii<=iendgp[icall]; ii++ ) {
      ibase = ndi[0][ii];
      isize = ndi[1][ii]-ibase+1;
      for( is=0; is<isize; is+=nblok1 ) {
        for( i=0; i<std::min(isize-is,nblok1); i++ ) {
          gxd[ibase+is+i] += gxig[iblok*nblok1+i];
          gyd[ibase+is+i] += gyig[iblok*nblok1+i];
          gzd[ibase+is+i] += gzig[iblok*nblok1+i];
        }
        iblok++;
      }
    }
  }
}
