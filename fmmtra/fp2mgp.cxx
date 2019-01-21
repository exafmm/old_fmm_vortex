#include "../misc/constants.h"

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj;
extern int **ndj,*nfj,*nc;
extern std::complex<double> (*bx)[mpsym],(*by)[mpsym],(*bz)[mpsym],*bnm,*bth;
extern float *fac;
extern int *istagp,*iendgp,*nvecd;
extern float *xjg,*yjg,*zjg,*gxjg,*gyjg,*gzjg;
extern float *brex,*brey,*brez,*bimx,*bimy,*bimz;

extern "C" void p2mgpu_(int*, double*, double*, double*, double*, double*, double*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*);

void p2m(int nmp, int mp, int lbj, double rb, int nlm) {
  int idev,mblok,nj,ncall,jj,icall,iblok,jc,jbase,j,jsize,jm;
  double op,rbd,epsd,xmind,ymind,zmind;

  idev = myrank%ngpu;
  mblok = 3;

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
        jc++;
      }
      jsize = jc-jbase;
      nvecd[iblok*mblok+10] = nfj[jj];
      nvecd[iblok*mblok+11] = jbase;
      nvecd[iblok*mblok+12] = jsize;
      op += nblok1*jsize;
      iblok++;
    }
    nvecd[0] = idev;
    nvecd[1] = iblok;
    nvecd[2] = mblok;
    nvecd[3] = jc;
    nvecd[4] = neq;
    nvecd[5] = myrank;
    nvecd[6] = mp;
    rbd = rb;
    epsd = eps;
    xmind = xmin;
    ymind = ymin;
    zmind = zmin;
    p2mgpu_(nvecd,&op,&rbd,&epsd,&xmind,&ymind,&zmind,xjg,yjg,zjg,gxjg,gyjg,gzjg,brex,brey,brez,bimx,bimy,bimz,fac);
    iblok = 0;
    for( jj=istagp[icall]; jj<=iendgp[icall]; jj++ ) {
      for( j=0; j<nmp; j++ ) {
        jm = iblok*nblok1+j;
        bx[jj+nlm][j] = std::complex<double>(brex[jm],bimx[jm]);
        by[jj+nlm][j] = std::complex<double>(brey[jm],bimy[jm]);
        bz[jj+nlm][j] = std::complex<double>(brez[jm],bimz[jm]);
      }
      iblok++;
    }
  }
}
