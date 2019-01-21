#include "../misc/constants.h"

extern std::complex<double> (*ax)[mpsym],(*axo)[mpsym],(*bx)[mpsym];
extern std::complex<double> (*ay)[mpsym],(*ayo)[mpsym],(*by)[mpsym];
extern std::complex<double> (*az)[mpsym],(*azo)[mpsym],(*bz)[mpsym];
extern std::complex<double> *bnm,*bth,*ynm,***dnm;
extern float *fac,*sr,*anm;
extern double (*px)[nspm],(*pxo)[nspm],(*gx)[nspm];
extern double (*py)[nspm],(*pyo)[nspm],(*gy)[nspm];
extern double (*pz)[nspm],(*pzo)[nspm],(*gz)[nspm];
extern double *xsp,*ysp,*zsp,**psm,**p2g;
extern int *ibase,*isize,*jbase,*jsize;
extern int **istamd,**iendmd,*jstamd,*jendmd,*nicall;
extern int *istagp,*iendgp,**jstagp,**jendgp,*njcall,*nvecd,*njj;
extern double *bxmd,*bymd,*bzmd,(*xmd)[3],(*ymd)[3],(*zmd)[3],(*pos)[3],(*pod)[3];
extern float *xig,*yig,*zig,*gxig,*gyig,*gzig,*vig;
extern float *xjg,*yjg,*zjg,*gxjg,*gyjg,*gzjg,*vjg,*sjg;
extern float *gxdg,*gydg,*gzdg,*vdg;
extern float *arex,*arey,*arez,*aimx,*aimy,*aimz,*brex,*brey,*brez,*bimx,*bimy,*bimz;

extern void memoryuse();
extern void memoryfree();

void alloc() {
  int i,j;
  if( nfmm%2 == 1 ) {
    px = new double [nbne][nspm];
    pxo = new double [nbne][nspm];
    gx = new double [nbnes][nspm];
    py = new double [nbne][nspm];
    pyo = new double [nbne][nspm];
    gy = new double [nbnes][nspm];
    pz = new double [nbne][nspm];
    pzo = new double [nbne][nspm];
    gz = new double [nbnes][nspm];
    xsp = new double [nspm];
    ysp = new double [nspm];
    zsp = new double [nspm];
    psm = new double* [8];
    for( i=0; i<8; i++ ) psm[i] = new double [nspm*nspm];
    p2g = new double* [nwmax];
    for( i=0; i<nwmax; i++ ) p2g[i] = new double [nwmax];
    mem = nspm*nbne*6*8+nspm*nbnes*3*8+nspm*3*8+nspm*nspm*8*8+nwmax*nwmax*8;
    memoryuse();
    m2lrp2p = 5000;
  } else {
    ax = new std::complex<double> [nbne][mpsym];
    axo = new std::complex<double> [nbne][mpsym];
    bx = new std::complex<double> [nbnes][mpsym];
    ay = new std::complex<double> [nbne][mpsym];
    ayo = new std::complex<double> [nbne][mpsym];
    by = new std::complex<double> [nbnes][mpsym];
    az = new std::complex<double> [nbne][mpsym];
    azo = new std::complex<double> [nbne][mpsym];
    bz = new std::complex<double> [nbnes][mpsym];
    bnm = new std::complex<double> [4*mpmax];
    bth = new std::complex<double> [mpmax];
    sr = new float [mpmax*mpmax];
    ynm = new std::complex<double> [4*mpmax];
    dnm = new std::complex<double>** [2*nrbm];
    for( i=0; i<2*nrbm; i++ ) {
      dnm[i] = new std::complex<double>* [mpcmp];
      for( j=0; j<mpcmp; j++ ) dnm[i][j] = new std::complex<double> [mpmax];
    }
    fac = new float [4*mpmax];
    anm = new float [4*mpmax];
    mem = mpsym*nbne*6*16+mpsym*nbnes*3*16+mpmax*17*16+mpmax*mpmax*4+mpmax*mpcmp*2*nrbm*16;
    memoryuse();
    m2lrp2p = 200;
  }

  if( ndev == 1 ) {
    ibase = new int [nbnp];
    isize = new int [nbnp];
    jbase = new int [nbnp];
    jsize = new int [nbnp];
    istamd = new int* [nebm];
    for( i=0; i<nebm; i++ ) istamd[i] = new int [nbnp];
    iendmd = new int* [nebm];
    for( i=0; i<nebm; i++ ) iendmd[i] = new int [nbnp];
    jstamd = new int [nbnp];
    jendmd = new int [nbnp];
    nicall = new int [nbnp];
    bxmd = new double [mimax];
    bymd = new double [mimax];
    bzmd = new double [mimax];
    xmd = new double [mimax][3];
    ymd = new double [mimax][3];
    zmd = new double [mimax][3];
    pos = new double [mimax][3];
    pod = new double [mimax][3];
    mem = nbnp*7*4+nebm*nbnp*2*4+mimax*15*8;
    memoryuse();
  } else if( ndev == 2 ) {
    ibase = new int [nbnp];
    isize = new int [nbnp];
    jbase = new int [nbnp];
    jsize = new int [nbnp];
    istagp = new int [nbne];
    iendgp = new int [nbne];
    jstagp = new int* [nebm];
    for( i=0; i<nebm; i++ ) jstagp[i] = new int [nbnp];
    jendgp = new int* [nebm];
    for( i=0; i<nebm; i++ ) jendgp[i] = new int [nbnp];
    njcall = new int [nbnp];
    nvecd = new int [npmax];
    njj = new int [nbnp];
    xig = new float [nimax];
    yig = new float [nimax];
    zig = new float [nimax];
    gxig = new float [nimax];
    gyig = new float [nimax];
    gzig = new float [nimax];
    vig = new float [nimax];
    xjg = new float [njmax];
    yjg = new float [njmax];
    zjg = new float [njmax];
    gxjg = new float [njmax];
    gyjg = new float [njmax];
    gzjg = new float [njmax];
    vjg = new float [njmax];
    sjg = new float [njmax];
    gxdg = new float [nimax];
    gydg = new float [nimax];
    gzdg = new float [nimax];
    vdg = new float [nimax];
    arex = new float [njmax];
    arey = new float [njmax];
    arez = new float [njmax];
    aimx = new float [njmax];
    aimy = new float [njmax];
    aimz = new float [njmax];
    brex = new float [njmax];
    brey = new float [njmax];
    brez = new float [njmax];
    bimx = new float [njmax];
    bimy = new float [njmax];
    bimz = new float [njmax];
    mem = nbnp*6*4+nbne*2*4+nebm*nbnp*2*4+nimax*18*4+njmax*20*4;
    memoryuse();
  }
}

void dealloc() {
  int i,j;
  if( nfmm%2 == 1 ) {
    delete[] px;
    delete[] pxo;
    delete[] gx;
    delete[] py;
    delete[] pyo;
    delete[] gy;
    delete[] pz;
    delete[] pzo;
    delete[] gz;
    delete[] xsp;
    delete[] ysp;
    delete[] zsp;
    for( i=0; i<8; i++ ) delete[] psm[i];
    delete[] psm;
    for( i=0; i<nwmax; i++ ) delete[] p2g[i];
    delete[] p2g;
    mem = nspm*nbne*6*8+nspm*nbnes*3*8+nspm*3*8+nspm*nspm*8*8+nwmax*nwmax*8;
    memoryfree();
  } else {
    delete[] ax;
    delete[] axo;
    delete[] bx;
    delete[] ay;
    delete[] ayo;
    delete[] by;
    delete[] az;
    delete[] azo;
    delete[] bz;
    delete[] bnm;
    delete[] bth;
    delete[] sr;
    delete[] ynm;
    for( i=0; i<2*nrbm; i++ ) {
      for( j=0; j<mpcmp; j++ ) delete[] dnm[i][j];
      delete[] dnm[i];
    }
    delete[] dnm;
    delete[] fac;
    delete[] anm;
    mem = mpsym*nbne*6*16+mpsym*nbnes*3*16+mpmax*17*16+mpmax*mpmax*4+mpmax*mpcmp*2*nrbm*16;
    memoryfree();
  }  

  if( ndev == 1 ) {
    delete[] ibase;
    delete[] isize;
    delete[] jbase;
    delete[] jsize;
    for( i=0; i<nebm; i++ ) delete[] istamd[i];
    delete[] istamd;
    for( i=0; i<nebm; i++ ) delete[] iendmd[i];
    delete[] iendmd;
    delete[] jstamd;
    delete[] jendmd;
    delete[] nicall;
    delete[] bxmd;
    delete[] bymd;
    delete[] bzmd;
    delete[] xmd;
    delete[] ymd;
    delete[] zmd;
    delete[] pos;
    delete[] pod;
    mem = nbnp*7*4+nebm*nbnp*2*4+mimax*15*8;
    memoryfree();
  } else if( ndev == 2 ) {
    delete[] ibase;
    delete[] isize;
    delete[] jbase;
    delete[] jsize;
    delete[] istagp;
    delete[] iendgp;
    for( i=0; i<nebm; i++ ) delete[] jstagp[i];
    delete[] jstagp;
    for( i=0; i<nebm; i++ ) delete[] jendgp[i];
    delete[] jendgp;
    delete[] njcall;
    delete[] nvecd;
    delete[] njj;
    delete[] xig;
    delete[] yig;
    delete[] zig;
    delete[] gxig;
    delete[] gyig;
    delete[] gzig;
    delete[] vig;
    delete[] xjg;
    delete[] yjg;
    delete[] zjg;
    delete[] gxjg;
    delete[] gyjg;
    delete[] gzjg;
    delete[] vjg;
    delete[] sjg;
    delete[] gxdg;
    delete[] gydg;
    delete[] gzdg;
    delete[] vdg;
    delete[] arex;
    delete[] arey;
    delete[] arez;
    delete[] aimx;
    delete[] aimy;
    delete[] aimz;
    delete[] brex;
    delete[] brey;
    delete[] brez;
    delete[] bimx;
    delete[] bimy;
    delete[] bimz;
    mem = nbnp*6*4+nbne*2*4+nebm*nbnp*2*4+nimax*18*4+njmax*20*4;
    memoryfree();
  }
}
