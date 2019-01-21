#include "../misc/constants.h"

extern int *nej,*nfo,*nc;
extern std::complex<double> (*bx)[mpsym],(*by)[mpsym],(*bz)[mpsym];
extern std::complex<double> *ynm,***dnm;
extern int *istagp,*iendgp,*nvecd;
extern float *arex,*arey,*arez,*aimx,*aimy,*aimz;
extern float *brex,*brey,*brez,*bimx,*bimy,*bimz;
extern float *ynmre,*ynmim,*dnmre,*dnmim;

extern void memoryuse();
extern void memoryfree();
extern void boxc(int, int, int*);
extern void boxn1(int*, int&, int);
extern "C" void m2lgpu_(int*, double*, double*, double*,
        float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*);

void m2m1(int nmp, int mp, int lev, int lbj, int lbjo, int* nlbj, double rb) {
  int idev,mblok,mpdnm,mmdnm,m,n,npm,nmm,je,k,nk,nmk,ii,ib,j,nj,ncall,jj,icall,iblok,jc,nfjc,jb,jbase,jsize,nfjp,jm;
  double op,rbd,epsd;

  idev = myrank%ngpu;
  mblok = 3;
  mpdnm = (4*mp*mp*mp-mp)/3;
  mmdnm = mpdnm*2*nrbm;

  ynmre = new float [4*mpmax];
  ynmim = new float [4*mpmax];
  dnmre = new float [mmdnm];
  dnmim = new float [mmdnm];
  mem = mpmax*8*4+mmdnm*2*4;
  memoryuse();

  for( m=0; m<2*mp; m++ ) {
    for( n=m; n<2*mp; n++ ) {
      npm = n*n+n+m;
      nmm = n*n+n-m;
      ynmre[npm] = std::real(ynm[npm]);
      ynmre[nmm] = std::real(ynm[nmm]);
      ynmim[npm] = std::imag(ynm[npm]);
      ynmim[nmm] = std::imag(ynm[nmm]);
    }
  }
  for( je=0; je<2*nrbm; je++ ) {
    for( n=0; n<mp; n++ ) {
      for( m=0; m<=n; m++ ) {
        for( k=-n; k<=n; k++ ) {
          nk = n*(n+1)+k;
          nmk = (4*n*n*n+6*n*n+5*n)/3+m*(2*n+1)+k+je*mpdnm;
          dnmre[nmk] = std::real(dnm[je][m][nk]);
          dnmim[nmk] = std::imag(dnm[je][m][nk]);
        }
      }
    }
  }

  if( lev==1 ) lbj = 8;
  for( ii=0; ii<lbj; ii++ ) {
    ib = ii+nlbj[lev-1];
    for( j=0; j<nmp; j++ ) {
      bx[ib][j] = 0;
      by[ib][j] = 0;
      bz[ib][j] = 0;
    }
  }

  nj = 0;
  ncall = 0;
  istagp[0] = 0;
  for( jj=0; jj<lbjo; jj++ ) {
    nj += nblok1;
    if( nj > njmax ) {
      iendgp[ncall] = jj-1;
      ncall++;
      istagp[ncall] = jj;
      nj = nblok1;
    }
  }
  iendgp[ncall] = lbjo-1;
  ncall++;

  for( icall=0; icall<ncall; icall++ ) {
    iblok = 0;
    jc = 0;
    op = 0;
    for( jj=istagp[icall]; jj<=iendgp[icall]; jj++ ) {
      nfjc = nfo[jj]%8;
      boxc(nfjc,3,nc);
      nc[0] = 4-nc[0]*2;
      nc[1] = 4-nc[1]*2;
      nc[2] = 4-nc[2]*2;
      boxn1(nc,je,3);
      jb = jj+nlbj[lev];
      jbase = jc;
      for( j=0; j<nmp; j++ ) {
        brex[jc] = std::real(bx[jb][j]);
        brey[jc] = std::real(by[jb][j]);
        brez[jc] = std::real(bz[jb][j]);
        bimx[jc] = std::imag(bx[jb][j]);
        bimy[jc] = std::imag(by[jb][j]);
        bimz[jc] = std::imag(bz[jb][j]);
        jc++;
      }
      jsize = jc-jbase;
      nvecd[iblok*mblok+10] = 1;
      nvecd[iblok*mblok+11] = jbase;
      nvecd[iblok*mblok+12] = je+1;
      op += nblok1*jsize;
      iblok++;
    }

    nvecd[0] = idev;
    nvecd[1] = iblok;
    nvecd[2] = mblok;
    nvecd[3] = jc;
    nvecd[4] = 0;
    nvecd[5] = myrank;
    nvecd[6] = mp;
    nvecd[7] = nrbm;
    rbd = rb;
    epsd = eps;
    m2lgpu_(nvecd,&op,&rbd,&epsd,arex,arey,arez,aimx,aimy,aimz,brex,brey,brez,bimx,bimy,bimz,
            ynmre,ynmim,dnmre,dnmim);
    iblok = 0;
    for( jj=istagp[icall]; jj<=iendgp[icall]; jj++ ) {
      nfjp = nfo[jj]/8;
      if( lev == 1 ) {
        ib = nfjp+nlbj[lev-1];
      } else {
        ib = nej[nfjp]+nlbj[lev-1];
      }
      for( j=0; j<nmp; j++ ) {
        jm = iblok*nblok1+j;
        bx[ib][j] += std::complex<double>(arex[jm],aimx[jm]);
        by[ib][j] += std::complex<double>(arey[jm],aimy[jm]);
        bz[ib][j] += std::complex<double>(arez[jm],aimz[jm]);
      }
      iblok++;
    }
  }
  delete[] ynmre;
  delete[] ynmim;
  delete[] dnmre;
  delete[] dnmim;
  mem= mpmax*8*4+mmdnm*2*4;
  memoryfree();
}
