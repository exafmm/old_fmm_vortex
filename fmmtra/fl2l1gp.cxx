#include "../misc/constants.h"

extern int *nfi,*neo,*nc;
extern std::complex<double> (*ax)[mpsym],(*axo)[mpsym],(*ay)[mpsym],(*ayo)[mpsym],(*az)[mpsym],(*azo)[mpsym];
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

void l2l1(int nmp, int mp, int lbi, double rb) {
  int idev,mblok,mpdnm,mmdnm,m,n,npm,nmm,je,k,nk,nmk,lbio,ii,i,ni,ncall,icall,iblok,ic,nfip,nfic,ib,jbase,jsize,im;
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

  lbio = lbi;
  if( lbio < 8 ) lbio = 8;
  for( ii=0; ii<lbio; ii++ ) {
    for( i=0; i<nmp; i++ ) {
      axo[ii][i] = ax[ii][i];
      ayo[ii][i] = ay[ii][i];
      azo[ii][i] = az[ii][i];
    }
  }

  ni = 0;
  ncall = 0;
  istagp[0] = 0;
  for( ii=0; ii<lbi; ii++ ) {
    ni += nblok1;
    if( ni > njmax ) {
      iendgp[ncall] = ii-1;
      ncall++;
      istagp[ncall] = ii;
      ni = nblok1;
    }
  }
  iendgp[ncall] = lbi-1;
  ncall++;

  for( icall=0; icall<ncall; icall++ ) {
    iblok = 0;
    ic = 0;
    op = 0;
    for( ii=istagp[icall]; ii<=iendgp[icall]; ii++ ) {
      nfip = nfi[ii]/8;
      nfic = nfi[ii]%8;
      boxc(nfic,3,nc);
      nc[0] = nc[0]*2+2;
      nc[1] = nc[1]*2+2;
      nc[2] = nc[2]*2+2;
      boxn1(nc,je,3);
      ib = neo[nfip];
      jbase = ic;
      for( i=0; i<nmp; i++ ) {
        brex[ic] = std::real(axo[ib][i]);
        brey[ic] = std::real(ayo[ib][i]);
        brez[ic] = std::real(azo[ib][i]);
        bimx[ic] = std::imag(axo[ib][i]);
        bimy[ic] = std::imag(ayo[ib][i]);
        bimz[ic] = std::imag(azo[ib][i]);
        ic++;
      }
      jsize = ic-jbase;
      nvecd[iblok*mblok+10] = 1;
      nvecd[iblok*mblok+11] = jbase;
      nvecd[iblok*mblok+12] = je+1;
      op += nblok1*jsize;
      iblok++;
    }

    nvecd[0] = idev;
    nvecd[1] = iblok;
    nvecd[2] = mblok;
    nvecd[3] = ic;
    nvecd[4] = 2;
    nvecd[5] = myrank;
    nvecd[6] = mp;
    nvecd[7] = nrbm;
    rbd = rb;
    epsd = eps;
    m2lgpu_(nvecd,&op,&rbd,&epsd,arex,arey,arez,aimx,aimy,aimz,brex,brey,brez,bimx,bimy,bimz,
            ynmre,ynmim,dnmre,dnmim);
    iblok = 0;
    for( ii=istagp[icall]; ii<=iendgp[icall]; ii++ ) {
      for( i=0; i<nmp; i++ ) {
        im = iblok*nblok1+i;
        ax[ii][i] = std::complex<double>(arex[im],aimx[im]);
        ay[ii][i] = std::complex<double>(arey[im],aimy[im]);
        az[ii][i] = std::complex<double>(arez[im],aimz[im]);
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
