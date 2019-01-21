#include "../misc/constants.h"

extern int *nfi,*nfj,*nc,**npx,**neij,*nij,*njb;
extern float *sr;
extern std::complex<double> (*ax)[mpsym],(*bx)[mpsym],(*ay)[mpsym],(*by)[mpsym],(*az)[mpsym],(*bz)[mpsym];
extern std::complex<double> *ynm,***dnm;

extern void boxc(int, int, int*);
extern void boxn1(int*, int&, int);
extern void spharot(std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, int, std::complex<double>**);

void m2l(int nmp, int mp, int lbi, int lbj, int lev, int ini, double rb) {
  int i,j,ii,ix,iy,iz,ij,jj,jb,jx,jy,jz,je,k,jk,jks,n,nk,nks,jkn,jnk;
  double xijc,yijc,zijc,rh,rhj,rhjk,rhjn;
  std::complex<double> aax[mpsym],aay[mpsym],aaz[mpsym],bbx[mpsym],bby[mpsym],bbz[mpsym];
  std::complex<double> aaxd[mpsym],aayd[mpsym],aazd[mpsym],bbxd[mpsym],bbyd[mpsym],bbzd[mpsym];
  std::complex<double> cnm,aaxdd,aaydd,aazdd;

  if( ini != 0 ) {
    for( i=0; i<ini; i++ ) {
      for( j=0; j<nmp; j++ ) {
        ax[i][j] = 0;
        ay[i][j] = 0;
        az[i][j] = 0;
      }
    }
  }
  for( ii=0; ii<lbi; ii++ ) {
    boxc(nfi[ii],3,nc);
    ix = nc[0];
    iy = nc[1];
    iz = nc[2];
    for( ij=0; ij<nij[ii]; ij++ ) {
      jj = neij[ij][ii];
      jb = njb[jj];
      for( j=0; j<nmp; j++ ) {
        bbx[j] = bx[jb][j];
        bby[j] = by[jb][j];
        bbz[j] = bz[jb][j];
      }
      boxc(nfj[jj],3,nc);
      jx = nc[0]+npx[0][jj]*pow(2,lev);
      jy = nc[1]+npx[1][jj]*pow(2,lev);
      jz = nc[2]+npx[2][jj]*pow(2,lev);
      xijc = (ix-jx)*rb;
      yijc = (iy-jy)*rb;
      zijc = (iz-jz)*rb;
      nc[0] = (ix-jx)+3;
      nc[1] = (iy-jy)+3;
      nc[2] = (iz-jz)+3;
      boxn1(nc,je,3);
      rh = sqrt(xijc*xijc+yijc*yijc+zijc*zijc)+eps;
      spharot(bbx,bby,bbz,bbxd,bbyd,bbzd,mp,dnm[je]);
      rhj = 1;
      for( j=0; j<mp; j++ ) {
        rhjk = rhj;
        rhj *= rh;
        for( k=0; k<=j; k++ ) {
          jk = j*j+j+k;
          jks = j*(j+1)/2+k;
          aaxdd = 0;
          aaydd = 0;
          aazdd = 0;
          rhjn = rhjk;
          rhjk *= rh;
          for( n=abs(k); n<mp; n++ ) {
            rhjn *= rh;
            nk = n*n+n+k;
            nks = n*(n+1)/2+k;
            jkn = jk*mp*mp+nk;
            jnk = (j+n)*(j+n)+j+n;
            cnm = sr[jkn]/rhjn*ynm[jnk];
            aaxdd += bbxd[nks]*cnm;
            aaydd += bbyd[nks]*cnm;
            aazdd += bbzd[nks]*cnm;
          }
          aaxd[jks] = aaxdd;
          aayd[jks] = aaydd;
          aazd[jks] = aazdd;
        }
      }
      spharot(aaxd,aayd,aazd,aax,aay,aaz,mp,dnm[je+nrbm]);
      for( j=0; j<nmp; j++ ) {
        ax[ii][j] += aax[j];
        ay[ii][j] += aay[j];
        az[ii][j] += aaz[j];
      }
    }
  }
  for( jj=0; jj<lbj; jj++ ) {
    jb = njb[jj];
    for( j=0; j<nmp; j++ ) {
      bx[jb][j] = 0;
      by[jb][j] = 0;
      bz[jb][j] = 0;
    }
  }
}
