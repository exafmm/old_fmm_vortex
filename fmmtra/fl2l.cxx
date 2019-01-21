#include "../misc/constants.h"

extern float *anm;
extern std::complex<double> *ynm,***dnm;

extern void spharot(std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, int, std::complex<double>**);

void fl2l(std::complex<double>* aax, std::complex<double>* aay, std::complex<double>* aaz, std::complex<double>* aaxd, std::complex<double>* aayd, std::complex<double>* aazd, double rh, int mp, int je) {
  int j,k,jk,jks,n,jnk,nk,nks;
  std::complex<double> cnm,aaxdd,aaydd,aazdd;

  spharot(aaxd,aayd,aazd,aax,aay,aaz,mp,dnm[je]);
  for( j=0; j<mp; j++ ) {
    for( k=0; k<=j; k++ ) {
      jk = j*j+j+k;
      jks = j*(j+1)/2+k;
      aaxdd = 0;
      aaydd = 0;
      aazdd = 0;
      for( n=j; n<mp; n++ ) {
        jnk = (n-j)*(n-j)+n-j;
        nk = n*n+n+k;
        nks = n*(n+1)/2+k;
        cnm = anm[jnk]*anm[jk]/anm[nk]*pow(rh,n-j)*ynm[jnk];
        aaxdd += aax[nks]*cnm;
        aaydd += aay[nks]*cnm;
        aazdd += aaz[nks]*cnm;
      }
      aaxd[jks] = aaxdd;
      aayd[jks] = aaydd;
      aazd[jks] = aazdd;
    }
  }
  spharot(aaxd,aayd,aazd,aax,aay,aaz,mp,dnm[je+nrbm]);
}
