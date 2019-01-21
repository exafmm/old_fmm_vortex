#include "../misc/constants.h"

extern float *anm;
extern std::complex<double> *ynm,***dnm;

extern void spharot(std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, int, std::complex<double>**);

void fm2m(std::complex<double>* bbx, std::complex<double>* bby, std::complex<double>* bbz, std::complex<double>* bbxd, std::complex<double>* bbyd, std::complex<double>* bbzd, double rh, int mp, int je) {
  int j,k,jk,jks,n,jnk,jnks,nm;
  std::complex<double> cnm,bbxdd,bbydd,bbzdd;

  spharot(bbxd,bbyd,bbzd,bbx,bby,bbz,mp,dnm[je]);
  for( j=0; j<mp; j++ ) {
    for( k=0; k<=j; k++ ) {
      jk = j*j+j+k;
      jks = j*(j+1)/2+k;
      bbxdd = 0;
      bbydd = 0;
      bbzdd = 0;
      for( n=0; n<=j-abs(k); n++ ) {
        jnk = (j-n)*(j-n)+j-n+k;
        jnks = (j-n)*(j-n+1)/2+k;
        nm = n*n+n;
        cnm = pow(-1,n)*anm[nm]*anm[jnk]/anm[jk]*pow(rh,n)*ynm[nm];
        bbxdd += bbx[jnks]*cnm;
        bbydd += bby[jnks]*cnm;
        bbzdd += bbz[jnks]*cnm;
      }
      bbxd[jks] = bbxdd;
      bbyd[jks] = bbydd;
      bbzd[jks] = bbzdd;
    }
  }
  spharot(bbxd,bbyd,bbzd,bbx,bby,bbz,mp,dnm[je+nrbm]);
}
