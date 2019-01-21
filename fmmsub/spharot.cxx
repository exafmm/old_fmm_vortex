#include "../misc/constants.h"

void spharot(std::complex<double>* yx, std::complex<double>* yy, std::complex<double>* yz, std::complex<double>* yxr, std::complex<double>* yyr, std::complex<double>* yzr, int mp, std::complex<double>** dnm) {
  int n,m,nms,k,nk,nks;
  std::complex<double> yxd,yyd,yzd;

  for( n=0; n<mp; n++ ) {
    for( m=0; m<=n; m++ ) {
      nms = n*(n+1)/2+m;
      yxd = 0;
      yyd = 0;
      yzd = 0;
      for( k=-n; k<=-1; k++ ) {
        nk = n*(n+1)+k;
        nks = n*(n+1)/2-k;
        yxd += dnm[m][nk]*conj(yx[nks]);
        yyd += dnm[m][nk]*conj(yy[nks]);
        yzd += dnm[m][nk]*conj(yz[nks]);
      }
      for( k=0; k<=n; k++ ) {
        nk = n*(n+1)+k;
        nks = n*(n+1)/2+k;
        yxd += dnm[m][nk]*yx[nks];
        yyd += dnm[m][nk]*yy[nks];
        yzd += dnm[m][nk]*yz[nks];
      }
      yxr[nms] = yxd;
      yyr[nms] = yyd;
      yzr[nms] = yzd;
    }
  }
}
