#include "../misc/constants.h"

extern double **psm;

void pm2m(double* ggx, double* ggy, double* ggz, double* ggxd, double* ggyd, double* ggzd, int nmp, int ni) {
  int j,k,njk;
  double ggxdd,ggydd,ggzdd;

  for( j=0; j<nmp; j++ ) {
    ggxdd = 0;
    ggydd = 0;
    ggzdd = 0;
    for( k=0; k<nmp; k++ ) {
      njk = k*nmp+j;
      ggxdd += ggxd[k]*psm[ni][njk];
      ggydd += ggyd[k]*psm[ni][njk];
      ggzdd += ggzd[k]*psm[ni][njk];
    }
    ggx[j] = ggxdd;
    ggy[j] = ggydd;
    ggz[j] = ggzdd;
  }
}
