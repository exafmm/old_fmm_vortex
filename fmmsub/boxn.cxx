#include "../misc/constants.h"

extern int **nxs,*nfn;

void boxn(int n, int ndi, int lev) {
  int i,j,nfx,nfy,nfz;
  for( i=0; i<n; i++ ) nfn[i] = 0;
  for( i=0; i<lev; i++ ) {
    for( j=0; j<n; j++ ) {
      nfx = nxs[0][j]%2;
      nxs[0][j] /= 2;
      nfn[j] += nfx*int(pow(2,ndi*i+1));
    }

    if( ndi >= 2 ) {
      for( j=0; j<n; j++ ) {
        nfy = nxs[1][j]%2;
        nxs[1][j] /= 2;
        nfn[j] += nfy*int(pow(2,ndi*i));
      }
    }

    if( ndi == 3 ) {
      for( j=0; j<n; j++ ) {
        nfz = nxs[2][j]%2;
        nxs[2][j] /= 2;
        nfn[j] += nfz*int(pow(2,ndi*i+2));
      }
    }

  }
}
