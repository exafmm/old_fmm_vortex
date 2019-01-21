#include <cmath>

void boxn1(int* nc, int& n, int lev) {
  int i,nfx,nfy,nfz;
  n = 0;
  for( i=0; i<lev; i++ ) {
    nfx = nc[0]%2;
    nc[0] /= 2;
    n += nfx*2*int(pow(8,i));

    nfy = nc[1]%2;
    nc[1] /= 2;
    n += nfy*int(pow(8,i));

    nfz = nc[2]%2;
    nc[2] /= 2;
    n += nfz*4*int(pow(8,i));
  }
}
