#include "../misc/constants.h"

extern float *gxi,*gyi,*gzi,*gxd,*gyd,*gzd,*vd;

extern void memoryuse();
extern void memoryfree();
extern void dirgp(int, int, int, int, int);

void dird(int n0, int n1, int n2, int n3) {
  int i;

  gxd = new float [npmax];
  gyd = new float [npmax];
  gzd = new float [npmax];
  vd  = new float [npmax];
  mem = npmax*4*4;
  memoryuse();
  if( neq == 5 ) {
    for( i=0; i<npmax; i++ ) gxi[i] = 0;
  } else if ( neq == 6 ) {
    for( i=0; i<npmax; i++ ) gyi[i] = 0;
  } else if ( neq == 7 ) {
    for( i=0; i<npmax; i++ ) gzi[i] = 0;
  }
  dirgp(n0,n1,n2,n3,neq+10);
  if( neq == 5 ) {
    for( i=n0; i<=n1; i++ ) gxi[i] = gxd[i];
  } else if ( neq == 6 ) {
    for( i=n0; i<=n1; i++ ) gyi[i] = gxd[i];
  } else if ( neq == 7 ) {
    for( i=n0; i<=n1; i++ ) gzi[i] = gxd[i];
  }
  delete[] gxd;
  delete[] gyd;
  delete[] gzd;
  delete[] vd;
  mem = npmax*4*4;
  memoryfree();
}
