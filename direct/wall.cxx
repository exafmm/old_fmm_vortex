#include "../misc/parameters.h"
#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern float *gxd,*gyd,*gzd;

extern void memoryuse();
extern void memoryfree();
extern double ierfc(double);

void dir(int n0, int n1, int n2, int n3) {
  int i,j;
  double vist,gxdd,gydd,gzdd,dxij,dyij,dzij,rij,sij,rsij,cutoff;

  gxd = new float [npmax];
  gyd = new float [npmax];
  gzd = new float [npmax];
  mem = npmax*3*4;
  memoryuse();

  vist = 1.0/sqrt(2*vis*dt);
  for( i=n0; i<=n1; i++ ) {
    gxdd = 0;
    gydd = 0;
    gzdd = 0;
    for( j=n2; j<=n3; j++ ){
      dxij = xi[i]-xj[j];
      dyij = yi[i]-yj[j];
      dzij = std::abs(zi[i]-zj[j]);
      rij = dxij*dxij+dyij*dyij+dzij*dzij+eps;
      sij = sj[j]*sj[j];
      if( rij < 50*sij ) {
        cutoff = (erfc((dzij-dz/2)*vist)-erfc((dxij+dz/2)*vist))*
                 (0.5/vist*(ierfc((dxij+dx)*vist)+ierfc((dxij-dx)*vist)-2*ierfc(dxij*vist)))*
                 (0.5/vist*(ierfc((dyij+dy)*vist)+ierfc((dyij-dy)*vist)-2*ierfc(dyij*vist)));
        gxdd += gxj[j]*cutoff;
        gydd += gyj[j]*cutoff;
        gzdd += gzj[j]*cutoff;
      }
    }
    gxd[i] = gxdd;
    gyd[i] = gydd;
    gzd[i] = gzdd;
  }
  for( i=0; i<npmax; i++ ) {
    gxi[i] = 0;
    gyi[i] = 0;
    gzi[i] = 0;
  }
  for( i=n0; i<=n1; i++ ) {
    gxi[i] = gxd[i];
    gyi[i] = gyd[i];
    gzi[i] = gzd[i];
  }

  delete[] gxd;
  delete[] gyd;
  delete[] gzd;
  mem = npmax*3*4;
  memoryfree();
}
