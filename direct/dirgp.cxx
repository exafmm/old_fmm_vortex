#include "../misc/parameters.h"
#include "../misc/constants.h"

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern float *gxd,*gyd,*gzd,*vd;
extern int *nvecd;
extern float *xig,*yig,*zig,*gxig,*gyig,*gzig,*vig;
extern float *xjg,*yjg,*zjg,*gxjg,*gyjg,*gzjg,*vjg,*sjg;
extern float *gxdg,*gydg,*gzdg,*vdg;

extern void memoryuse();
extern void memoryfree();
extern "C" void p2pgpu_(int*, double*, double*, double*, double*, double*, double*, double*,
        float*, float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*);

void dirgp(int n0, int n1, int n2, int n3, int neqd) {
  int idev,i,nicall,njcall,icall,iwork1,iwork2,ista,iend,ibase,isize,iblok,is,mblok;
  int jcall,jwork1,jwork2,jsta,jend,jbase,jsize,nj;
  double op,visd,epsd,dxd,dyd,dzd,dxyzd;

  nvecd = new int [nimax];
  xig  = new float [nimax];
  yig  = new float [nimax];
  zig  = new float [nimax];
  gxig = new float [nimax];
  gyig = new float [nimax];
  gzig = new float [nimax];
  vig  = new float [nimax];
  xjg  = new float [njmax];
  yjg  = new float [njmax];
  zjg  = new float [njmax];
  gxjg = new float [njmax];
  gyjg = new float [njmax];
  gzjg = new float [njmax];
  vjg  = new float [njmax];
  sjg  = new float [njmax];
  gxdg = new float [nimax];
  gydg = new float [nimax];
  gzdg = new float [nimax];
  vdg  = new float [nimax];
  mem = nimax*12*4+njmax*8*4;
  memoryuse();

  idev = myrank%2;
  for( i=n0; i<=n1; i++ ) {
    gxd[i] = 0;
    gyd[i] = 0;
    gzd[i] = 0;
  }
  nicall = (n1-n0+1)/nimax+1;
  njcall = (n3-n2+1)/njmax+1;
  for( icall=0; icall<nicall; icall++ ) {
    iwork1 = (n1-n0+1)/nicall;
    iwork2 = (n1-n0+1)%nicall;
    ista = icall*iwork1+n0+std::min(icall,iwork2);
    iend = ista+iwork1-1;
    if( iwork2 > icall ) iend++;
    ibase = ista;
    isize = iend-ibase+1;
    iblok = 0;
    for( is=0; is<isize; is+=nblok0 ) {
      for( i=0; i<std::min(isize-is,nblok0); i++ ) {
        xig[iblok*nblok0+i] = xi[ibase+is+i];
        yig[iblok*nblok0+i] = yi[ibase+is+i];
        zig[iblok*nblok0+i] = zi[ibase+is+i];
        gxig[iblok*nblok0+i] = gxi[ibase+is+i];
        gyig[iblok*nblok0+i] = gyi[ibase+is+i];
        gzig[iblok*nblok0+i] = gzi[ibase+is+i];
        vig[iblok*nblok0+i] = vi[ibase+is+i];
      }
      for( i=isize-is; i<nblok0; i++ ) {
        xig[iblok*nblok0+i] = 0;
        yig[iblok*nblok0+i] = 0;
        zig[iblok*nblok0+i] = 0;
        gxig[iblok*nblok0+i] = 0;
        gyig[iblok*nblok0+i] = 0;
        gzig[iblok*nblok0+i] = 0;
        vig[iblok*nblok0+i] = 0;
      }
      iblok++;
    }
    mblok = 3;
    for( jcall=0; jcall<njcall; jcall++ ) {
      jwork1 = (n3-n2+1)/njcall;
      jwork2 = (n3-n2+1)%njcall;
      jsta = jcall*jwork1+n2+std::min(jcall,jwork2);
      jend = jsta+jwork1;
      if( jwork2 > jcall ) jend++;
      jbase = jsta;
      jsize = jend-jbase;
      for( i=0; i<iblok; i++ ) {
        nvecd[i*mblok+10] = 1;
        nvecd[i*mblok+11] = 0;
        nvecd[i*mblok+12] = jsize;
      }
      for( i=jsta; i<jend; i++ ) {
        nj = i-jsta;
        xjg[nj] = xj[i];
        yjg[nj] = yj[i];
        zjg[nj] = zj[i];
        gxjg[nj] = gxj[i];
        gyjg[nj] = gyj[i];
        gzjg[nj] = gzj[i];
        vjg[nj] = vj[i];
        sjg[nj] = sj[i];
      }
      nj++;
      op = (double) isize*jsize;
      nvecd[0] = idev;
      nvecd[1] = iblok;
      nvecd[2] = mblok;
      nvecd[3] = nj;
      nvecd[4] = neqd;
      nvecd[5] = myrank;
      visd = vis;
      epsd = eps;
      dxd = dx;
      dyd = dy;
      dzd = dz;
      dxyzd = dxyz;
      p2pgpu_(nvecd,&op,&visd,&epsd,&dxd,&dyd,&dzd,&dxyzd,xig,yig,zig,gxig,gyig,gzig,vig,
              xjg,yjg,zjg,gxjg,gyjg,gzjg,vjg,sjg,gxdg,gydg,gzdg,vdg);
      iblok = 0;
      for( is=0; is<isize; is+=nblok0 ) {
        for( i=0; i<std::min(isize-is,nblok0); i++ ) {
          gxd[ibase+is+i] += gxdg[iblok*nblok0+i];
          gyd[ibase+is+i] += gydg[iblok*nblok0+i];
          gzd[ibase+is+i] += gzdg[iblok*nblok0+i];
          vd[ibase+is+i] += vdg[iblok*nblok0+i];
        }
        iblok++;
      }
    }
  }
  delete[] nvecd;
  delete[] xig;
  delete[] yig;
  delete[] zig;
  delete[] gxig;
  delete[] gyig;
  delete[] gzig;
  delete[] vig;
  delete[] xjg;
  delete[] yjg;
  delete[] zjg;
  delete[] gxjg;
  delete[] gyjg;
  delete[] gzjg;
  delete[] vjg;
  delete[] sjg;
  delete[] gxdg;
  delete[] gydg;
  delete[] gzdg;
  delete[] vdg;
  mem = nimax*12*4+njmax*8*4;
  memoryfree();
}
