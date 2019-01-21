#define MAIN
#include "../misc/parameters.h"
#include "../misc/constants.h"
#include "../misc/arrays.h"
#undef MAIN
extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*vi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern float *xo,*yo,*zo,*xe,*ye,*ze,*gxe,*gye,*gze,*se;

extern void memoryuse();
extern void memoryfree();
extern void fmm(int, int, int, double*);
extern void nlevel(int);

int main(){
  char *udir,*umpi,*udev,*ufmm,*ugeo,*usmt,*uequ;
  char ufile[32];
  int j,nbyte,i,nbb,ksta,kend,ksum,lev,it,np,ni,nj,lmaxd;
  const int nmax = 2.1e6;
  double s,sv,tic,toc,t1,t2;
  std::fstream fid1,fid2,fid4;

  ncheck = 0;
  umem = 0;
  nprocs = 1;
  myrank = 0;
  xi = new float [npmax];
  yi = new float [npmax];
  zi = new float [npmax];
  gxi = new float [npmax];
  gyi = new float [npmax];
  gzi = new float [npmax];
  vi = new float [npmax];
  xj = new float [npmax];
  yj = new float [npmax];
  zj = new float [npmax];
  gxj = new float [npmax];
  gyj = new float [npmax];
  gzj = new float [npmax];
  vj = new float [npmax];
  sj = new float [npmax];
  xo = new float [npmax];
  yo = new float [npmax];
  zo = new float [npmax];
  xe = new float [npmax];
  ye = new float [npmax];
  ze = new float [npmax];
  gxe = new float [npmax];
  gye = new float [npmax];
  gze = new float [npmax];
  se = new float [npmax];
  mem = npmax*25*4;
  memoryuse();

  udir = getenv("Z0_DIR");
  umpi = getenv("Z1_MPI");
  udev = getenv("Z2_DEV");
  ufmm = getenv("Z3_FMM");
  ugeo = getenv("Z4_GEO");
  usmt = getenv("Z5_SMT");
  uequ = getenv("Z6_EQU");
  ndir = atoi(udir);
  nmpi = atoi(umpi);
  ndev = atoi(udev);
  nfmm = atoi(ufmm);
  ngeo = atoi(ugeo);
  nsmt = atoi(usmt);
  nequ = atoi(uequ);
  printf("||-----------------------------------------------------------------||\n");
  printf("||   fmm :  %s     geo :  %s     smt :  %s     equ :  %s          run  ||\n",
    ufmm,ugeo,usmt,uequ);
  printf("||-----------------------------------------------------------------||\n");
  sprintf(ufile,"../../dat/%s%s%s",ugeo,usmt,uequ);

  fid1.open("otest.dat",std::ios::out|std::ios::binary|std::ios::ate);
  fid2.open("otime.dat",std::ios::out|std::ios::binary|std::ios::ate);
  fid4.open("../../dat/coordinates",std::ios::in|std::ios::binary);
  for( j=0; j<npmax-nmax+1; j+=nmax ) {
    fid4.read((char *)(&nbyte),sizeof(int));
    binary_read<float>(fid4,xi+j,nmax);
    binary_read<float>(fid4,yi+j,nmax);
    binary_read<float>(fid4,zi+j,nmax);
    binary_read<float>(fid4,xe+j,nmax);
    binary_read<float>(fid4,ye+j,nmax);
    binary_read<float>(fid4,ze+j,nmax);
    binary_read<float>(fid4,gxe+j,nmax);
    binary_read<float>(fid4,gye+j,nmax);
    binary_read<float>(fid4,gze+j,nmax);
  }
  fid4.close();

  if( ngeo == 0 || nequ > 3 ) {
    npb = 0;
  } else {
    npb = 3;
  } 
  nbb = pow(2,npb); 
  stx = 0;
  s = 0;
  ksta = 1;
  kend = 2*nbb+1;
  ksum = 6*nbb+3;
  if( ngeo == 2 ) { 
    stx = pi/2;
    s = 1;
  } else if( ngeo == 3 ) {
    ksta = nbb+1;
    kend = nbb+1;
    ksum = 5*nbb+3; 
  } 
  xmin = -pi;
  ymin = -pi;
  zmin = -pi;
  rd = 2*pi;

  lev = 2;
  for( it=0; it<25; it++ ) {
    np = pow(10,(it+24)/8.0);
    ni = np;
    nj = np;
    printf("N = %d\n",np);
    dx = 2*pi*pow(np,-1.0/3);
    dy = dx;
    dz = dx;
    dxyz = dx;

    for( i=0; i<nj; i++ ) {
      se[i] = dx;
    }
    sv = s*pow(dx,3);

    for( i=0; i<nj; i++ ) {
      xj[i] = xe[i];
      yj[i] = ye[i];
      zj[i] = ze[i];
      gxj[i] = gxe[i]/nj;
      gyj[i] = gye[i]/nj;
      gzj[i] = gze[i]/nj;
      sj[i] = se[i];
      vj[i] = pow(dx,3);
    }

    for( i=0; i<ni; i++ ) {
      gxi[i] = gxe[i]/ni;
      gyi[i] = gye[i]/ni-sv;
      gzi[i] = gze[i]/ni;
      vi[i] = pow(dx,3);
    }

    for( lmaxd=-lev; lmaxd=-lev-1; lmaxd-- ) {
      lmax = lmaxd;
      for( i=0; i<9; i++ ) tfmm[i] = 0;
      tic = get_time();
      fmm(ni,nj,nequ,tfmm);
      toc = get_time();
      if( lmaxd == -lev ) {
        t1 = toc-tic;
        printf("%g\n",t1);
      } else {
        t2 = toc-tic;
        if( t2 < t1 ) {
          nlevel[nj];
          lev++;
        }
        printf("%g\n",t2);
      }
      for( i=0; i<9; i++ ) {
        nbyte = sizeof(double);
        fid2.write((char *)(&nbyte),sizeof(int));
        fid2.write((char *)(&tfmm[i]),sizeof(double));
        fid2.write((char *)(&nbyte),sizeof(int));
      }
    }
    nbyte = sizeof(int)*2+sizeof(double)*2;
    fid2.write((char *)(&nbyte),sizeof(int));
    fid2.write((char *)(&np),sizeof(int));
    fid2.write((char *)(&lev),sizeof(int));
    fid2.write((char *)(&t1),sizeof(double));
    fid2.write((char *)(&t2),sizeof(double));
    fid2.write((char *)(&nbyte),sizeof(int));
  }

  fid1.close();
  fid2.close();
  delete[] xi;
  delete[] yi;
  delete[] zi;
  delete[] gxi;
  delete[] gyi;
  delete[] gzi;
  delete[] vi;
  delete[] xj;
  delete[] yj;
  delete[] zj;
  delete[] gxj;
  delete[] gyj;
  delete[] gzj;
  delete[] vj;
  delete[] sj;
  delete[] xo;
  delete[] yo;
  delete[] zo;
  delete[] xe;
  delete[] ye;
  delete[] ze;
  delete[] gxe;
  delete[] gye;
  delete[] gze;
  delete[] se;
  mem = npmax*25*4;
  memoryfree();
}
