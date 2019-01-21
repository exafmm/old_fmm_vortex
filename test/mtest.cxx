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
extern void dir(int, int, int, int);
extern void dird(int, int, int, int);

int main(){
  char *udir,*umpi,*udev,*ufmm,*ugeo,*usmt,*uequ;
  char ufile[32];
  int j,nbyte,i,it,np,ni,nj;
  const int nmax = 2100000;
  double td,tm,tic,toc,errn,errd,erre;
  std::fstream fid1,fid3,fid4;

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

  fid1.open("mtest.dat",std::ios::out|std::ios::binary|std::ios::ate);
  if( ndir == 0 ) fid3.open(ufile,std::ios::in|std::ios::binary);
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
  neq = nequ;

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

    if( ndir == 0 ) {
      fid3.read((char *)(&nbyte),sizeof(int));
      fid3.read((char *)(&td),sizeof(double));
      fid3.read((char *)(&nbyte),sizeof(int));
    } else {
      for( i=0; i<ni; i++ ) {
        gxi[i] = gxe[i]/ni;
        gyi[i] = gye[i]/ni;
        gzi[i] = gze[i]/ni;
        vi[i] = pow(dx,3);
      }
      tic = get_time();
      dird(0,ni-1,0,nj-1);
      toc = get_time();
      td = toc-tic;
      for( i=0; i<ni; i++ ) {
        xo[i] = gxi[i];
        yo[i] = gyi[i];
        zo[i] = gzi[i];
      } 
    }
    printf("host   : %g\n",td);

    for( i=0; i<ni; i++ ) {
      gxi[i] = gxe[i]/ni;
      gyi[i] = gye[i]/ni;
      gzi[i] = gze[i]/ni;
      vi[i] = pow(dx,3);
    }
    tic = get_time();
    dir(0,ni-1,0,nj-1);
    toc = get_time();
    tm = toc-tic;
    if( ndev == 1 ) {
      printf("mdg3   : %g\n",tm);
    } else if( ndev == 2 ) {
      printf("gpu    : %g\n",tm);
    } else {
      printf("Z2_DEV must be 1 or 2\n");
    }
    errn = 0;
    for( i=0; i<ni; i++ ) {
      if( ndir == 0 ) {
        fid3.read((char *)(&nbyte),sizeof(int));
        fid3.read((char *)(&xo[i]),sizeof(float));
        fid3.read((char *)(&yo[i]),sizeof(float));
        fid3.read((char *)(&zo[i]),sizeof(float));
        fid3.read((char *)(&nbyte),sizeof(int));
      }
      errd = pow(gxi[i]-xo[i],2)+pow(gyi[i]-yo[i],2)+pow(gzi[i]-zo[i],2);
      erre = pow(xo[i],2)+pow(yo[i],2)+pow(zo[i],2);
      errn += errd/erre/ni;
    }
    errn = sqrt(errn);
    printf("error  : %g\n\n",errn);
    nbyte = sizeof(int)+sizeof(double)*3;
    fid1.write((char *)(&nbyte),sizeof(int));
    fid1.write((char *)(&np),sizeof(int));
    fid1.write((char *)(&td),sizeof(double));
    fid1.write((char *)(&tm),sizeof(double));
    fid1.write((char *)(&errn),sizeof(double));
    fid1.write((char *)(&nbyte),sizeof(int));
  }
  fid1.close();
  if( ndir == 0 ) {
    fid3.close();
  }
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
