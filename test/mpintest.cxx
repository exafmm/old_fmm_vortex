#include "mpi.h"
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
extern void mpirange(int, int, int&, int&);
extern void mpiallreduce(int, int, float*);
extern void mpiprei(int, int&);
extern void mpiprej(int, int&);
extern void mpiposti(int, int&);
extern void mpipostj(int, int&);

int main(int argc, char **argv){
  char *udir,*umpi,*udev,*ufmm,*ugeo,*usmt,*uequ;
  char ufile[32];
  int ncheck,mi,mj,neqd,nbyte,j,i,nbb,ksta,kend,ksum,it,np,ni,nj,k,l,id,jsta,jend,mprocs;
  const int nmax=2100000;
  double tfmm[9],s,sv,cut,tic,xb,yb,zb,toc,td,tf,errn,errd,erre;
  double (*tall)[9] = new double [mimax][9];
  std::fstream fid1,fid2,fid3,fid4;

  ncheck = 0;
  umem = 0;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

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

  if( myrank == 0 ) {
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
  
    if( ndir != 1 ) {
      fid1.open("ntest.dat",std::ios::out|std::ios::binary);
      fid2.open("ntime.dat",std::ios::out|std::ios::binary);
    }
    if( ndir < 1 ) {
      fid3.open(ufile,std::ios::in|std::ios::binary);
    } else if( ndir == 1 ) {
      fid3.open(ufile,std::ios::out|std::ios::binary|std::ios::ate);
    }
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
  }
  MPI_Bcast(&ndir,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&nmpi,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&ndev,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&nfmm,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&ngeo,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&nsmt,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&nequ,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(xi,npmax,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(yi,npmax,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(zi,npmax,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(xe,npmax,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(ye,npmax,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(ze,npmax,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(gxe,npmax,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(gye,npmax,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(gze,npmax,MPI_FLOAT,0,MPI_COMM_WORLD);

  if( ngeo == 0 || nequ > 3 ) {
    npb = 0;
  } else {
    npb = 3;
  }
  nbb = int(pow(2,npb));
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
  neq = nequ;

  for( it=0; it<25; it++ ) {
    np = int(pow(10,(it+24)/8.0));
    ni = np;
    nj = np;
    if( myrank == 0 ) printf("N = %d\n",np);
    dx = 2*pi*pow(np,-1.0/3);
    dy = dx;
    dz = dx;
    dxyz = dx;

    for( i=0; i<nj; i++ ) {
      se[i] = dx;
    }
    sv = s*pow(dx,3);
    if( ndir != 0 ) {
      cut = 2*nbb*pi;
      mprocs = 1;
      if ( 5 <= nequ && nequ <= 7 ) mprocs = nprocs;
      for( i=0; i<ni; i++ ) {
        xo[i] = 0;
        yo[i] = 0;
        zo[i] = 0;
        gxi[i] = gxe[i]/ni/mprocs;
        gyi[i] = (gye[i]/ni-sv)/mprocs;
        gzi[i] = gze[i]/ni/mprocs;
        vi[i] = pow(dx,3);
      }

      MPI_Barrier(MPI_COMM_WORLD);
      tic = MPI_Wtime();
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
      mpirange(0,nj-1,jsta,jend);
      dir(0,ni-1,jsta,jend);
      mpiallreduce(0,ni-1,gxi);
      mpiallreduce(0,ni-1,gyi);
      mpiallreduce(0,ni-1,gzi);
      for( i=0; i<ni; i++ ) {
        xo[i] += gxi[i];
        yo[i] += gyi[i];
        zo[i] += gzi[i];
        gxi[i] = gxe[i]/ni;
        gyi[i] = gye[i]/ni-sv;
        gzi[i] = gze[i]/ni;
        vi[i] = pow(dx,3);
      }
      if( 5 <= nequ && nequ <= 7 ) {
        for( i=0; i<ni; i++ ) {
          gxi[i] = 0;
          gyi[i] = 0;
          gzi[i] = 0;
        }
      }
      if( ngeo != 0 ) {
        id = 0;
        for( i=1; i<=2*nbb+1; i++ ) {
          for( j=1; j<=2*nbb+1; j++ ) {
            for( k=ksta; k<=kend; k++ ) {
              if( i != nbb+1 || j != nbb+1 || k != nbb+1 ) {
                for( l=0; l<nj; l++ ) {
                  xb = xe[l]+(k-nbb-1)*stx;
                  while( xb < -pi ) {
                    xb += 2*pi;
                  }
                  while( xb > pi ) {
                    xb -= 2*pi;
                  }
                  xb += (i-nbb-1)*2*pi;
                  yb = ye[l]+(j-nbb-1)*2*pi;
                  zb = ze[l]+(k-nbb-1)*2*pi;
                  if( -cut < xb && xb < cut &&
                      -cut < yb && yb < cut &&
                      -cut < zb && zb < cut ) {
                    xj[id] = xb;
                    yj[id] = yb;
                    zj[id] = zb;
                    gxj[id] = gxe[l]/nj;
                    gyj[id] = gye[l]/nj;
                    gzj[id] = gze[l]/nj;
                    sj[id] = se[l];
                    vj[id] = pow(dx,3);
                    id++;
                  }
                }
              }
              if( id < npmax-nj || i+j+k == ksum ) {
                mpirange(0,id-1,jsta,jend);
                dird(0,ni-1,jsta,jend);
                mpiallreduce(0,ni-1,gxi);
                mpiallreduce(0,ni-1,gyi);
                mpiallreduce(0,ni-1,gzi);
                for( l=0; l<ni; l++ ) {
                  xo[l] += gxi[l];
                  yo[l] += gyi[l];
                  zo[l] += gzi[l];
                  gxi[l] = gxe[l]/ni;
                  gyi[l] = gye[l]/ni-sv;
                  gzi[l] = gze[l]/ni;
                  vi[l] = pow(dx,3);
                }
                if( 5 <= nequ && nequ <= 7 ) {
                  for( l=0; l<ni; l++ ) {
                    gxi[l] = 0;
                    gyi[l] = 0;
                    gzi[l] = 0;
                  }
                }
                id = 0;
              }
            }
          }
        }
      }
      MPI_Barrier(MPI_COMM_WORLD);
      toc = MPI_Wtime();
      td = toc-tic;
    } else {
      if( myrank == 0 ) {
        fid3.read((char *)(&nbyte),sizeof(int));
        fid3.read((char *)(&td),sizeof(double));
        fid3.read((char *)(&nbyte),sizeof(int));
      }
    }
    if( myrank == 0 ) printf("direct : %g\n",td);

    if( ndir == 1 ) {
      if( myrank == 0 ) {
        nbyte = sizeof(double);
        fid3.write((char *)(&nbyte),sizeof(int));
        fid3.write((char *)(&td),sizeof(double));
        fid3.write((char *)(&nbyte),sizeof(int));
        for ( i=0; i<ni; i++ ) {
          nbyte = sizeof(float)*3;
          fid3.write((char *)(&nbyte),sizeof(int));
          fid3.write((char *)(&xo[i]),sizeof(float));
          fid3.write((char *)(&yo[i]),sizeof(float));
          fid3.write((char *)(&zo[i]),sizeof(float));
          fid3.write((char *)(&nbyte),sizeof(int));
        }
      }
    } else {
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

      for( i=0; i<9; i++ ) tfmm[i] = 0;
      for( i=0; i<ni; i++ ) {
        gxi[i] = gxe[i]/ni;
        gyi[i] = gye[i]/ni-sv;
        gzi[i] = gze[i]/ni;
        vi[i] = pow(dx,3);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      tic = MPI_Wtime();
      mpiprei(ni,mi);
      mpiprej(nj,mj);
      fmm(mi,mj,nequ,tfmm);
      mpiposti(mi,ni);
      mpipostj(mj,nj);
      MPI_Barrier(MPI_COMM_WORLD);
      toc = MPI_Wtime();
      tf = toc-tic;
      tf = 0;
      for ( i=0; i<9; i++ ) tf += tfmm[i];
      MPI_Allgather(tfmm,9,MPI_DOUBLE,tall,9,MPI_DOUBLE,MPI_COMM_WORLD);
      if( myrank == 0 ) {
        printf("fmm    : %g\n",tf);
        for( j=0; j<nprocs; j++ ) {
          for( i=0; i<9; i++ ) {
            nbyte = sizeof(double);
            fid2.write((char *)(&nbyte),sizeof(int)); 
            fid2.write((char *)(&tall[j][i]),sizeof(double));
            fid2.write((char *)(&nbyte),sizeof(int));
          }
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
        fid1.write((char *)(&tf),sizeof(double));
        fid1.write((char *)(&errn),sizeof(double));
        fid1.write((char *)(&nbyte),sizeof(int));
      }
    }
  }

  if( myrank == 0 ) {
    if( ndir != 1 ) {
      fid1.close();
      fid2.close();
    }
    if( ndir <= 1 ) {
      fid3.close();
    }
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
  delete[] tall;
  mem = npmax*25*4;
  memoryfree();
  MPI_Finalize();
};
