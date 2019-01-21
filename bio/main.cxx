#define MAIN
#include "parameters.h"
#include "../misc/constants.h"
#include "../misc/arrays.h"
#include "fmm.h"
#undef MAIN

extern void memoryuse();
extern void memoryfree();
extern void geometry(int, int&, int, double);
extern void gauss_trgl();
extern void elm_geom();
extern void velocity(int, int, int, int, int);

int main(int argc, char *argv[]) {
  char *udir,*umpi,*udev,*ufmm,*ugeo,*usmt,*uequ;
  Options options;
  PetscInt i,np,ng,npo,nwc,nen,ijk,j,k,nbyte,ng3,ij;
  PetscInt nij[6*nwmax];
  PetscScalar dxg=0,dyg=0,dzg=0,xgmin=0,ygmin=0,zgmin=0,uwd,vwd,wwd;
  std::fstream fid;
  Mat A;
  Vec b,x;
  KSP ksp;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &argv, (char *) 0, (char *) 0);CHKERRQ(ierr);
  umem = 0;
  xi = new float [npmax];
  yi = new float [npmax];
  zi = new float [npmax];
  xw = new float [npmax];
  yw = new float [npmax];
  zw = new float [npmax];
  uw = new float [npmax];
  vw = new float [npmax];
  ww = new float [npmax];
  si = new double [nwmax];
  gt = new double [nwmax];
  gs = new double [nwmax];
  un = new double [nwmax];
  mem = npmax*13*4;
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

  np = 200;
  ng = 21;
  if( nge == 5 ) {
    dxg = 2;
    dyg = 2;
    dzg = 2;
    xgmin = -1;
    ygmin = -1;
    zgmin = -1;
  } else if( nge == 6 ) {
    dxg = 0.09;
    dyg = 0.09;
    dzg = 0.09;
    xgmin = 0.09;
    ygmin = -0.45;
    zgmin = -0.45;
  }

  npo = np;

  geometry(ndv,nwc,nge,rwb);

  gauss_trgl();
  elm_geom();

  for( i=0; i<nwn; i++ ) {
    xw[i] = ptr[0][i];
    yw[i] = ptr[1][i];
    zw[i] = ptr[2][i];

    un[i] = uin*vna[0][i]+vin*vna[1][i]+win*vna[2][i];
    si[i] = 0;
  }
  nen = 6;

// Start here
  int idx[nwmax];
  for( i=0; i<nwn; i++ ) idx[i] = i;
  MatCreateShell(PETSC_COMM_WORLD,nwn,nwn,PETSC_DECIDE,PETSC_DECIDE,&options,&A);
  MatShellSetOperation(A,MATOP_MULT, (void (*)(void)) mymatmult);
  VecCreate(PETSC_COMM_WORLD,&x);
  VecSetSizes(x,PETSC_DECIDE,nwn);
  VecSetFromOptions(x);
  VecDuplicate(x,&b);
  VecSetValues(x,nwn,idx,si,INSERT_VALUES);
  VecAssemblyBegin(x);
  VecAssemblyEnd(x);
  VecSetValues(b,nwn,idx,un,INSERT_VALUES);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);
  KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
  KSPSetFromOptions(ksp);
  KSPSolve(ksp,b,x);
  VecGetValues(x,nwn,idx,si);
  KSPDestroy(ksp);
  VecDestroy(x);
  VecDestroy(b);
  MatDestroy(A);
// End here

  ijk = 0;
  for( i=0; i<ng; i++ ) {
    for( j=0; j<ng; j++ ) {
      for( k=0; k<ng; k++ ) {
        xi[ijk] = xgmin+dxg*i/(ng-1);
        yi[ijk] = ygmin+dyg*j/(ng-1);
        zi[ijk] = zgmin+dzg*k/(ng-1);
        ijk++;
      }
    }
  }
  velocity(0,nwe,1,ng*ng*ng,nge);

  fid.open("bio.dat",std::ios::out | std::ios::binary);

  ng3 = ng*ng*ng;
  nbyte = sizeof(int)*7;
  fid.write((char *)(&nbyte),sizeof(int)); 
  fid.write((char *)(&nwe),sizeof(int));
  fid.write((char *)(&nwn),sizeof(int));
  fid.write((char *)(&nen),sizeof(int));
  fid.write((char *)(&np),sizeof(int));
  fid.write((char *)(&npo),sizeof(int));
  fid.write((char *)(&ng3),sizeof(int));
  fid.write((char *)(&nge),sizeof(int));
  fid.write((char *)(&nbyte),sizeof(int)); 
  nbyte = sizeof(float)*nwn*3+sizeof(double)*nwn*9;
  fid.write((char *)(&nbyte),sizeof(int)); 
  for( i=0; i<nwn; i++ ) fid.write((char *)(&xw[i]),sizeof(float));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&yw[i]),sizeof(float));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&zw[i]),sizeof(float));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vna[0][i]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vna[1][i]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vna[2][i]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vta[0][i]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vta[1][i]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vta[2][i]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vsa[0][i]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vsa[1][i]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&vsa[2][i]),sizeof(double));
  fid.write((char *)(&nbyte),sizeof(int));
  nbyte = sizeof(float)*ng3*6;
  fid.write((char *)(&nbyte),sizeof(int)); 
  for( i=0; i<ng3; i++ ) fid.write((char *)(&xi[i]),sizeof(float));
  for( i=0; i<ng3; i++ ) fid.write((char *)(&yi[i]),sizeof(float));
  for( i=0; i<ng3; i++ ) fid.write((char *)(&zi[i]),sizeof(float));
  for( i=0; i<ng3; i++ ) {
    uwd = uw[i]+uin;
    fid.write((char *)(&uwd),sizeof(double));
  }
  for( i=0; i<ng3; i++ ) {
    vwd = vw[i]+vin;
    fid.write((char *)(&vwd),sizeof(double));
  }
  for( i=0; i<ng3; i++ ) {
    wwd = ww[i]+win;
    fid.write((char *)(&wwd),sizeof(double));
  }
  fid.write((char *)(&nbyte),sizeof(int)); 

  for( i=0; i<nwe; i++ ) {
    for( j=0; j<6; j++ ) {
      ij = i*6+j;
      nij[ij] = ntr[j][i];
    }
  }
  nwn = nwe*6;
  nbyte = sizeof(double)*nwn*3;
  fid.write((char *)(&nbyte),sizeof(int));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&ptr[0][nij[i]]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&ptr[1][nij[i]]),sizeof(double));
  for( i=0; i<nwn; i++ ) fid.write((char *)(&ptr[2][nij[i]]),sizeof(double));
  fid.write((char *)(&nbyte),sizeof(int));

  fid.close();
  ierr = PetscFinalize();CHKERRQ(ierr);

  delete[] xi;
  delete[] yi;
  delete[] zi;
  delete[] xw;
  delete[] yw;
  delete[] zw;
  delete[] uw;
  delete[] vw;
  delete[] ww;
  delete[] si;
  delete[] gt;
  delete[] gs;
  delete[] un;
  mem = npmax*13*4;
  memoryfree();

  PetscFunctionReturn(0);
}
