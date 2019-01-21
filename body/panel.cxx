#define MAIN
#include "parameters.h"
#include "../misc/constants.h"
#include "../misc/arrays.h"
#include "matmult.h"
#undef MAIN

extern void memoryuse();
extern void geometry(int, int&, int, double);
extern void gauss_trgl();
extern void elm_geom();
extern void dimhij(int, double(*)[nwmax]);
extern void bicgstab(double*, double(*)[nwmax], double*, int);
extern void bicgstab2(double*, double(*)[nwmax], double*, int, int, int, int);
extern void velocity(int, int, int, int, int);
extern void intrude(int&, int&, int, int);
extern void savedata(int, int, int);

int main(int argc, char *argv[]) {
  char *udir,*umpi,*udev,*ufmm,*ugeo,*usmt,*uequ;
  Options options;
  PetscInt i,np,ng,npo,nwc,nen,ijk,j,k,nbyte,ng3;
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
  ut = new double [nwmax];
  us = new double [nwmax];
  dimh = new double [nwmax][nwmax];
  ust = new double [nwmax];
  gts = new double [nwmax];
  mem = npmax*20*4+nwmax*8*8+nwmax*nwmax*8;
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

  fid.open("panel.dat",std::ios::out | std::ios::binary);

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

  for( i=0; i<npmax; i++ ) {
    xj[i] = (double) rand()/RAND_MAX*dxg+xgmin;
    yj[i] = (double) rand()/RAND_MAX*dyg+ygmin;
    zj[i] = (double) rand()/RAND_MAX*dzg+zgmin;
    xo[i] = xj[i];
    yo[i] = yj[i];
    zo[i] = zj[i];
    gxj[i] = (double) rand()/RAND_MAX/np;
    gyj[i] = (double) rand()/RAND_MAX/np;
    gzj[i] = (double) rand()/RAND_MAX/np;
    vj[i] = pow(0.3,3);
    sj[i] = 0.3;
  }
  npo = np;

  geometry(ndv,nwc,nge,rwb);

  gauss_trgl();
  elm_geom();
  dimhij(nsv,dimh);

  for( i=0; i<nwn; i++ ) {
    xw[i] = ptr[0][i];
    yw[i] = ptr[1][i];
    zw[i] = ptr[2][i];

    un[i] = uin*vna[0][i]+vin*vna[1][i]+win*vna[2][i];
    ut[i] = uin*vta[0][i]+vin*vta[1][i]+win*vta[2][i];
    us[i] = uin*vsa[0][i]+vin*vsa[1][i]+win*vsa[2][i];
    ust[i] = us[i];
    ust[i+nwn] = -ut[i];
    si[i] = 0;
    gts[i] = 0;
    gts[i+nwn] = 0;
  }
  nen = 6;

// Start here
  int idx[nwmax];
  options.nsv = nsv;
  if( nsv == 1 ) {
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
  } else {
    for( i=0; i<2*nwn; i++ ) idx[i] = i;
    MatCreateShell(PETSC_COMM_WORLD,2*nwn,2*nwn,PETSC_DECIDE,PETSC_DECIDE,&options,&A);
    MatShellSetOperation(A,MATOP_MULT, (void (*)(void)) mymatmult);
    VecCreate(PETSC_COMM_WORLD,&x);
    VecSetSizes(x,PETSC_DECIDE,2*nwn);
    VecSetFromOptions(x);
    VecDuplicate(x,&b);
    VecSetValues(x,2*nwn,idx,gts,INSERT_VALUES);
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    VecSetValues(b,2*nwn,idx,ust,INSERT_VALUES);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
  }
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);
  KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
  KSPSetFromOptions(ksp);
  KSPSolve(ksp,b,x);
  if( nsv == 1 ) {
    VecGetValues(x,nwn,idx,si);
  } else {
    VecGetValues(x,2*nwn,idx,gts);
  }
  KSPDestroy(ksp);
  VecDestroy(x);
  VecDestroy(b);
  MatDestroy(A);
// End here

  for( i=0; i<nwn; i++ ) {
    si[i] *= 3;
    gt[i] = gts[i];
    gs[i] = gts[i+nwn];
  }

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
  velocity(0,nwe,nsv,ng*ng*ng,nge);

  intrude(np,npo,nge,2);

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
  nbyte = sizeof(float)*npo*3;
  fid.write((char *)(&nbyte),sizeof(int)); 
  for( i=0; i<npo; i++ ) fid.write((char *)(&xo[i]),sizeof(float));
  for( i=0; i<npo; i++ ) fid.write((char *)(&yo[i]),sizeof(float));
  for( i=0; i<npo; i++ ) fid.write((char *)(&zo[i]),sizeof(float));
  fid.write((char *)(&nbyte),sizeof(int)); 

  savedata(nwe,np,nge);

  fid.close();
  ierr = PetscFinalize();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
