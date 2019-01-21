typedef struct {
} Options;

extern void memoryuse();
extern void memoryfree();
extern void dir(int, int, int, int);
extern void fmm(int, int, int, double*);

PetscErrorCode mymatmult(Mat A,Vec X,Vec Y)
{
  Options *options;
  PetscInt ii=0,jc=0,ni,nj;
  PetscScalar x0,y0,z0,ptl;
  PetscScalar x,y,z,vnx,vny,vnz,qint;
  PetscScalar dxij,dyij,dzij,r,cf,Gn;
  PetscScalar ph[6],tfmm[9];
  PetscScalar *xx,*yy;
  PetscErrorCode ierr;
  char *udir,*umpi,*udev,*ufmm,*ugeo,*usmt,*uequ;
  float xmax,ymax,zmax,rdx,rdy,rdz;

  PetscFunctionBegin;
  ierr = MatShellGetContext(A, (void **) &options);CHKERRQ(ierr);
  ierr = VecSet(Y, 0.0);CHKERRQ(ierr);
  ierr = VecGetArray(X, &xx);CHKERRQ(ierr);
  ierr = VecGetArray(Y, &yy);CHKERRQ(ierr);

  ncheck = 0;
  umem = 0;
  nprocs = 1;
  myrank = 0;
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
  mem = npmax*15*4;
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
  neq = nequ;
  stx = 0;
  for(PetscInt i = 0; i < 9; i++) tfmm[i] = 0;

  xmin = 1e6;
  xmax = -1e6;
  ymin = 1e6;
  ymax = -1e6;
  zmin = 1e6;
  zmax = -1e6;
  // Prepare target data
  for(PetscInt l = 0; l < nwn; ++l) {
    xi[l] = ptr[0][l];
    yi[l] = ptr[1][l];
    zi[l] = ptr[2][l];
    xmin = std::min(xi[l],xmin);
    xmax = std::max(xi[l],xmin);
    ymin = std::min(yi[l],ymin);
    ymax = std::max(yi[l],ymin);
    zmin = std::min(zi[l],zmin);
    zmax = std::max(zi[l],zmin);
  }
  ni = nwn;

  // Prepare source data for far field
  // Run over target nodes
  for(PetscInt m = 0; m < nwn; ++m) {
    // Collocation point
    xj[m] = ptr[0][m];
    yj[m] = ptr[1][m];
    zj[m] = ptr[2][m];
    xmin = std::min(xj[m],xmin);
    xmax = std::max(xj[m],xmin);
    ymin = std::min(yj[m],ymin);
    ymax = std::max(yj[m],ymin);
    zmin = std::min(zj[m],zmin);
    zmax = std::max(zj[m],zmin);
    // Normal at collocation point
    gxj[m] = vna[0][m];
    gyj[m] = vna[1][m];
    gzj[m] = vna[2][m];
    // Run over source elements
    ptl = 0;
    for(PetscInt kk = 1; kk <= nne[0][m]; ++kk) {
      PetscInt k = nne[kk][m];
      // Run over quadratures
      for(PetscInt j = 0; j < ngl; ++j) {
        // Quadrature weights
        ph[0] = p1[j][k];
        ph[1] = p2[j][k];
        ph[2] = p3[j][k];
        ph[3] = p4[j][k];
        ph[4] = p5[j][k];
        ph[5] = p6[j][k];
        // Compute the Green's function derivative and apply the triangle quadrature
        for(PetscInt i = 0; i < 6; ++i) {
          if( ntr[i][k] == m ) {
            ptl += 0.5*hsg[j][k]*ph[i];
          }
        }
      }
    }
    vj[m] = ptl*xx[m];
  }
  nj = nwn;
  rdx = (xmax-xmin)*(1+1e-5);
  rdy = (ymax-ymin)*(1+1e-5);
  rdz = (zmax-zmin)*(1+1e-5);
  rd = std::max(rdx,rdy);
  rd = std::max(rd,rdz);

// Stage 1. Run for far field
//  dir(0,ni-1,0,nj-1);
//  for(PetscInt l = 0; l < ni; ++l) yy[l] = vi[l];
  fmm(ni,nj,-2,tfmm);
//  for(PetscInt l = 0; l < ni; ++l) std::cout << l << " " << vi[l] << " " << yy[l] << std::endl;

  jc = 0;
  // Prepare source data for principal integral
  // Run over target nodes
  for(PetscInt k = 0; k < nwe; ++k) {
    // Run over quadratures
    for(PetscInt j = 0; j < ngl; ++j) {
      // Quadrature weights
      ph[0] = p1[j][k];
      ph[1] = p2[j][k];
      ph[2] = p3[j][k];
      ph[3] = p4[j][k];
      ph[4] = p5[j][k];
      ph[5] = p6[j][k];
      x    = 0;
      y    = 0;
      z    = 0;
      vnx  = 0;
      vny  = 0;
      vnz  = 0;
      // Run over node contributions
      for(PetscInt i = 0; i < 6; ++i) {
        ii = ntr[i][k];
        // Quadrature point coordinates
        x += ptr[0][ii]*ph[i];
        y += ptr[1][ii]*ph[i];
        z += ptr[2][ii]*ph[i];
        // Normal at quadrature point
        vnx += vna[0][ii]*ph[i];
        vny += vna[1][ii]*ph[i];
        vnz += vna[2][ii]*ph[i];
      }
      xj[jc] = x;
      yj[jc] = y;
      zj[jc] = z;
      gxj[jc] = vnx;
      gyj[jc] = vny;
      gzj[jc] = vnz;
      xmin = std::min(xj[jc],xmin);
      xmax = std::max(xj[jc],xmin);
      ymin = std::min(yj[jc],ymin);
      ymax = std::max(yj[jc],ymin);
      zmin = std::min(zj[jc],zmin);
      zmax = std::max(zj[jc],zmin);
      vj[jc] = -0.5*hsg[j][k];
      jc++;
    }
  }
  nj = jc;
  rdx = (xmax-xmin)*(1+1e-5);
  rdy = (ymax-ymin)*(1+1e-5);
  rdz = (zmax-zmin)*(1+1e-5);
  rd = std::max(rdx,rdy);
  rd = std::max(rd,rdz);

// Stage 2. Run for principal integral
//  dir(0,ni-1,0,nj-1);
  fmm(ni,nj,-2,tfmm);
  for(PetscInt l = 0; l < ni; ++l) yy[l] += vi[l]*xx[l];

// Stage 3. Run for self contribution
  // Run over target nodes
  for(PetscInt l = 0; l < nwn; ++l) {
    x0 = ptr[0][l];
    y0 = ptr[1][l];
    z0 = ptr[2][l];
    ptl = 0;
    // Run over target elements
    for(PetscInt kk = 1; kk <= nne[0][l]; ++kk) {
      PetscInt k = nne[kk][l];
      // Run over quadratures
      for(PetscInt j = 0; j < ngl; ++j) {
        // Quadrature weights
        ph[0] = p1[j][k];
        ph[1] = p2[j][k];
        ph[2] = p3[j][k];
        ph[3] = p4[j][k];
        ph[4] = p5[j][k];
        ph[5] = p6[j][k];
        x    = 0;
        y    = 0;
        z    = 0;
        vnx  = 0;
        vny  = 0;
        vnz  = 0;
        qint = 0;
        // Run over node contributions
        for(PetscInt i = 0; i < 6; ++i) {
          ii = ntr[i][k];
          // Quadrature point coordinates
          x += ptr[0][ii]*ph[i];
          y += ptr[1][ii]*ph[i];
          z += ptr[2][ii]*ph[i];
          // Normal at quadrature point
          vnx += vna[0][ii]*ph[i];
          vny += vna[1][ii]*ph[i];
          vnz += vna[2][ii]*ph[i];
          if( ii == l ) qint = ph[i];
        }
        // (integral area) * (quadrature weight)
        cf = qint*0.5*hsg[j][k];
        // Compute the Green's function derivative and apply the triangle quadrature
        dxij = x0-x;
        dyij = y0-y;
        dzij = z0-z;
        r = sqrt(dxij*dxij+dyij*dyij+dzij*dzij);
        Gn = (vnx*dxij+vny*dyij+vnz*dzij)/(4*pi*r*r*r);
        ptl += Gn*cf*xx[l];
      }
    }
    yy[l] += ptl-xx[l];
  }

  ierr = VecRestoreArray(X,&xx);CHKERRQ(ierr);
  ierr = VecRestoreArray(Y,&yy);CHKERRQ(ierr);

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
  mem = npmax*15*4;
  memoryfree();
  PetscFunctionReturn(0);
}
