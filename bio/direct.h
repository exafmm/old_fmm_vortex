typedef struct {
} Options;

PetscErrorCode mymatmult(Mat A,Vec X,Vec Y)
{
  Options *options;
  PetscInt ii=0;
  PetscScalar x0,y0,z0,ptl;
  PetscScalar x,y,z,vnx,vny,vnz,qint;
  PetscScalar dxij,dyij,dzij,r,cf,Gn;
  PetscScalar ph[6];
  PetscScalar *xx,*yy;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatShellGetContext(A, (void **) &options);CHKERRQ(ierr);
  ierr = VecSet(Y, 0.0);CHKERRQ(ierr);
  ierr = VecGetArray(X, &xx);CHKERRQ(ierr);
  ierr = VecGetArray(Y, &yy);CHKERRQ(ierr);

  // Run over target nodes
  for(PetscInt l = 0; l < nwn; ++l) {
    x0 = ptr[0][l];
    y0 = ptr[1][l];
    z0 = ptr[2][l];
    ptl = 0;
    // Run over source nodes
    for(PetscInt m = 0; m < nwn; ++m) {
      // Run over source elements
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
          // Collocation point
          x = ptr[0][m];
          y = ptr[1][m];
          z = ptr[2][m];
          // Normal at collocation point
          vnx = vna[0][m];
          vny = vna[1][m];
          vnz = vna[2][m];
          // Compute the Green's function derivative and apply the triangle quadrature
          dxij = x0-x;
          dyij = y0-y;
          dzij = z0-z;
          r = sqrt(dxij*dxij+dyij*dyij+dzij*dzij)+eps;
          Gn = (vnx*dxij+vny*dyij+vnz*dzij)/(4*pi*r*r*r);
          for(PetscInt i = 0; i < 6; ++i) if( ntr[i][k] == m ) ptl += Gn*0.5*hsg[j][k]*ph[i]*xx[m];

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
          // (integral area) * (quadrature weight)
          cf = -0.5*hsg[j][k]/6;
          // Compute the Green's function derivative and apply the triangle quadrature
          dxij = x0-x;
          dyij = y0-y;
          dzij = z0-z;
          r = sqrt(dxij*dxij+dyij*dyij+dzij*dzij);
          Gn = (vnx*dxij+vny*dyij+vnz*dzij)/(4*pi*r*r*r);
          ptl += Gn*cf*xx[l];
        }
      }
    }
    yy[l] = ptl;
  }
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
    yy[l] += ptl;
  }
  for(PetscInt i = 0; i < nwn; i++) {
    yy[i] -= xx[i];
  }

  ierr = VecRestoreArray(X,&xx);CHKERRQ(ierr);
  ierr = VecRestoreArray(Y,&yy);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
