typedef struct {
  PetscInt nsv;
} Options;

PetscErrorCode mymatmult(Mat A,Vec X,Vec Y)
{
  Options *options;
  PetscInt nsv,i0,i1,i2,i3,i4,i5,i6,ii,itest;
  PetscInt im[nwmax];
  PetscScalar x0,y0,z0,ptl,ptl1,ptl2,ptl3,ptl4;
  PetscScalar x,y,z,vnx,vny,vnz,vtx,vty,vtz,vsx,vsy,vsz,qint;
  PetscScalar dxij,dyij,dzij,r,den,cf,Gn,Gtt,Gts,Gst,Gss;
  PetscScalar ph[6];
  PetscScalar *xx,*yy;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatShellGetContext(A, (void **) &options);CHKERRQ(ierr);
  ierr = VecSet(Y, 0.0);CHKERRQ(ierr);
  ierr = VecGetArray(X, &xx);CHKERRQ(ierr);
  ierr = VecGetArray(Y, &yy);CHKERRQ(ierr);
  nsv = options->nsv;

  for(PetscInt i = 0; i < nwn; ++i) im[i] = 0;
  // Run over source nodes
  for(PetscInt m = 0; m < nwn; ++m) {
    // Impulse
    im[m] = 1.0;
    // Run over target nodes
    for(PetscInt l = 0; l < nwn; ++l) {
      x0 = ptr[0][l];
      y0 = ptr[1][l];
      z0 = ptr[2][l];
      i0 = im[l];
      ptl = 0;
      ptl1 = 0;
      ptl2 = 0;
      ptl3 = 0;
      ptl4 = 0;
      // Run over source elements
      for(PetscInt k = 0; k < nwe; ++k) {
        i1 = ntr[0][k];
        i2 = ntr[1][k];
        i3 = ntr[2][k];
        i4 = ntr[3][k];
        i5 = ntr[4][k];
        i6 = ntr[5][k];
        // Does this source element contain the current source node?
        itest = im[i1]+im[i2]+im[i3]+im[i4]+im[i5]+im[i6]+i0;
        if( itest != 0 ) {
          // Run over quadratures
          for(PetscInt j = 0; j < ngl; ++j) {
            x    = 0;
            y    = 0;
            z    = 0;
            vnx  = 0;
            vny  = 0;
            vnz  = 0;
            vtx  = 0;
            vty  = 0;
            vtz  = 0;
            vsx  = 0;
            vsy  = 0;
            vsz  = 0;
            qint = 0;

            // This could be calculated on the fly from a reference element
            ph[0] = p1[j][k];
            ph[1] = p2[j][k];
            ph[2] = p3[j][k];
            ph[3] = p4[j][k];
            ph[4] = p5[j][k];
            ph[5] = p6[j][k];
            // Run over node contributions
            for(PetscInt i = 0; i < 6; ++i) {
              ii = ntr[i][k];
              // Quad point coordinates
              x += ptr[0][ii]*ph[i];
              y += ptr[1][ii]*ph[i];
              z += ptr[2][ii]*ph[i];
              // Normal at quad point
              vnx += vna[0][ii]*ph[i];
              vny += vna[1][ii]*ph[i];
              vnz += vna[2][ii]*ph[i];
              // Tangent at quad point
              vtx += vta[0][ii]*ph[i];
              vty += vta[1][ii]*ph[i];
              vtz += vta[2][ii]*ph[i];
              // Binormal at quad point
              vsx += vsa[0][ii]*ph[i];
              vsy += vsa[1][ii]*ph[i];
              vsz += vsa[2][ii]*ph[i];
              // Basis func value for active source point
              qint += im[ii]*ph[i];
            }

            // Compute the Green's function derivative and apply the triangle quadrature
            dxij = x0-x;
            dyij = y0-y;
            dzij = z0-z;
            r = sqrt(dxij*dxij+dyij*dyij+dzij*dzij);
            den = 4*pi*r*r*r;
            // Quad weight
            cf = 0.5*hsg[j][k];

            if (nsv == 1) {
              Gn = vnx*dxij+vny*dyij+vnz*dzij;
              // This is the Green function contribution at the quad point, and then subtracting out the singular part
              ptl += (qint - i0)*Gn/den*cf;
            } else if (nsv == 2) {
              Gtt = vta[0][l]*(vty*dzij-vtz*dyij)+vta[1][l]*(vtz*dxij-vtx*dzij)+vta[2][l]*(vtx*dyij-vty*dxij);
              Gts = vta[0][l]*(vsy*dzij-vsz*dyij)+vta[1][l]*(vsz*dxij-vsx*dzij)+vta[2][l]*(vsx*dyij-vsy*dxij);
              Gst = vsa[0][l]*(vty*dzij-vtz*dyij)+vsa[1][l]*(vtz*dxij-vtx*dzij)+vsa[2][l]*(vtx*dyij-vty*dxij);
              Gss = vsa[0][l]*(vsy*dzij-vsz*dyij)+vsa[1][l]*(vsz*dxij-vsx*dzij)+vsa[2][l]*(vsx*dyij-vsy*dxij);
              ptl1 += (qint-i0)*Gtt/den*cf;
              ptl2 += (qint-i0)*Gts/den*cf;
              ptl3 += (qint-i0)*Gst/den*cf;
              ptl4 += (qint-i0)*Gss/den*cf;
            }

          }
        }
      }
      if (nsv == 1) {
        yy[l] += (ptl - 0.5*i0)*xx[m];
      } else if (nsv == 2) {
        yy[l]     += ( ptl3 - 0.5*i0)*xx[m];
        yy[l]     += (-ptl4 + 0.5*i0)*xx[m+nwn];
        yy[l+nwn] += ( ptl1 - 0.5*i0)*xx[m];
        yy[l+nwn] += (-ptl2 + 0.5*i0)*xx[m+nwn];
      }
    }
    im[m] = 0.0; // Reset impulse
  }
  if (nsv == 1) {
    for(PetscInt i = 0; i < nwn; i++) {
      yy[i] -= 0.5*xx[i];
    }
  } else if (nsv == 2) {
    for(PetscInt i = nwn; i < 2*nwn; i++) {
      yy[i] -= 0.5*xx[i];
    }
  }

  ierr = VecRestoreArray(X,&xx);CHKERRQ(ierr);
  ierr = VecRestoreArray(Y,&yy);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

