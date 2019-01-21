#include "../misc/constants.h"
extern "C" {
#include "mdgrape3.h"
}

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi;
extern int **ndi,*nfj,*nc,**npx,**neij,*nij,*njb;
extern double (*gx)[nspm],(*gy)[nspm],(*gz)[nspm],*xsp,*ysp,*zsp,**p2g;
extern float *gxd,*gyd,*gzd;
extern double *bxmd,*bymd,*bzmd,(*xmd)[3],(*ymd)[3],(*zmd)[3],(*pos)[3],(*pod)[3];
extern int *ibase,*isize,*jbase,*jsize,*jstamd,*jendmd;

extern void boxc(int, int, int*);

void stm2p(int nmp, int mp, int lbi, int lbj, int lev, int ipb, double rb) {
  int imd,jmd,ncall,jj,ij,ii,icall,nmd,jb,jx,jy,jz,j,njsize,nmdd,i;
  double rsp,xjc,yjc,zjc,xgj,ygj,zgj,gxp[mjmax],gyp[mjmax],gzp[mjmax];
  M3_UNIT *n_unit;
  M3_CELL cell[nbnp];

  rsp = rb*sqrt(3.0)*0.5;
  imd = 0;
  jmd = 0;
  ncall = 0;
  jstamd[0] = 0;
  for( jj=0; jj<lbj; jj++ ) {
    if( nij[jj] != 0 ) {
      for( ij=0; ij<nij[jj]; ij++ ) {
        ii = neij[ij][jj];
        imd += ndi[1][ii]-ndi[0][ii]+1;
      }
      jmd += nmp;
      if( imd > mimax || jmd > mjmax ) {
        jendmd[ncall] = jj-1;
        ncall++;
        jstamd[ncall] = jj;
        imd = 0;
        for( ij=0; ij<nij[jj]; ij++ ) {
          ii = neij[ij][jj];
          imd += ndi[1][ii]-ndi[0][ii]+1;
        }
        jmd = nmp;
      }
    }
  }
  jendmd[ncall] = lbj-1;
  ncall++;

  for( icall=0; icall<ncall; icall++ ) {
    nmd = 0;
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      if( nij[jj] != 0 ) {
        jb = njb[jj];
        cell[jj].base = nmd;
        cell[jj].size = nmp;
        boxc(nfj[jj],3,nc);
        jx = nc[0]+npx[0][jj]*int(pow(2,lev));
        jy = nc[1]+npx[1][jj]*int(pow(2,lev));
        jz = nc[2]+npx[2][jj]*int(pow(2,lev));
        xjc = xmin+(jx-0.5+pow(0.5,ipb))*rb;
        yjc = ymin+(jy-0.5+pow(0.5,ipb))*rb;
        zjc = zmin+(jz-0.5+pow(0.5,ipb))*rb;
        for( j=0; j<nmp; j++ ) {
          pod[nmd][0] = xjc+xsp[j]*rsp;
          pod[nmd][1] = yjc+ysp[j]*rsp;
          pod[nmd][2] = zjc+zsp[j]*rsp;
          gxp[nmd] = gx[jb][j];
          gyp[nmd] = gy[jb][j];
          gzp[nmd] = gz[jb][j];
          nmd++;
        }
      }
    }
    njsize = 1;

    nmdd = 0;
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      if( nij[jj] != 0 ) {
        ibase[jj] = nmdd;
        isize[jj] = 0;
        for( ij=0; ij<nij[jj]; ij++ ) {
          ii = neij[ij][jj];
          isize[jj] += ndi[1][ii]-ndi[0][ii]+1;
          for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
            pos[nmdd][0] = xi[i];
            pos[nmdd][1] = yi[i];
            pos[nmdd][2] = zi[i];
            nmdd++;
          }
        }
      }
    }

    nmd = 0;
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      if( nij[jj] != 0 ) {
        for( j=0; j<nmp; j++ ) {
          bxmd[nmd] = pod[nmd][1]*gzp[nmd]-pod[nmd][2]*gyp[nmd];
          bymd[nmd] = gyp[nmd];
          bzmd[nmd] = gzp[nmd];
          nmd++;
        }
      }
    }

    n_unit = m3_allocate_unit("../../dat/stpp.tblmd3",M3_FORCE,xmins,xmaxs,0);
    m3_set_positions(n_unit,pod,nmd);

    m3_set_charges(n_unit,bxmd,nmd);
    m3_setup_overlap(n_unit);
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      m3_set_cells(n_unit,&cell[jj],njsize);
      m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],xmd+ibase[jj]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    m3_set_charges(n_unit,bymd,nmd);
    m3_setup_overlap(n_unit);
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      m3_set_cells(n_unit,&cell[jj],njsize);
      m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],ymd+ibase[jj]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    m3_set_charges(n_unit,bzmd,nmd);
    m3_setup_overlap(n_unit);
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      m3_set_cells(n_unit,&cell[jj],njsize);
      m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],zmd+ibase[jj]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    nmdd = 0;
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      for( ij=0; ij<nij[jj]; ij++ ) {
        ii = neij[ij][jj];
        for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
          xgj = gxi[i]*xmd[nmdd][0]+gyi[i]*xmd[nmdd][1]+gzi[i]*xmd[nmdd][2];
          ygj = gxi[i]*ymd[nmdd][0]+gyi[i]*ymd[nmdd][1]+gzi[i]*ymd[nmdd][2];
          zgj = gxi[i]*zmd[nmdd][0]+gyi[i]*zmd[nmdd][1]+gzi[i]*zmd[nmdd][2];
          gxd[i] += 0.375/pi*(yi[i]*zgj-zi[i]*ygj-xgj);
          gxd[i] += 0.375/pi*gxi[i]*(yi[i]*zmd[nmdd][0]-zi[i]*ymd[nmdd][0]-xmd[nmdd][0]);
          gyd[i] += 0.375/pi*gxi[i]*(yi[i]*zmd[nmdd][1]-zi[i]*ymd[nmdd][1]-xmd[nmdd][1]);
          gzd[i] += 0.375/pi*gxi[i]*(yi[i]*zmd[nmdd][2]-zi[i]*ymd[nmdd][2]-xmd[nmdd][2]);
          nmdd++;
        }
      }
    }

    nmd = 0;
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      if( nij[jj] != 0 ) {
        for( j=0; j<nmp; j++ ) {
          bxmd[nmd] = gxp[nmd];
          bymd[nmd] = pod[nmd][2]*gxp[nmd]-pod[nmd][0]*gzp[nmd];
          bzmd[nmd] = gzp[nmd];
          nmd++;
        }
      }
    }

    m3_set_charges(n_unit,bxmd,nmd);
    m3_setup_overlap(n_unit);
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      m3_set_cells(n_unit,&cell[jj],njsize);
      m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],xmd+ibase[jj]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    m3_set_charges(n_unit,bymd,nmd);
    m3_setup_overlap(n_unit);
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      m3_set_cells(n_unit,&cell[jj],njsize);
      m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],ymd+ibase[jj]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    m3_set_charges(n_unit,bzmd,nmd);
    m3_setup_overlap(n_unit);
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      m3_set_cells(n_unit,&cell[jj],njsize);
      m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],zmd+ibase[jj]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    nmdd = 0;
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      for( ij=0; ij<nij[jj]; ij++ ) {
        ii = neij[ij][jj];
        for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
          xgj = gxi[i]*xmd[nmdd][0]+gyi[i]*xmd[nmdd][1]+gzi[i]*xmd[nmdd][2];
          ygj = gxi[i]*ymd[nmdd][0]+gyi[i]*ymd[nmdd][1]+gzi[i]*ymd[nmdd][2];
          zgj = gxi[i]*zmd[nmdd][0]+gyi[i]*zmd[nmdd][1]+gzi[i]*zmd[nmdd][2];
          gyd[i] += 0.375/pi*(zi[i]*xgj-xi[i]*zgj-ygj);
          gxd[i] += 0.375/pi*gyi[i]*(zi[i]*xmd[nmdd][0]-xi[i]*zmd[nmdd][0]-ymd[nmdd][0]);
          gyd[i] += 0.375/pi*gyi[i]*(zi[i]*xmd[nmdd][1]-xi[i]*zmd[nmdd][1]-ymd[nmdd][1]);
          gzd[i] += 0.375/pi*gyi[i]*(zi[i]*xmd[nmdd][2]-xi[i]*zmd[nmdd][2]-ymd[nmdd][2]);
          nmdd++;
        }
      }
    }

    nmd = 0;
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      if( nij[jj] != 0 ) {
        for( j=0; j<nmp; j++ ) {
          bxmd[nmd] = gxp[nmd];
          bymd[nmd] = gyp[nmd];
          bzmd[nmd] = pod[nmd][0]*gyp[nmd]-pod[nmd][1]*gxp[nmd];
          nmd++;
        }
      }
    }

    m3_set_charges(n_unit,bxmd,nmd);
    m3_setup_overlap(n_unit);
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      m3_set_cells(n_unit,&cell[jj],njsize);
      m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],xmd+ibase[jj]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    m3_set_charges(n_unit,bymd,nmd);
    m3_setup_overlap(n_unit);
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      m3_set_cells(n_unit,&cell[jj],njsize);
      m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],ymd+ibase[jj]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    m3_set_charges(n_unit,bzmd,nmd);
    m3_setup_overlap(n_unit);
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      m3_set_cells(n_unit,&cell[jj],njsize);
      m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],zmd+ibase[jj]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    nmdd = 0;
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      for( ij=0; ij<nij[jj]; ij++ ) {
        ii = neij[ij][jj];
        for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
          xgj = gxi[i]*xmd[nmdd][0]+gyi[i]*xmd[nmdd][1]+gzi[i]*xmd[nmdd][2];
          ygj = gxi[i]*ymd[nmdd][0]+gyi[i]*ymd[nmdd][1]+gzi[i]*ymd[nmdd][2];
          zgj = gxi[i]*zmd[nmdd][0]+gyi[i]*zmd[nmdd][1]+gzi[i]*zmd[nmdd][2];
          gzd[i] += 0.375/pi*(xi[i]*ygj-yi[i]*xgj-zgj);
          gxd[i] += 0.375/pi*gzi[i]*(xi[i]*ymd[nmdd][0]-yi[i]*xmd[nmdd][0]-zmd[nmdd][0]);
          gyd[i] += 0.375/pi*gzi[i]*(xi[i]*ymd[nmdd][1]-yi[i]*xmd[nmdd][1]-zmd[nmdd][1]);
          gzd[i] += 0.375/pi*gzi[i]*(xi[i]*ymd[nmdd][2]-yi[i]*xmd[nmdd][2]-zmd[nmdd][2]);
          nmdd++;
        }
      }
    }

  }

}
