#include "../misc/constants.h"
extern "C" {
#include "mdgrape3.h"
}

extern float *xi,*yi,*zi;
extern int **ndi,*nfj,*nc,**npx,**neij,*nij,*njb;
extern double (*gx)[nspm],(*gy)[nspm],(*gz)[nspm],*xsp,*ysp,*zsp,**p2g;
extern float *gxd,*gyd,*gzd;
extern double *bxmd,*bymd,*bzmd,(*xmd)[3],(*ymd)[3],(*zmd)[3],(*pos)[3];
extern int *ibase,*isize,*jbase,*jsize,*jstamd,*jendmd;

extern void boxc(int, int, int*);

void bsm2p(int nmp, int mp, int lbi, int lbj, int lev, int ipb, double rb) {
  int imd,jmd,ncall,jj,ij,ii,icall,nmd,jb,jx,jy,jz,j,njsize,nmdd,i;
  double rsp,xjc,yjc,zjc;
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
  jendmd[ncall] = lbi-1;
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
          pos[nmd][0] = xjc+xsp[j]*rsp;
          pos[nmd][1] = yjc+ysp[j]*rsp;
          pos[nmd][2] = zjc+zsp[j]*rsp;
          bxmd[nmd] = gx[jb][j];
          bymd[nmd] = gy[jb][j];
          bzmd[nmd] = gz[jb][j];
          nmd++;
        }
      }
    }
    njsize = 1;

    n_unit = m3_allocate_unit("../../dat/bspp.tblmd3",M3_FORCE,xmins,xmaxs,0);
    m3_set_positions(n_unit,pos,nmd);

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

    m3_set_charges(n_unit,bxmd,nmd);
    m3_setup_overlap(n_unit);
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      if( nij[jj] != 0 ) {
        m3_set_cells(n_unit,&cell[jj],njsize);
        m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],xmd+ibase[jj]);
      }
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    m3_set_charges(n_unit,bymd,nmd);
    m3_setup_overlap(n_unit);
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      if( nij[jj] != 0 ) {
        m3_set_cells(n_unit,&cell[jj],njsize);
        m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],ymd+ibase[jj]);
      }
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    m3_set_charges(n_unit,bzmd,nmd);
    m3_setup_overlap(n_unit);
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      if( nij[jj] != 0 ) {
        m3_set_cells(n_unit,&cell[jj],njsize);
        m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],zmd+ibase[jj]);
      }
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    nmdd = 0;
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      for( ij=0; ij<nij[jj]; ij++ ) {
        ii = neij[ij][jj];
        for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
          gxd[i] += 0.25/pi*(ymd[nmdd][2]-zmd[nmdd][1]);
          gyd[i] += 0.25/pi*(zmd[nmdd][0]-xmd[nmdd][2]);
          gzd[i] += 0.25/pi*(xmd[nmdd][1]-ymd[nmdd][0]);
          nmdd++;
        }
      }
    }
    m3_free_unit(n_unit);

  }

}
