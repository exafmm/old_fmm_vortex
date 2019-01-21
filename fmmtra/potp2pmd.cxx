#include "../misc/constants.h"
extern "C" {
#include "mdgrape3.h"
}

extern float *xi,*yi,*zi,*xj,*yj,*zj,*gxj;
extern int **ndj,**ndi,*nej,*nij,**neij;
extern float *gxd,*gyd,*gzd,*vd;
extern double *bxmd,(*xmd)[3],(*pos)[3];
extern int *ibase,*isize,*jbase,*jsize,**istamd,**iendmd,*jstamd,*jendmd,*nicall;

void potp2p(int lbi, int lbj) {
  int imd,jmd,njcall,ic,jj,imdd,ij,ii,jcall,icall,nmd,j,njsize,nmdd,i;
  M3_UNIT *n_unit;
  M3_CELL cell[nbnp];

  imd = 0;
  jmd = 0;
  njcall = 0;
  jstamd[0] = 0;
  ic = 0;
  for( jj=0; jj<lbj; jj++ ) {
    if( nij[jj] != 0 ) {
      imdd = 0;
      ic = 0;
      istamd[0][jj] = 0;
      for( ij=0; ij<nij[jj]; ij++ ) {
        ii = neij[ij][jj];
        imd += ndi[1][ii]-ndi[0][ii]+1;
        imdd += ndi[1][ii]-ndi[0][ii]+1;
        if( imdd > mimax ) {
          iendmd[ic][jj] = ij-1;
          ic++;
          istamd[ic][jj] = ij;
          imdd = ndi[1][ii]-ndi[0][ii]+1;
        }
      }
      iendmd[ic][jj] = nij[jj]-1;
      jmd += std::max(ndj[1][jj]-ndj[0][jj]+1,12);
      if( ic != 0 ) {
        if( jj > jstamd[njcall] ) {
          nicall[njcall] = 1;
          jendmd[njcall] = jj-1;
          njcall++;
          jstamd[njcall] = jj;
        }     
        if( jj != lbj ) {
          nicall[njcall] = ic+1;
          jendmd[njcall] = jj;
          njcall++;
          jstamd[njcall] = jj+1;
          imd = 0;
          jmd = 0;
        }
      } else if ( imd > mimax || jmd > mjmax ) {
        nicall[njcall] = ic+1;
        jendmd[njcall] = jj-1;
        njcall++;
        jstamd[njcall] = jj;
        imd = 0;
        for( ij=0; ij<nij[jj]; ij++ ) {
          ii = neij[ij][jj];
          imd += ndi[1][ii]-ndi[0][ii]+1;
        }
        jmd = std::max(ndj[1][jj]-ndj[0][jj]+1,12);
      }
    }
  }
  nicall[njcall] = ic+1;
  jendmd[njcall] = lbj-1;
  njcall++;

  for( jcall=0; jcall<njcall; jcall++ ) {
    for( icall=0; icall<nicall[jcall]; icall++ ) {

      nmd = 0;
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          cell[jj].base = nmd;
          cell[jj].size = ndj[1][jj]-ndj[0][jj]+1;
          for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
            pos[nmd][0] = xj[j];
            pos[nmd][1] = yj[j];
            pos[nmd][2] = zj[j];
            bxmd[nmd] = gxj[j];
            nmd++;
          }
          if( cell[jj].size < 12 ) {
            for( j=0; j<12-cell[jj].size; j++ ) {
              pos[nmd][0] = 0;
              pos[nmd][1] = 0;
              pos[nmd][2] = 0;
              bxmd[nmd] = 0;
              nmd++;
            }
            cell[jj].size = 12;
          }
        }
      }
      njsize = 1;

      n_unit = m3_allocate_unit("../../dat/pp.tblmd3",M3_POTENTIAL,xmins,xmaxs,0);
      m3_set_positions(n_unit,pos,nmd);

      nmdd = 0;
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          ibase[jj] = nmdd;
          for( ij=istamd[icall][jj]; ij<=iendmd[icall][jj]; ij++ ) {
            ii = neij[ij][jj];
            for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
              pos[nmdd][0] = xi[i];
              pos[nmdd][1] = yi[i];
              pos[nmdd][2] = zi[i];
              nmdd++;
            }
          }
          isize[jj] = nmdd-ibase[jj];
        }
      }

      m3_set_charges(n_unit,bxmd,nmd);
      m3_setup_overlap(n_unit);
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          m3_set_cells(n_unit,&cell[jj],njsize);
          m3_calculate_potentials(n_unit,pos+ibase[jj],isize[jj],bxmd+ibase[jj]);
        }
      }
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);
      nmdd = 0;
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          for( ij=istamd[icall][jj]; ij<=iendmd[icall][jj]; ij++ ) {
            ii = neij[ij][jj];
            for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
              vd[i] += 0.25/pi*bxmd[nmdd];
              nmdd++;
            }
          }
        }
      }

      m3_free_unit(n_unit);

      nmd = 0;
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          cell[jj].base = nmd;
          cell[jj].size = ndj[1][jj]-ndj[0][jj]+1;
          for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
            pos[nmd][0] = xj[j];
            pos[nmd][1] = yj[j];
            pos[nmd][2] = zj[j];
            bxmd[nmd] = gxj[j];
            nmd++;
          }
          if( cell[jj].size < 12 ) {
            for( j=0; j<12-cell[jj].size; j++ ) {
              pos[nmd][0] = 0;
              pos[nmd][1] = 0;
              pos[nmd][2] = 0;
              bxmd[nmd] = 0;
              nmd++;
            }
            cell[jj].size = 12;
          }
        }
      }
      njsize = 1;

      n_unit = m3_allocate_unit("../../dat/force.tblmd3",M3_FORCE,xmins,xmaxs,0);
      m3_set_positions(n_unit,pos,nmd);
      
      nmdd = 0;
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) { 
          ibase[jj] = nmdd;
          for( ij=istamd[icall][jj]; ij<=iendmd[icall][jj]; ij++ ) {
            ii = neij[ij][jj];
            for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
              pos[nmdd][0] = xi[i];
              pos[nmdd][1] = yi[i];
              pos[nmdd][2] = zi[i];
              nmdd++;
            }
          }
          isize[jj] = nmdd-ibase[jj];
        }
      }

      m3_set_charges(n_unit,bxmd,nmd);
      m3_setup_overlap(n_unit);
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          m3_set_cells(n_unit,&cell[jj],njsize);
          m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],xmd+ibase[jj]);
        }
      }
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);

      nmdd = 0;
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          for( ij=istamd[icall][jj]; ij<=iendmd[icall][jj]; ij++ ) {
            ii = neij[ij][jj];
            for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
              gxd[i] += 0.25/pi*xmd[nmdd][0];
              gyd[i] += 0.25/pi*xmd[nmdd][1];
              gzd[i] += 0.25/pi*xmd[nmdd][2];
              nmdd++;
            }
          }
        }
      }
      m3_free_unit(n_unit);

    }
  }

}
