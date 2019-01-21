#include "../misc/constants.h"
extern "C" {
#include "mdgrape3.h"
}

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi,*xj,*yj,*zj,*gxj,*gyj,*gzj,*vj,*sj;
extern int **ndj,**ndi,*nej,*nij,**neij;
extern float *gxd,*gyd,*gzd;
extern double *bxmd,*bymd,*bzmd,(*xmd)[3],(*ymd)[3],(*zmd)[3],(*pos)[3];
extern int *ibase,*isize,*jbase,*jsize,**istamd,**iendmd,*jstamd,*jendmd,*nicall;

void stp2p(int lbi, int lbj) {
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
            pos[nmd][0] = xj[j]/sj[0];
            pos[nmd][1] = yj[j]/sj[0];
            pos[nmd][2] = zj[j]/sj[0];
            bxmd[nmd] = gxj[j]/pow(sj[0],3);
            bymd[nmd] = gyj[j]/pow(sj[0],3);
            bzmd[nmd] = gzj[j]/pow(sj[0],3);
            nmd++;
          }
          if( cell[jj].size < 12 ) {
            for( j=0; j<12-cell[jj].size; j++ ) {
              pos[nmd][0] = 0;
              pos[nmd][1] = 0;
              pos[nmd][2] = 0;
              bxmd[nmd] = 0;
              bymd[nmd] = 0;
              bzmd[nmd] = 0;
              nmd++;
            }
            cell[jj].size = 12;
          }
        }
      }
      njsize = 1;

      n_unit = m3_allocate_unit("../../dat/bsg.tblmd3",M3_POTENTIAL,xminn,xmaxn,0);
      m3_set_positions(n_unit,pos,nmd);

      nmdd = 0;
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          ibase[jj] = nmdd;
          for( ij=istamd[icall][jj]; ij<=iendmd[icall][jj]; ij++ ) {
            ii = neij[ij][jj];
            for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
              pos[nmdd][0] = xi[i]/sj[0];
              pos[nmdd][1] = yi[i]/sj[0];
              pos[nmdd][2] = zi[i]/sj[0];
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
      m3_set_charges(n_unit,bymd,nmd);
      m3_setup_overlap(n_unit);
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          m3_set_cells(n_unit,&cell[jj],njsize);
          m3_calculate_potentials(n_unit,pos+ibase[jj],isize[jj],bymd+ibase[jj]);
        }
      }
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);
      m3_set_charges(n_unit,bzmd,nmd);
      m3_setup_overlap(n_unit);
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          m3_set_cells(n_unit,&cell[jj],njsize);
          m3_calculate_potentials(n_unit,pos+ibase[jj],isize[jj],bzmd+ibase[jj]);
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
              gxd[i] -= 0.25/pi*(gyi[i]*bzmd[nmdd]-gzi[i]*bymd[nmdd]);
              gyd[i] -= 0.25/pi*(gzi[i]*bxmd[nmdd]-gxi[i]*bzmd[nmdd]);
              gzd[i] -= 0.25/pi*(gxi[i]*bymd[nmdd]-gyi[i]*bxmd[nmdd]);
              nmdd++;
            }
          }
        }
      }

      m3_free_unit(n_unit);

      nmd = 0;
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          cell[jj].size = ndj[1][jj]-ndj[0][jj]+1;
          for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
            pos[nmd][0] = xj[j]/sj[0];
            pos[nmd][1] = yj[j]/sj[0];
            pos[nmd][2] = zj[j]/sj[0];
            nmd++;
          }
          if( cell[jj].size < 12 ) {
            for( j=0; j<12-cell[jj].size; j++ ) {
              pos[nmd][0] = 0;
              pos[nmd][1] = 0;
              pos[nmd][2] = 0;
              nmd++;
            }
            cell[jj].size = 12;
          }
        }
      }

      n_unit = m3_allocate_unit("../../dat/stg.tblmd3",M3_FORCE,xminn,xmaxn,0);
      m3_set_positions(n_unit,pos,nmd);

      nmdd = 0;
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          for( ij=istamd[icall][jj]; ij<=iendmd[icall][jj]; ij++ ) {
            ii = neij[ij][jj];
            for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
              pos[nmdd][0] = xi[i]/sj[0];
              pos[nmdd][1] = yi[i]/sj[0];
              pos[nmdd][2] = zi[i]/sj[0];
              nmdd++;
            }
          }
        }
      }

      nmd = 0;
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          cell[jj].size = ndj[1][jj]-ndj[0][jj]+1;
          for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
            bxmd[nmd] = (yj[j]*gzj[j]-zj[j]*gyj[j])/pow(sj[0],4);
            bymd[nmd] = gyj[j]/pow(sj[0],4);
            bzmd[nmd] = gzj[j]/pow(sj[0],4);
            nmd++;
          }
          if( cell[jj].size < 12 ) {
            for( j=0; j<12-cell[jj].size; j++ ) {
              bxmd[nmd] = 0;
              bymd[nmd] = 0;
              bzmd[nmd] = 0;
              nmd++;
            }
            cell[jj].size = 12;
          }
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
      m3_set_charges(n_unit,bymd,nmd);
      m3_setup_overlap(n_unit);
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          m3_set_cells(n_unit,&cell[jj],njsize);
          m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],ymd+ibase[jj]);
        }
      }
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);
      m3_set_charges(n_unit,bzmd,nmd);
      m3_setup_overlap(n_unit);
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          m3_set_cells(n_unit,&cell[jj],njsize);
          m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],zmd+ibase[jj]);
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
              gxd[i] += 0.75/pi*gxi[i]*(yi[i]*zmd[nmdd][0]-zi[i]*ymd[nmdd][0]-xmd[nmdd][0]);
              gyd[i] += 0.75/pi*gxi[i]*(yi[i]*zmd[nmdd][1]-zi[i]*ymd[nmdd][1]-xmd[nmdd][1]);
              gzd[i] += 0.75/pi*gxi[i]*(yi[i]*zmd[nmdd][2]-zi[i]*ymd[nmdd][2]-xmd[nmdd][2]);
              nmdd++;
            }
          }
        }
      }

      nmd = 0;
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          cell[jj].size = ndj[1][jj]-ndj[0][jj]+1;
          for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
            bxmd[nmd] = gxj[j]/pow(sj[0],4);
            bymd[nmd] = (zj[j]*gxj[j]-xj[j]*gzj[j])/pow(sj[0],4);
            bzmd[nmd] = gzj[j]/pow(sj[0],4);
            nmd++;
          }
          if( cell[jj].size < 12 ) {
            for( j=0; j<12-cell[jj].size; j++ ) {
              bxmd[nmd] = 0;
              bymd[nmd] = 0;
              bzmd[nmd] = 0;
              nmd++;
            } 
            cell[jj].size = 12;
          }   
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
      m3_set_charges(n_unit,bymd,nmd);
      m3_setup_overlap(n_unit);
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          m3_set_cells(n_unit,&cell[jj],njsize);
          m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],ymd+ibase[jj]);
        }   
      }       
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);
      m3_set_charges(n_unit,bzmd,nmd);
      m3_setup_overlap(n_unit);
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          m3_set_cells(n_unit,&cell[jj],njsize);
          m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],zmd+ibase[jj]);
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
              gxd[i] += 0.75/pi*gyi[i]*(zi[i]*xmd[nmdd][0]-xi[i]*zmd[nmdd][0]-ymd[nmdd][0]);
              gyd[i] += 0.75/pi*gyi[i]*(zi[i]*xmd[nmdd][1]-xi[i]*zmd[nmdd][1]-ymd[nmdd][1]);
              gzd[i] += 0.75/pi*gyi[i]*(zi[i]*xmd[nmdd][2]-xi[i]*zmd[nmdd][2]-ymd[nmdd][2]);
              nmdd++;
            }
          }
        }
      }

      nmd = 0;
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          cell[jj].size = ndj[1][jj]-ndj[0][jj]+1;
          for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
            bxmd[nmd] = gxj[j]/pow(sj[0],4);
            bymd[nmd] = gyj[j]/pow(sj[0],4);
            bzmd[nmd] = (xj[j]*gyj[j]-yj[j]*gxj[j])/pow(sj[0],4);
            nmd++;
          }
          if( cell[jj].size < 12 ) {
            for( j=0; j<12-cell[jj].size; j++ ) {
              bxmd[nmd] = 0;
              bymd[nmd] = 0;
              bzmd[nmd] = 0;
              nmd++;
            } 
            cell[jj].size = 12;
          }   
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
      m3_set_charges(n_unit,bymd,nmd);
      m3_setup_overlap(n_unit);
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          m3_set_cells(n_unit,&cell[jj],njsize);
          m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],ymd+ibase[jj]);
        }   
      }       
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);
      m3_set_charges(n_unit,bzmd,nmd);
      m3_setup_overlap(n_unit);
      for( jj=jstamd[jcall]; jj<=jendmd[jcall]; jj++ ) {
        if( nij[jj] != 0 ) {
          m3_set_cells(n_unit,&cell[jj],njsize);
          m3_calculate_forces(n_unit,pos+ibase[jj],isize[jj],zmd+ibase[jj]);
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
              gxd[i] += 0.75/pi*gzi[i]*(xi[i]*ymd[nmdd][0]-yi[i]*xmd[nmdd][0]-zmd[nmdd][0]);
              gyd[i] += 0.75/pi*gzi[i]*(xi[i]*ymd[nmdd][1]-yi[i]*xmd[nmdd][1]-zmd[nmdd][1]);
              gzd[i] += 0.75/pi*gzi[i]*(xi[i]*ymd[nmdd][2]-yi[i]*xmd[nmdd][2]-zmd[nmdd][2]);
              nmdd++;
            }
          }
        }
      }
      m3_free_unit(n_unit);

    }
  }

}
