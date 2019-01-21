#include "../misc/constants.h"
extern "C" {
#include "mdgrape3.h"
}

extern float *xj,*yj,*zj,*gxj,*gyj,*gzj;
extern int **ndj,*nfj,*nc;
extern double (*gx)[nspm],(*gy)[nspm],(*gz)[nspm],*xsp,*ysp,*zsp,**p2g;
extern double *bxmd,*bymd,*bzmd,(*xmd)[3],(*ymd)[3],(*zmd)[3],(*pos)[3];
extern int *ibase,*isize,*jbase,*jsize,*jstamd,*jendmd;

extern void boxc(int, int, int*);

void p2m(int nmp, int mp, int lbj, double rb, int nlm) {
  int jmd,ncall,jj,icall,nmd,j,njsize,nmdd,k;
  double rsp,xjc,yjc,zjc,gxdd,gydd,gzdd;
  M3_UNIT *n_unit;
  M3_CELL cell[nbnp];

  rsp = rb*sqrt(3.0)*0.5;

  jmd = 0;
  ncall = 0;
  jstamd[0] = 0;
  for( jj=0; jj<lbj; jj++ ) {
    jmd += ndj[1][jj]-ndj[0][jj]+1;
    if( jmd > mjmax ) {
      jendmd[ncall] = jj-1;
      ncall++;
      jstamd[ncall] = jj;
      jmd = ndj[1][jj]-ndj[0][jj]+1;
    }
  }
  jendmd[ncall] = lbj-1;
  ncall++;

  for( icall=0; icall<ncall; icall++ ) {

    nmd = 0;
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      cell[jj].base = nmd;
      cell[jj].size = ndj[1][jj]-ndj[0][jj]+1;
      for( j=ndj[0][jj]; j<=ndj[1][jj]; j++ ) {
        pos[nmd][0] = xj[j];
        pos[nmd][1] = yj[j];
        pos[nmd][2] = zj[j];
        bxmd[nmd] = gxj[j];
        bymd[nmd] = gyj[j];
        bzmd[nmd] = gzj[j];
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
    njsize = 1;

    n_unit = m3_allocate_unit("../../dat/pp.tblmd3",M3_POTENTIAL,xmins,xmaxs,0);
    m3_set_positions(n_unit,pos,nmd);

    nmdd = 0;
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      ibase[jj] = (jj-jstamd[icall])*nmp;
      isize[jj] = nmp;
      boxc(nfj[jj],3,nc);
      xjc = xmin+(nc[0]+0.5)*rb;
      yjc = ymin+(nc[1]+0.5)*rb;
      zjc = zmin+(nc[2]+0.5)*rb;
      for( j=0; j<nmp; j++ ) {
        pos[nmdd][0] = xjc+xsp[j]*rsp*2;
        pos[nmdd][1] = yjc+ysp[j]*rsp*2;
        pos[nmdd][2] = zjc+zsp[j]*rsp*2;
        nmdd++;
      }
    }

    m3_set_charges(n_unit,bxmd,nmd);
    m3_setup_overlap(n_unit);
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      m3_set_cells(n_unit,&cell[jj],njsize);
      m3_calculate_potentials(n_unit,pos+ibase[jj],isize[jj],bxmd+ibase[jj]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    m3_set_charges(n_unit,bymd,nmd);
    m3_setup_overlap(n_unit);
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      m3_set_cells(n_unit,&cell[jj],njsize);
      m3_calculate_potentials(n_unit,pos+ibase[jj],isize[jj],bymd+ibase[jj]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    m3_set_charges(n_unit,bzmd,nmd);
    m3_setup_overlap(n_unit);
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      m3_set_cells(n_unit,&cell[jj],njsize);
      m3_calculate_potentials(n_unit,pos+ibase[jj],isize[jj],bzmd+ibase[jj]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);

    nmd = 0;
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      for( j=0; j<nmp; j++ ) {
        gxdd = 0;
        gydd = 0;
        gzdd = 0;
        for( k=0; k<nmp; k++ ) {
          nmd = (jj-jstamd[icall])*nmp+k;
          gxdd += p2g[k][j]*0.25/pi*bxmd[nmd]*rsp;
          gydd += p2g[k][j]*0.25/pi*bymd[nmd]*rsp;
          gzdd += p2g[k][j]*0.25/pi*bzmd[nmd]*rsp;
        }
        gx[jj+nlm][j] = gxdd;
        gy[jj+nlm][j] = gydd;
        gz[jj+nlm][j] = gzdd;
      }
    }
    m3_free_unit(n_unit);

  }

}
