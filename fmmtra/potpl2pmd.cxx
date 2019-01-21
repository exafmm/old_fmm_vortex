#include "../misc/constants.h"
extern "C" {
#include "mdgrape3.h"
}

extern float *xi,*yi,*zi;
extern int **ndi,*nfi,*nc;
extern double (*px)[nspm],(*py)[nspm],(*pz)[nspm],*xsp,*ysp,*zsp,**p2g;
extern float *gxd,*gyd,*gzd,*vd;
extern double *bxmd,(*xmd)[3],(*pos)[3];
extern int *ibase,*isize,*jbase,*jsize,*jstamd,*jendmd;

extern void boxc(int, int, int*);

void potl2p(int nmp, int mp, int lbi, double rb) {
  int imd,jmd,ncall,ii,icall,nmd,j,k,njsize,nmdd,i;
  double rsp,xic,yic,zic,gxdd;
  M3_UNIT *n_unit;
  M3_CELL cell[nbnp];

  rsp = rb*sqrt(3.0)*0.5;
  imd = 0;
  jmd = 0;
  ncall = 0;
  jstamd[0] = 0;
  for( ii=0; ii<lbi; ii++ ) {
    imd += ndi[1][ii]-ndi[0][ii]+1;
    jmd += nmp;
    if( imd > mimax || jmd > mjmax ) {
      jendmd[ncall] = ii-1;
      ncall++;
      jstamd[ncall] = ii;
      imd = ndi[1][ii]-ndi[0][ii]+1;
      jmd = nmp;
    }
  }
  jendmd[ncall] = lbi-1;
  ncall++;

  for( icall=0; icall<ncall; icall++ ) {
    nmd = 0;
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      cell[ii].base = nmd;
      cell[ii].size = nmp;
      boxc(nfi[ii],3,nc);
      xic = xmin+(nc[0]+0.5)*rb;
      yic = ymin+(nc[1]+0.5)*rb;
      zic = zmin+(nc[2]+0.5)*rb;
      for( j=0; j<nmp; j++ ) {
        gxdd = 0;
        for( k=0; k<nmp; k++ ) {
          gxdd += p2g[k][j]*px[k][ii]*rsp;
        }
        pos[nmd][0] = xic+xsp[j]*rsp*2;
        pos[nmd][1] = yic+ysp[j]*rsp*2;
        pos[nmd][2] = zic+zsp[j]*rsp*2;
        bxmd[nmd] = gxdd;
        nmd++;
      }
    }
    njsize = 1;

    n_unit = m3_allocate_unit("../../dat/pp.tblmd3",M3_POTENTIAL,xmins,xmaxs,0);
    m3_set_positions(n_unit,pos,nmd);

    nmdd = 0;
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      ibase[ii] = nmdd;
      isize[ii] = ndi[1][ii]-ndi[0][ii]+1;
      for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
        pos[nmdd][0] = xi[i];
        pos[nmdd][1] = yi[i];
        pos[nmdd][2] = zi[i];
        nmdd++;
      }
    }

    m3_set_charges(n_unit,bxmd,nmd);
    m3_setup_overlap(n_unit);
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      m3_set_cells(n_unit,&cell[ii],njsize);
      m3_calculate_potentials(n_unit,pos+ibase[ii],isize[ii],bxmd+ibase[ii]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    nmdd = 0;
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
        vd[i] += 0.25/pi*bxmd[nmdd];
        nmdd++;
      }
    }
    m3_free_unit(n_unit);

    nmd = 0;
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      cell[ii].base = nmd;
      cell[ii].size = nmp;
      boxc(nfi[ii],3,nc);
      xic = xmin+(nc[0]+0.5)*rb;
      yic = ymin+(nc[1]+0.5)*rb;
      zic = zmin+(nc[2]+0.5)*rb;
      for( j=0; j<nmp; j++ ) {
        gxdd = 0;
        for( k=0; k<nmp; k++ ) {
          gxdd += p2g[k][j]*px[k][ii]*rsp;
        }
        pos[nmd][0] = xic+xsp[j]*rsp*2;
        pos[nmd][1] = yic+ysp[j]*rsp*2;
        pos[nmd][2] = zic+zsp[j]*rsp*2;
        bxmd[nmd] = gxdd;
        nmd++;
      }
    }
    njsize = 1;

    n_unit = m3_allocate_unit("../../dat/force.tblmd3",M3_FORCE,xmins,xmaxs,0);
    m3_set_positions(n_unit,pos,nmd);
  
    nmdd = 0;
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      ibase[ii] = nmdd;
      isize[ii] = ndi[1][ii]-ndi[0][ii]+1;
      for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
        pos[nmdd][0] = xi[i];
        pos[nmdd][1] = yi[i];
        pos[nmdd][2] = zi[i];
        nmdd++;
      }
    }

    m3_set_charges(n_unit,bxmd,nmd); 
    m3_setup_overlap(n_unit);
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      m3_set_cells(n_unit,&cell[ii],njsize);
      m3_calculate_forces(n_unit,pos+ibase[ii],isize[ii],xmd+ibase[ii]);
    } 
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    nmdd = 0;
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      for( i=ndi[0][ii]; i<=ndi[1][ii]; i++ ) {
        gxd[i] += 0.25/pi*xmd[nmdd][0];
        gyd[i] += 0.25/pi*xmd[nmdd][1];
        gzd[i] += 0.25/pi*xmd[nmdd][2];
        nmdd++;
      }
    } 
    m3_free_unit(n_unit);

  }

}
