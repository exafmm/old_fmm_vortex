#include "../misc/constants.h"
extern "C" {
#include "mdgrape3.h"
}

extern int *nfi,*nfj,*nlbj,*nc,**npx,**neij,*nij,*njb;
extern double (*px)[nspm],(*gx)[nspm],(*py)[nspm],(*gy)[nspm],(*pz)[nspm],(*gz)[nspm],*xsp,*ysp,*zsp;
extern double *bxmd,*bymd,*bzmd,(*pos)[3];
extern int *ibase,*isize,*jbase,*jsize,*jstamd,*jendmd;

extern void boxc(int, int, int*);

void m2l(int nmp, int mp, int lbi, int lbj, int lev, int ini, double rb) {
  int i,j,imd,jmd,ncall,jj,icall,nmd,jb,njsize,nmdd,ij,ii,ib;
  double rsp,xjc,yjc,zjc,xic,yic,zic;
  M3_UNIT *n_unit;
  M3_CELL cell[nbnp];

  rsp = rb*sqrt(3.0)*0.5;
  if( ini != 0 ) {
    for( i=0; i<ini; i++ ) {
      for( j=0; j<nmp; j++ ) {
        px[i][j] = 0;
        py[i][j] = 0;
        pz[i][j] = 0;
      }
    }
  }

  imd = 0;
  jmd = 0;
  ncall = 0;
  jstamd[0] = 0;
  for( jj=0; jj<lbj; jj++ ) {
    imd += nij[jj]*nmp;
    jmd += nmp;
    if( imd > mimax || jmd > mjmax ) {
      jendmd[ncall] = jj-1;
      ncall++;
      jstamd[ncall] = jj;
      imd = nij[jj]*nmp;
      jmd = nmp;
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
        xjc = (nc[0]+npx[0][jj]*pow(2,lev)+0.5)*rb;
        yjc = (nc[1]+npx[1][jj]*pow(2,lev)+0.5)*rb;
        zjc = (nc[2]+npx[2][jj]*pow(2,lev)+0.5)*rb;
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

    n_unit = m3_allocate_unit("../../dat/ppf.tblmd3",M3_POTENTIAL,xminn,xmaxn,0);
    m3_set_positions(n_unit,pos,nmd);

    nmdd = 0;
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      if( nij[jj] != 0 ) {
        ibase[jj] = nmdd;
        isize[jj] = nij[jj]*nmp;
        for( ij=0; ij<nij[jj]; ij++ ) {
          ii = neij[ij][jj];
          boxc(nfi[ii],3,nc);
          xic = (nc[0]+0.5)*rb;
          yic = (nc[1]+0.5)*rb;
          zic = (nc[2]+0.5)*rb;
          for( j=0; j<nmp; j++ ) {
            pos[nmdd][0] = xic+xsp[j]*rsp;
            pos[nmdd][1] = yic+ysp[j]*rsp;
            pos[nmdd][2] = zic+zsp[j]*rsp;
            nmdd++;
          }
        }
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

    nmdd = 0;
    for( jj=jstamd[icall]; jj<=jendmd[icall]; jj++ ) {
      for( ij=0; ij<nij[jj]; ij++ ) {
        ib = neij[ij][jj];
        for( j=0; j<nmp; j++ ) {
          px[ib][j] += 0.25/pi*bxmd[nmdd];
          py[ib][j] += 0.25/pi*bymd[nmdd];
          pz[ib][j] += 0.25/pi*bzmd[nmdd];
          nmdd++;
        }
      }
    }
    m3_free_unit(n_unit);

  }

  for( jj=0; jj<lbj; jj++ ) {
    jb = njb[jj];
    for( j=0; j<nmp; j++ ) {
      gx[jb][j] = 0;
      gy[jb][j] = 0;
      gz[jb][j] = 0;
    }
  }
}
