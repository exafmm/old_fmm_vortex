#include "../misc/constants.h"
extern "C" {
#include "mdgrape3.h"
}

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi;
extern int **ndi,*nfi,*nc;
extern double (*px)[nspm],(*py)[nspm],(*pz)[nspm],*xsp,*ysp,*zsp,**p2g;
extern float *gxd,*gyd,*gzd;
extern double *bxmd,*bymd,*bzmd,(*xmd)[3],(*ymd)[3],(*zmd)[3],(*pos)[3],(*pod)[3];
extern int *ibase,*isize,*jbase,*jsize,*jstamd,*jendmd;

extern void boxc(int, int, int*);

void stl2p(int nmp, int mp, int lbi, double rb) {
  int imd,jmd,ncall,ii,icall,nmd,j,k,njsize,nmdd,i;
  double rsp,xic,yic,zic,gxdd,gydd,gzdd,xgj,ygj,zgj,gxp[mjmax],gyp[mjmax],gzp[mjmax];
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
        gydd = 0;
        gzdd = 0;
        for( k=0; k<nmp; k++ ) {
          gxdd += p2g[k][j]*px[ii][k]*rsp;
          gydd += p2g[k][j]*py[ii][k]*rsp;
          gzdd += p2g[k][j]*pz[ii][k]*rsp;
        }
        pod[nmd][0] = xic+xsp[j]*rsp*2;
        pod[nmd][1] = yic+ysp[j]*rsp*2;
        pod[nmd][2] = zic+zsp[j]*rsp*2;
        gxp[nmd] = gxdd;
        gyp[nmd] = gydd;
        gzp[nmd] = gzdd;
        nmd++;
      }
    }
    njsize = 1;

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

    nmd = 0;
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      for( j=0; j<nmp; j++ ) {
        bxmd[nmd] = pod[nmd][1]*gzp[nmd]-pod[nmd][2]*gyp[nmd];
        bymd[nmd] = gyp[nmd];
        bzmd[nmd] = gzp[nmd];
        nmd++;
      }
    }

    n_unit = m3_allocate_unit("../../dat/stpp.tblmd3",M3_FORCE,xmins,xmaxs,0);
    m3_set_positions(n_unit,pod,nmd);

    m3_set_charges(n_unit,bxmd,nmd);
    m3_setup_overlap(n_unit);
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      m3_set_cells(n_unit,&cell[ii],njsize);
      m3_calculate_forces(n_unit,pos+ibase[ii],isize[ii],xmd+ibase[ii]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    m3_set_charges(n_unit,bymd,nmd);
    m3_setup_overlap(n_unit);
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      m3_set_cells(n_unit,&cell[ii],njsize);
      m3_calculate_forces(n_unit,pos+ibase[ii],isize[ii],ymd+ibase[ii]);
    } 
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    m3_set_charges(n_unit,bzmd,nmd);
    m3_setup_overlap(n_unit);
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      m3_set_cells(n_unit,&cell[ii],njsize);
      m3_calculate_forces(n_unit,pos+ibase[ii],isize[ii],zmd+ibase[ii]);
    }   
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    nmdd = 0;
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
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

    nmd = 0;
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      for( j=0; j<nmp; j++ ) {
        bxmd[nmd] = gxp[nmd];
        bymd[nmd] = pod[nmd][2]*gxp[nmd]-pod[nmd][0]*gzp[nmd];
        bzmd[nmd] = gzp[nmd];
        nmd++;
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
    m3_set_charges(n_unit,bymd,nmd);
    m3_setup_overlap(n_unit);
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      m3_set_cells(n_unit,&cell[ii],njsize);
      m3_calculate_forces(n_unit,pos+ibase[ii],isize[ii],ymd+ibase[ii]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    m3_set_charges(n_unit,bzmd,nmd);
    m3_setup_overlap(n_unit);
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      m3_set_cells(n_unit,&cell[ii],njsize);
      m3_calculate_forces(n_unit,pos+ibase[ii],isize[ii],zmd+ibase[ii]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    nmdd = 0;
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
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

    nmd = 0;
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      for( j=0; j<nmp; j++ ) {
        bxmd[nmd] = gxp[nmd];
        bymd[nmd] = gyp[nmd];
        bzmd[nmd] = pod[nmd][0]*gyp[nmd]-pod[nmd][1]*gxp[nmd];
        nmd++;
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
    m3_set_charges(n_unit,bymd,nmd);
    m3_setup_overlap(n_unit);
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      m3_set_cells(n_unit,&cell[ii],njsize);
      m3_calculate_forces(n_unit,pos+ibase[ii],isize[ii],ymd+ibase[ii]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    m3_set_charges(n_unit,bzmd,nmd);
    m3_setup_overlap(n_unit);
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
      m3_set_cells(n_unit,&cell[ii],njsize);
      m3_calculate_forces(n_unit,pos+ibase[ii],isize[ii],zmd+ibase[ii]);
    }
    m3_start_overlap_calculation(n_unit);
    m3_wait_overlap_calculation(n_unit);
    nmdd = 0;
    for( ii=jstamd[icall]; ii<=jendmd[icall]; ii++ ) {
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
    m3_free_unit(n_unit);

  }

}
