#include "../misc/constants.h"
extern "C" {
#include "mdgrape3.h"
}

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*sj;
extern float *gxd,*gyd,*gzd;
extern double *bxmd,*bymd,*bzmd,(*xmd)[3],(*ymd)[3],(*zmd)[3],(*pos)[3];

extern void memoryuse();
extern void memoryfree();

void dird(int n0, int n1, int n2, int n3) {
  int i,nicall,njcall,icall,iwork1,iwork2,ista,iend,jcall,jwork1,jwork2,jsta,jend,nmd,nmdd;
  double xgj,ygj,zgj;
  M3_UNIT *n_unit;

  gxd = new float [npmax];
  gyd = new float [npmax];
  gzd = new float [npmax];
  bxmd = new double [mimax];
  bymd = new double [mimax];
  bzmd = new double [mimax];
  xmd = new double [mimax][3];
  ymd = new double [mimax][3];
  zmd = new double [mimax][3];
  pos = new double [mimax][3];
  mem = npmax*3*4+mimax*15*8;
  memoryuse();

  for( i=n0; i<=n1; i++ ) {
    gxd[i] = 0;
    gyd[i] = 0;
    gzd[i] = 0;
  }

  nicall = (n1-n0+1)/mimax+1;
  njcall = (n3-n2+1)/mjmax+1;
  for( icall=0; icall<nicall; icall++ ) {
    iwork1 = (n1-n0+1)/nicall;
    iwork2 = (n1-n0+1)%nicall;
    ista = icall*iwork1+n0+std::min(icall,iwork2);
    iend = ista+iwork1-1;
    if( iwork2 > icall ) iend++;

    for( jcall=0; jcall<njcall; jcall++ ) {
      jwork1 = (n3-n2+1)/njcall;
      jwork2 = (n3-n2+1)%njcall;
      jsta = jcall*jwork1+n2+std::min(jcall,jwork2);
      jend = jsta+jwork1;
      if( jwork2 > jcall ) jend++;

      for( i=jsta; i<jend; i++ ) {
        nmd = i-jsta;
        pos[nmd][0] = xj[i]/sj[0];
        pos[nmd][1] = yj[i]/sj[0];
        pos[nmd][2] = zj[i]/sj[0];
        bxmd[nmd] = gxj[i]/pow(sj[0],3);
        bymd[nmd] = gyj[i]/pow(sj[0],3);
        bzmd[nmd] = gzj[i]/pow(sj[0],3);
      }
      nmd++;

      n_unit = m3_allocate_unit("../../dat/bsaf.tblmd3",M3_POTENTIAL,xminn,xmaxn,0);
      m3_set_positions(n_unit,pos,nmd);
 
      for( i=ista; i<=iend; i++ ) {
        nmdd = i-ista;
        pos[nmdd][0] = xi[i]/sj[0];
        pos[nmdd][1] = yi[i]/sj[0];
        pos[nmdd][2] = zi[i]/sj[0];
      }
      nmdd++;

      m3_set_charges(n_unit,bxmd,nmd);
      m3_setup_overlap(n_unit);
      m3_calculate_potentials(n_unit,pos,nmdd,bxmd);
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);
      m3_set_charges(n_unit,bymd,nmd);
      m3_setup_overlap(n_unit);
      m3_calculate_potentials(n_unit,pos,nmdd,bymd);
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);
      m3_set_charges(n_unit,bzmd,nmd);
      m3_setup_overlap(n_unit);
      m3_calculate_potentials(n_unit,pos,nmdd,bzmd);
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);
        
      m3_free_unit(n_unit);

      for( i=ista; i<=iend; i++ ) {
        nmd = i-ista;
        gxd[i] += 0.25/pi*(gyi[i]*bzmd[nmd]-gzi[i]*bymd[nmd]);
        gyd[i] += 0.25/pi*(gzi[i]*bxmd[nmd]-gxi[i]*bzmd[nmd]);
        gzd[i] += 0.25/pi*(gxi[i]*bymd[nmd]-gyi[i]*bxmd[nmd]);
      }

      for( i=jsta; i<jend; i++ ) {
        nmd = i-jsta;
        pos[nmd][0] = xj[i]/sj[0];
        pos[nmd][1] = yj[i]/sj[0];
        pos[nmd][2] = zj[i]/sj[0];
        bxmd[nmd] = (yj[i]*gzj[i]-zj[i]*gyj[i])/pow(sj[0],4);
        bymd[nmd] = gyj[i]/pow(sj[0],4);
        bzmd[nmd] = gzj[i]/pow(sj[0],4);
      }
      nmd++;

      n_unit = m3_allocate_unit("../../dat/staf.tblmd3",M3_FORCE,xminn,xmaxn,0);
      m3_set_positions(n_unit,pos,nmd);

      for( i=ista; i<=iend; i++ ) {
        nmdd = i-ista;
        pos[nmdd][0] = xi[i]/sj[0];
        pos[nmdd][1] = yi[i]/sj[0];
        pos[nmdd][2] = zi[i]/sj[0];
      }
      nmdd++;

      m3_set_charges(n_unit,bxmd,nmd);
      m3_setup_overlap(n_unit);
      m3_calculate_forces(n_unit,pos,nmdd,xmd);
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);
      m3_set_charges(n_unit,bymd,nmd);
      m3_setup_overlap(n_unit);
      m3_calculate_forces(n_unit,pos,nmdd,ymd);
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);
      m3_set_charges(n_unit,bzmd,nmd);
      m3_setup_overlap(n_unit);
      m3_calculate_forces(n_unit,pos,nmdd,zmd);
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);

      for( i=ista; i<=iend; i++ ) {
        nmd = i-ista;
        xgj = gxi[i]*xmd[nmd][0]+gyi[i]*xmd[nmd][1]+gzi[i]*xmd[nmd][2];
        ygj = gxi[i]*ymd[nmd][0]+gyi[i]*ymd[nmd][1]+gzi[i]*ymd[nmd][2];
        zgj = gxi[i]*zmd[nmd][0]+gyi[i]*zmd[nmd][1]+gzi[i]*zmd[nmd][2];
        gxd[i] += 0.75/pi*(yi[i]*zgj-zi[i]*ygj-xgj);
      }

      for( i=jsta; i<jend; i++ ) {
        nmd = i-jsta;
        bxmd[nmd] = gxj[i]/pow(sj[0],4);
        bymd[nmd] = (zj[i]*gxj[i]-xj[i]*gzj[i])/pow(sj[0],4);
        bzmd[nmd] = gzj[i]/pow(sj[0],4);
      }
      nmd++;

      m3_set_charges(n_unit,bxmd,nmd);
      m3_setup_overlap(n_unit);
      m3_calculate_forces(n_unit,pos,nmdd,xmd);
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);
      m3_set_charges(n_unit,bymd,nmd);
      m3_setup_overlap(n_unit);
      m3_calculate_forces(n_unit,pos,nmdd,ymd);
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);
      m3_set_charges(n_unit,bzmd,nmd);
      m3_setup_overlap(n_unit);
      m3_calculate_forces(n_unit,pos,nmdd,zmd);
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);

      for( i=ista; i<=iend; i++ ) {
        nmd = i-ista;
        xgj = gxi[i]*xmd[nmd][0]+gyi[i]*xmd[nmd][1]+gzi[i]*xmd[nmd][2];
        ygj = gxi[i]*ymd[nmd][0]+gyi[i]*ymd[nmd][1]+gzi[i]*ymd[nmd][2];
        zgj = gxi[i]*zmd[nmd][0]+gyi[i]*zmd[nmd][1]+gzi[i]*zmd[nmd][2];
        gyd[i] += 0.75/pi*(zi[i]*xgj-xi[i]*zgj-ygj);
      }

      for( i=jsta; i<jend; i++ ) {
        nmd = i-jsta;
        bxmd[nmd] = gxj[i]/pow(sj[0],4);
        bymd[nmd] = gyj[i]/pow(sj[0],4);
        bzmd[nmd] = (xj[i]*gyj[i]-yj[i]*gxj[i])/pow(sj[0],4);
      }
      nmd++;
  
      m3_set_charges(n_unit,bxmd,nmd);
      m3_setup_overlap(n_unit);
      m3_calculate_forces(n_unit,pos,nmdd,xmd);
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);
      m3_set_charges(n_unit,bymd,nmd);
      m3_setup_overlap(n_unit);
      m3_calculate_forces(n_unit,pos,nmdd,ymd);
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);
      m3_set_charges(n_unit,bzmd,nmd);
      m3_setup_overlap(n_unit);
      m3_calculate_forces(n_unit,pos,nmdd,zmd);
      m3_start_overlap_calculation(n_unit);
      m3_wait_overlap_calculation(n_unit);

      for( i=ista; i<=iend; i++ ) {
        nmd = i-ista;
        xgj = gxi[i]*xmd[nmd][0]+gyi[i]*xmd[nmd][1]+gzi[i]*xmd[nmd][2];
        ygj = gxi[i]*ymd[nmd][0]+gyi[i]*ymd[nmd][1]+gzi[i]*ymd[nmd][2];
        zgj = gxi[i]*zmd[nmd][0]+gyi[i]*zmd[nmd][1]+gzi[i]*zmd[nmd][2];
        gzd[i] += 0.75/pi*(xi[i]*ygj-yi[i]*xgj-zgj);
      }

      m3_free_unit(n_unit);
    }
  }

  for( i=0; i<npmax; i++ ) {
    gxi[i] = 0;
    gyi[i] = 0;
    gzi[i] = 0;
  }
  for( i=n0; i<=n1; i++ ) {
    gxi[i] = gxd[i];
    gyi[i] = gyd[i];
    gzi[i] = gzd[i];
  }

  delete[] gxd;
  delete[] gyd;
  delete[] gzd;
  delete[] bxmd;
  delete[] bymd;
  delete[] bzmd;
  delete[] xmd;
  delete[] ymd;
  delete[] zmd;
  delete[] pos;
  mem = npmax*3*4+mimax*15*8;
  memoryfree();
}
