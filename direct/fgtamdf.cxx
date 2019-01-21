#include "../misc/constants.h"
extern "C" {
#include "mdgrape3.h"
}

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*sj;
extern float *gxd,*gyd;
extern double *bxmd,(*pos)[3];

extern void memoryuse();
extern void memoryfree();

void dird(int n0, int n1, int n2, int n3) {
  int i,nicall,njcall,icall,iwork1,iwork2,ista,iend,jcall,jwork1,jwork2,jsta,jend,nmd,nmdd;
  M3_UNIT *n_unit;

  gxd = new float [npmax];
  gyd = new float [npmax];
  bxmd = new double [mimax];
  pos = new double [mimax][3];
  mem = npmax*2*4+mimax*4*8;
  memoryuse();

  if( neq==5 ) {
    for( i=n2; i<=n3; i++ ) gyd[i] = gxj[i];
  } else if( neq==6 ) {
    for( i=n2; i<=n3; i++ ) gyd[i] = gyj[i];
  } else if( neq==7 ) {
    for( i=n2; i<=n3; i++ ) gyd[i] = gzj[i];
  }
  for( i=n0; i<=n1; i++ ) {
    gxd[i] = 0;
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
        bxmd[nmd] = gyd[i]/pow(sj[0],3);
      }
      nmd++;

      n_unit = m3_allocate_unit("../../dat/fgta.tblmd3",M3_POTENTIAL,xmins,xmaxs,0);
      m3_set_softening(n_unit,0.002);
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
        
      m3_free_unit(n_unit);

      for( i=ista; i<=iend; i++ ) {
        nmd = i-ista;
        gxd[i] += 0.25/pi*bxmd[nmd];
      }
    }
  }

  if( neq==5 ) {
    for( i=0; i<npmax; i++ ) gxi[i] = 0;
    for( i=n0; i<=n1; i++ ) gxi[i] = gxd[i];
  } else if( neq==6 ) {
    for( i=0; i<npmax; i++ ) gyi[i] = 0;
    for( i=n0; i<=n1; i++ ) gyi[i] = gxd[i];
  } else if( neq==7 ) {
    for( i=0; i<npmax; i++ ) gzi[i] = 0;
    for( i=n0; i<=n1; i++ ) gzi[i] = gxd[i];
  }

  delete[] gxd;
  delete[] gyd;
  delete[] bxmd;
  delete[] pos;
  mem = npmax*2*4+mimax*4*8;
  memoryfree();
}
