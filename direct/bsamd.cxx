#include "../misc/constants.h"
extern "C" {
#include "mdgrape3.h"
}

extern float *xi,*yi,*zi,*gxi,*gyi,*gzi;
extern float *xj,*yj,*zj,*gxj,*gyj,*gzj,*sj;
extern double *bxmd,*bymd,*bzmd,(*xmd)[3],(*ymd)[3],(*zmd)[3],(*pos)[3];

extern void memoryuse();
extern void memoryfree();

void dir(int n0, int n1, int n2, int n3) {
  int i,nicall,njcall,icall,iwork1,iwork2,ista,iend,jcall,jwork1,jwork2,jsta,jend,nmd,nmdd;
  M3_UNIT *n_unit;

  bxmd = new double [mimax];
  bymd = new double [mimax];
  bzmd = new double [mimax];
  xmd = new double [mimax][3];
  ymd = new double [mimax][3];
  zmd = new double [mimax][3];
  pos = new double [mimax][3];
  mem = mimax*15*8;
  memoryuse();

  for( i=0; i<npmax; i++ ) {
    gxi[i] = 0;
    gyi[i] = 0;
    gzi[i] = 0;
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
        bxmd[nmd] = gxj[i]/sj[0]/sj[0];
        bymd[nmd] = gyj[i]/sj[0]/sj[0];
        bzmd[nmd] = gzj[i]/sj[0]/sj[0];
      }
      nmd++;

      n_unit = m3_allocate_unit("../../dat/bsa.tblmd3",M3_FORCE,xminn,xmaxn,0);
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
        
      m3_free_unit(n_unit);

      for( i=ista; i<=iend; i++ ) {
        nmd = i-ista;
        gxi[i] += 0.25/pi*(ymd[nmd][2]-zmd[nmd][1]);
        gyi[i] += 0.25/pi*(zmd[nmd][0]-xmd[nmd][2]);
        gzi[i] += 0.25/pi*(xmd[nmd][1]-ymd[nmd][0]);
      }
    }
  }

  delete[] bxmd;
  delete[] bymd;
  delete[] bzmd;
  delete[] xmd;
  delete[] ymd;
  delete[] zmd;
  delete[] pos;
  mem = mimax*15*8;
  memoryfree();
}
