#include "mpi.h"
#include "../misc/constants.h"

void mpiallreduce(int n0, int n1, float* var) {
  int ncall,icall,iwork1,iwork2,ista,iend,nbuf,i;
  float send[mpim],recv[mpim];

  ncall = (n1-n0+1)/mpim+1;
  for( icall=0; icall<ncall; icall++ ) {
    iwork1 = (n1-n0+1)/ncall;
    iwork2 = (n1-n0+1)%ncall;
    ista = icall*iwork1+n0+std::min(icall,iwork2);
    iend = ista+iwork1-1;
    if( iwork2 > icall ) iend++;
    nbuf = iend-ista+1;
    for( i=ista; i<=iend; i++ ) send[i-ista] = var[i];
    if( nprocs == 1 ) {
      for( i=0; i<nbuf; i++ ) recv[i] = send[i];
    } else {
      MPI_Allreduce(send,recv,nbuf,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    }
    for( i=ista; i<=iend; i++ ) var[i] = recv[i-ista];
  }
}
