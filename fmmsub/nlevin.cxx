#include <fstream>

int main() {
  int nbyte,i,mpi,ndf,nge,nsmt,neq,ienv,nlev[1000000];
  std::fstream fid;

  fid.open("../../dat/nlevel",std::ios::in|std::ios::binary);
  fid.read((char *)(&nbyte),sizeof(int));
  fid.read((char *)(&nlev),sizeof(nlev));
  fid.close();

//  for( i=0; i<1000000; i++ ) nlev[i] = 100000000;

  for( mpi=0; mpi<10; mpi++ ) {
    for( ndf=0; ndf<6; ndf++ ) {
      for( nge=0; nge<10; nge++ ) {
        for( nsmt=0; nsmt<2; nsmt++ ) {
          for( neq=0; neq<10; neq++ ) {
            ienv = 100000*mpi+10000*ndf+1000*nge+100*nsmt+10*neq;
            if( ndf == 0 ) {
              if( neq == 0 ) {
                if( nge == 0 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 70000;
                  nlev[ienv+2] = 700000;
                  nlev[ienv+3] = 5000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 1 ) {
                  nlev[ienv+0] = 15000;
                  nlev[ienv+1] = 150000;
                  nlev[ienv+2] = 3000000;
                  nlev[ienv+3] = 10000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 2 ) {
                  nlev[ienv+0] = 20000;
                  nlev[ienv+1] = 200000;
                  nlev[ienv+2] = 2000000;
                  nlev[ienv+3] = 10000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 3 ) {
                  nlev[ienv+0] = 80000;
                  nlev[ienv+1] = 500000;
                  nlev[ienv+2] = 3000000;
                  nlev[ienv+3] = 10000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 4 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 50000;
                  nlev[ienv+2] = 100000;
                  nlev[ienv+3] = 200000;
                  nlev[ienv+4] = 1000000;
                  nlev[ienv+5] = 5000000;
                } else if( nge == 5 ) {
                  nlev[ienv+0] = 8000;
                  nlev[ienv+1] = 15000;
                  nlev[ienv+2] = 50000;
                  nlev[ienv+3] = 200000;
                  nlev[ienv+4] = 1000000;
                  nlev[ienv+5] = 5000000;
                } else if( nge == 6 ) {
                  nlev[ienv+0] = 8000;
                  nlev[ienv+1] = 15000;
                  nlev[ienv+2] = 50000;
                  nlev[ienv+3] = 200000;
                  nlev[ienv+4] = 1000000;
                  nlev[ienv+5] = 5000000;
                } else if( nge == 7 ) {
                  nlev[ienv+0] = 8000;
                  nlev[ienv+1] = 15000;
                  nlev[ienv+2] = 50000;
                  nlev[ienv+3] = 200000;
                  nlev[ienv+4] = 1000000;
                  nlev[ienv+5] = 5000000;
                } else if( nge == 8 ) {
                  nlev[ienv+0] = 7000;
                  nlev[ienv+1] = 15000;
                  nlev[ienv+2] = 50000;
                  nlev[ienv+3] = 200000;
                  nlev[ienv+4] = 1000000;
                  nlev[ienv+5] = 5000000;
                }
              } else if( 1 <= neq && neq <= 3 ) {
                if( nge == 0 ) {
                  nlev[ienv+0] = 7000;
                  nlev[ienv+1] = 70000;
                  nlev[ienv+2] = 1000000;
                  nlev[ienv+3] = 5000000;
                  nlev[ienv+4] = 10000000;
                  nlev[ienv+5] = 30000000;
                } else if( nge == 1 ) {
                  nlev[ienv+0] = 15000;
                  nlev[ienv+1] = 150000;
                  nlev[ienv+2] = 1500000;
                  nlev[ienv+3] = 15000000;
                  nlev[ienv+4] = 150000000;
                  nlev[ienv+5] = 300000000;
                } else if( nge == 2 ) {
                  nlev[ienv+0] = 15000;
                  nlev[ienv+1] = 150000;
                  nlev[ienv+2] = 1500000;
                  nlev[ienv+3] = 15000000;
                  nlev[ienv+4] = 150000000;
                  nlev[ienv+5] = 300000000;
                } else if( nge == 3 ) {
                  nlev[ienv+0] = 15000;
                  nlev[ienv+1] = 150000;
                  nlev[ienv+2] = 1500000;
                  nlev[ienv+3] = 15000000;
                  nlev[ienv+4] = 150000000;
                  nlev[ienv+5] = 300000000;
                } else if( nge == 4 ) {
                  nlev[ienv+0] = 3000;
                  nlev[ienv+1] = 7000;
                  nlev[ienv+2] = 20000;
                  nlev[ienv+3] = 80000;
                  nlev[ienv+4] = 200000;
                  nlev[ienv+5] = 1000000;
                } else if( nge == 5 ) {
                  nlev[ienv+0] = 6000;
                  nlev[ienv+1] = 10000;
                  nlev[ienv+2] = 20000;
                  nlev[ienv+3] = 80000;
                  nlev[ienv+4] = 200000;
                  nlev[ienv+5] = 1000000;
                } else if( nge == 6 ) {
                  nlev[ienv+0] = 3000;
                  nlev[ienv+1] = 5000;
                  nlev[ienv+2] = 16000;
                  nlev[ienv+3] = 50000;
                  nlev[ienv+4] = 200000;
                  nlev[ienv+5] = 1000000;
                } else if( nge == 7 ) {
                  nlev[ienv+0] = 6000;
                  nlev[ienv+1] = 10000;
                  nlev[ienv+2] = 20000;
                  nlev[ienv+3] = 80000;
                  nlev[ienv+4] = 200000;
                  nlev[ienv+5] = 1000000;
                } else if( nge == 8 ) {
                  nlev[ienv+0] = 3000;
                  nlev[ienv+1] = 20000;
                  nlev[ienv+2] = 40000;
                  nlev[ienv+3] = 80000;
                  nlev[ienv+4] = 200000;
                  nlev[ienv+5] = 1000000;
                }
              } else if( neq == 4 ) {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              } else if( neq == 8 ) {
                nlev[ienv+0] = 5000;
                nlev[ienv+1] = 50000000;
                nlev[ienv+2] = 50000000;
                nlev[ienv+3] = 50000000;
                nlev[ienv+4] = 50000000;
                nlev[ienv+5] = 50000000;
              } else if( neq == 9 ) {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              } else {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              }
            } else if( ndf == 1 ) {
              if( neq == 0 ) {
                if( nge == 0 ) {
                  nlev[ienv+0] = 50000;
                  nlev[ienv+1] = 500000;
                  nlev[ienv+2] = 5000000;
                  nlev[ienv+3] = 50000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 1 ) {
                  nlev[ienv+0] = 80000;
                  nlev[ienv+1] = 500000;
                  nlev[ienv+2] = 3000000;
                  nlev[ienv+3] = 10000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 2 ) {
                  nlev[ienv+0] = 80000;
                  nlev[ienv+1] = 500000;
                  nlev[ienv+2] = 3000000;
                  nlev[ienv+3] = 10000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 3 ) {
                  nlev[ienv+0] = 80000;
                  nlev[ienv+1] = 500000;
                  nlev[ienv+2] = 3000000;
                  nlev[ienv+3] = 10000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 4 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 50000;
                  nlev[ienv+2] = 100000;
                  nlev[ienv+3] = 200000;
                  nlev[ienv+4] = 1000000;
                  nlev[ienv+5] = 5000000;
                } else if( nge == 5 ) {
                  nlev[ienv+0] = 8000;
                  nlev[ienv+1] = 15000;
                  nlev[ienv+2] = 50000;
                  nlev[ienv+3] = 200000;
                  nlev[ienv+4] = 1000000;
                  nlev[ienv+5] = 5000000;
                } else if( nge == 6 ) {
                  nlev[ienv+0] = 8000;
                  nlev[ienv+1] = 15000;
                  nlev[ienv+2] = 50000;
                  nlev[ienv+3] = 200000;
                  nlev[ienv+4] = 1000000;
                  nlev[ienv+5] = 5000000;
                } else if( nge == 7 ) {
                  nlev[ienv+0] = 8000;
                  nlev[ienv+1] = 15000;
                  nlev[ienv+2] = 50000;
                  nlev[ienv+3] = 200000;
                  nlev[ienv+4] = 1000000;
                  nlev[ienv+5] = 5000000;
                } else if( nge == 8 ) {
                  nlev[ienv+0] = 7000;
                  nlev[ienv+1] = 15000;
                  nlev[ienv+2] = 50000;
                  nlev[ienv+3] = 200000;
                  nlev[ienv+4] = 1000000;
                  nlev[ienv+5] = 5000000;
                }
              } else if( 1 <= neq && neq <= 3 ) {
                if( nge == 0 ) {
                  nlev[ienv+0] = 70000;
                  nlev[ienv+1] = 700000;
                  nlev[ienv+2] = 1500000;
                  nlev[ienv+3] = 5000000;
                  nlev[ienv+4] = 10000000;
                  nlev[ienv+5] = 30000000;
                } else if( nge == 1 ) {
                  nlev[ienv+0] = 15000;
                  nlev[ienv+1] = 150000;
                  nlev[ienv+2] = 1500000;
                  nlev[ienv+3] = 15000000;
                  nlev[ienv+4] = 150000000;
                  nlev[ienv+5] = 300000000;
                } else if( nge == 2 ) {
                  nlev[ienv+0] = 15000;
                  nlev[ienv+1] = 150000;
                  nlev[ienv+2] = 1500000;
                  nlev[ienv+3] = 15000000;
                  nlev[ienv+4] = 150000000;
                  nlev[ienv+5] = 300000000;
                } else if( nge == 3 ) {
                  nlev[ienv+0] = 15000;
                  nlev[ienv+1] = 150000;
                  nlev[ienv+2] = 1500000;
                  nlev[ienv+3] = 15000000;
                  nlev[ienv+4] = 150000000;
                  nlev[ienv+5] = 300000000;
                } else if( nge == 4 ) {
                  nlev[ienv+0] = 3000;
                  nlev[ienv+1] = 7000;
                  nlev[ienv+2] = 20000;
                  nlev[ienv+3] = 80000;
                  nlev[ienv+4] = 200000;
                  nlev[ienv+5] = 1000000;
                } else if( nge == 5 ) {
                  nlev[ienv+0] = 6000;
                  nlev[ienv+1] = 10000;
                  nlev[ienv+2] = 20000;
                  nlev[ienv+3] = 80000;
                  nlev[ienv+4] = 200000;
                  nlev[ienv+5] = 1000000;
                } else if( nge == 6 ) {
                  nlev[ienv+0] = 3000;
                  nlev[ienv+1] = 5000;
                  nlev[ienv+2] = 16000;
                  nlev[ienv+3] = 50000;
                  nlev[ienv+4] = 200000;
                  nlev[ienv+5] = 1000000;
                } else if( nge == 7 ) {
                  nlev[ienv+0] = 6000;
                  nlev[ienv+1] = 10000;
                  nlev[ienv+2] = 20000;
                  nlev[ienv+3] = 80000;
                  nlev[ienv+4] = 200000;
                  nlev[ienv+5] = 1000000;
                } else if( nge == 8 ) {
                  nlev[ienv+0] = 3000;
                  nlev[ienv+1] = 20000;
                  nlev[ienv+2] = 40000;
                  nlev[ienv+3] = 80000;
                  nlev[ienv+4] = 200000;
                  nlev[ienv+5] = 1000000;
                }
              } else if( neq == 4 ) {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              } else if( neq == 8 ) {
                nlev[ienv+0] = 5000;
                nlev[ienv+1] = 50000000;
                nlev[ienv+2] = 50000000;
                nlev[ienv+3] = 50000000;
                nlev[ienv+4] = 50000000;
                nlev[ienv+5] = 50000000;
              } else if( neq == 9 ) {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              } else {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              }
            } else if( ndf == 2 ) {
              if( neq == 0 ) {
                if( nge == 0 ) {
                  nlev[ienv+0] = 100000;
                  nlev[ienv+1] = 1500000;
                  nlev[ienv+2] = 20000000;
                  nlev[ienv+3] = 200000000;
                  nlev[ienv+4] = 300000000;
                  nlev[ienv+5] = 400000000;
                } else if( nge == 1 ) {
                  nlev[ienv+0] = 800000;
                  nlev[ienv+1] = 6000000;
                  nlev[ienv+2] = 100000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 2 ) {
                  nlev[ienv+0] = 800000;
                  nlev[ienv+1] = 6000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 3 ) {
                  nlev[ienv+0] = 1000000;
                  nlev[ienv+1] = 6000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 4 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 100000;
                  nlev[ienv+2] = 1000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 5 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 5000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 6 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 5000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 7 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 5000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 8 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 5000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                }
              } else if( 1 <= neq && neq <= 3 ) {
                if( nge == 0 ) {
                  nlev[ienv+0] = 150000;
                  nlev[ienv+1] = 2000000;
                  nlev[ienv+2] = 100000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 1 ) {
                  nlev[ienv+0] = 1000000;
                  nlev[ienv+1] = 2000000;
                  nlev[ienv+2] = 100000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 2 ) {
                  nlev[ienv+0] = 1000000;
                  nlev[ienv+1] = 2000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 3 ) {
                  nlev[ienv+0] = 1000000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 4 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 100000;
                  nlev[ienv+2] = 1000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 5 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 6 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 7 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 8 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                }
              } else if( neq == 4 ) {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              } else if( neq == 8 ) {
                nlev[ienv+0] = 5000;
                nlev[ienv+1] = 50000000;
                nlev[ienv+2] = 50000000;
                nlev[ienv+3] = 50000000;
                nlev[ienv+4] = 50000000;
                nlev[ienv+5] = 50000000;
              } else if( neq == 9 ) {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              } else {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              }
            } else if( ndf == 3 ) {
              if( neq == 0 ) {
                if( nge == 0 ) {
                  nlev[ienv+0] = 300000;
                  nlev[ienv+1] = 1000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 1 ) {
                  nlev[ienv+0] = 800000;
                  nlev[ienv+1] = 6000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 2 ) {
                  nlev[ienv+0] = 800000;
                  nlev[ienv+1] = 6000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 3 ) {
                  nlev[ienv+0] = 50000;
                  nlev[ienv+1] = 5000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 4 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 100000;
                  nlev[ienv+2] = 1000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 5 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 5000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 6 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 5000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 7 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 5000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 8 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 5000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                }
              } else if( 1 <= neq && neq <= 3 ) {
                if( nge == 0 ) {
                  nlev[ienv+0] = 500000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 1 ) {
                  nlev[ienv+0] = 500000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 2 ) {
                  nlev[ienv+0] = 500000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 3 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 4 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 100000;
                  nlev[ienv+2] = 1000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 5 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 6 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 7 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 8 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                }
              } else if( neq == 4 ) {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              } else if( neq == 8 ) {
                nlev[ienv+0] = 5000;
                nlev[ienv+1] = 50000000;
                nlev[ienv+2] = 50000000;
                nlev[ienv+3] = 50000000;
                nlev[ienv+4] = 50000000;
                nlev[ienv+5] = 50000000;
              } else if( neq == 9 ) {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              } else {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              }
            } else if( ndf == 4 ) {
              if( neq == 0 ) {
                if( nge == 0 ) {
                  nlev[ienv+0] = 70000;
                  nlev[ienv+1] = 500000;
                  nlev[ienv+2] = 4000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 200000000;
                  nlev[ienv+5] = 1000000000;
                } else if( nge == 1 ) {
                  nlev[ienv+0] = 70000;
                  nlev[ienv+1] = 500000;
                  nlev[ienv+2] = 4000000;
                  nlev[ienv+3] = 200000000;
                  nlev[ienv+4] = 300000000;
                  nlev[ienv+5] = 400000000;
                } else if( nge == 2 ) {
                  nlev[ienv+0] = 70000;
                  nlev[ienv+1] = 500000;
                  nlev[ienv+2] = 4000000;
                  nlev[ienv+3] = 200000000;
                  nlev[ienv+4] = 300000000;
                  nlev[ienv+5] = 400000000;
                } else if( nge == 3 ) {
                  nlev[ienv+0] = 70000;
                  nlev[ienv+1] = 500000;
                  nlev[ienv+2] = 4000000;
                  nlev[ienv+3] = 200000000;
                  nlev[ienv+4] = 300000000;
                  nlev[ienv+5] = 400000000;
                } else if( nge == 4 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 100000;
                  nlev[ienv+2] = 1000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 5 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 5000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 6 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 5000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 7 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 5000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 8 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 5000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                }
              } else if( 1 <= neq && neq <= 3 ) {
                if( nge == 0 ) {
                  nlev[ienv+0] = 70000;
                  nlev[ienv+1] = 500000;
                  nlev[ienv+2] = 4000000;
                  nlev[ienv+3] = 30000000;
                  nlev[ienv+4] = 200000000;
                  nlev[ienv+5] = 300000000;
                } else if( nge == 1 ) {
                  nlev[ienv+0] = 70000;
                  nlev[ienv+1] = 500000;
                  nlev[ienv+2] = 4000000;
                  nlev[ienv+3] = 30000000;
                  nlev[ienv+4] = 200000000;
                  nlev[ienv+5] = 300000000;
                } else if( nge == 2 ) {
                  nlev[ienv+0] = 70000;
                  nlev[ienv+1] = 500000;
                  nlev[ienv+2] = 4000000;
                  nlev[ienv+3] = 30000000;
                  nlev[ienv+4] = 200000000;
                  nlev[ienv+5] = 300000000;
                } else if( nge == 3 ) {
                  nlev[ienv+0] = 70000;
                  nlev[ienv+1] = 500000;
                  nlev[ienv+2] = 4000000;
                  nlev[ienv+3] = 30000000;
                  nlev[ienv+4] = 200000000;
                  nlev[ienv+5] = 300000000;
                } else if( nge == 4 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 100000;
                  nlev[ienv+2] = 1000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 5 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 6 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 7 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 8 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                }
              } else if( neq == 4 ) {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              } else if( neq == 8 ) {
                nlev[ienv+0] = 5000;
                nlev[ienv+1] = 50000000;
                nlev[ienv+2] = 50000000;
                nlev[ienv+3] = 50000000;
                nlev[ienv+4] = 50000000;
                nlev[ienv+5] = 50000000;
              } else if( neq == 9 ) {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              } else {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              }
            } else if( ndf == 5 ) {
              if( neq == 0 ) {
                if( nge == 0 ) {
                  nlev[ienv+0] = 70000;
                  nlev[ienv+1] = 600000;
                  nlev[ienv+2] = 5000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 1 ) {
                  nlev[ienv+0] = 70000;
                  nlev[ienv+1] = 600000;
                  nlev[ienv+2] = 5000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 2 ) {
                  nlev[ienv+0] = 70000;
                  nlev[ienv+1] = 600000;
                  nlev[ienv+2] = 5000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 3 ) {
                  nlev[ienv+0] = 70000;
                  nlev[ienv+1] = 600000;
                  nlev[ienv+2] = 5000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 4 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 100000;
                  nlev[ienv+2] = 1000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 5 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 5000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 6 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 5000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 7 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 5000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 8 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 5000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                }
              } else if( 1 <= neq && neq <= 3 ) {
                if( nge == 0 ) {
                  nlev[ienv+0] = 100000;
                  nlev[ienv+1] = 700000;
                  nlev[ienv+2] = 5000000;
                  nlev[ienv+3] = 30000000;
                  nlev[ienv+4] = 40000000;
                  nlev[ienv+5] = 50000000;
                } else if( nge == 1 ) {
                  nlev[ienv+0] = 100000;
                  nlev[ienv+1] = 700000;
                  nlev[ienv+2] = 5000000;
                  nlev[ienv+3] = 30000000;
                  nlev[ienv+4] = 40000000;
                  nlev[ienv+5] = 50000000;
                } else if( nge == 2 ) {
                  nlev[ienv+0] = 100000;
                  nlev[ienv+1] = 700000;
                  nlev[ienv+2] = 5000000;
                  nlev[ienv+3] = 30000000;
                  nlev[ienv+4] = 40000000;
                  nlev[ienv+5] = 50000000;
                } else if( nge == 3 ) {
                  nlev[ienv+0] = 100000;
                  nlev[ienv+1] = 700000;
                  nlev[ienv+2] = 5000000;
                  nlev[ienv+3] = 30000000;
                  nlev[ienv+4] = 40000000;
                  nlev[ienv+5] = 50000000;
                } else if( nge == 4 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 100000;
                  nlev[ienv+2] = 1000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 5 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 6 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 7 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                } else if( nge == 8 ) {
                  nlev[ienv+0] = 10000;
                  nlev[ienv+1] = 3000000;
                  nlev[ienv+2] = 10000000;
                  nlev[ienv+3] = 20000000;
                  nlev[ienv+4] = 30000000;
                  nlev[ienv+5] = 40000000;
                }
              } else if( neq == 4 ) {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              } else if( neq == 8 ) {
                nlev[ienv+0] = 5000;
                nlev[ienv+1] = 50000000;
                nlev[ienv+2] = 50000000;
                nlev[ienv+3] = 50000000;
                nlev[ienv+4] = 50000000;
                nlev[ienv+5] = 50000000;
              } else if( neq == 9 ) {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              } else {
                nlev[ienv+0] = 10000;
                nlev[ienv+1] = 100000;
                nlev[ienv+2] = 1000000;
                nlev[ienv+3] = 20000000;
                nlev[ienv+4] = 30000000;
                nlev[ienv+5] = 40000000;
              }
            }
          }
        }
      }
    }
  }
  fid.open("../../dat/nlevel",std::ios::out|std::ios::binary);
  nbyte = sizeof(nlev);
  fid.write((char *)(&nbyte),sizeof(int));
  fid.write((char *)(&nlev),sizeof(nlev));
  fid.write((char *)(&nbyte),sizeof(int));
  fid.close();
}
