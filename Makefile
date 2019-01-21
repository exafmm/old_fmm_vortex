.SUFFIXES: .cxx .o

CCS = mpicxx -O3 -g -ffast-math
CCM = $(CCS)
EXS = ./a.out
EXM = mpirun -np 4 ./a.out
Z1_MPI=0
Z2_DEV=0
Z3_FMM=0
Z4_GEO=0
Z5_SMT=0
Z6_EQU=0


ifeq ($(Z1_MPI),0)
CC = $(CCS)
EX = $(EXS)
else
CC = $(CCM)
EX = $(EXM)
endif

FMD = ../fmm
TRA = ../fmmtra
SUB = ../fmmsub
DID = ../direct
SORT = ../sort
MISC = ../misc
MPID = ../mpi

# Box structure
FMM += $(SUB)/alloc.o $(SUB)/boxdatai.o $(SUB)/boxdataj.o\
	$(SUB)/boxdatak.o $(SUB)/boxn.o $(SUB)/boxn1.o $(SUB)/boxc.o\
	$(SUB)/jpbox.o $(SUB)/jsbox.o $(SUB)/jcbox.o $(SUB)/jpsub.o\
	$(SORT)/sortvar.o $(SORT)/sort.o
ifeq ($(Z2_DEV),1)
  FMM += $(SUB)/jibox.o
else
  FMM += $(SUB)/ijbox.o
endif

# FMM
ifeq ($(Z3_FMM),0)
  FMM += $(TRA)/fm2m.o $(TRA)/fm2m2.o $(TRA)/fl2l.o $(TRA)/fl2l2.o\
	$(SUB)/ftrans.o $(SUB)/spharot.o $(SUB)/multipole.o\
	$(SUB)/multipoled.o $(MISC)/factorial.o
# FMM-GPU
  ifeq ($(Z2_DEV),2)
    FMM += $(TRA)/fp2mgp.o $(TRA)/fm2m1gp.o $(TRA)/fm2lgp.o\
	$(TRA)/fl2l1gp.o $(TRA)/bsfl2pgp.o $(TRA)/potfl2pgp.o $(TRA)/potbemfl2p.o $(TRA)/p2pgp.o
    PAR += $(TRA)/p2pgp.o
    ifeq ($(Z6_EQU),1)
      FMM += $(TRA)/stfl2pgp.o
    else
      ifeq ($(Z6_EQU),2)
        FMM += $(TRA)/trfl2pgp.o
      else
        FMM += $(TRA)/mifl2pgp.o
      endif
    endif
# FMM-HOST
  else
    FMM += $(TRA)/fp2m.o $(TRA)/fm2m1.o $(TRA)/fl2l1.o $(TRA)/bsfl2p.o $(TRA)/potfl2p.o\
	$(TRA)/potbemfl2p.o
    ifeq ($(Z2_DEV),1)
      FMM += $(TRA)/fm2lmd.o
    else
      FMM += $(TRA)/fm2l.o
    endif
    ifeq ($(Z6_EQU),1)
      FMM += $(TRA)/stfl2p.o
    else
      ifeq ($(Z6_EQU),2)
        FMM += $(TRA)/trfl2p.o
      else
        FMM += $(TRA)/mifl2p.o
      endif
    endif
  endif
# PPM
else
  ifeq ($(Z3_FMM),1)
    FMM += $(TRA)/pm2m.o $(TRA)/pm2m1.o $(TRA)/pm2m2.o $(TRA)/pl2l.o\
	$(TRA)/pl2l1.o $(TRA)/pl2l2.o $(SUB)/ptrans.o $(MISC)/inv.o
# PPM-HOST
    ifeq ($(Z2_DEV),0)
      FMM += $(TRA)/pp2m.o $(TRA)/pm2l.o $(TRA)/bspl2p.o $(TRA)/potpl2p.o
      ifeq ($(Z6_EQU),1)
        FMM += $(TRA)/stpl2p.o
      else
        ifeq ($(Z6_EQU),2)
          FMM += $(TRA)/trpl2p.o
        else
          FMM += $(TRA)/mipl2p.o
        endif
      endif
    else
# PPM-MDG
      ifeq ($(Z2_DEV),1)
        MDG += $(TRA)/pp2mmd.o $(TRA)/pm2lmd.o $(TRA)/bspl2pmd.o $(TRA)/potpl2pmd.o
        ifeq ($(Z6_EQU),1)
          MDG += $(TRA)/stpl2pmd.o
        else
          ifeq ($(Z6_EQU),2)
            MDG += $(TRA)/trpl2pmd.o
          else
            MDG += $(TRA)/mipl2pmd.o
          endif
        endif
# PPM-GPU
      else
        FMM += $(TRA)/pp2mgp.o $(TRA)/pm2lgp.o $(TRA)/bspl2pgp.o $(TRA)/potpl2pgp.o\
	$(TRA)/p2pgp.o
        PAR += $(TRA)/p2pgp.o
        ifeq ($(Z6_EQU),1)
          FMM += $(TRA)/stpl2pgp.o
        else
          ifeq ($(Z6_EQU),2)
            FMM += $(TRA)/trpl2pgp.o
          else
            FMM += $(TRA)/mipl2pgp.o
          endif
        endif
      endif
    endif
  else
# TREE
    ifeq ($(Z3_FMM),2)
      FMM += $(TRA)/fm2m.o $(TRA)/fm2m2.o\
        $(SUB)/ftrans.o $(SUB)/spharot.o $(SUB)/multipole.o\
        $(SUB)/multipoled.o $(MISC)/factorial.o
# TREE-HOST
      ifeq ($(Z2_DEV),0)
        FMM += $(TRA)/fp2m.o $(TRA)/fm2m1.o $(TRA)/bsfm2p.o $(TRA)/potfm2p.o
        ifeq ($(Z6_EQU),1)
          FMM += $(TRA)/stfm2p.o
        else
          ifeq ($(Z6_EQU),2)
            FMM += $(TRA)/trfm2p.o
          else
            FMM += $(TRA)/mifm2p.o
          endif
        endif
      else
# TREE-MDG
        ifeq ($(Z2_DEV),1)
          FMM += $(TRA)/fp2m.o $(TRA)/fm2m1.o $(TRA)/bsfm2pmd.o $(TRA)/potfm2pmd.o
          ifeq ($(Z6_EQU),1)
            FMM += $(TRA)/stfm2pmd.o
          else
            ifeq ($(Z6_EQU),2)
              FMM += $(TRA)/trfm2pmd.o
            else
              FMM += $(TRA)/mifm2pmd.o
            endif
          endif
# TREE-GPU
        else
          FMM += $(TRA)/fp2mgp.o $(TRA)/fm2m1gp.o\
        $(TRA)/bsfm2pgp.o $(TRA)/potfm2pgp.o $(TRA)/p2pgp.o
          PAR += $(TRA)/p2pgp.o
          ifeq ($(Z6_EQU),1)
            FMM += $(TRA)/stfm2pgp.o
          else
            ifeq ($(Z6_EQU),2)
              FMM += $(TRA)/trfm2pgp.o
            else
              FMM += $(TRA)/mifm2pgp.o
            endif
          endif
        endif
      endif
    else
# TREE-PPM
      ifeq ($(Z3_FMM),3)
        FMM += $(TRA)/pm2m.o $(TRA)/pm2m1.o $(TRA)/pm2m2.o\
        $(SUB)/ptrans.o $(MISC)/inv.o
# TREE-PPM-HOST
        ifeq ($(Z2_DEV),0)
          FMM += $(TRA)/pp2m.o $(TRA)/bspm2p.o $(TRA)/potpm2p.o
          ifeq ($(Z6_EQU),1)
            FMM += $(TRA)/stpm2p.o
          else
            ifeq ($(Z6_EQU),2)
              FMM += $(TRA)/trpm2p.o
            else
              FMM += $(TRA)/mipm2p.o
            endif
          endif
        else
# TREE-PPM-MDG
          ifeq ($(Z2_DEV),1)
            MDG += $(TRA)/pp2mmd.o $(TRA)/bspm2pmd.o $(TRA)/potpm2pmd.o
            ifeq ($(Z6_EQU),1)
              MDG += $(TRA)/stpm2pmd.o
            else
              ifeq ($(Z6_EQU),2)
                MDG += $(TRA)/trpm2pmd.o
              else
                MDG += $(TRA)/mipm2pmd.o
              endif
            endif
# TREE-PPM-GPU
          else
            FMM += $(TRA)/pp2mgp.o $(TRA)/bspm2pgp.o $(TRA)/potpm2pgp.o\
        $(TRA)/p2pgp.o
            PAR += $(TRA)/p2pgp.o
            ifeq ($(Z6_EQU),1)
              FMM += $(TRA)/stpm2pgp.o
            else
              ifeq ($(Z6_EQU),2)
                FMM += $(TRA)/trpm2pgp.o
              else
                FMM += $(TRA)/mipm2pgp.o
              endif
            endif
          endif
        endif
      endif
    endif
  endif
endif
# P2P-HOST
ifeq ($(Z2_DEV),0)
  FMM += $(TRA)/wallp2p.o $(MISC)/ierfc.o
  ifeq ($(Z5_SMT),0)
    FMM += $(TRA)/bsgp2p.o $(TRA)/psegp2p.o $(TRA)/fgtgp2p.o
    PAR += $(TRA)/psegp2p.o
  else
    FMM += $(TRA)/bsap2p.o $(TRA)/pseap2p.o $(TRA)/fgtap2p.o
    PAR += $(TRA)/pseap2p.o
  endif
  ifeq ($(Z6_EQU),3)
    ifeq ($(Z5_SMT),0)
      FMM += $(TRA)/migp2p.o
    else
      FMM += $(TRA)/miap2p.o
    endif
  else
    ifeq ($(Z6_EQU),2)
      ifeq ($(Z5_SMT),0)
        FMM += $(TRA)/trgp2p.o
      else
        FMM += $(TRA)/trap2p.o
      endif
    else
      ifeq ($(Z5_SMT),0)
        FMM += $(TRA)/stgp2p.o
      else
        FMM += $(TRA)/stap2p.o
      endif
    endif
  endif
  FMM += $(TRA)/remp2p.o $(TRA)/potp2p.o $(TRA)/potbemp2p.o
# P2P-MDG
else
  ifeq ($(Z2_DEV),1)
    MDG += $(TRA)/remp2pmd.o $(TRA)/potp2pmd.o $(TRA)/potbemp2p.o
    ifeq ($(Z5_SMT),0)
      MDG += $(TRA)/bsgp2pmd.o $(TRA)/psegp2pmd.o $(TRA)/fgtgp2pmd.o
      PAR += $(TRA)/psegp2pmd.o
    else
      MDG += $(TRA)/bsap2pmd.o $(TRA)/pseap2pmd.o $(TRA)/fgtap2pmd.o
      PAR += $(TRA)/pseap2pmd.o
    endif
    ifeq ($(Z6_EQU),3)
      ifeq ($(Z5_SMT),0)
        MDG += $(TRA)/migp2pmd.o
      else
        MDG += $(TRA)/miap2pmd.o
      endif
    else
      ifeq ($(Z6_EQU),2)
        ifeq ($(Z5_SMT),0)
          MDG += $(TRA)/trgp2pmd.o
        else
          MDG += $(TRA)/trap2pmd.o
        endif
      else
        ifeq ($(Z5_SMT),0)
          MDG += $(TRA)/stgp2pmd.o
        else
          MDG += $(TRA)/stap2pmd.o
        endif
      endif
    endif
    FMM += $(MDG) $(TRA)/wallp2p.o $(MISC)/ierfc.o
# P2P-GPU
  else
    FMM += $(TRA)/wallp2p.o $(MISC)/ierfc.o
    ifeq ($(Z5_SMT),0)
      FMM += $(TRA)/bsgp2pgp.o $(TRA)/psegp2pgp.o $(TRA)/fgtgp2pgp.o
      PAR += $(TRA)/psegp2pgp.o
    else
      FMM += $(TRA)/bsap2pgp.o $(TRA)/pseap2pgp.o $(TRA)/fgtap2pgp.o
      PAR += $(TRA)/pseap2pgp.o
    endif
    ifeq ($(Z6_EQU),3)
      ifeq ($(Z5_SMT),0)
        FMM += $(TRA)/migp2pgp.o
      else
        FMM += $(TRA)/miap2pgp.o
      endif
    else
      ifeq ($(Z6_EQU),2)
        ifeq ($(Z5_SMT),0)
          FMM += $(TRA)/trgp2pgp.o
        else
          FMM += $(TRA)/trap2pgp.o
        endif
      else
        ifeq ($(Z5_SMT),0)
          FMM += $(TRA)/stgp2pgp.o
        else
          FMM += $(TRA)/stap2pgp.o
        endif
      endif
    endif
    FMM += $(TRA)/remp2pgp.o $(TRA)/potp2pgp.o $(TRA)/potbemp2p.o
  endif
endif
ifeq ($(Z1_MPI),0)
  FMM += $(SUB)/nlevel.o $(SUB)/boxallocate.o $(SORT)/sorti.o $(SORT)/sortj.o\
	$(SUB)/boxpart.o
  PAR += $(SUB)/boxallocate.o
else
  MPI = $(MPID)/mpiallreduce.o $(MPID)/mpibcast.o\
	$(MPID)/mpisendrecv.o $(MPID)/mpialltoallv.o\
	$(MPID)/mpinlevel.o $(MPID)/mpiboxallocate.o\
	$(MPID)/mpiprei.o $(MPID)/mpiprej.o\
	$(MPID)/mpiposti.o $(MPID)/mpipostj.o\
	$(SORT)/mpisorti.o $(SORT)/mpisortj.o\
	$(SORT)/mpisortvar.o $(MPID)/mpisendp2p.o\
        $(MPID)/mpiboxl2g.o $(MPID)/mpiboxg2l.o\
        $(SUB)/setvert.o $(SUB)/setedge.o $(MPID)/mpiboxpart.o

  ifeq ($(Z2_DEV),1)
    MPI += $(MPID)/mpijicnt.o
  else
    MPI += $(MPID)/mpiijcnt.o
  endif
  ifeq ($(Z3_FMM),0)
    MPI += $(MPID)/mpireducemf.o
  else
    MPI += $(MPID)/mpireducemp.o
  endif
  PAR += $(MPID)/mpiboxallocate.o
endif

par:
	@$(RM) $(DEF) $(PAR) *.o *.mod
purge:
	@$(RM) $(DEF) $(PAR) *.o *.out *.mod $(FMD)/*.o $(TRA)/*.o $(SUB)/*.o $(DID)/*.o $(SORT)/*.o $(MISC)/*.o $(MPID)/*.o
cleanall:
	@$(RM) body/*.dat body/*.o body/*.out body/*.mod
	@$(RM) channel/*.dat channel/*.o channel/*.out channel/*.mod
	@$(RM) direct/*.dat direct/*.o direct/*.out direct/*.mod
	@$(RM) fmm/*.dat fmm/*.o fmm/*.out fmm/*.mod
	@$(RM) fmmsub/*.dat fmmsub/*.o fmmsub/*.out fmmsub/*.mod
	@$(RM) fmmtra/*.dat fmmtra/*.o fmmtra/*.out fmmtra/*.mod
	@$(RM) isotropic/*.dat isotropic/*.o isotropic/*.out isotropic/*.mod
	@$(RM) misc/*.dat misc/*.o misc/*.out misc/*.mod
	@$(RM) mpi/*.dat mpi/*.o mpi*/out mpi/*.mod
	@$(RM) ring/*.dat ring/*.o ring/*.out ring/*.mod
	@$(RM) shear/*.dat shear/*.o shear/*.out shear/*.mod
	@$(RM) sort/*.dat sort/*.o sort/*.out sort/*.mod
	@$(RM) test/*.dat test/*.o test/*.out test/*.mod
cleansvn:
	@$(RM) -r body/.svn
	@$(RM) -r channel/.svn
	@$(RM) -r direct/.svn
	@$(RM) -r fmm/.svn
	@$(RM) -r fmmsub/.svn
	@$(RM) -r fmmtra/.svn
	@$(RM) -r isotropic/.svn
	@$(RM) -r manual/.svn
	@$(RM) -r md/.svn
	@$(RM) -r misc/.svn
	@$(RM) -r mpi/.svn
	@$(RM) -r ring/.svn
	@$(RM) -r shear/.svn
	@$(RM) -r sort/.svn
	@$(RM) -r test/.svn
	@$(RM) -r .svn

#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules
