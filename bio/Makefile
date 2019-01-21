include ../Makefile

ifeq ($(Z6_EQU),-2)
  DIR = $(DID)/potbem.o $(DID)/potbemd.o
else
  ifeq ($(Z6_EQU),-1)
    DIR = $(DID)/pot.o $(DID)/potd.o
  else
    DIR = $(DID)/bsg.o $(DID)/bsgd.o
  endif
endif
DIR += $(MISC)/memory.o
FMO = $(FMD)/fmmf.o
ifeq ($(Z1_MPI),0)
  FMO += $(SUB)/ifbox.o
else
  MPI += $(MPID)/mpiifbox.o
  MPI += $(MPID)/mpisendm2lf.o
  MPI += $(MPID)/mpirange.o
endif

BEM = $(MISC)/geometry.o $(MISC)/bem3d.o $(MISC)/velocity.o $(MISC)/memory.o

ifeq ($(Z1_MPI),0)
  OBJN = $(DEF) ntest.o $(DIR) $(FMO) $(FMM)
else
  OBJN = $(DEF) mpintest.o $(DIR) $(FMO) $(FMM) $(MPI)
endif
OBJSP = $(DEF) main.o $(BEM) $(DIR) $(FMO) $(FMM)

n: $(OBJN)
	$(CLINKER) $? $(PETSC_SNES_LIB) $(CCL) $(MDGL)

m: $(OBJSP)
	$(CLINKER) $? $(PETSC_SNES_LIB) $(MDGL)
	@$(RM) *.o

mpintest.o: mpintest.cxx
	$(CCM) -c $?
$(MPI): $(MPI:.o=.cxx)
	$(CCM) $(MDGI) -c $*.cxx -o $@
