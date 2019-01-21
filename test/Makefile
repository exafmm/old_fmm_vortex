include ../Makefile

ifeq ($(Z6_EQU),0)
  ifeq ($(Z2_DEV),1)
    ifeq ($(Z5_SMT),0)
      MDIR = $(DID)/bsgd.o $(DID)/bsgmd.o
      DIR = $(DID)/bsgmd.o $(DID)/bsgmdf.o
    else
      MDIR = $(DID)/bsad.o $(DID)/bsamd.o
      DIR = $(DID)/bsamd.o $(DID)/bsamdf.o
    endif
  else
    ifeq ($(Z2_DEV),2)
      ifeq ($(Z5_SMT),0)
        MDIR = $(DID)/bsgd.o $(DID)/bsggp.o $(DID)/dirgp.o
        DIR = $(DID)/bsggp.o $(DID)/bsggpd.o $(DID)/dirgp.o
      else
        MDIR = $(DID)/bsad.o $(DID)/bsagp.o $(DID)/dirgp.o
        DIR = $(DID)/bsagp.o $(DID)/bsagpd.o $(DID)/dirgp.o
      endif
    else
      ifeq ($(Z5_SMT),0)
        DIR = $(DID)/bsg.o $(DID)/bsgd.o
      else
        DIR = $(DID)/bsa.o $(DID)/bsad.o
      endif
    endif
  endif
else

  ifeq ($(Z6_EQU),1)
    ifeq ($(Z2_DEV),1)
      ifeq ($(Z5_SMT),0)
        MDIR = $(DID)/stgd.o $(DID)/stgmd.o
        DIR = $(DID)/stgmd.o $(DID)/stgmdf.o
      else
        MDIR = $(DID)/stad.o $(DID)/stamd.o
        DIR = $(DID)/stamd.o $(DID)/stamdf.o
      endif
    else
      ifeq ($(Z2_DEV),2)
        ifeq ($(Z5_SMT),0)
          MDIR = $(DID)/stgd.o $(DID)/stggp.o $(DID)/dirgp.o
          DIR = $(DID)/stggp.o $(DID)/stggpd.o $(DID)/dirgp.o
        else
          MDIR = $(DID)/stad.o $(DID)/stagp.o $(DID)/dirgp.o
          DIR = $(DID)/stagp.o $(DID)/stagpd.o $(DID)/dirgp.o
        endif
      else
        ifeq ($(Z5_SMT),0)
          DIR = $(DID)/stg.o $(DID)/stgd.o
        else
          DIR = $(DID)/sta.o $(DID)/stad.o
        endif
      endif
    endif
  else

    ifeq ($(Z6_EQU),2)
      ifeq ($(Z2_DEV),1)
        ifeq ($(Z5_SMT),0)
          MDIR = $(DID)/trgd.o $(DID)/trgmd.o
          DIR = $(DID)/trgmd.o $(DID)/trgmdf.o
        else
          MDIR = $(DID)/trad.o $(DID)/tramd.o
          DIR = $(DID)/tramd.o $(DID)/tramdf.o
        endif
      else
        ifeq ($(Z2_DEV),2)
          ifeq ($(Z5_SMT),0)
            MDIR = $(DID)/trgd.o $(DID)/trggp.o $(DID)/dirgp.o
            DIR = $(DID)/trggp.o $(DID)/trggpd.o $(DID)/dirgp.o
          else
            MDIR = $(DID)/trad.o $(DID)/tragp.o $(DID)/dirgp.o
            DIR = $(DID)/tragp.o $(DID)/tragpd.o $(DID)/dirgp.o
          endif
        else
          ifeq ($(Z5_SMT),0)
            DIR = $(DID)/trg.o $(DID)/trgd.o
          else
            DIR = $(DID)/tra.o $(DID)/trad.o
          endif
        endif
      endif
    else

      ifeq ($(Z6_EQU),3)
        ifeq ($(Z2_DEV),1)
          ifeq ($(Z5_SMT),0)
            MDIR = $(DID)/migd.o $(DID)/migmd.o
            DIR = $(DID)/migmd.o $(DID)/migmdf.o
          else
            MDIR = $(DID)/miad.o $(DID)/miamd.o
            DIR = $(DID)/miamd.o $(DID)/miamdf.o
          endif
        else
          ifeq ($(Z2_DEV),2)
            ifeq ($(Z5_SMT),0)
              MDIR = $(DID)/migd.o $(DID)/miggp.o $(DID)/dirgp.o
              DIR = $(DID)/miggp.o $(DID)/miggpd.o $(DID)/dirgp.o
            else
              MDIR = $(DID)/miad.o $(DID)/miagp.o $(DID)/dirgp.o
              DIR = $(DID)/miagp.o $(DID)/miagpd.o $(DID)/dirgp.o
            endif
          else
            ifeq ($(Z5_SMT),0)
              DIR = $(DID)/mig.o $(DID)/migd.o
            else
              DIR = $(DID)/mia.o $(DID)/miad.o
            endif
          endif
        endif
      else

        ifeq ($(Z6_EQU),4)
          ifeq ($(Z2_DEV),1)
            ifeq ($(Z5_SMT),0)
              MDIR = $(DID)/psegd.o $(DID)/psegmd.o
              DIR = $(DID)/psegmd.o $(DID)/psegmdf.o
            else
              MDIR = $(DID)/psead.o $(DID)/pseamd.o
              DIR = $(DID)/pseamd.o $(DID)/pseamdf.o
            endif
          else
            ifeq ($(Z2_DEV),2)
              ifeq ($(Z5_SMT),0)
                MDIR = $(DID)/psegd.o $(DID)/pseggp.o $(DID)/dirgp.o
                DIR = $(DID)/pseggp.o $(DID)/pseggpd.o $(DID)/dirgp.o
              else
                MDIR = $(DID)/psead.o $(DID)/pseagp.o $(DID)/dirgp.o
                DIR = $(DID)/pseagp.o $(DID)/pseagpd.o $(DID)/dirgp.o
              endif
            else
              ifeq ($(Z5_SMT),0)
                DIR = $(DID)/pseg.o $(DID)/psegd.o
              else
                DIR = $(DID)/psea.o $(DID)/psead.o
              endif
            endif
          endif
        else

          ifeq ($(Z6_EQU),8)
            ifeq ($(Z2_DEV),1)
              MDIR = $(DID)/rem.o $(DID)/remd.o
              DIR = $(DID)/rem.o $(DID)/remd.o
            else
              ifeq ($(Z2_DEV),2)
                MDIR = $(DID)/remd.o $(DID)/remgp.o $(DID)/dirgp.o
                DIR = $(DID)/remgp.o $(DID)/remgpd.o $(DID)/dirgp.o
              else
                DIR = $(DID)/rem.o $(DID)/remd.o
              endif
            endif
          else

            ifeq ($(Z6_EQU),9)
              MDIR = $(DID)/wall.o $(DID)/walld.o
              DIR = $(DID)/wall.o $(DID)/walld.o
            else

              ifeq ($(Z6_EQU),-1)
                ifeq ($(Z2_DEV),1)
                  MDIR = $(DID)/potd.o $(DID)/potmd.o
                  DIR = $(DID)/potmd.o $(DID)/potmdf.o
                else
                  ifeq ($(Z2_DEV),2)
                    MDIR = $(DID)/potd.o $(DID)/potgp.o $(DID)/dirgp.o
                    DIR = $(DID)/potgp.o $(DID)/potgpd.o $(DID)/dirgp.o
                  else
                    DIR = $(DID)/pot.o $(DID)/potd.o
                  endif
                endif
              else

                ifeq ($(Z2_DEV),1)
                  ifeq ($(Z5_SMT),0)
                    MDIR = $(DID)/fgtgd.o $(DID)/fgtgmd.o
                    DIR = $(DID)/fgtgmd.o $(DID)/fgtgmdf.o
                  else
                    MDIR = $(DID)/fgtad.o $(DID)/fgtamd.o
                    DIR = $(DID)/fgtamd.o $(DID)/fgtamdf.o
                  endif
                else
                  ifeq ($(Z2_DEV),2)
                    ifeq ($(Z5_SMT),0)
                      MDIR = $(DID)/fgtgd.o $(DID)/fgtggp.o $(DID)/dirgp.o
                      DIR = $(DID)/fgtggp.o $(DID)/fgtggpd.o $(DID)/dirgp.o
                    else
                      MDIR = $(DID)/fgtad.o $(DID)/fgtagp.o $(DID)/dirgp.o
                      DIR = $(DID)/fgtagp.o $(DID)/fgtagpd.o $(DID)/dirgp.o
                    endif
                  else
                    ifeq ($(Z5_SMT),0)
                      DIR = $(DID)/fgtg.o $(DID)/fgtgd.o
                    else
                      DIR = $(DID)/fgta.o $(DID)/fgtad.o
                    endif
                  endif
                endif
              endif
            endif
          endif
        endif
      endif
    endif
  endif
endif
ifeq ($(Z4_GEO),0)
#  DEF += parameterf.o
  ifeq ($(Z3_FMM),0)
    FMO = $(FMD)/fmmf.o
  else
    ifeq ($(Z3_FMM),1)
      FMO = $(FMD)/fmmf.o
    else
      FMO = $(FMD)/treef.o
    endif
  endif
  ifeq ($(Z1_MPI),0)
    FMO += $(SUB)/ifbox.o
  else
    MPI += $(MPID)/mpiifbox.o
  endif
else
  ifeq ($(Z4_GEO),1)
#    DEF += parameterp.o
    ifeq ($(Z3_FMM),0)
      FMO = $(FMD)/fmmp.o
    else
      ifeq ($(Z3_FMM),1)
        FMO = $(FMD)/fmmp.o
      else
        FMO = $(FMD)/treep.o
      endif
    endif
    ifeq ($(Z1_MPI),0)
      FMO += $(SUB)/ipbox.o
    else
      MPI += $(MPID)/mpiipbox.o
    endif
  else
    ifeq ($(Z4_GEO),2)
#      DEF += parameters.o
      ifeq ($(Z3_FMM),0)
        FMO = $(FMD)/fmms.o
      else
        ifeq ($(Z3_FMM),1)
          FMO = $(FMD)/fmms.o
        else
          FMO = $(FMD)/trees.o
        endif
      endif
      ifeq ($(Z1_MPI),0)
        FMO += $(SUB)/isbox.o
      else
        MPI += $(MPID)/mpiisbox.o
      endif
    else
#      DEF += parameterc.o
      ifeq ($(Z3_FMM),0)
        FMO = $(FMD)/fmmc.o
      else
        ifeq ($(Z3_FMM),1)
          FMO = $(FMD)/fmmc.o
        else
          FMO = $(FMD)/treec.o
        endif
      endif
      ifeq ($(Z1_MPI),0)
        FMO += $(SUB)/icbox.o
      else
        MPI += $(MPID)/mpiicbox.o
      endif
    endif
  endif
endif
ifeq ($(Z4_GEO),2)
  ifeq ($(Z3_FMM),0)
    MPI += $(MPID)/mpisendm2lfs.o
  else
    ifeq ($(Z3_FMM),2)
      MPI += $(MPID)/mpisendm2lfs.o
    else
      MPI += $(MPID)/mpisendm2lps.o
    endif
  endif
else
  ifeq ($(Z3_FMM),0)
    MPI += $(MPID)/mpisendm2lf.o
  else
    ifeq ($(Z3_FMM),2)
      MPI += $(MPID)/mpisendm2lf.o
    else
      MPI += $(MPID)/mpisendm2lp.o
    endif
  endif
endif
MPI += $(MPID)/mpirange.o
DIR += $(MISC)/memory.o
MDIR += $(MISC)/memory.o
VM = $(MISC)/statistics.o
ifeq ($(Z1_MPI),0)
  VM += $(MISC)/rbf.o
  PAR += $(MISC)/rbf.o
  OBJM = $(DEF) mtest.o $(MDIR)
  OBJN = $(DEF) ntest.o $(DIR) $(FMO) $(FMM)
  OBJN2 = $(DEF) ntest2.o $(DIR) $(FMO) $(FMM)
  OBJO = $(DEF) otest.o $(DIR) $(FMO) $(FMM)
  OBJP = $(DEF) ptest.o $(DIR) $(FMO) $(FMM)
  OBJR = $(DEF) rtest.o $(DIR) $(FMO) $(FMM) $(VM)
else
  MVM = $(MPID)/mpirbf.o
  PAR += $(MPID)/mpirbf.o
  OBJM = $(DEF) mpimtest.o $(MDIR) $(MPID)/mpirange.o $(MPID)/mpiallreduce.o
  OBJN = $(DEF) mpintest.o $(DIR) $(FMO) $(FMM) $(MPI)
  OBJN2 = $(DEF) mpintest2.o $(DIR) $(FMO) $(FMM) $(MPI)
  OBJO = $(DEF) mpiotest.o $(DIR) $(FMO) $(FMM) $(MPI)
  OBJP = $(DEF) mpiptest.o $(DIR) $(FMO) $(FMM) $(MPI)
  OBJR = $(DEF) mpirtest.o $(DIR) $(FMO) $(FMM) $(MPI) $(VM) $(MVM)
endif
OBJQ = $(DEF) qtest.o $(DIR)
OBJT = $(DEF) ttest.o $(DIR)
PAR += $(DIR) $(MDIR) $(FMO)
.cxx.o:
	$(CCS) -c $< -o $@
m: $(OBJM)
	$(CC) $(OBJM) $(CCL) $(MDGL)
n: $(OBJN)
	$(CC) $(OBJN) $(CCL) $(MDGL)
n2: $(OBJN2)
	$(CC) $(OBJN2) $(CCL) $(MDGL)
o: $(OBJO)
	$(CC) $(OBJO) $(CCL) $(MDGL)
p: $(OBJP)
	$(CC) $(OBJP) $(CCL) $(MDGL)
q: $(OBJQ)
	$(CC) $(OBJQ) $(CCL) $(MDGL)
r: $(OBJR)
	$(CC) $(OBJR) $(CCL) $(MDGL) -lfftw3
t: $(OBJT)
	$(CC) $(OBJT) $(CCL) $(MDGL)

mpimtest.o: mpimtest.cxx
	$(CCM) -c $?
mpintest.o: mpintest.cxx
	$(CCM) -c $?
mpintest2.o: mpintest2.cxx
	$(CCM) -c $?
mpiotest.o: mpiotest.cxx
	$(CCM) -c $?
mpiptest.o: mpiptest.cxx
	$(CCM) -c $?
mpirtest.o: mpirtest.cxx
	$(CCM) -c $?
$(MPI): $(MPI:.o=.cxx)
	$(CCM) $(MDGI) -c $*.cxx -o $@
$(MVM): $(MVM:.o=.cxx)
	$(CCM) $(MDGI) -c $*.cxx -o $@
$(DID)/bsgmd.o: $(DID)/bsgmd.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/bsgmdf.o: $(DID)/bsgmdf.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/bsamd.o: $(DID)/bsamd.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/bsamdf.o: $(DID)/bsamdf.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/stgmd.o: $(DID)/stgmd.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/stgmdf.o: $(DID)/stgmdf.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/stamd.o: $(DID)/stamd.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/stamdf.o: $(DID)/stamdf.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/trgmd.o: $(DID)/trgmd.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/trgmdf.o: $(DID)/trgmdf.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/tramd.o: $(DID)/tramd.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/tramdf.o: $(DID)/tramdf.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/migmd.o: $(DID)/migmd.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/migmdf.o: $(DID)/migmdf.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/miamd.o: $(DID)/miamd.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/miamdf.o: $(DID)/miamdf.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/psegmd.o: $(DID)/psegmd.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/psegmdf.o: $(DID)/psegmdf.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/pseamd.o: $(DID)/pseamd.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/pseamdf.o: $(DID)/pseamdf.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/fgtgmd.o: $(DID)/fgtgmd.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/fgtgmdf.o: $(DID)/fgtgmdf.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/fgtamd.o: $(DID)/fgtamd.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/fgtamdf.o: $(DID)/fgtamdf.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/potmd.o: $(DID)/potmd.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(DID)/potmdf.o: $(DID)/potmdf.cxx
	$(CCS) $(MDGI) -c $? -o $@
$(MDG): $(MDG:.o=.cxx)
	$(CCS) $(MDGI) -c $*.cxx -o $@

ex_m:
	make m; $(EX)
ex_n:
	make n; $(EX)
ex_o:
	make o; $(EX)
ex_p:
	make p; $(EX)
ex_q:
	make q; $(EX)
ex_r:
	make r; $(EX)
ex_t:
	make t; $(EX)
all:
	./mtest.sh
	./ntest.sh
	./otest.sh
	./ptest.sh
	./qtest.sh
	./rtest.sh
	./ttest.sh

