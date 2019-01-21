include ../Makefile

BEM = $(MISC)/geometry.o $(MISC)/bem3d.o $(MISC)/velocity.o\
	$(MISC)/bicgstab.o $(MISC)/bicgstab2.o $(MISC)/intrude.o\
	$(MISC)/savedata.o $(MISC)/memory.o

OBJSP = panel.o $(BEM)

panels: $(OBJSP)
	$(CLINKER) $? $(PETSC_SNES_LIB) $(MDGL)
