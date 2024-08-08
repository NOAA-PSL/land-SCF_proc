SHELL=	/bin/ksh
FC      =$(FCMP)
LIBS    =$(LIBSM) 
OBJS	=driver_fSCF.o SCFaggregate_mod.o

CMD     =calcfSCF
$(CMD):	$(OBJS)
	${FC} $(FFLAGS) $(OBJS) $(LIBS) -o $(CMD)

SCFaggregate_mod.o: SCFaggregate_mod.f90
	$(FC)  $(FFLAGS) -c SCFaggregate_mod.f90

driver_fSCF.o: SCFaggregate_mod.o driver_fSCF.f90
	$(FC)  $(FFLAGS) -c driver_fSCF.f90

install:
	-cp -p $(CMD) ../exec/.
clean:
	-rm -f $(OBJS) *.mod $(CMD)
