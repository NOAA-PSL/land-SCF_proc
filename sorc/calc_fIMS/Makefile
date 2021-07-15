SHELL=	/bin/ksh
FC      =$(FCMP)
LIBS    =$(LIBSM) 
OBJS	=driver_fIMS.o IMSaggregate_mod.o

CMD     =calcfIMS
$(CMD):	$(OBJS)
	${FC} $(FFLAGS) $(LDFLG) $(OBJS) $(LIBS) -o $(CMD)

IMSaggregate_mod.o: IMSaggregate_mod.f90
	$(FC)  $(FFLAGS) -c IMSaggregate_mod.f90

driver_fIMS.o: IMSaggregate_mod.o driver_fIMS.f90
	$(FC)  $(FFLAGS) -c driver_fIMS.f90

install:
	-cp -p $(CMD) ../../exec/.
clean:
	-rm -f $(OBJS) *.mod $(CMD)
