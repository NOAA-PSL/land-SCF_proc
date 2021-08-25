SHELL=	/bin/ksh
FC      =$(FCMP)
LIBS    =$(LIBSM) 
OBJS	=driver_fIMS_serial.o IMSaggregate_serial_mod.o

CMD     =calcfIMS
$(CMD):	$(OBJS)
	${FC} $(FFLAGS) $(LDFLG) $(OBJS) $(LIBS) -o $(CMD)

IMSaggregate_serial_mod.o: IMSaggregate_serial_mod.f90
	$(FC)  $(FFLAGS) -c IMSaggregate_serial_mod.f90

driver_fIMS_serial.o: IMSaggregate_serial_mod.o driver_fIMS_serial.f90
	$(FC)  $(FFLAGS) -c driver_fIMS_serial.f90

install:
	-cp -p $(CMD) ../../exec/.
clean:
	-rm -f $(OBJS) *.mod $(CMD)
