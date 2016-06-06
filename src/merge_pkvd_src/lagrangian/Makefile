#----------------------------------------------------------------------
# OPTIONS FOR RUNNING ON SCINET-GPC MACHINES
#----------------------------------------------------------------------

F90 = mpif90
OPTIMIZE =  -O4 -w -mcmodel=large -shared-intel -openmp

CC = mpicc
CCOPTIMIZE = -O2 -w 
#----------------------------------------------------------------------
# DO NOT MODIFY THIS FILE
#----------------------------------------------------------------------

OPTIONS =  $(OPTIMIZE) $(OPT) -module .
CCOPTIONS = $(CCOPTIMIZE) $(CCOPT)

exmddir = ../../external_modules

EXEC = merge_pkvd_lag
mods = sort2.o exclusion.o mpivars.o textlib.o timing_diagnostics.o 

OBJS = $(mods) merge_pkvd_lag.o

.SUFFIXES: .o .f .f90 .F90 .c

$(exmddir)/%.o: $(exmddir)/%.f90 
	$(F90) $(OPTIONS) -c $< -o $@
$(exmddir)/%.o: $(exmddir)/%.F90 
	$(F90) $(OPTIONS) -c $< -o $@
$(exmddir)/%.o: $(exmddir)/%.c
	$(CC) $(CCOPTIONS) -c $< -o $@
%.o: %.f
	$(F90) $(OPTIONS) -c $< 
%.o: %.f90
	$(F90) $(OPTIONS) -c $< 
%.o: %.F90
	$(F90) $(OPTIONS) -c $< 
%.o: %.c
	$(CC) $(CCOPTIONS) -c $< 

$(EXEC): $(OBJS) 
	$(F90) $(OPTIONS) $(OBJS) -o  $(EXEC)  

clean:
	rm -f *~ *.o *.mod $(exmddir)/*.o $(exmddir)/*.mod $(EXEC)

