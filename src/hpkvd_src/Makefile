#----------------------------------------------------------------------
# OPTIONS FOR RUNNING ON SCINET-GPC MACHINES
#----------------------------------------------------------------------

F90 = mpif90
F77 = mpif77
OPTIMIZE =  -O4 -w -mcmodel=large -shared-intel -openmp
FFTWOMP = -lfftw3f_omp

CC = mpicc
CCOPTIMIZE = -O2 -w 

#----------------------------------------------------------------------
# DO NOT MODIFY THIS FILE
#----------------------------------------------------------------------

FFTLIB = -L$(FFTW_PATH)/lib $(FFTWOMP) -lfftw3f_mpi -lfftw3f
FFTINC = -I$(FFTW_PATH)/include 

OPTIONS =  $(OPTIMIZE) $(FFTLIB) $(FFTINC)
CCOPTIONS = $(CDEFS) $(CCOPTIMIZE)

exmddir = ./external_modules
healpixmods = \
	    $(exmddir)/healpix_types.o\
	    $(exmddir)/cgetEnvironment.o\
	    $(exmddir)/extension.o\
	    $(exmddir)/long_intrinsic.o\
	    $(exmddir)/misc_utils.o \
	    $(exmddir)/num_rec.o\
	    $(exmddir)/bit_manipulation.o\
	    $(exmddir)/indmed.o\
	    $(exmddir)/statistics.o\
	    $(exmddir)/pix_tools.o $(exmddir)/fitstools.o\
	    $(exmddir)/head_fits.o \

rfdir   = ./RandomField
rf_mods  = \
	$(rfdir)/globalvars.o\
	$(rfdir)/grid.o\
	$(rfdir)/fftw_interface.o\
	$(rfdir)/tiles.o\
	$(rfdir)/cosmology.o\
	$(rfdir)/random.o\
	$(rfdir)/pktable.o\
	$(rfdir)/gaussian_field.o\
	$(rfdir)/time.o\
	$(rfdir)/collapse.o\
	$(rfdir)/growth.o\
	$(rfdir)/RandomFieldWrapper.o

s2cdir = ./SlabToCube
s2c_mods = $(s2cdir)/SlabToCube.o

EXEC_s = stkbound
EXEC_b = binit_g
EXEC_e = run_hom_ellipse_tab
EXEC_h = hpkvd
EXEC_m = merge_pkvd
EXEC_t = make_maptable
EXEC_c = make_cltable
EXEC_p = pks2map
EXEC_h2f = hp2fits
EXEC_ts = test
hpkvd_mods = openmpvars.o myio.o memorytracking.o mpivars.o intreal_types.o memory_management.o arrays.o input_parameters.o textlib.o timing_diagnostics.o Type_Kinds.o Endian_Utility.o \
	   Dlin_params_commod.o Dlin_tab_commod.o EdeS_commod.o  \
	   bound_commod.o bvalues_commod.o ccoordstrat_commod.o \
	   count_commod.o cran3_commod.o \
	   cran3b_commod.o cubic_params_commod.o decay_commod.o \
	   equil_commod.o equil_sph_commod.o etafiles1_commod.o \
	   etafiles2_commod.o etafiles_commod.o evalues_commod.o \
	   evalues_sph_commod.o evpvtable_commod.o fac_commod.o \
	   h0_commod.o hard_commod.o \
	   ipart_commod.o lad_commod.o lists_commod.o ncorr_commod.o \
	   nodes_commod.o params_commod.o params_sph_commod.o \
	   particle_commod.o powerc1_commod.o powerc2_commod.o \
	   qtime_commod.o quant_commod.o rd_commod.o rpart_commod.o \
	   seq_name_commod.o sig_commod.o table_commod.o tabvals_sph_commod.o\
	   tabvec_commod.o tau_commod.o univer2_commod.o univer_commod.o \
	   univercurv_commod.o wantevpvtab_commod.o io.o

merge_mods = intreal_types.o array_bounds.o ipart_commod.o rpart_commod.o \
	     mparam.o nodes_commod.o bound_commod.o lists_commod.o mnodes.o 

OBJS_s = stkbound.o
OBJS_b = intreal_types.o binit_g.o 
OBJS_e = intreal_types.o homeosubs.o psubs_Dlinear.o hom_sphere_tab.o \
	 run_hom_ellipse_tab.o 
OBJS_h = $(hpkvd_mods)  $(rf_mods) $(s2c_mods)\
	 peakvoidsubs.o fftsubs.o dens_power2.o homeosubs.o \
	 psubs_Dlinear.o hom_sphere_tab.o hpkvd.o
OBJS_m = $(hpkvd_mods) mparam.o\
	 BitSubs.o mnodes.o \
	 mbound.o mlists.o psubs_Dlinear.o merge_subs.o merge_pkvd.o
OBJS_t = cosmology.o textlib.o bbps_profile.o integrate_profiles.o maptable.o \
	 make_maptable.o
OBJS_c = cosmology.o textlib.o bbps_profile.o integrate_profiles.o cltable.o \
	make_cltable.o
OBJS_p = textlib.o mpivars.o healpixvars.o flatskyvars.o fitsvars.o\
	 $(healpixmods) cosmology.o bbps_profile.o integrate_profiles.o \
	 maptable.o haloproject.o pksc.o pks2map.o 	 
OBJS_h2f = healpixvars.o fitsvars.o $(healpixmods) hp2fits.o	 

EXEC = $(EXEC_b) $(EXEC_e) $(EXEC_h) $(EXEC_m) $(EXEC_t) $(EXEC_p)
OBJS = $(OBJS_b) $(OBJS_e) $(OBJS_h) $(OBJS_m) $(OBJS_t) $(OBJS_p)

.SUFFIXES: .o .f .f90 .F90 .c .C

$(s2cdir)/%.o: $(s2cdir)/%.f90 
	$(F90) $(OPTIONS) -c $< -o $@

$(rfdir)/%.o: $(rfdir)/%.f90 
	$(F90) $(OPTIONS) -c $< -o $@
$(rfdir)/%.o: $(rfdir)/%.F90 
	$(F90) $(OPTIONS) -c $< -o $@
$(rfdir)/%.o: $(rfdir)/%.F
	$(F77) $(OPTIONS) -c $< -o $@
$(rfdir)/%.o: $(rfdir)/%.f
	$(F77) $(OPTIONS) -c $< -o $@
$(rfdir)/%.o: $(rfdir)/%.C
	$(CC) $(OPTIONS) -c $< -o $@
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
%.o: %.C
	$(CC) $(CCOPTIONS) -c $< 

$(EXEC_s): $(OBJS_s) 
	$(F90) $(OPTIONS) $(OBJS_s)     -o  ../bin/$(EXEC_s)  
$(EXEC_b): $(OBJS_b) 
	$(F90) $(OPTIONS) $(OBJS_b)     -o  ../bin/$(EXEC_b)  
$(EXEC_e): $(OBJS_e) 
	$(F90) $(OPTIONS) $(OBJS_e)     -o  ../bin/$(EXEC_e)  
$(EXEC_h): $(OBJS_h) 
	$(F90) $(OPTIONS) $(OBJS_h) $(FFTLIB) -o  ../bin/$(EXEC_h)  
$(EXEC_m): $(OBJS_m) 
	$(F90) $(OPTIONS) $(OBJS_m) -lm -o  ../bin/$(EXEC_m)  
$(EXEC_t): $(OBJS_t) 
	$(F90) $(OPTIONS) $(OBJS_t) -lm -o  ../bin/$(EXEC_t)  
$(EXEC_c): $(OBJS_c) 
	$(F90) $(OPTIONS) $(OBJS_c) -lm -o  ../bin/$(EXEC_c)  
$(EXEC_p): $(OBJS_p) 
	$(F90) $(OPTIONS) $(OBJS_p) -L$(CFITSIO_LIBDIR) -lm -lcfitsio -o  ../bin/$(EXEC_p)  
$(EXEC_h2f): $(OBJS_h2f) 
	$(F90) $(OPTIONS) $(OBJS_h2f) -L$(CFITSIO_LIBDIR) -lm -lcfitsio -o  ../bin/$(EXEC_h2f)  
$(EXEC_ts): test.o mpivars.o
	$(F90) $(OPTIONS) test.o mpivars.o -o ../bin/$(EXEC_ts)
clean:
	rm -f *~ *.o *.mod   $(rfdir)/*.o   $(rfdir)/*.mod\
                           $(exmddir)/*.o $(exmddir)/*.mod\
                            $(s2cdir)/*.o  $(s2cdir)/*.mod

