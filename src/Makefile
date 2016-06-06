#----------------------------------------------------------------------
# OPTIONS FOR RUNNING ON SCINET-GPC MACHINES
#----------------------------------------------------------------------

F90 = mpif90
F77 = mpif77
OPTIMIZE =  -O4 -w -mcmodel=large -shared-intel 
FFTWOMP = -lfftw3f_omp
MODFLAG = -module 
OMPLIB = -openmp

CC = mpicc
C++ = mpiCC
CCOPTIMIZE = -w 

#----------------------------------------------------------------------
# DO NOT MODIFY THIS FILE
#----------------------------------------------------------------------

FFTLIB = -L$(FFTW_PATH)/lib $(FFTWOMP) -lfftw3f_mpi -lfftw3f
FFTINC = -I$(FFTW_PATH)/include 

CCOPTIONS = $(CDEFS) $(CCOPTIMIZE)

# EXTERNAL MODULES
exdir = ./modules/External
hpx_mods = \
	    $(exdir)/healpix_types.o\
	    $(exdir)/cgetEnvironment.o\
	    $(exdir)/extension.o\
	    $(exdir)/long_intrinsic.o\
	    $(exdir)/misc_utils.o \
	    $(exdir)/num_rec.o\
	    $(exdir)/bit_manipulation.o\
	    $(exdir)/indmed.o\
	    $(exdir)/statistics.o\
	    $(exdir)/pix_tools.o\
	    $(exdir)/fitstools.o\
	    $(exdir)/head_fits.o 

ex_mods  = $(exdir)/intreal_types.o\
	 $(exdir)/openmpvars.o\
	 $(exdir)/mpivars.o\
	 $(exdir)/textlib.o\
       	 $(exdir)/Type_Kinds.o\
       	 $(exdir)/Endian_Utility.o\
	 $(exdir)/timing_diagnostics.o\
       	 $(exdir)/myio.o\
       	 $(exdir)/memorytracking.o \
       	 $(exdir)/memory_management.o

# IO MODULES
iodir   = ./modules/IO
io_mods = $(iodir)/io.o

# RANDOM FIELD MODULES
rfdir   = ./modules/RandomField
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
	$(rfdir)/chi2zeta.o\
	$(rfdir)/RandomFieldWrapper.o

# SLAB TO CUBE FILES
scdir = ./modules/SlabToCube
sc_mods = $(scdir)/SlabToCube.o

# HPKVD MODULE AND OBJECT FILES
hpdir  = ./hpkvd_src/
hpkvd_mods = \
       $(hpdir)/arrays.o \
       $(hpdir)/input_parameters.o \
       $(hpdir)/Dlin_params_commod.o \
       $(hpdir)/Dlin_tab_commod.o \
       $(hpdir)/EdeS_commod.o \
       $(hpdir)/bound_commod.o \
       $(hpdir)/bvalues_commod.o \
       $(hpdir)/ccoordstrat_commod.o \
       $(hpdir)/count_commod.o \
       $(hpdir)/cran3_commod.o \
       $(hpdir)/cran3b_commod.o \
       $(hpdir)/cubic_params_commod.o \
       $(hpdir)/decay_commod.o \
       $(hpdir)/equil_commod.o \
       $(hpdir)/equil_sph_commod.o \
       $(hpdir)/etafiles1_commod.o \
       $(hpdir)/etafiles2_commod.o \
       $(hpdir)/etafiles_commod.o \
       $(hpdir)/evalues_commod.o \
       $(hpdir)/evalues_sph_commod.o \
       $(hpdir)/evpvtable_commod.o \
       $(hpdir)/fac_commod.o \
       $(hpdir)/h0_commod.o \
       $(hpdir)/hard_commod.o \
       $(hpdir)/ipart_commod.o \
       $(hpdir)/lad_commod.o \
       $(hpdir)/lists_commod.o \
       $(hpdir)/ncorr_commod.o \
       $(hpdir)/nodes_commod.o \
       $(hpdir)/params_commod.o \
       $(hpdir)/params_sph_commod.o \
       $(hpdir)/particle_commod.o \
       $(hpdir)/powerc1_commod.o \
       $(hpdir)/powerc2_commod.o \
       $(hpdir)/qtime_commod.o \
       $(hpdir)/quant_commod.o \
       $(hpdir)/rd_commod.o \
       $(hpdir)/rpart_commod.o \
       $(hpdir)/seq_name_commod.o \
       $(hpdir)/sig_commod.o \
       $(hpdir)/table_commod.o \
       $(hpdir)/tabvals_sph_commod.o \
       $(hpdir)/tabvec_commod.o \
       $(hpdir)/tau_commod.o \
       $(hpdir)/univer2_commod.o \
       $(hpdir)/univer_commod.o \
       $(hpdir)/univercurv_commod.o \
       $(hpdir)/wantevpvtab_commod.o \
       $(hpdir)/io.o 

hpkvd_objs = \
       $(hpdir)/peakvoidsubs.o \
       $(hpdir)/fftsubs.o \
       $(hpdir)/dens_power2.o \
       $(hpdir)/homeosubs.o \
       $(hpdir)/psubs_Dlinear.o \
       $(hpdir)/hom_sphere_tab.o \
       $(hpdir)/hpkvd.o

etab_objs = $(hpdir)/run_hom_ellipse_tab.o 
	 
# MERGE_PKVD MODULE AND OBJECT FILES
mgdir  = ./merge_pkvd_src/
merge_mods = \
       $(mgdir)/bound_commod.o\
       $(mgdir)/intreal_types.o\
       $(mgdir)/array_bounds.o\
       $(mgdir)/rpart_commod.o\
       $(mgdir)/mparam.o\
       $(mgdir)/ipart_commod.o\
       $(mgdir)/lists_commod.o\
       $(mgdir)/nodes_commod.o\
       $(mgdir)/mnodes.o\
       $(mgdir)/mbound.o\
       $(mgdir)/mlists.o\
       $(mgdir)/merge_subs.o\
       $(mgdir)/mpivars.o\
       $(mgdir)/textlib.o\
       $(mgdir)/timing_diagnostics.o\
       $(mgdir)/decay_commod.o\
       $(mgdir)/Dlin_tab_commod.o\
       $(mgdir)/Dlin_params_commod.o\
       $(mgdir)/input_parameters.o\
       $(mgdir)/psubs_mods.o\
       $(mgdir)/psubs_Dlinear.o

merge_objs = $(mgdir)/BitSubs.o\
	     $(mgdir)/merge_pkvd.o

# PKS2MAP MODULE AND OBJECT FILES
pmdir  = ./pks2map_src/
p2m_mods = \
	 $(exdir)/mpivars.o\
	 $(exdir)/textlib.o\
	 $(rfdir)/random.o\
	 $(pmdir)/healpixvars.o\
	 $(pmdir)/flatskyvars.o\
	 $(pmdir)/fitsvars.o\
	 $(pmdir)/cosmology.o\
	 $(pmdir)/bbps_profile.o\
	 $(pmdir)/integrate_profiles.o\
	 $(pmdir)/maptable.o\
	 $(pmdir)/haloproject.o\
	 $(pmdir)/pksc.o

p2m_objs = $(pmdir)/pks2map.o 	 

mmt_mods = \
	 $(exdir)/textlib.o\
	 $(pmdir)/cosmology.o\
	 $(pmdir)/bbps_profile.o\
	 $(pmdir)/integrate_profiles.o\
	 $(pmdir)/maptable.o

mmt_objs = $(pmdir)/make_maptable.o

moddir = ./modules/
bindir = ../bin

OPTIONS       =  $(OPTIMIZE) $(FFTLIB) $(FFTINC) $(MODFLAG)$(moddir) $(OMPLIB)
OPTIONS_MERGE =  $(OPTIMIZE) $(MODFLAG)$(moddir)

EXEC_h = hpkvd
OBJS_h = $(ex_mods) $(rf_mods) $(sc_mods) $(hpkvd_mods) $(hpkvd_objs) 

EXEC_m = merge_pkvd
OBJS_m = $(merge_mods) $(merge_objs)

EXEC_t = make_maptable
OBJS_t = $(mmt_mods) $(mmt_objs)

EXEC_p = pks2map
OBJS_p = $(hpx_mods) $(p2m_mods) $(p2m_objs)

.SUFFIXES: .o .f .f90 .F90 .c .C

$(scdir)/%.o: $(scdir)/%.f90 
	$(F90) $(OPTIONS) -c $< -o $@

$(rfdir)/%.o: $(rfdir)/%.f
	$(F77) $(OPTIONS) -c $< -o $@
$(rfdir)/%.o: $(rfdir)/%.f90
	$(F90) $(OPTIONS) -c $< -o $@

$(mgdir)/%.o: $(mgdir)/%.f90
	$(F90) $(OPTIONS) -c $< -o $@
$(mgdir)/%.o: $(mgdir)/%.F90
	$(F90) $(OPTIONS) -c $< -o $@

$(exdir)/%.o: $(exdir)/%.F90 
	$(F90) $(OPTIONS) -c $< -o $@
$(exdir)/%.o: $(exdir)/%.f90 
	$(F90) $(OPTIONS) -c $< -o $@
$(exdir)/%.o: $(exdir)/%.c
	$(CC)  $(CCOPTIONS) -c $< -o $@
$(exdir)/%.o: $(exdir)/%.C
	$(C++)  $(CCOPTIONS) -c $< -o $@

$(hpdir)/%.o: $(hpdir)/%.f90 
	$(F90) $(OPTIONS) -c $< -o $@
$(hpdir)/%.o: $(hpdir)/%.C 
	$(C++)  $(OPTIONS) -c $< -o $@

$(pmdir)/%.o: $(pmdir)/%.f90 
	$(F90) $(OPTIONS) -c $< -o $@
$(pmdir)/%.o: $(pmdir)/%.C 
	$(C++)  $(OPTIONS) -c $< -o $@

$(mgdir)/%.o: $(mgdir)/%.f90 
	$(F90) $(OPTIONS_MERGE) -c $< -o $@
$(mgdir)/%.o: $(mgdir)/%.C 
	$(C++) $(OPTIONS_MERGE) -c $< -o $@

$(EXEC_h): $(OBJS_h) 
	$(F90) $(OPTIONS) $(OBJS_h) $(FFTLIB) -o  $(bindir)/$(EXEC_h)  
$(EXEC_m): $(OBJS_m) 
	$(F90) $(OPTIONS) $(OBJS_m) $(FFTLIB) -o  $(bindir)/$(EXEC_m)
$(EXEC_t): $(OBJS_t) 
	$(F90) $(OPTIONS) $(OBJS_t) -lm -o  $(bindir)/$(EXEC_t)  
$(EXEC_p): $(OBJS_p) 
	$(F90) $(OPTIONS) $(OBJS_p) -L$(CFITSIO_LIBDIR) -lm -lcfitsio -o\
	$(bindir)/$(EXEC_p)  

EXEC = $(bindir)/$(EXEC_h) $(bindir)/$(EXEC_m) \
       $(bindir)/$(EXEC_t) $(bindir)/$(EXEC_p)
OBJS = $(OBJS_h) $(OBJS_m) $(OBJS_t) $(OBJS_p)

clean:
	rm -f $(EXEC) $(OBJS) $(moddir)/*.mod $(rfdir)/*.o $(hpdir)/*.o\
	      $(mgdir)/*.o $(exdir)/*.o $(scdir)/*.o     

