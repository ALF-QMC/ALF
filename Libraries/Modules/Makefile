.PHONY : lib clean

LIB=modules_90.a

SRCS = Mat_subroutines_mod.F90 \
	alf_hdf5_mod.F90 \
	check_omp_num_threads_mod.F90 \
	entanglement_mod.F90 \
	errors_mod.F90 \
	files_mod.F90 \
	fourier_mod.F90 \
	histograms_mod.F90 \
	histograms_v2_mod.F90 \
	iso8601_datetime_mod.F90 \
	lattices_interface_mod.F90 \
	lattices_v3_mod.F90 \
	log_mesh_mod.F90 \
	matrix_mod.F90 \
	maxent_mod.F90 \
	maxent_stoch_mod.F90 \
	mpi_shared_mem_mod_v2.F90 \
	mymats_mod.F90 \
	natural_constants_mod.F90 \
	precdef_mod.F90 \
	random_wrap_mod.F90 \
	runtime_error_mod.F90

# Dependency helper variables
DEPSFILE := .deps.mk
DEPSGEN := ../../Prog/gen_deps.py

lib: $(LIB)

# OBJS and MODS are auto-generated into $(DEPSFILE) by gen_deps.py
-include $(DEPSFILE)

$(LIB): $(DEPSFILE) $(OBJS)
	ar -r $(LIB) $(OBJS)

lattices_x86-64.so: runtime_error_mod.F90 lattices_interface_mod.F90 matrix_mod.F90 lattices_v3_mod.F90
	(make clean)
	gfortran -O3 -march=x86-64 -c runtime_error_mod.F90
	gfortran -O3 -march=x86-64 -c matrix_mod.F90
	gfortran -O3 -march=x86-64 -c lattices_v3_mod.F90
	gfortran -O3 -march=x86-64 -shared -fPIC lattices_interface_mod.F90 -o lattices_x86-64.so \
			runtime_error_mod.o matrix_mod.o lattices_v3_mod.o

lattices_armv8.4-a.so: runtime_error_mod.F90 lattices_interface_mod.F90 matrix_mod.F90 lattices_v3_mod.F90
	(make clean)
	gfortran -O3 -march=armv8.4-a -c runtime_error_mod.F90
	gfortran -O3 -march=armv8.4-a -c matrix_mod.F90
	gfortran -O3 -march=armv8.4-a -c lattices_v3_mod.F90
	gfortran -O3 -march=armv8.4-a -shared -fPIC lattices_interface_mod.F90 -o lattices_armv8.4-a.so \
			runtime_error_mod.o matrix_mod.o lattices_v3_mod.o

%.o: %.F90
	$(ALF_FC) -c $(ALF_FLAGS_MODULES) $<

$(DEPSFILE): $(SRCS) $(DEPSGEN)
	python3 $(DEPSGEN) $(SRCS) > $(DEPSFILE)

clean: $(DEPSFILE)
	rm -f $(OBJS) $(MODS) $(LIB) $(DEPSFILE)

#unused modules: histograms histograms_v2 log_mesh natural_constants
#unused, except for in Prog/FFA_Orginals/: precdef
