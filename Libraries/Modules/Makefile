LIB=modules_90.a
OBJS=runtime_error_mod.o errors_mod.o files_mod.o   fourier_mod.o histograms_mod.o histograms_v2_mod.o lattices_v3_mod.o \
	log_mesh_mod.o mymats_mod.o matrix_mod.o maxent_mod.o   maxent_stoch_mod.o   mpi_shared_mem_mod_v2.o \
	natural_constants_mod.o precdef_mod.o random_wrap_mod.o Mat_subroutines_mod.o entanglement_mod.o   alf_hdf5_mod.o \
	check_omp_num_threads_mod.o
MODS=runtime_error_mod.mod errors.mod   files_mod.mod fourier.mod   histograms.mod   histograms_v2.mod   lattices_v3.mod   \
	log_mesh.mod   mymats.mod   matrix.mod   maxent_mod.mod maxent_stoch_mod.mod mpi_shared_memory.mod \
	natural_constants.mod   precdef.mod   random_wrap.mod   mat_subroutines.mod   entanglement_mod.mod alf_hdf5.mod \
	check_omp_num_threads_mod.mod

$(LIB): $(OBJS)
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

alf_hdf5_mod.o: runtime_error_mod.o
errors_mod.o: runtime_error_mod.o mymats_mod.o random_wrap_mod.o
fourier_mod.o: runtime_error_mod.o maxent_mod.o maxent_stoch_mod.o matrix_mod.o
lattices_v3_mod.o: runtime_error_mod.o mymats_mod.o 
maxent_mod.o: mymats_mod.o errors_mod.o runtime_error_mod.o
histograms_v2_mod.o: log_mesh_mod.o
log_mesh_mod.o: runtime_error_mod.o
mat_subroutines.o: runtime_error_mod.o
mpi_shared_mem_mod_v2.o: runtime_error_mod.o
mymats.o: runtime_error_mod.o

clean:
	rm -f $(OBJS) $(MODS) $(LIB)

#unused modules: histograms histograms_v2 log_mesh natural_constants
#unused, except for in Prog/FFA_Orginals/: precdef
