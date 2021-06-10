LIB=modules_90.a
OBJS=errors_mod.o files_mod.o   fourier_mod.o histograms_mod.o histograms_v2_mod.o lattices_v3_mod.o \
	log_mesh_mod.o mymats_mod.o matrix_mod.o maxent_mod.o   maxent_stoch_mod.o  mpi_shared_mem_mod_v2.o  \
	natural_constants_mod.o precdef_mod.o random_wrap_mod.o Mat_subroutines.o entanglement_mod.o
MODS=errors.mod   files_mod.mod fourier.mod   histograms.mod   histograms_v2.mod   lattices_v3.mod   \
	log_mesh.mod   mymats.mod   matrix.mod   maxent_mod.mod maxent_stoch_mod.mod mpi_shared_memory.mod \
	natural_constants.mod   precdef.mod   random_wrap.mod  entanglement_mod.mod

$(LIB): $(OBJS)
	ar -r $(LIB) $(OBJS)

%.o: %.F90
	$(ALF_FC) -c $(ALF_FLAGS_MODULES) $<

errors_mod.o: mymats_mod.o random_wrap_mod.o
fourier_mod.o: maxent_mod.o maxent_stoch_mod.o matrix_mod.o
lattices_v3_mod.o: mymats_mod.o 
maxent_mod.o: mymats_mod.o errors_mod.o
histograms_v2_mod.o: log_mesh_mod.o

clean:
	rm -f $(OBJS) $(MODS) $(LIB)

#unused modules: histograms histograms_v2 log_mesh natural_constants
#unused, except for in Prog/FFA_Orginals/: precdef
