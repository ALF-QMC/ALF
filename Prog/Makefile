.PHONY : Compile tidy clean program parse_parameters echo_compile_mods

OHAM = `./parse_ham.py --print_obj_list`

OBJS = Hamiltonians/LRC_mod.o control_mod.o Fields_mod.o Operator_mod.o WaveFunction_mod.o observables_mod.o \
	ContainerElementBase_mod.o DynamicMatrixArray_mod.o OpTTypes_mod.o OpT_time_dependent.o\
	Predefined_Int_mod.o Predefined_Obs_mod.o Predefined_Latt_mod.o Predefined_Hop_mod.o  Predefined_Trial_mod.o \
	Hamiltonian_main_mod.o QDRP_decompose_mod.o udv_state_mod.o Hop_mod.o UDV_WRAP_mod.o \
	wrapul_mod.o   cgr1_mod.o   wrapur_mod.o   cgr2_2_mod.o   upgrade_mod.o   Set_random_mod.o \
	Global_mod.o Langevin_HMC_mod.o Wrapgr_mod.o tau_m_mod.o tau_p_mod.o \
	main.o

MODS = control.mod fields_mod.mod global_mod.mod hop_mod.mod lrc_mod.mod observables.mod \
	operator_mod.mod predefined_int.mod predefined_lattices.mod predefined_obs.mod \
	predefined_hoppings.mod  predefined_trial.mod qdrp_mod.mod tau_m_mod.mod tau_p_mod.mod \
	udv_state_mod.mod udv_wrap_mod.mod wavefunction_mod.mod wrapgr_mod.mod hamiltonian.mod \
	opt_time_dependent.mod containerelementbase_mod.mod opttypes_mod.mod dynamicmatrixarray_mod.mod langevin_hmc_mod.mod \
	wrapul_mod.mod cgr1_mod.mod wrapur_mod.mod cgr2_2_mod.mod upgrade_mod.mod set_random.mod

program: git.h parse_parameters echo_compile_mods $(OBJS)
	make $(OHAM)
	@echo "Link program"
	$(ALF_FC) -o "ALF.out" $(OBJS) $(OHAM) $(ALF_LIB)

Hamiltonian_main_mod.o: Hamiltonians_interface.h Hamiltonians_case.h

git.h:
	./Git_config.sh

parse_parameters:
	@echo "Parsing Hamiltonian parameters"
	./parse_ham.py --create_read_write_par

Hamiltonians_interface.h:
	./parse_ham.py --create_hamiltonians_interface

Hamiltonians_case.h:
	./parse_ham.py --create_hamiltonians_case

%.o: %.F90
	$(ALF_FC) -c -o $@ $(ALF_FLAGS_PROG) $<

tidy:
	rm -f $(OBJS) $(OHAM) $(MODS) hamiltonian*.mod hamiltonian*.smod
	rm -f git.h Hamiltonians_interface.h Hamiltonians_case.h Hamiltonians/Hamiltonian_*_read_write_parameters.F90

clean: tidy
	rm -f "ALF.out"

echo_compile_mods:
	@echo "Compiling program modules"
