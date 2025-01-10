.PHONY : Compile tidy clean program parse_parameters echo_compile_mods

OHAM = `while read -r ham; do \
	echo "Hamiltonians/Hamiltonian_$${ham}_smod.o"; \
done < Hamiltonians.list`

OBJS = control_mod.o Fields_mod.o Operator_mod.o WaveFunction_mod.o observables_mod.o \
	ContainerElementBase_mod.o DynamicMatrixArray_mod.o OpTTypes_mod.o \
	Predefined_Int_mod.o Predefined_Obs_mod.o Predefined_Latt_mod.o Predefined_Hop_mod.o  Predefined_Trial_mod.o \
	Hamiltonian_main_mod.o QDRP_decompose_mod.o udv_state_mod.o Hop_mod.o \
	gfun_mod.o upgrade_mod.o Set_random_mod.o Stepwlk.o main.o brw_init.o

MODS = control.mod fields_mod.mod hop_mod.mod observables.mod \
	operator_mod.mod predefined_int.mod predefined_lattices.mod predefined_obs.mod \
	predefined_hoppings.mod  predefined_trial.mod qdrp_mod.mod \
	udv_state_mod.mod  wavefunction_mod.mod hamiltonian.mod \
	containerelementbase_mod.mod opttypes_mod.mod dynamicmatrixarray_mod.mod \
	gfun_mod.mod upgrade_mod.mod set_random.mod stepwlk_mod.mod brw_init.mod

program: git.h parse_parameters Hamiltonians_interface.h Hamiltonians_case.h echo_compile_mods $(OBJS)
	make $(OHAM)
	@echo "Link program"
	$(ALF_FC) -o "CPMC.out" $(OBJS) $(OHAM) $(ALF_LIB)

git.h:
	./Git_config.sh

parse_parameters:
	@echo "Parsing Hamiltonian parameters"
	./parse_ham.py --create_read_write_par

Hamiltonians_interface.h: Hamiltonians.list
	@echo "interface" > "Hamiltonians_interface.h"
	@while read -r ham; do \
		echo "  module subroutine Ham_Alloc_$$ham()" >> "Hamiltonians_interface.h"; \
		echo "  end subroutine Ham_Alloc_$$ham" >> "Hamiltonians_interface.h"; \
	done < Hamiltonians.list 
	@echo "end interface" >> "Hamiltonians_interface.h"

Hamiltonians_case.h: Hamiltonians.list
	@rm -f "Hamiltonians_case.h"
	@while read -r ham; do \
		echo "Case ('$$ham')" >> "Hamiltonians_case.h"; \
		echo "  call Ham_Alloc_$$ham()" >> "Hamiltonians_case.h"; \
	done < Hamiltonians.list 

%.o: %.F90
	$(ALF_FC) -c -o $@ $(ALF_FLAGS_PROG) $<

tidy:
	rm -f $(OBJS) $(OHAM) $(MODS) hamiltonian*.mod hamiltonian*.smod
	rm -f git.h Hamiltonians_interface.h Hamiltonians_case.h Hamiltonians/Hamiltonian_*_read_write_parameters.F90

clean: tidy
	rm -f "CPMC.out"

echo_compile_mods:
	@echo "Compiling program modules"
