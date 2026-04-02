.PHONY : Compile tidy clean program

# Source files (single source of truth; OBJS and MODS are generated into DEPSFILE)
SRCS = Hamiltonians/LRC_mod.F90 control_mod.F90 Fields_mod.F90 Operator_mod.F90 WaveFunction_mod.F90 observables_mod.F90 \
	ContainerElementBase_mod.F90 DynamicMatrixArray_mod.F90 OpTTypes_mod.F90 OpT_time_dependent.F90 \
	Hamiltonian_main_mod.F90 QDRP_decompose_mod.F90 udv_state_mod.F90 Hop_mod.F90 UDV_WRAP_mod.F90 \
	Predefined_Int_mod.F90 Predefined_Obs_mod.F90 Predefined_Latt_mod.F90 Predefined_Hop_mod.F90 Predefined_Trial_mod.F90 \
	wrapul_mod.F90 cgr1_mod.F90 wrapur_mod.F90 cgr2_2_mod.F90 upgrade_mod.F90 Set_random_mod.F90 \
	Global_mod.F90 Langevin_HMC_mod.F90 Wrapgr_mod.F90 tau_m_mod.F90 tau_p_mod.F90 QMC_runtime_var_mod.F90 \
	main.F90

# Dependency helper variables
DEPSFILE := .deps.mk
DEPSGEN := ./gen_deps.py
HAM_LIST_D := $(wildcard Hamiltonians.list.d/*)

all: program

# OBJS and MODS are auto-generated into $(DEPSFILE) by gen_deps.py
-include $(DEPSFILE)

program: ALF.out

ALF.out: $(DEPSFILE) $(OBJS)
	@echo "Link program"
	$(ALF_FC) -o "ALF.out" $(OBJS) $(ALF_LIB)

OHAM := $(shell ./parse_ham.py --print_obj_list 2>/dev/null)

Hamiltonian_main_mod.o: Hamiltonians_interface.h Hamiltonians_case.h

$(OHAM): .parse_parameters Hamiltonian_main_mod.o

main.o: git.h

# Regenerate git.h whenever the commit hash, branch, or dirty state changes.
# .git/HEAD changes on branch switch; .git/index changes on commit/stage.
# $(SRCS) catches unstaged modifications that flip the -dirty suffix.
git.h: $(wildcard ../.git/HEAD ../.git/index) $(SRCS)
	./Git_config.sh

.parse_parameters: Hamiltonians.list $(HAM_LIST_D)
	@echo "Parsing Hamiltonian parameters"
	./parse_ham.py --create_read_write_par
	@touch $@

Hamiltonians_interface.h: Hamiltonians.list $(HAM_LIST_D)
	./parse_ham.py --create_hamiltonians_interface

Hamiltonians_case.h: Hamiltonians.list $(HAM_LIST_D)
	./parse_ham.py --create_hamiltonians_case

# The order-only prerequisite (after |) guarantees .echo_compile_mods is built
# before any .o file under parallel make (-j N), so the banner appears first.
# Because it is order-only, Make never timestamps .echo_compile_mods against
# the .o files, so it never triggers recompilation.
%.o: %.F90 | .echo_compile_mods
	$(ALF_FC) -c $(ALF_FLAGS_PROG) -o $@ $<

# Generate dependency file that maps object files to the object(s) providing
# the Fortran modules they USE. This enables correct parallel builds.
$(DEPSFILE): $(SRCS) $(DEPSGEN) Hamiltonians.list $(HAM_LIST_D) $(shell ./parse_ham.py --print_src_list)
	python3 $(DEPSGEN) $(SRCS) $(shell ./parse_ham.py --print_src_list) > $(DEPSFILE)

tidy: $(DEPSFILE)
	rm -f $(OBJS) $(OHAM) $(MODS) hamiltonian*.mod hamiltonian*.smod $(DEPSFILE) .echo_compile_mods
	rm -f git.h .parse_parameters Hamiltonians_interface.h Hamiltonians_case.h Hamiltonians/Hamiltonian_*_read_write_parameters.F90

clean: tidy
	rm -f "ALF.out"

# Stamp file: printed once when the dependency file is (re-)generated, then
# touched so Make can track whether the echo has already fired.  Deleted by
# tidy/clean so the banner reappears on the next fresh build.
.echo_compile_mods: $(DEPSFILE)
	@echo "Compiling program modules"
	@touch $@
