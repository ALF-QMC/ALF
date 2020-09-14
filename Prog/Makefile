.PHONY : Compile tidy clean Examples Z2_Slave Z2_Matter Hub_Can Kondo Hubbard tV Hubbard_Plain_Vanilla LRC

OBJS= Hamiltonians/LRC_mod.o  Set_random.o control_mod.o Fields_mod.o Operator_mod.o WaveFunction_mod.o observables_mod.o \
      Predefined_Int_mod.o Predefined_Obs_mod.o Predefined_Latt_mod.o Predefined_Hop_mod.o  Predefined_Trial_mod.o \
      $(HAMILTONIAN) QDRP_decompose_mod.o udv_state_mod.o Hop_mod.o UDV_WRAP_mod.o \
      Global_mod.o Wrapgr_mod.o tau_m_mod.o tau_p_mod.o main.o wrapul.o cgr1.o wrapur.o cgr2_2.o upgrade.o 
MODS= control.mod fields_mod.mod global_mod.mod hop_mod.mod lrc_mod.mod observables.mod \
      operator_mod.mod predefined_int.mod predefined_lattices.mod predefined_obs.mod \
      predefined_hoppings.mod  predefined_trial.mod qdrp_mod.mod tau_m_mod.mod tau_p_mod.mod \
      udv_state_mod.mod udv_wrap_mod.mod wavefunction_mod.mod wrapgr_mod.mod hamiltonian.mod
OHAM= Hamiltonian_Examples_mod.o Hamiltonians/Hamiltonian_Z2_slave_spins_mod.o Hamiltonians/Hamiltonian_Z2_Matter_mod.o \
      Hamiltonians/Hamiltonian_Hub_Canonical_mod.o Hamiltonians/Hamiltonian_Kondo_mod.o Hamiltonians/Hamiltonian_KN_Kondo_chain.o \
      Hamiltonians/Hamiltonian_Hubbard_mod.o Hamiltonians/Hamiltonian_Hubbard_Plain_Vanilla_mod.o  Hamiltonians/Hamiltonian_tV_mod.o \
      Hamiltonians/Hamiltonian_LRC_mod.o
BINS= Examples.out Z2_slave_spins.out Z2_Matter.out Hubb_Can.out Kondo.out KN_Kondo.out Hubbard.out Hubbard_Plain_Vanilla.out tV.out \
      LRC.out

      
all:
	(make tidy Examples)
	(make tidy Hubbard)
	(make tidy Hubbard_Plain_Vanilla)
	(make tidy tV)
	(make tidy Kondo)
	(make tidy LRC)
	(make tidy Z2_Slave)
	(make tidy Z2_Matter)
	(make tidy Hub_Can)
	(make tidy Kondo_KN)

Examples:
	./Git_config.sh
	(make Compile TARGET="Examples.out" HAMILTONIAN="Hamiltonian_Examples_mod.o")

Hubbard:
	./Git_config.sh
	(make Compile TARGET="Hubbard.out" HAMILTONIAN="Hamiltonians/Hamiltonian_Hubbard_mod.o")

Hubbard_Plain_Vanilla:
	./Git_config.sh
	(make Compile TARGET="Hubbard_Plain_Vanilla.out" HAMILTONIAN="Hamiltonians/Hamiltonian_Hubbard_Plain_Vanilla_mod.o")

tV:
	./Git_config.sh
	(make Compile TARGET="tV.out" HAMILTONIAN="Hamiltonians/Hamiltonian_tV_mod.o")
Kondo:
	./Git_config.sh
	(make Compile TARGET="Kondo.out" HAMILTONIAN="Hamiltonians/Hamiltonian_Kondo_mod.o")

LRC:
	./Git_config.sh
	(make Compile TARGET="LRC.out" HAMILTONIAN="Hamiltonians/Hamiltonian_LRC_mod.o")

Kondo_KN:
	./Git_config.sh
	(make Compile TARGET="KN_Kondo.out" HAMILTONIAN="Hamiltonians/Hamiltonian_KN_Kondo_chain.o")

Z2_Slave:
	./Git_config.sh
	(make Compile TARGET="Z2_slave_spins.out" HAMILTONIAN="Hamiltonians/Hamiltonian_Z2_slave_spins_mod.o")

Z2_Matter:
	./Git_config.sh
	(make Compile TARGET="Z2_Matter.out" HAMILTONIAN="Hamiltonians/Hamiltonian_Z2_Matter_mod.o")

Hub_Can:
	./Git_config.sh
	(make Compile TARGET="Hubb_Can.out" HAMILTONIAN="Hamiltonians/Hamiltonian_Hub_Canonical_mod.o")



Compile: $(OBJS)
	$(ALF_FC) -o $(TARGET) $(OBJS) $(ALF_LIB)

.SUFFIXES: .F90 .f
.f.o .F90.o:
	$(ALF_FC) -c -o $@ $(ALF_FLAGS_PROG) $<

tidy:
	rm -f $(OBJS) $(MODS) $(OHAM)

clean: tidy
	rm -f $(BINS)
