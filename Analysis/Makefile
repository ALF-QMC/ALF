.PHONY : all tidy clean
BINS= cov_scal.out cov_tau.out cov_tau_ph.out Max_SAC.out cov_eq.out cov_mut.out ana.out
OBJS= cov_scal.o   cov_tau.o   cov_tau_ph.o   Max_SAC.o   cov_eq.o   cov_mut.o   ana.o
OBJS1=Predefined_Latt_mod.o ana_mod.o
MODS=predefined_lattices.mod ana_mod.mod


all: $(BINS)

Predefined_Latt_mod.o: ../Prog/Predefined_Latt_mod.F90
			$(ALF_FC) -c -o $@ $(ALF_FLAGS_ANA) $<

cov_tau_ph.o: cov_tau.F90
	$(ALF_FC) -c -o $@ $(ALF_FLAGS_ANA) -DPartHole $<

%.o: %.F90
	$(ALF_FC) -c -o $@ $(ALF_FLAGS_ANA) $<

%.out: $(OBJS1) %.o
	$(ALF_FC) -o $@ $^ $(ALF_LIB)

tidy:
	rm -f $(OBJS) $(OBJS1) $(MODS) cov_tau_ph.F90

clean: tidy
	rm -f $(BINS)
