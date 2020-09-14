.PHONY : all tidy clean
BINS= cov_scal.out cov_tau.out cov_tau_ph.out Max_SAC.out   cov_eq.out
OBJS= cov_scal.o   cov_tau.o   cov_tau_ph.o   Max_SAC.o   cov_eq.o
OBJS1 = Predefined_Latt_mod.o 
MODS=predefined_lattices.mod

all: $(OBJS1) $(BINS) 

cov_tau_ph.F90 : cov_tau.F90
	cpp -w -P -DPartHole cov_tau.F90 -o cov_tau_ph.F90

%.o: %.F90
	$(ALF_FC) -c $*.F90 $(ALF_FLAGS_ANA) 

%.out: %.F90
	$(ALF_FC) -c $(ALF_FLAGS_ANA) $<
	$(ALF_FC) -o $@ $*.o $(OBJS1) $(ALF_LIB)

tidy:
	rm -f $(OBJS) $(OBJS1) $(MODS) cov_tau_ph.F90

clean: tidy
	rm -f $(BINS)
