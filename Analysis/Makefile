.PHONY : all tidy clean
BINS=cov_eq.out cov_scal.out cov_tau.out cov_tau_ph.out Max_SAC.out
OBJS=cov_eq.o   cov_scal.o   cov_tau.o   cov_tau_ph.o   Max_SAC.o  weights_corr_mod.o

all: $(BINS)

cov_tau_ph.F90 : cov_tau.F90
	cpp -w -P -DPartHole cov_tau.F90 -o cov_tau_ph.F90

weights_corr_mod.o : weights_corr_mod.F90
	$(ALF_FC) -c $(ALF_FLAGS_ANA) $<

%.out: %.F90 weights_corr_mod.o
	$(ALF_FC) -c $(ALF_FLAGS_ANA) $<
	$(ALF_FC) -o $@ $*.o $(ALF_LIB) weights_corr_mod.o

tidy:
	rm -f $(OBJS) cov_tau_ph.F90

clean: tidy
	rm -f $(BINS)
