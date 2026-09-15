.PHONY : all tidy clean

SRCS = Max_SAC.F90 \
	ana.F90 \
	ana_hdf5.F90 \
	ana_mod.F90 \
	convert_latt.F90 \
	convert_local.F90 \
	convert_scal.F90 \
	cov_eq.F90 \
	cov_mut.F90 \
	cov_scal.F90 \
	cov_tau.F90 \
	maxent_wrapper_mod.F90

BINS= cov_scal.out cov_tau.out cov_tau_ph.out Max_SAC.out cov_eq.out cov_mut.out ana.out
BINS_H5=ana_hdf5.out convert_scal.out convert_latt.out convert_local.out
OBJS1=Predefined_Latt_mod.o ana_mod.o maxent_wrapper_mod.o

# Dependency helper variables
DEPSFILE := .deps.mk
DEPSGEN := ../Prog/gen_deps.py

all: $(DEPSFILE) $(BINS) $(if $(findstring DHDF5,$(ALF_FLAGS_ANA)),$(BINS_H5))

# OBJS and MODS are auto-generated into $(DEPSFILE) by gen_deps.py
-include $(DEPSFILE)
ana_mod.o: Predefined_Latt_mod.o
cov_eq.o: Predefined_Latt_mod.o

cov_tau_ph.o: cov_tau.F90
	$(ALF_FC) -c -o $@ $(ALF_FLAGS_ANA) -DPartHole $<

Predefined_Latt_mod.o: ../Prog/Predefined_Latt_mod.F90
	$(ALF_FC) -c -o $@ $(ALF_FLAGS_ANA) $<

%.o: %.F90
	$(ALF_FC) -c -o $@ $(ALF_FLAGS_ANA) $<

%.out: $(OBJS1) %.o
	$(ALF_FC) -o $@ $^ $(ALF_LIB)

$(DEPSFILE): $(SRCS) $(DEPSGEN)
	python3 $(DEPSGEN) $(SRCS) > $(DEPSFILE)

tidy: $(DEPSFILE)
	rm -f $(OBJS) cov_tau_ph.o $(MODS) $(DEPSFILE) Predefined_Latt_mod.o predefined_lattices.mod

clean: tidy
	rm -f $(BINS) $(BINS_H5)
