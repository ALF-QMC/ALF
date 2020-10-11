.PHONY : all lib ana Z2_Matter Kondo Hubbard tV Hubbard_Plain_Vanilla LRC
# To add a new Hamiltonain named My_New_Hamiltonian, add "My_New_Hamiltonian"  to the above list

all: lib ana
	cd Prog && $(MAKE) all

lib:
	cd Libraries && $(MAKE)
ana: lib
	cd Analysis && $(MAKE)
Hubbard: lib
	cd Prog && $(MAKE) Hubbard
Hubbard_Plain_Vanilla: lib
	cd Prog && $(MAKE) Hubbard_Plain_Vanilla
tV: lib
	cd Prog && $(MAKE) tV
Kondo: lib
	cd Prog && $(MAKE) Kondo
LRC: lib
	cd Prog && $(MAKE) LRC
Z2_Matter: lib
	cd Prog && $(MAKE) Z2_Matter

# To add a new Hamiltonain named My_New_Hamiltonian  uncomment the lines below.
# My_New_Hamiltonian: lib
#	cd Prog && $(MAKE) My_New_Hamiltonian

.PHONY : clean cleanall cleanprog cleanlib cleanana tidy tidyprog tidyana help
clean: cleanall
cleanall: cleanprog cleanlib cleanana
tidy: cleanlib tidyana tidyprog
cleanprog:
	cd Prog && $(MAKE) clean 
cleanlib:
	cd Libraries && $(MAKE) clean
cleanana:
	cd Analysis && $(MAKE) clean
tidyana:
	cd Analysis && $(MAKE) tidy
tidyprog:
	cd Prog && $(MAKE) tidy

help:
	@echo "The following are some of the valid targets of this Makefile"
	@echo "lib, ana, clean, cleanall, cleanprog, cleanlib, cleanana"
	@echo "Z2_Matter, Kondo,  Hubbard, Hubbard_Plain_Vanilla, tV"
