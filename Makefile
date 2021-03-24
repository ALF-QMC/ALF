.PHONY : all lib ana program all
# To add a new Hamiltonain named My_New_Hamiltonian, add "My_New_Hamiltonian"  to the above list

all: lib ana program

lib:
	cd Libraries && $(MAKE)
ana: lib
	cd Analysis && $(MAKE)
program: lib
	cd Prog && $(MAKE) program

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
