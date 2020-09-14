.PHONY : all lib ana Examples Z2_Slave Z2_Matter Hub_Can Kondo Hubbard tV Hubbard_Plain_Vanilla LRC

all: lib ana
	cd Prog && $(MAKE) all

lib:
	cd Libraries && $(MAKE)
ana: lib
	cd Analysis && $(MAKE)
Examples: lib
	cd Prog && $(MAKE) Examples
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
Z2_Slave: lib
	cd Prog && $(MAKE) Z2_Slave
Z2_Matter: lib
	cd Prog && $(MAKE) Z2_Matter
Hub_Can: lib
	cd Prog && $(MAKE) Hub_Can

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
	@echo "Examples, Z2_Slave, Z2_Matter, Hub_Can, Kondo,  Hubbard, Hubbard_Plain_Vanilla, tV"
