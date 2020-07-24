.PHONY : all lib ana Examples Z2_Slave Z2_Matter Hub_Can Kondo_Honey Hubbard

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
Z2_Slave: lib
	cd Prog && $(MAKE) Z2_Slave
Z2_Matter: lib
	cd Prog && $(MAKE) Z2_Matter
Hub_Can: lib
	cd Prog && $(MAKE) Hub_Can
Kondo_Honey: lib
	cd Prog && $(MAKE) Kondo_Honey

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
	@echo "Examples, Z2_Slave, Z2_Matter, Hub_Can, Kondo_Honey,  Hubbard"
