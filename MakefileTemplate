.PHONY : all lib ana program test_compiler_set


all: lib ana program

lib: test_compiler_set
	@echo "Compiling Libraries"
	cd Libraries && $(MAKE)
ana: lib
	@echo "Compiling Analysis"
	cd Analysis && $(MAKE)
program: lib
	@echo "Compiling Program"
	cd Prog && $(MAKE) program

.PHONY : clean cleanall cleanprog cleanlib cleanana tidy tidyprog tidyana help
clean: cleanall
cleanall: cleanprog cleanlib cleanana
tidy: cleanlib tidyana tidyprog
cleanprog:
	@echo "Cleaning up Prog/"
	cd Prog && $(MAKE) clean 
cleanlib:
	@echo "Cleaning up Libraries/"
	cd Libraries && $(MAKE) clean
cleanana:
	@echo "Cleaning up Analysis/"
	cd Analysis && $(MAKE) clean
tidyana:
	cd Analysis && $(MAKE) tidy
tidyprog:
	cd Prog && $(MAKE) tidy

help:
	@echo "The following are some of the valid targets of this Makefile"
	@echo "all, lib, ana, program, clean, cleanall, cleanprog, cleanlib, cleanana"


# Test whether enviroment variable ALF_FC is set and is a valid command.
test_compiler_set:
	@if [ -z ${ALF_FC} ]; then \
	  printf "\n\033[0;31m Environment variable ALF_FC not set.\n"; \
	  printf " Please source configure.sh before compilation.\033[0m\n\n"; \
	  exit 1; \
	fi
	@if [ ! `command -v ${ALF_FC}` ]; then \
	  printf "\n\033[0;31m Compiler '${ALF_FC}' does not exist.\n"; \
	  printf " Please install or choose other in configure.sh.\033[0m\n"; \
	  exit 1; \
	fi
