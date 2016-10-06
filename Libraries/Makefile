FC=$(f90)
FLAGS=$(FL)

all:
	(cd Modules;make FC=$(FC) FLAGS="$(FLAGS)") ;\
	(cd MyNag;make FC=$(FC) FLAGS="$(FLAGS)") 

cleanmodules:
	(cd Modules; make clean)

cleannag:
	(cd MyNag; make clean)

clean: cleanmodules cleannag

