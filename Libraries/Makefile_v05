FC=$(f90)
FLAGS=$(FL)

all:
	(cd Modules;make -f Makefile_v05 FC=$(FC) FLAGS="$(FLAGS)")
	(cd libqrref; make FC=$(FC) FLAGS="$(FLAGS)")

cleanmodules:
	(cd Modules; make -f Makefile_v05 clean; cd ..)

cleanqrref:
	(cd libqrref; make clean)

clean: cleanmodules cleanqrref
