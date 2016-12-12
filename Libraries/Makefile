FC=$(f90)
FLAGS=$(FL)

all:
	(cd Modules;make FC=$(FC) FLAGS="$(FLAGS)")
	(cd libqrref; make FC=$(FC) FLAGS="$(FLAGS)")
#	(cd MyNag;make FC=$(FC) FLAGS="$(FLAGS)") 

cleanmodules:
	(cd Modules; make clean; cd ..)

#cleannag:
#	(cd MyNag; make clean)

cleanqrref:
	(cd libqrref; make clean)

clean: cleanmodules cleanqrref

