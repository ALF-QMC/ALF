.PHONY : all cleanmodules cleanqrref clean

all:
	(cd Modules;  make )
	(cd libqrref; make )
	(cd libmscbdecomp; make )

cleanmodules:
	cd Modules; make clean

cleanmscbdecomp:
	cd libmscbdecomp; make clean

cleanqrref:
	cd libqrref; make clean

clean: cleanmodules cleanqrref cleanmscbdecomp
