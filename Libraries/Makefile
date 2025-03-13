.PHONY : all cleanmodules cleanqrref cleanpfapack clean

all:
	(cd Modules;  make )
	(cd libqrref; make )
	(cd pfapack; make )

cleanmodules:
	cd Modules; make clean

cleanqrref:
	cd libqrref; make clean

cleanpfapack:
	cd pfapack; make clean

clean: cleanmodules cleanqrref cleanpfapack
