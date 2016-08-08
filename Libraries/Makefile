FC=$(f90)
FLAGS=$(FL)

all:
	(cd Modules;make FC=$(FC) FLAGS="$(FLAGS)") ;\
	(cd MyNag;make FC=$(FC) FLAGS="$(FLAGS)") ;\
	(cd MyLin;make FC=$(FC) FLAGS="$(FLAGS)")

clean:
	(cd Modules;make clean);\
	(cd MyNag;make clean);\
	(cd MyLin;make clean)

