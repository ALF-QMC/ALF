FLAGS= -c   -O3 -std=f2003
SUFFIX= -qsuffix=f=f90
LF= 
LIB=modules_90.a

all: 
	(make -f Compile FC="$(FC)" FLAGS="$(FLAGS)"  LIB="$(LIB)" )

clean:
	(make -f Compile  clean ) ;\
 	rm -f *.mod *~ \#* 
