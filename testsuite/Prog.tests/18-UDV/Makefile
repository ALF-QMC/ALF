FC= $(mpif90)
FC= $(f90)
FLAGS= $(FL) -std=f2003
LF = $(Lflags) -std=f2003
LIBS= $(Libs)/Modules/modules_90.a \
      $(Libs)/MyNag/libnag.a \
      $(LIB_BLAS_LAPACK)

all:
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) 

clean: 	
	(make -f Compile   clean );\
	rm *.mod *~ \#* 
