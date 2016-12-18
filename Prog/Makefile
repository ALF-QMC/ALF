FC= $(mpif90)
FC= $(f90)
FLAGS= $(FL) -std=f2003
LF = $(Lflags) -std=f2003
LIBS= $(Libs)/Modules/modules_90.a \
      $(LIB_BLAS_LAPACK) \
      $(Libs)/libqrref/libqrref.a

Examples:
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Examples  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) 

clean: 	
	(make -f Compile_Examples   clean );\
	rm *.mod *~ \#* 
