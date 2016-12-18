FC= $(mpif90)
FC= $(f90)
FLAGS= $(FL) 
LIBS= $(Libs)/Modules/modules_90.a \
      $(Libs)/libqrref/libqrref.a  \
      $(LIB_BLAS_LAPACK)

all:
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_cov  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) ;\
	(make -f Compile_scal FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) ;\
	(make -f Compile_eq  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) 

clean: 	
	(make -f Compile_cov  clean );\
	(make -f Compile_scal clean );\
	(make -f Compile_eq  clean );\
	rm *.mod *~ \#* *.out
