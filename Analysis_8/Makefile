FC= $(mpif90)
FC= $(f90)
FLAGS= $(FL) 
LIBS= $(Libs)/Modules/modules_90.a \
      $(Libs)/MyNag/libnag.a \
      $(Libs)/MyLin/liblin.a \
      $(LIB_BLAS_LAPACK)

all:
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_en  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) ;\
	(make -f Compile_cov  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) ;\
	(make -f Compile_Max  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) ;\
	(make -f Compile_Max_ph  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) ;\
	(make -f Compile_eq  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) 
clean: 	
	(make -f Compile_eq  clean );\
	(make -f Compile_cov  clean );\
	(make -f Compile_Max  clean );\
	(make -f Compile_Max_ph  clean );\
	(make -f Compile_en  clean );\
	rm *.mod *~ \#* *.out
