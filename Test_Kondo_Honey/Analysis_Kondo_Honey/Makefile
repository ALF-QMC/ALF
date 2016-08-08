FC= $(mpif90)
FC= $(f90)
FLAGS= $(FL) 
DIR1=$(HOME)
LIBS= $(DIR1)/lib_90/Modules/modules_90.a \
      $(DIR1)/lib_90/MyEis/libeis.a \
      $(DIR1)/lib_90/MyNag/libnag.a \
      $(DIR1)/lib_90/MyLin/liblin.a \
      $(DIR1)/lib_90/Ran/libran.a \
      $(LIB_BLAS_LAPACK)

all:
	cp $(DIR1)/lib_90/Modules/*.mod . ;\
	(make -f Compile_en  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) ;\
	(make -f Compile_eq  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) 
clean: 	
	(make -f Compile_eq  clean );\
	(make -f Compile_en  clean );\
	rm *.mod *~ \#* *.out
