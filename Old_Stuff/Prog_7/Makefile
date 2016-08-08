FC= $(mpif90)
FC= $(f90)
FLAGS= $(FL) 
LF = $(Lflags)
LIBS= $(Libs)/Modules/modules_90.a \
      $(Libs)/MyEis/libeis.a \
      $(Libs)/MyNag/libnag.a \
      $(Libs)/MyLin/liblin.a \
      $(LIB_BLAS_LAPACK)

Hub:
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Hub  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) 

SPT:
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_SPT  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" LF="$(LF)" ) 

clean: 	
	(make -f Compile_Hub  clean );\
	(make -f Compile_SPT  clean );\
	rm *.mod *~ \#* 
