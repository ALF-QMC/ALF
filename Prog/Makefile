FC= $(f90)
FLAGS= $(FL)
LF = $(Lflags)
LIBS= $(Libs)/Modules/modules_90.a \
      $(LIB_BLAS_LAPACK) \
      $(Libs)/libqrref/libqrref.a

Examples:
	./Git_config.sh; \
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Examples  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" F90USEFULFLAGS="$(F90USEFULFLAGS)")

clean: 	
	(make -f Compile_Examples   clean );\
	rm -f *.mod *~ \#* 
