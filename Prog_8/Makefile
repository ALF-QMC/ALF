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

Hub_Ising:
	./Git_config.sh; \
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Hub_Ising  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" F90USEFULFLAGS="$(F90USEFULFLAGS)")

Hub:
	./Git_config.sh; \
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Hub  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" F90USEFULFLAGS="$(F90USEFULFLAGS)")

Hub_Can:
	./Git_config.sh; \
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Hub_Can  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" F90USEFULFLAGS="$(F90USEFULFLAGS)")

SPT:
	./Git_config.sh; \
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_SPT  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" LF="$(LF)" F90USEFULFLAGS="$(F90USEFULFLAGS)")

Ising:
	./Git_config.sh; \
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Ising  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" LF="$(LF)" F90USEFULFLAGS="$(F90USEFULFLAGS)")

Kondo_Honey:
	./Git_config.sh; \
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Kondo_Honey  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" LF="$(LF)" F90USEFULFLAGS="$(F90USEFULFLAGS)")

clean: 	
	(make -f Compile_Examples   clean );\
	(make -f Compile_Hub_Can   clean );\
	(make -f Compile_Hub_Ising   clean );\
	(make -f Compile_Hub   clean );\
	(make -f Compile_SPT   clean );\
	(make -f Compile_Ising clean );\
	(make -f Compile_Kondo_Honey clean );\
	rm *.mod *~ \#* 
