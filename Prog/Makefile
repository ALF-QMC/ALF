FC= $(f90)
FLAGS= $(FL)
LF = $(Lflags)
LIBS= $(Libs)/Modules/modules_90.a \
      $(LIB_BLAS_LAPACK) \
      $(Libs)/libqrref/libqrref.a \
      -lz

Examples:
	./Git_config.sh; \
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Examples  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" F90USEFULFLAGS="$(F90USEFULFLAGS)")

Z2_Slave:
	./Git_config.sh; \
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Z2_slave  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" F90USEFULFLAGS="$(F90USEFULFLAGS)")

Hub_Ising:
	./Git_config.sh; \
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Hub_Ising  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" F90USEFULFLAGS="$(F90USEFULFLAGS)")

Hub_Can:
	./Git_config.sh; \
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Hub_Can  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" F90USEFULFLAGS="$(F90USEFULFLAGS)")

Nematic_Dirac:
	./Git_config.sh; \
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_NematicDirac  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" F90USEFULFLAGS="$(F90USEFULFLAGS)")

SPT:
	./Git_config.sh; \
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_SPT  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" LF="$(LF)" F90USEFULFLAGS="$(F90USEFULFLAGS)")

Z2-FLstar:
	./Git_config.sh; \
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Z2-FLstar  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" LF="$(LF)" F90USEFULFLAGS="$(F90USEFULFLAGS)")

Deconf-Kondo:
	./Git_config.sh; \
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_deconf-Kondo  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" LF="$(LF)" F90USEFULFLAGS="$(F90USEFULFLAGS)")

QSH:
	./Git_config.sh; \
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_QSH  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" LF="$(LF)" F90USEFULFLAGS="$(F90USEFULFLAGS)")


Kondo_Honey:
	./Git_config.sh; \
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Kondo_Honey  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" LF="$(LF)" F90USEFULFLAGS="$(F90USEFULFLAGS)")

clean: 	
	(make -f Compile_Examples   clean );\
	(make -f Compile_Hub_Can   clean );\
	(make -f Compile_Hub_Ising   clean );\
	(make -f Compile_SPT   clean );\
	(make -f Compile_Z2-FLstar   clean );\
	(make -f Compile_deconf-Kondo   clean );\
	(make -f Compile_Kondo_Honey clean );\
	(make -f Compile_Z2_slave clean );\
	(make -f Compile_NematicDirac clean );\
	(make -f Compile_QSH clean );\
	rm -f *.mod *~ \#* 
