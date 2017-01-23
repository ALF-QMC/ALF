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

Hub_Ising:
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Hub_Ising  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) 

Hub:
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Hub  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) 

Hub_Can:
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Hub_Can  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) 

SPT:
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_SPT  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" LF="$(LF)" ) 

Ising:
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Ising  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" LF="$(LF)" ) 

Kondo_Honey:
	cp $(Libs)/Modules/*.mod . ;\
	(make -f Compile_Kondo_Honey  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" LF="$(LF)" ) 

clean: 	
	(make -f Compile_Examples   clean );\
	(make -f Compile_Hub_Ising   clean );\
	(make -f Compile_Hub   clean );\
	(make -f Compile_SPT   clean );\
	(make -f Compile_Ising clean );\
	(make -f Compile_Kondo_Honey clean );\
	rm *.mod *~ \#* 
