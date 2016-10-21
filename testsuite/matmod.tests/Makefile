FC=gfortran
INCLUDES=-I ../../Libraries/Modules/
LIBS=../../Libraries/Modules/modules_90.a -llapack -lblas ../../Libraries/MyNag/libnag.a
FLAGS=-L ../../Libraries/Modules/
OBJS=1-diag-real-matrix.bin 2-inv-r0.bin 3-inv-R-variable.bin 4-inv-R-variable-1.bin 5-inv-R1.bin 6-inv_C.bin 7-inv_C_variable.bin \
     8-inv-C1.bin 9-inv-C1_test2.bin 10-initd_C.bin

%.bin : %.f90
	$(FC) $(FLAGS) $(INCLUDES)  $< $(LIBS) -o $@

all: $(OBJS)

clean:
	rm $(OBJS)
