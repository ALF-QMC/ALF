FC=gfortran
INCLUDES=-I ../../Libraries/Modules/ -I ../../Prog_8/
LIBS=../../Libraries/Modules/modules_90.a -llapack -lblas ../../Libraries/MyNag/libnag.a
FLAGS=-L ../../Libraries/Modules/
OBJS=10-Op-mmultL.bin 12-get-matrix-blocks.bin 14-Op-Wrapdo.bin 1-copy-sel-rows.bin 3-opmult.bin 5-FillExpOps.bin \
 6-opexpmult-2.bin 7-opexpmultct-2.bin 8-scalematrix-2.bin 11-Op-mmultR.bin  13-Op-Wrapup.bin 15-cgr.bin 2-copy-sel-columns.bin \
 4-opmultct.bin 6-opexpmult-1.bin 7-opexpmultct-1.bin 8-scalematrix-1.bin 9-Op-Phase.bin

%.bin : %.f90
	$(FC) $(FLAGS) $(INCLUDES)  $< $(LIBS) -o $@

15-cgr.bin: 15-cgr.f90
	$(FC) $(FLAGS) $(INCLUDES)  $< ../../Prog_8/cgr1.o ../../Prog_8/UDV_WRAP.o $(LIBS) -o $@

all: $(OBJS)

default: all

clean:
	rm $(OBJS)
