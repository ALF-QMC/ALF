.PHONY : all clean

all: gpu_linalg_bench.out

gpu_linalg_bench.out: gpu_linalg_bench.F90
	$(ALF_FC) $(ALF_FLAGS_PROG) -o $@ $< $(ALF_LIB)

clean:
	rm -f gpu_linalg_bench.out *.o *.mod
