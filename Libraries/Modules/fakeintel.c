// This file is a fake implementation of the MKL function mkl_serv_get_cpu_true,
// which is used by MKL to determine if the CPU is an Intel CPU.
// This is used to trick MKL into thinking that it is running on an Intel CPU,
// which allows it to use the optimized code paths for Intel CPUs, even when running on an AMD CPU.
// This has shown to improve performance on AMD CPUs significantly, especially on AMD Zen CPUs.
// See: https://danieldk.eu/Software/Misc/Intel-MKL-on-AMD-Zen

int mkl_serv_intel_cpu_true(void) {
  return 1;
}
 
typedef int (*fakeintel_fptr)(void);
 
fakeintel_fptr mkl_serv_get_cpu_true(void) {
  return &mkl_serv_intel_cpu_true;
}
