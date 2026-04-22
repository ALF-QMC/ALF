int mkl_serv_intel_cpu_true(void) {
  return 1;
}
 
typedef int (*fakeintel_fptr)(void);
 
fakeintel_fptr mkl_serv_get_cpu_true(void) {
  return &mkl_serv_intel_cpu_true;
}
