//compile using:
/*
 * gcc -std=c99 -O -I /brokendisk/lib/clBLAS/include/ -c clzhemm.c  -L /brokendisk/lib/clBLAS/lib64/ -lclBLAS -lOpenCL
 * $ g++ -O -I /brokendisk/lib/clBLAS/include/ example_zhemm.cpp clzhemm.o -L /brokendisk/lib/clBLAS/lib64/ -lclBLAS -lOpenCL -o out.exe
 * */



#include <sys/types.h>
#include <stdio.h>
#include <string.h>

#include <CL/cl.h>

/* This example uses predefined matrices and their characteristics for
 * simplicity purpose.
 */

#define M  4
#define N  3

static const cl_double2 alpha = {{10, 10}};

static const cl_double2 A[M*M] = {
    {{11, 12}}, {{-1, -1}}, {{-1, -1}}, {{-1, -1}},
    {{21, 22}}, {{22, 23}}, {{-1, -1}}, {{-1, -1}},
    {{31, 32}}, {{32, 33}}, {{33, 34}}, {{-1, -1}},
    {{41, 61}}, {{42, 62}}, {{43, 73}}, {{44, 23}}
};
static const size_t lda = M;

static const cl_double2 B[M*N] = {
    {{11, -21}},  {{-12, 23}}, {{13, 33}},
    {{21, 12}},   {{22, -10}}, {{23, 5}},
    {{31, 1}},    {{-32, 65}}, {{33, -1}},
    {{1, 41}},    {{-33, 42}}, {{12, 43}},
};
static const size_t ldb = N;

static const cl_double2 beta = {{20, 20}};

static cl_double2 C[M*N] = {
    {{11, 11}},  {{-12, 12}}, {{13, 33}},
    {{21, -32}}, {{22,  -1}}, {{23, 0}},
    {{31, 13}},  {{-32, 78}}, {{33, 45}},
    {{41, 14}},  {{0,   42}}, {{43, -1}},
};
static const size_t ldc = N;

static void
printResult(void)
{
    size_t i, j, nrows;

    printf("Result:\n");

    nrows = (sizeof(C) / sizeof(cl_double2)) / ldc;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ldc; j++) {
            printf("<%9.2f, %-9.2f> ", *((double*)(C + i * ldc + j)), *(((double*)(C + i * ldc + j)) + 1));
        }
        printf("\n");
    }
}

extern "C" {
void initopenclandclblas(int32_t* info);
void clalfzhemm(char* side, char* uplo, int32_t* m, int32_t* n, double* alpha, double* A, int32_t* lda, double* B, int32_t* ldb, double* beta, double* C, int32_t* ldc, int32_t* info);
void teardown(int32_t* t);
}

int main(void)
{
    int32_t res;
    initopenclandclblas(&res);
    char side = 'L';
    char uplo = 'L';
    int m = M;
    int n = N;
    clalfzhemm(&side, &uplo, &m, &n, (double*)&alpha, (double*)A, &m, (double*)B, &n, (double*)&beta, (double*)C, &n, &res);
            printResult();
teardown(&res);
}
