/*  Copyright (C) 2017 The ALF project
 
  This file is part of the ALF project.
 
     The ALF project is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.
 
     The ALF project is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
 
     You should have received a copy of the GNU General Public License
     along with Foobar.  If not, see http://www.gnu.org/licenses/.
     
     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
     
     - It is our hope that this program makes a contribution to the scientific community. Being
       part of that community we feel that it is reasonable to require you to give an attribution
       back to the original authors if you have benefitted from this program.
       Guidelines for a proper citation can be found on the project's homepage
       http://alf.physik.uni-wuerzburg.de .
       
     - We require the preservation of the above copyright notice and this license in all original files.
     
     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
 
     - If you make substantial changes to the program we require you to either consider contributing
       to the ALF project or to mark your material in a reasonable way as different from the original version.
*/

#include <inttypes.h>
#include <stdio.h>
#include <clBLAS.h>

//some global state. funny things will happen in multithreading
cl_context ctx = 0;
cl_command_queue queue = 0;

/*
 * CLErrString --
 *
 *      Utility function that converts an OpenCL status into a human
 *      readable string.
 *
 * Results:
 *      const char * pointer to a static string.
 */

static const char *
CLErrString(cl_int status) {
   static struct { cl_int code; const char *msg; } error_table[] = {
      { CL_SUCCESS, "success" },
      { CL_DEVICE_NOT_FOUND, "device not found", },
      { CL_DEVICE_NOT_AVAILABLE, "device not available", },
      { CL_COMPILER_NOT_AVAILABLE, "compiler not available", },
      { CL_MEM_OBJECT_ALLOCATION_FAILURE, "mem object allocation failure", },
      { CL_OUT_OF_RESOURCES, "out of resources", },
      { CL_OUT_OF_HOST_MEMORY, "out of host memory", },
      { CL_PROFILING_INFO_NOT_AVAILABLE, "profiling not available", },
      { CL_MEM_COPY_OVERLAP, "memcopy overlaps", },
      { CL_IMAGE_FORMAT_MISMATCH, "image format mismatch", },
      { CL_IMAGE_FORMAT_NOT_SUPPORTED, "image format not supported", },
      { CL_BUILD_PROGRAM_FAILURE, "build program failed", },
      { CL_MAP_FAILURE, "map failed", },
      { CL_INVALID_VALUE, "invalid value", },
      { CL_INVALID_DEVICE_TYPE, "invalid device type", },
      { 0, NULL },
   };
   static char unknown[25];
   int ii;

   for (ii = 0; error_table[ii].msg != NULL; ii++) {
      if (error_table[ii].code == status) {
         return error_table[ii].msg;
      }
   }

   snprintf(unknown, sizeof unknown, "unknown error %d", status);
   return unknown;
}

/*
 * PrintPlatform --
 *
 *      Dumps everything about the given platform ID.
 *
 * Results:
 *      void.
 */

static void
PrintPlatform(cl_platform_id platform) {
   static struct { cl_platform_info param; const char *name; } props[] = {
      { CL_PLATFORM_PROFILE, "profile" },
      { CL_PLATFORM_VERSION, "version" },
      { CL_PLATFORM_NAME, "name" },
      { CL_PLATFORM_VENDOR, "vendor" },
      { CL_PLATFORM_EXTENSIONS, "extensions" },
      { 0, NULL },
   };
   cl_device_id* deviceList;
   cl_uint numDevices;
   cl_int status;
   char buf[65536];
   size_t size;
   int ii;

   for (ii = 0; props[ii].name != NULL; ii++) {
      status = clGetPlatformInfo(platform, props[ii].param, sizeof buf, buf, &size);
      if (status != CL_SUCCESS) {
         fprintf(stderr, "platform[%p]: Unable to get %s: %s\n",
                 platform, props[ii].name, CLErrString(status));
         continue;
      }
      if (size > sizeof buf) {
         fprintf(stderr, "platform[%p]: Huge %s (%d bytes)!  Truncating to %d\n",
                 platform, props[ii].name, size, sizeof buf);
      }
      printf("platform[%p]: %s: %s\n", platform, props[ii].name, buf);
   }

   if ((status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL,
                                0, NULL, &numDevices)) != CL_SUCCESS) {
      fprintf(stderr, "platform[%p]: Unable to query the number of devices: %s\n",
              platform, CLErrString(status));
      return;
   }
   printf("platform[%p]: Found %d device(s).\n", platform, numDevices);

   deviceList = (cl_device_id*) malloc(numDevices * sizeof(cl_device_id));
   if ((status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL,
                                numDevices, deviceList, NULL)) != CL_SUCCESS) {
      fprintf(stderr, "platform[%p]: Unable to enumerate the devices: %s\n",
              platform, CLErrString(status));
      free(deviceList);
      return;
   }

//    for (ii = 0; ii < numDevices; ii++) {
//       PrintDevice(deviceList[ii]);
//    }

   free(deviceList);
}


/** Initialize OpenCL and clBLAS.
 * We also Initialize some global state variables.
 */

void initOpenCLandclBlas(int32_t* info)
{
    cl_int err;
    cl_uint numPlatforms;
    cl_platform_id platform = 0;
    cl_device_id device = 0;
    cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };
    /*Query the number of available platforms*/
    err = clGetPlatformIDs(0, NULL, &numPlatforms);
    printf("Found %d platform(s).\n", numPlatforms);

    cl_platform_id* platformList = (cl_platform_id*) malloc(sizeof(cl_platform_id) * numPlatforms);
    if ((err = clGetPlatformIDs(numPlatforms, platformList, NULL)) != CL_SUCCESS) {
       fprintf(stderr, "Unable to enumerate the platforms: %s\n",
               CLErrString(err));
       exit(1);
    }

    for (int ii = 0; ii < numPlatforms; ii++) {
       PrintPlatform(platformList[ii]);
    }

    free(platformList);
    /* Setup OpenCL environment. */
    err = clGetPlatformIDs(1, &platform, NULL);
    if (err != CL_SUCCESS) {
        printf( "clGetPlatformIDs() failed with %d\n", err );
        *info = 1;
        return;
    }

    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 1, &device, NULL);
    if (err != CL_SUCCESS) {
        printf( "clGetDeviceIDs() failed with %d\n", err );
        *info = 1;
        return;
    }

    props[1] = (cl_context_properties)platform;
    ctx = clCreateContext(props, 1, &device, NULL, NULL, &err);
    if (err != CL_SUCCESS) {
        printf( "clCreateContext() failed with %d\n", err );
        *info = 1;
        return;
    }
    queue = clCreateCommandQueue(ctx, device, 0, &err);
    if (err != CL_SUCCESS) {
        printf( "clCreateCommandQueue() failed with %d\n", err );
        clReleaseContext(ctx);
        *info = 1;
        return;
    }
    
    /* Setup clblas. */
    err = clblasSetup();
    if (err != CL_SUCCESS) {
        printf("clblasSetup() failed with %d\n", err);
        clReleaseCommandQueue(queue);
        clReleaseContext(ctx);
        *info = 1;
        return;
    }
}

void clalfzhemm(char* side, char* uplo, int32_t* m, int32_t* n, double* alpha, double* A, int32_t* lda, double* B, int32_t* ldb, double* beta, double* C, int32_t* ldc, int32_t* info)
{
    *info = 0;
    cl_int err;
    //map some arguments to clBLAS types
    const clblasOrder order = clblasColumnMajor;//FIXME: should be standard Fortran order
    const clblasSide zhemmside = (*side == 'R' ? clblasRight:clblasLeft);
    const clblasUplo zhemmuplo = (*uplo == 'U' ? clblasUpper:clblasLower);
    int ka = (*side == 'R' ? *n : *m);
    //FIXME: allow non-standard sizes of matrices...
    cl_event event = NULL;
/* Prepare OpenCL memory objects and place matrices inside them. */
    cl_mem bufA = clCreateBuffer(ctx, CL_MEM_READ_ONLY, 2*ka * ka * sizeof(*A), NULL, &err);
    cl_mem bufB = clCreateBuffer(ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 2 * *m * *n * sizeof(*B), B, &err);
    cl_mem bufC = clCreateBuffer(ctx, CL_MEM_READ_WRITE, 2 * *m * *n * sizeof(*C), NULL, &err);
    cl_double2 alp = {alpha[0], alpha[1]};
    cl_double2 bet = {beta[0], beta[1]};
    err = clEnqueueWriteBuffer(queue, bufA, CL_TRUE, 0,
        2 * ka * ka * sizeof(*A), A, 0, NULL, NULL);
//     err = clEnqueueWriteBuffer(queue, bufB, CL_TRUE, 0,
//         2 * *m * *n * sizeof(*B), B, 0, NULL, NULL);
    err = clEnqueueWriteBuffer(queue, bufC, CL_TRUE, 0,
        2 * *m * *n * sizeof(*C), C, 0, NULL, NULL);

    /* Call clblas function. */
    err = clblasZhemm(order, zhemmside, zhemmuplo, *m, *n, alp, bufA,
                         0, *lda, bufB, 0, *ldb, bet, bufC, 0, *ldc, 1, &queue,
                         0, NULL, &event);
    if (err != CL_SUCCESS) {
        printf("clblasSsymm() failed with %d\n", err);
        *info = 1;
    }
    else {
        /* Wait for calculations to be finished. */
        err = clWaitForEvents(1, &event);

        /* Fetch results of calculations from GPU memory. */
        err = clEnqueueReadBuffer(queue, bufC, CL_TRUE, 0, 2* *m * *n * sizeof(*C),
                                  C, 0, NULL, NULL);

        /* At this point you will get the result of SYMM placed in C array. */
    }
    *info += err;
}

/** Tidy and release our objects
 */
void teardown(int32_t* t)
{
    /* Finalize work with clblas. */
    clblasTeardown();

    /* Release OpenCL working objects. */
    clReleaseCommandQueue(queue);
    clReleaseContext(ctx);
}
