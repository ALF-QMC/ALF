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

/** Initialize OpenCL and clBLAS.
 * We also Initialize some global state variables.
 */

void initOpenCLandclBLas(int32_t* info)
{
    cl_int err;
    cl_platform_id platform = 0;
    cl_device_id device = 0;
    cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };
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

void clalfzhemm(int32_t* info)
{
}

/** Tidy and release our objects
 */
void teardown(int32_t* info)
{
    /* Finalize work with clblas. */
    clblasTeardown();

    /* Release OpenCL working objects. */
    clReleaseCommandQueue(queue);
    clReleaseContext(ctx);
}
