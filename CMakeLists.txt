# CMake project file for General QMCT

##################################################
# Define the project and the depencies that it has
##################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(General_QMCT Fortran)

# Set the General QMCT version
SET(VERSION 1.0)

# Add our local modlues to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Uncomment if it is required that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
#OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING.  You should  review this file and make sure the flags 
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
# Locate and set parallelization libraries.  There are some CMake peculiarities
# taken care of here, such as the fact that the FindOpenMP routine doesn't know
# about Fortran.
INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)
# Setup the LAPACK libraries.  This also takes care of peculiarities, such as
# the fact the searching for MKL requires a C compiler, and that the results
# are not stored in the cache. 
INCLUDE(${CMAKE_MODULE_PATH}/SetUpLAPACK.cmake)

# There is an error in CMAKE with this flag for pgf90.  Unset it
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
IF(FCNAME STREQUAL "pgf90")
    UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")

############################################################
# Define the actual files and folders that make up the build
############################################################

# Define the executable name
SET(PROGEXE qmct)
SET(BLA_STATIC TRUE)

#SET(CMAKE_AR gcc-ar)
#SET(CMAKE_RANLIB gcc-ranlib)

# Define some directories
SET(SRC ${CMAKE_SOURCE_DIR}/src)
SET(LIB ${CMAKE_SOURCE_DIR}/lib)
SET(BIN ${CMAKE_SOURCE_DIR}/bin)

# Define the library names
SET(MODULES modules)
SET(MYEIS myeis)
SET(MYNAG mynag)
SET(MYLIN mylin)

SET(SRCMODULES ${SRC}/Modules)
SET(SRCEIS ${SRC}/MyEis)
SET(SRCNAG ${SRC}/MyNag)
SET(SRCLIN ${SRC}/MyLin)
SET(SRCPROG ${SRC}/Prog)

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LIB})

# The sources for the libraries and have it placed in the lib folder
ADD_SUBDIRECTORY(${SRCMODULES} ${LIB}/Modules)
ADD_SUBDIRECTORY(${SRCEIS} ${LIB}/MyEis)

ADD_SUBDIRECTORY(${SRCNAG} ${LIB}/MyNAG)
ADD_SUBDIRECTORY(${SRCLIN} ${LIB}/MyLin)

# The source for the FOO binary and have it placed in the bin folder
ADD_SUBDIRECTORY(${SRCPROG} ${BIN})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)
