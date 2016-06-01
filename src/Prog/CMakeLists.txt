########################################
# Set up how to compile the source files
########################################

# Add the source files
SET(Common_src ${SRCPROG}/Hop_mod.f90
            ${SRCPROG}/Operator.f90
            ${SRCPROG}/UDV_WRAP.F90
            ${SRCPROG}/cgr1.f90
            ${SRCPROG}/cgr2.f90
            ${SRCPROG}/cgr2_1.f90
            ${SRCPROG}/cgr2_2.f90
            ${SRCPROG}/control_mod.F90
            ${SRCPROG}/gperp.F90
            ${SRCPROG}/inconfc.F90
            ${SRCPROG}/main.F90
            ${SRCPROG}/nranf.f90
            ${SRCPROG}/outconfc.F90
            ${SRCPROG}/print_bin_mod.F90
            ${SRCPROG}/tau_m.f90
            ${SRCPROG}/upgrade.f90
            ${SRCPROG}/wrapgrdo.f90
            ${SRCPROG}/wrapgrup.f90
            ${SRCPROG}/wrapul.f90
            ${SRCPROG}/wrapur.f90
)

#####################################
# Tell how to install this executable
#####################################

IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
#INSTALL(TARGETS ${PROGEXE} RUNTIME DESTINATION bin)
            
# Define the executables in terms of the source files
#Add in this loop further models
FOREACH(MODELTYPE Hub SPT Ising)
SET(${MODELTYPE}_src ${SRCPROG}/Hamiltonian_${MODELTYPE}.F90 ${Common_src})
ADD_EXECUTABLE(${PROGEXE}_${MODELTYPE} ${${MODELTYPE}_src})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################
TARGET_LINK_LIBRARIES(${PROGEXE}_${MODELTYPE} ${MODULES} ${MYEIS} ${MYNAG} ${MYLIN})

# Uncomment if you need to link to BLAS and LAPACK
TARGET_LINK_LIBRARIES(${PROGEXE}_${MODELTYPE} ${BLAS_LIBRARIES}
                                ${LAPACK_LIBRARIES}
                                ${CMAKE_THREAD_LIBS_INIT})

# Uncomment if you have parallization
IF(USE_OPENMP)
    SET_TARGET_PROPERTIES(${PROGEXE}_${MODELTYPE} PROPERTIES
                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ELSEIF(USE_MPI)
    SET_TARGET_PROPERTIES(${PROGEXE}_${MODELTYPE} PROPERTIES
                          COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
                          LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${PROGEXE}_${MODELTYPE} ${MPI_Fortran_LIBRARIES})
ENDIF(USE_OPENMP)
INSTALL(TARGETS ${PROGEXE}_${MODELTYPE} RUNTIME DESTINATION bin)
ENDFOREACH()
