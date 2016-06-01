#Modules library
SET(Modules_src ${SRCMODULES}/Files_mod.f90 ${SRCMODULES}/Histogram.f90 ${SRCMODULES}/Histogram_v2.f90 ${SRCMODULES}/Natural_constants.f90 ${SRCMODULES}/Random_Wrap.f90 ${SRCMODULES}/errors.f90 ${SRCMODULES}/fourier.f90 ${SRCMODULES}/lattices_v3.f90
${SRCMODULES}/log_mesh.f90 ${SRCMODULES}/mat_mod.f90 ${SRCMODULES}/matrix.f90 ${SRCMODULES}/maxent.f90 ${SRCMODULES}/maxent_stoch.f90 ${SRCMODULES}/maxent_stoch.f90 ${SRCMODULES}/precdef.mod.f90 ${SRCMODULES}/smooth_stoch.f90 ${SRCMODULES}/tmp.f90)
ADD_LIBRARY(${MODULES} STATIC ${Modules_src})
