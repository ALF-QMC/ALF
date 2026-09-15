cmake_minimum_required(VERSION 3.5...4.2)
Project(Analysis.tests C Fortran)

add_executable(1-maxent-kernel-overflow 1-maxent-kernel-overflow.F90)

set(analysistests 1-maxent-kernel-overflow)

foreach(iter ${analysistests})
  foreach(iflag ${FLAG_LIST})
    target_compile_options(${iter} PRIVATE ${iflag})
  endforeach(iflag)
  target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Analysis/maxent_wrapper_mod.o)
  target_link_libraries(${iter} ${LIB_LIST})
  target_include_directories(${iter} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../Analysis/)
  target_include_directories(${iter} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules/)
endforeach(iter)

enable_testing()
add_test(1-maxent-kernel-overflow 1-maxent-kernel-overflow)
