cmake_minimum_required(VERSION 2.8.12)
Project(libmscbdecomp.tests C Fortran)
find_package(LAPACK)

add_executable(1-ZeroDiag 1-ZeroDiag.F90)
add_executable(2-HomogeneousDiag 2-HomogeneousDiag.F90)
add_executable(3-TracelessDiag 3-TracelessDiag.F90)
add_executable(4-GeneralDiag 4-GeneralDiag.F90)
add_executable(5-EulerExp 5-EulerExponential.F90)
add_executable(6-FullExp 6-FullExponential.F90)

set (libmscbdecomptests 1-ZeroDiag 2-HomogeneousDiag 3-TracelessDiag 4-GeneralDiag 5-EulerExp 6-FullExp)

set_target_properties(${libmscbdecomptests} PROPERTIES COMPILE_FLAGS "-I ${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/libmscbdecomp/ -L ${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/libmscbdecomp/")

foreach(iter ${libmscbdecomptests})
target_link_libraries(${iter} ${LAPACK_LIBRARIES})
target_link_libraries(${iter} ${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/libmscbdecomp/libmscbdecomp.a)

endforeach(iter)

enable_testing()
add_test(1-ZeroDiag 1-ZeroDiag)
add_test(2-HomogeneousDiag 2-HomogeneousDiag)
add_test(3-TracelessDiag 3-TracelessDiag)
add_test(4-GeneralDiag 4-GeneralDiag)
add_test(5-EulerExp 5-EulerExp)
add_test(6-FullExp 6-FullExp)
