cmake_minimum_required(VERSION 2.8.12)
Project(matmod.tests C Fortran)

add_subdirectory(matmod.tests)
add_subdirectory(Prog.tests)
add_subdirectory(libmscbdecomp_mod)
enable_testing()
