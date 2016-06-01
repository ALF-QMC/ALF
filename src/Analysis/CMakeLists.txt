add_executable(jackv5 jackv5.f90)
add_executable(cov_eq cov_eq.F90)
add_executable(cov_tau cov_tau.f90)

foreach(analysisbinary jackv5 cov_eq cov_tau)
target_link_libraries(${analysisbinary} ${MODULES} ${MYEIS} ${MYNAG} ${MYLIN})
target_link_libraries(${analysisbinary} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endforeach()
