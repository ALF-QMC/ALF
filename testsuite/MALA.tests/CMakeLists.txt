# Reuse the objects from the normal Make build, excluding the main program.
# Read its generated object list so registered Hamiltonians are linked as well.
set(PROG_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../Prog")
file(READ "${PROG_DIR}/.deps.mk" dependencies)
string(FIND "${dependencies}" "\n\n" objects_end)
string(SUBSTRING "${dependencies}" 0 ${objects_end} object_list)
string(REGEX MATCHALL "[A-Za-z0-9_./+-]+\\.o" objects "${object_list}")
list(REMOVE_ITEM objects main.o)
set(prog_objects)
foreach(object ${objects})
  list(APPEND prog_objects "${PROG_DIR}/${object}")
endforeach()

add_executable(mala-measurement-cache measurement_cache.F90)
target_compile_options(mala-measurement-cache PRIVATE ${FLAG_LIST})
target_include_directories(mala-measurement-cache PRIVATE
  "${PROG_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/../../Libraries/Modules")
target_link_libraries(mala-measurement-cache ${prog_objects} ${LIB_LIST})

find_program(PYTHON_EXECUTABLE NAMES python3)
if(NOT PYTHON_EXECUTABLE)
  message(FATAL_ERROR "MALA integration tests require Python 3")
endif()
enable_testing()
add_test(NAME mala-measurement-cache COMMAND mala-measurement-cache)
foreach(scenario measurements restart)
  add_test(NAME mala-${scenario}
    COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/integration.py"
      "${PROG_DIR}/ALF.out" ${scenario})
endforeach()
set_tests_properties(mala-measurement-cache mala-measurements mala-restart
  PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 120 LABELS "MALA")
