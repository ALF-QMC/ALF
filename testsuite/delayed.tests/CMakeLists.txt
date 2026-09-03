cmake_minimum_required(VERSION 3.5...4.2)
Project(delayed.tests NONE)

enable_testing()

# Unlike the tests under Prog.tests this one drives the finished binary rather
# than a module, so it needs Prog/ALF.out and a build whose run directory looks
# like the plain one -- a TEMPERING binary expects Temp_* and would fail here
# for reasons of its own. It is therefore opt-in: the CI job that wants it sets
# ALF_TEST_DELAYED, and the ordinary matrix is left as it was.
if(DEFINED ENV{ALF_TEST_DELAYED})
  add_test(NAME 39-delayed-vs-immediate
           COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/delayed_vs_immediate.sh
                      ${CMAKE_CURRENT_SOURCE_DIR}/../../Prog/ALF.out
                      ${CMAKE_CURRENT_SOURCE_DIR}
                      ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(39-delayed-vs-immediate PROPERTIES SKIP_RETURN_CODE 77)
endif()