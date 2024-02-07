cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

project(Run_through)

# add_executable(run_through src/run_through_3_diag_matrix.cpp)

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
 
# Link runTests with what we want to test and the GTest and pthread library
add_executable(runTests src/tests.cpp)
target_link_libraries(runTests ${GTEST_LIBRARIES} pthread)