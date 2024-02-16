cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

project(Task2)

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
 
# Link runTests with what we want to test and the GTest and pthread library
add_executable(runTests test.cpp dense_CSR.cpp)
target_link_libraries(runTests ${GTEST_LIBRARIES} pthread)