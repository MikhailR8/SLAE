cmake_minimum_required(VERSION 3.22)
project(libSLAE)

add_subdirectory(src)

if (WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()