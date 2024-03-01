cmake_minimum_required(VERSION 3.22)
project(test_speed_dense_CSR)

add_executable(speed dense_CSR.cpp simple_iter.cpp testspeed.cpp)
add_executable(iter dense_CSR.cpp simple_iter.cpp testiter.cpp)
