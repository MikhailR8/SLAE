cmake_minimum_required(VERSION 3.22)
project(test)

add_executable(test ../src/libSLAE/dense_CSR.cpp ../src/libSLAE/Iter_solvers.cpp 
../src/libSLAE/Test_speed_tools.cpp test.cpp)