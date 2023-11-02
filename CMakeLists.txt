cmake_minimum_required(VERSION 3.0)

project(tke Fortran CXX C CUDA)
add_subdirectory(src)
add_subdirectory(examples)
if(DEFINED ENABLE_TESTS)
add_subdirectory(tests)
endif(DEFINED ENABLE_TESTS)
