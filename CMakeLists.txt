cmake_minimum_required(VERSION 3.0)

if(ENABLE_CUDA)
project(tke Fortran CXX C CUDA)
else()
project(tke Fortran CXX C)
endif()

add_subdirectory(src)
add_subdirectory(examples)
if(ENABLE_TESTS)
add_subdirectory(tests)
endif()
