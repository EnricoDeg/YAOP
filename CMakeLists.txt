cmake_minimum_required(VERSION 3.0)

set(LANG_CXX CXX)
set(LANG_C C)

if(ENABLE_FORTRAN)
set(LANG_FORTRAN Fortran)
endif()

if(ENABLE_CUDA)
set(LANG_CUDA CUDA)
endif()

project(tke ${LANG_CXX} ${LANG_C} ${LANG_FORTRAN} ${LANG_CUDA})

add_subdirectory(src)
add_subdirectory(examples)
if(ENABLE_TESTS)
add_subdirectory(tests)
endif()
