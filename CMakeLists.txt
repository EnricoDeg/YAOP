cmake_minimum_required(VERSION 3.0)

set(LANG_CXX CXX)

if(ENABLE_C)
    set(LANG_C C)
endif()

if(ENABLE_FORTRAN)
    set(LANG_C C)
    set(LANG_FORTRAN Fortran)
endif()

if(ENABLE_CUDA)
    set(LANG_CUDA CUDA)
endif()

if(ENABLE_HIP)
    message(FATAL_ERROR "HIP cmake configuration is work in progress")
    set(LANG_HIP HIP)
endif()

if(ENABLE_CUDA AND ENABLE_HIP)
    message(FATAL_ERROR "CUDA and HIP are both enabled: not allowed")
endif()

project(yaop ${LANG_CXX} ${LANG_C} ${LANG_FORTRAN} ${LANG_CUDA} ${LANG_HIP})

if(ENABLE_CUDA)
    find_package(CUDA QUIET REQUIRED)
endif()

add_subdirectory(src)

if(ENABLE_EXAMPLES)
    add_subdirectory(examples)
endif()

if(ENABLE_TESTS)
    add_subdirectory(tests)
endif()
