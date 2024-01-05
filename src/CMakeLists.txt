if(ENABLE_CUDA)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 80)
    endif()
    SET(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -lineinfo -DCUDA ")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCUDA")
    SET(SOURCE_CUDA backends/GPU/CUDA/cuda_check.cu)
endif()

if(ENABLE_HIP)
    # add architectures
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHIP")
    SET(SOURCE_HIP backends/GPU/HIP/cuda_check.cpp)
endif()

include_directories(${CMAKE_BINARY_DIR})

if(ENABLE_FORTRAN)
    set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
    set(SOURCE_FORTRAN bindings/Fortran/mod_YAOP.f90)
    set(SOURCE_C bindings/C/YAOP_capi.cpp)
endif()

if(ENABLE_C)
    set(SOURCE_C bindings/C/YAOP_capi.cpp)
endif()

if(ENABLE_CUDA)
    set_source_files_properties(YAOP.cpp PROPERTIES LANGUAGE CUDA)
elseif(ENABLE_HIP)
    set_source_files_properties(YAOP.cpp PROPERTIES LANGUAGE HIP)
endif()

set(SOURCE_EXE  YAOP.cpp
                ${SOURCE_C}
                ${SOURCE_FORTRAN}
                ${SOURCE_CUDA}
                ${SOURCE_HIP}
                shared/utils.cpp)

add_library(yaop SHARED ${SOURCE_EXE})
target_include_directories(yaop PRIVATE ${PROJECT_SOURCE_DIR})
target_include_directories(yaop PRIVATE ${PROJECT_SOURCE_DIR}/externals/mdspan/include)
set_property(TARGET yaop PROPERTY CXX_STANDARD 17)
if(ENABLE_CUDA)
    set_property(TARGET yaop PROPERTY CUDA_SEPARABLE_COMPILATION ON)
endif()

# Install the library
install (
  TARGETS yaop
  LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/lib
  ARCHIVE DESTINATION ${PROJECT_SOURCE_DIR}/lib)

# Install the library Fortran module
if(ENABLE_FORTRAN)
    install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
            DESTINATION ${PROJECT_SOURCE_DIR})
endif()

# Header Files
set(HEADERS_CXX YAOP.hpp)
install(FILES ${HEADERS_CXX} DESTINATION ${PROJECT_SOURCE_DIR}/include)
if(ENABLE_C)
    set(HEADERS_C bindings/C/YAOP.h)
    install(FILES ${HEADERS_C} DESTINATION ${PROJECT_SOURCE_DIR}/include)
endif()
