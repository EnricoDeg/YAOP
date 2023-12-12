
if(ENABLE_CUDA)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 80)
    endif()
    SET(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -lineinfo ")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCUDA")
endif()

include_directories(${CMAKE_BINARY_DIR})

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

if(ENABLE_FORTRAN)
    set(SOURCE_FORTRAN bindings/Fortran/mod_TKE.f90)
endif()

if(ENABLE_CUDA)
    set(SOURCE_EXE  TKE.cpp
                    bindings/C/TKE_capi.cpp
                    ${SOURCE_FORTRAN}
                    backends/TKE_backend.cpp
                    backends/GPU/CUDA/cuda_check.cu
                    backends/GPU/gpu_kernels.cpp
                    backends/GPU/TKE_gpu.cpp
                    backends/CPU/TKE_cpu.cpp
                    shared/utils.cpp
                    shared/interface/data_struct.cpp)
    set_source_files_properties(backends/GPU/gpu_kernels.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(backends/GPU/TKE_gpu.cpp PROPERTIES LANGUAGE CUDA)
else()
    set(SOURCE_EXE  TKE.cpp
                    bindings/C/TKE_capi.cpp
                    ${SOURCE_FORTRAN}
                    backends/TKE_backend.cpp
                    backends/CPU/TKE_cpu.cpp
                    shared/utils.cpp
                    shared/interface/data_struct.cpp)
endif()

add_library(tke SHARED ${SOURCE_EXE})
target_include_directories(tke PRIVATE ${PROJECT_SOURCE_DIR})
target_include_directories(tke PRIVATE ${PROJECT_SOURCE_DIR}/externals/mdspan/include)
set_property(TARGET tke PROPERTY CXX_STANDARD 17)

# Install the library
install (
  TARGETS tke
  LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/lib
  ARCHIVE DESTINATION ${PROJECT_SOURCE_DIR}/lib)

# Install the library Fortran module
if(ENABLE_FORTRAN)
    install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
            DESTINATION ${PROJECT_SOURCE_DIR})
endif()

# ------------ Header Files
set(HEADERS_CXX TKE.hpp)
set(HEADERS_C bindings/C/TKE.h)
install(FILES ${HEADERS_CXX} DESTINATION ${PROJECT_SOURCE_DIR}/include)
install(FILES ${HEADERS_C} DESTINATION ${PROJECT_SOURCE_DIR}/include)
