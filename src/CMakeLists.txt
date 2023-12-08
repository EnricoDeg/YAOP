
if(ENABLE_CUDA)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 80)
    endif()
    SET(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -lineinfo ")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCUDA")
endif()

include_directories(${CMAKE_BINARY_DIR})

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

if(ENABLE_CUDA)
    set(SOURCE_EXE  TKE.cpp
                    TKE_capi.cpp
                    mod_TKE.f90
                    TKE_backend.cpp
                    cuda_check.cu
                    cuda_kernels.cu
                    TKE_cuda.cu
                    TKE_cpu.cpp
                    utils.cpp
                    data_struct.cpp)
else()
    set(SOURCE_EXE  TKE.cpp
                    TKE_capi.cpp
                    mod_TKE.f90
                    TKE_backend.cpp
                    TKE_cpu.cpp
                    utils.cpp
                    data_struct.cpp)
endif()

add_library(tke SHARED ${SOURCE_EXE})
target_include_directories(tke PRIVATE ${PROJECT_SOURCE_DIR})
target_include_directories(tke PRIVATE ${PROJECT_SOURCE_DIR}/externals/mdspan/include)

# Install the library
install (
  TARGETS tke
  LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/lib
  ARCHIVE DESTINATION ${PROJECT_SOURCE_DIR}/lib)

# Install the library Fortran module
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
        DESTINATION ${PROJECT_SOURCE_DIR})

# ------------ Header Files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/
  DESTINATION ${PROJECT_SOURCE_DIR}/include
  FILES_MATCHING
  PATTERN "TKE.h*"
  )
