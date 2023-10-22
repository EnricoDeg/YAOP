
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 80)
endif()

include_directories(${CMAKE_BINARY_DIR})

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

set(SOURCE_EXE  TKE.cpp
                TKE_capi.cpp
                mod_TKE.f90
                TKE_backend.cpp
                cuda_check.cu
                TKE_cuda.cu
                utils.cpp)

add_library(tke SHARED ${SOURCE_EXE})

#add_executable(testf test.f90)
#set_property(TARGET testf PROPERTY LINKER_LANGUAGE Fortran)
#target_link_libraries (testf tke)

#add_executable(testcxx test.cpp)
#target_link_libraries (testcxx tke)

# Install the library
install (
  TARGETS tke
  LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/lib
  ARCHIVE DESTINATION ${PROJECT_SOURCE_DIR}/lib)

# Install the library Fortran module
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} 
        DESTINATION ${PROJECT_SOURCE_DIR})

# Install executables.
#install (TARGETS
#  testf testcxx # executables
#  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)

# ------------ Header Files	
install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/
  DESTINATION ${PROJECT_SOURCE_DIR}/include
  FILES_MATCHING
  PATTERN "*.h*"
  )
