
if(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    if(ENABLE_CUDA)
        SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mallocatable=03 -acc=verystrict -Minfo=accel,inline -gpu=cc80")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -acc=verystrict -Minfo=accel,inline -gpu=cc80")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc=verystrict -Minfo=accel,inline -gpu=cc80")
    endif()
endif()

add_executable(example_basic_cxx example_basic.cpp)
target_link_libraries (example_basic_cxx tke)
target_include_directories(example_basic_cxx PRIVATE ${PROJECT_SOURCE_DIR})
if(ENABLE_CUDA)
    set_property(TARGET example_basic_cxx PROPERTY CUDA_SEPARABLE_COMPILATION ON)
endif()

if(ENABLE_C)
    add_executable(example_basic_c example_basic.c)
    target_link_libraries (example_basic_c tke)
    target_include_directories(example_basic_c PRIVATE ${PROJECT_SOURCE_DIR})
endif()

if(ENABLE_FORTRAN)
    set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
    add_executable(example_basic_f example_basic.f90)
    set_property(TARGET example_basic_f PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries (example_basic_f tke)
endif()


install (TARGETS
  example_basic_cxx # executables
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/examples)

if(ENABLE_C)
    install (TARGETS
      example_basic_c # executables
      RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/examples)
endif()

if(ENABLE_FORTRAN)
    install (TARGETS
      example_basic_f # executables
      RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/examples)
endif()
