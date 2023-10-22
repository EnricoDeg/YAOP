
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mallocatable=03 -acc=verystrict -Minfo=accel,inline -gpu=cc80,cuda12.2")

add_executable(example_basic_cxx example_basic.cpp)
target_link_libraries (example_basic_cxx tke)
target_include_directories(example_basic_cxx PRIVATE ${PROJECT_SOURCE_DIR}/src)

add_executable(example_basic_c example_basic.c)
target_link_libraries (example_basic_c tke)
target_include_directories(example_basic_c PRIVATE ${PROJECT_SOURCE_DIR}/src)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
add_executable(example_basic_f example_basic.f90)
set_property(TARGET example_basic_f PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries (example_basic_f tke)
#target_include_directories(example_basic_f PRIVATE ${CMAKE_BINARY_DIR})


install (TARGETS
  example_basic_cxx # executables
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/examples)

install (TARGETS
  example_basic_c # executables
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/examples)

install (TARGETS
  example_basic_f # executables
  RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/examples)

