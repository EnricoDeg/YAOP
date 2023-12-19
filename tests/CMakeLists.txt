# GoogleTest requires at least C++14 / TKE requires C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# Do NOT install gtest with the project
option(INSTALL_GTEST "Enable installation of googletest." OFF)
FetchContent_MakeAvailable(googletest)

enable_testing()

# get_index_range
add_executable(
  get_index_range
  get_index_range.cpp
)

target_include_directories(get_index_range PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries (get_index_range yaop)

target_link_libraries(
  get_index_range
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(get_index_range)

if(ENABLE_CUDA OR ENABLE_HIP)

else()
    # calc_mxl_2
    add_executable(
      calc_mxl_2
      calc_mxl_2.cpp
    )
    target_include_directories(calc_mxl_2 PRIVATE ${PROJECT_SOURCE_DIR})
    target_include_directories(calc_mxl_2 PRIVATE ${PROJECT_SOURCE_DIR}/externals/mdspan/include)
    target_link_libraries (calc_mxl_2 yaop)
    target_link_libraries(
      calc_mxl_2
      GTest::gtest_main
    )
    include(GoogleTest)
    gtest_discover_tests(calc_mxl_2)

    # calc_diffusivity
    add_executable(
      calc_diffusivity
      calc_diffusivity.cpp
    )
    target_include_directories(calc_diffusivity PRIVATE ${PROJECT_SOURCE_DIR})
    target_include_directories(calc_diffusivity PRIVATE ${PROJECT_SOURCE_DIR}/externals/mdspan/include)
    target_link_libraries (calc_diffusivity yaop)
    target_link_libraries(
      calc_diffusivity
      GTest::gtest_main
    )
    include(GoogleTest)
    gtest_discover_tests(calc_diffusivity)

    # cpu_mdspan_impl
    add_executable(
      cpu_mdspan_impl
      cpu_mdspan_impl.cpp
    )
    target_include_directories(cpu_mdspan_impl PRIVATE ${PROJECT_SOURCE_DIR})
    target_include_directories(cpu_mdspan_impl PRIVATE ${PROJECT_SOURCE_DIR}/externals/mdspan/include)
    target_link_libraries (cpu_mdspan_impl yaop)
    target_link_libraries(
      cpu_mdspan_impl
      GTest::gtest_main
    )
    include(GoogleTest)
    gtest_discover_tests(cpu_mdspan_impl)

endif()
