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
    # calculate_mxl_2
    add_executable(
      calculate_mxl_2
      calculate_mxl_2.cpp
    )

    target_include_directories(calculate_mxl_2 PRIVATE ${PROJECT_SOURCE_DIR})
    target_include_directories(calculate_mxl_2 PRIVATE ${PROJECT_SOURCE_DIR}/externals/mdspan/include)
    target_link_libraries (calculate_mxl_2 yaop)

    target_link_libraries(
      calculate_mxl_2
      GTest::gtest_main
    )

    include(GoogleTest)
    gtest_discover_tests(calculate_mxl_2)
endif()
