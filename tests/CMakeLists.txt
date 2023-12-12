# GoogleTest requires at least C++14
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

add_executable(
  get_index_range
  get_index_range.cpp
)

target_include_directories(get_index_range PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries (get_index_range tke)

target_link_libraries(
  get_index_range
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(get_index_range)
