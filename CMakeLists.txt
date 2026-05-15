cmake_minimum_required(VERSION 3.12)
project(chromo_standalone CXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Vendored Armadillo (header-only path in include/).
# Built without LAPACK/BLAS: the solver only uses element-wise ops, sqrt, abs,
# min/max, and sub2ind, which are all header-only.
add_compile_definitions(ARMA_DONT_USE_LAPACK ARMA_DONT_USE_BLAS)

add_library(chromosphere STATIC
    chromosphere.cpp
    chromo_init.cpp
    util/chromo_util.cpp
)
target_include_directories(chromosphere PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_executable(chromo_main chromo_main.cpp)
target_link_libraries(chromo_main PRIVATE chromosphere)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(chromosphere PUBLIC OpenMP::OpenMP_CXX)
endif()

# Tests
enable_testing()
add_executable(chromo_tests tests/chromo_tests.cpp)
target_link_libraries(chromo_tests PRIVATE chromosphere)
add_test(NAME chromo_tests COMMAND chromo_tests)
