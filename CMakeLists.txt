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
    src/profiling.cpp
    src/eos.cpp
    src/grid.cpp
    src/state.cpp
    src/flux.cpp
    src/rhs.cpp
    src/integrators.cpp
    scenarios/scenario.cpp
    scenarios/mesh.cpp
    scenarios/model_c7.cpp
    scenarios/model_column.cpp
    scenarios/model_flare.cpp
    scenarios/analytic_canopy.cpp
    scenarios/data_file_parser.cpp
    scenarios/pfss_field_line.cpp
)
target_include_directories(chromosphere PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_executable(chromo_main chromo_main.cpp)
target_link_libraries(chromo_main PRIVATE chromosphere)

# Stage-1 offline cross-check: compare every CRASH zAv row with the actual
# double/log-domain C++ Saha implementation (not a duplicated Python formula).
add_executable(eos_validate_saha util/eos/validate_saha.cpp)
target_link_libraries(eos_validate_saha PRIVATE chromosphere)

add_executable(eos_validate_production_table util/eos/validate_production_table.cpp)
target_link_libraries(eos_validate_production_table PRIVATE chromosphere)

add_executable(eos_validate_refined_table util/eos/validate_refined_table.cpp)
target_link_libraries(eos_validate_refined_table PRIVATE chromosphere)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(chromosphere PUBLIC OpenMP::OpenMP_CXX)
endif()

# Tests
enable_testing()
add_executable(chromo_tests tests/chromo_tests.cpp)
target_link_libraries(chromo_tests PRIVATE chromosphere)
add_test(NAME chromo_tests COMMAND chromo_tests)
set(EOS_PRODUCTION_TABLE ${CMAKE_CURRENT_SOURCE_DIR}/data/eos/gamma1_hydrogen_v1.dat)
add_test(NAME eos_production_table
         COMMAND eos_validate_production_table ${EOS_PRODUCTION_TABLE})
add_test(NAME eos_production_table_checksum
         COMMAND ${CMAKE_COMMAND}
                 -DROOT=${CMAKE_CURRENT_SOURCE_DIR}
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_eos_checksum.cmake)

foreach(GAMMA_GUARD IN ITEMS enable_te single_fluid ionization)
    add_test(NAME gamma_reject_${GAMMA_GUARD}
             COMMAND ${CMAKE_COMMAND}
                     -DEXE=$<TARGET_FILE:chromo_main>
                     -DTABLE=${EOS_PRODUCTION_TABLE}
                     -DGUARD=${GAMMA_GUARD}
                     -DOUT=${CMAKE_CURRENT_BINARY_DIR}/gamma_guard_unused.dat
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_gamma_guard.cmake)
endforeach()

add_test(NAME gamma_output_sidecar
         COMMAND ${CMAKE_COMMAND}
                 -DEXE=$<TARGET_FILE:chromo_main>
                 -DTABLE=${EOS_PRODUCTION_TABLE}
                 -DOUT=${CMAKE_CURRENT_BINARY_DIR}/gamma_output_test.dat
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_gamma_output.cmake)

set(PROJECT_PYTHON "${CMAKE_CURRENT_SOURCE_DIR}/.venv/bin/python")
if(EXISTS "${PROJECT_PYTHON}")
    add_test(NAME conductive_flux_python
             COMMAND ${PROJECT_PYTHON}
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_conductive_flux.py)
endif()
