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

option(CHROMO_ENABLE_OPENMP "Enable OpenMP parallel execution" OFF)

# The release stores the conserved state (and the mesh metrics that share
# chromosphere::Vec) in DOUBLE. That is the default and the only supported
# production configuration; it matches SWMF/BATS-R-US, whose shipped builds all
# set PRECISION = ${DOUBLEPREC}. See docs/studies/numerics/state_precision_release_cutover.md.
#
# CHROMO_STATE_FLOAT32 reverts to the pre-cutover float storage. It exists ONLY
# to reproduce the float32 round-off sensitivity as a controlled diagnostic; it
# is not exercised by scripts/release_validation.sh and no release number may be
# quoted from it.
option(CHROMO_STATE_FLOAT32 "DIAGNOSTIC: store the conserved state in float (non-production)" OFF)
if(CHROMO_STATE_FLOAT32)
    add_compile_definitions(CHROMO_STATE_FLOAT32=1)
    message(WARNING
        "CHROMO_STATE_FLOAT32=ON — DIAGNOSTIC build. The release stores the "
        "conserved state in double; do not quote release numbers from this tree.")
endif()

# Retired option name. Double is now the default, so silently accepting the old
# flag would produce a build that looks configured and is not.
if(DEFINED CHROMO_STATE_FLOAT64)
    message(FATAL_ERROR
        "CHROMO_STATE_FLOAT64 was retired: double is the release default. "
        "Delete this build tree and reconfigure with no precision flag, or use "
        "-DCHROMO_STATE_FLOAT32=ON for the diagnostic float build.")
endif()

# Vendored Armadillo (header-only path in include/).
# Built without LAPACK/BLAS: the solver only uses element-wise ops, sqrt, abs,
# min/max, and sub2ind, which are all header-only.
add_compile_definitions(ARMA_DONT_USE_LAPACK ARMA_DONT_USE_BLAS)

# Two solvers, two source directories, one shared infrastructure layer.
#   src/*.cpp             shared: profiling, EOS, Grid geometry/broadcast
#   src/single_fluid/     RELEASE: U = (rho, rho u, E) equilibrium mixture
#   src/two_fluid/        HISTORICAL research: seven-row carrier state
# Neither solver directory depends on the other.
add_library(chromosphere STATIC
    src/profiling.cpp
    src/eos.cpp
    src/grid.cpp
    src/single_fluid/mixture.cpp
src/single_fluid/exact_rs.cpp
    src/single_fluid/integrator.cpp
    src/two_fluid/state.cpp
    src/two_fluid/flux.cpp
    src/two_fluid/rhs.cpp
    src/two_fluid/integrators.cpp
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
    ${CMAKE_CURRENT_SOURCE_DIR}/src
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

if(CHROMO_ENABLE_OPENMP)
    find_package(OpenMP QUIET)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(chromosphere PUBLIC OpenMP::OpenMP_CXX)
    elseif(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang")
        set(LIBOMP_ROOT "" CACHE PATH "Root of an LLVM/OpenMP installation")
        if(NOT LIBOMP_ROOT AND DEFINED ENV{LIBOMP_ROOT})
            set(LIBOMP_ROOT "$ENV{LIBOMP_ROOT}")
        endif()
        set(_libomp_roots ${LIBOMP_ROOT} /opt/homebrew/opt/libomp /usr/local/opt/libomp)
        find_path(LIBOMP_INCLUDE_DIR omp.h HINTS ${_libomp_roots} PATH_SUFFIXES include)
        find_library(LIBOMP_LIBRARY NAMES omp libomp HINTS ${_libomp_roots} PATH_SUFFIXES lib)
        if(NOT LIBOMP_INCLUDE_DIR OR NOT LIBOMP_LIBRARY)
            message(FATAL_ERROR
                "CHROMO_ENABLE_OPENMP=ON, but OpenMP was not found. Install libomp or set LIBOMP_ROOT.")
        endif()
        target_compile_options(chromosphere PUBLIC -Xpreprocessor -fopenmp)
        target_include_directories(chromosphere PUBLIC ${LIBOMP_INCLUDE_DIR})
        target_link_libraries(chromosphere PUBLIC ${LIBOMP_LIBRARY})
    else()
        message(FATAL_ERROR
            "CHROMO_ENABLE_OPENMP=ON, but a working OpenMP C++ runtime was not found")
    endif()
    target_compile_definitions(chromosphere PUBLIC CHROMO_USE_OPENMP=1)
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

foreach(RELEASE_GUARD IN ITEMS enable_te single_fluid two_fluid ionization)
    add_test(NAME release_reject_${RELEASE_GUARD}
             COMMAND ${CMAKE_COMMAND}
                     -DEXE=$<TARGET_FILE:chromo_main>
                     -DTABLE=${EOS_PRODUCTION_TABLE}
                     -DGUARD=${RELEASE_GUARD}
                     -DOUT=${CMAKE_CURRENT_BINARY_DIR}/release_guard_unused.dat
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_release_guard.cmake)
endforeach()

add_test(NAME release_output_sidecar
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
