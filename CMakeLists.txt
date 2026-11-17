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
# The diagnostic float32 build option was removed on 2026-11-16, once the
# precision question was closed.

# Vendored Armadillo (header-only path in include/).
# Built without LAPACK/BLAS: the solver only uses element-wise ops, sqrt, abs,
# min/max, and sub2ind, which are all header-only.
add_compile_definitions(ARMA_DONT_USE_LAPACK ARMA_DONT_USE_BLAS)

# Solvers by directory, over one shared infrastructure layer.
#   src/core/             shared: Grid, EOS, profiling, OpenMP helpers, physics formulas
#   src/single_fluid/     RELEASE: U = (rho, rho u, E) equilibrium mixture
#   scenarios/            model_column IC+BC, mesh, C7 profile
#   legacy/two_fluid/     HISTORICAL research: seven-row carrier state  (legacy only)
#   legacy/two_temp/      closed test study: U = (rho, rho u, E, E_e), separate
#                         electron and heavy-particle temperatures. Reuses only
#                         single_fluid/exact_rs.hpp (the standalone SWMF
#                         ideal-gas Riemann solver port)              (legacy only)
#   legacy/scenarios/     model_c7, model_flare, model_gentle, analytic_canopy,
#                         pfss_field_line (two-fluid), model_column_2t (legacy only)
# No solver directory depends on another's solver code.
add_library(chromosphere STATIC
    src/core/profiling.cpp
    src/core/eos.cpp
    src/core/grid.cpp
    src/single_fluid/mixture.cpp
    src/single_fluid/exact_rs.cpp
    src/single_fluid/integrator.cpp
    scenarios/scenario.cpp
    scenarios/mesh.cpp
    scenarios/c7_profile.cpp
    scenarios/model_column.cpp
    scenarios/data_file_parser.cpp
)

# The frozen two-fluid solver and its scenarios (legacy/README.md). Off by
# default, so a default build neither compiles nor links any of it; asking a
# default build for a legacy scenario fails immediately with a message naming
# this option.
option(CHROMO_ENABLE_LEGACY "Build the frozen two-fluid and two-temperature solvers and their scenarios" OFF)
if(CHROMO_ENABLE_LEGACY)
    target_sources(chromosphere PRIVATE
        legacy/two_fluid/state.cpp
        legacy/two_fluid/flux.cpp
        legacy/two_fluid/rhs.cpp
        legacy/two_fluid/integrators.cpp
        legacy/scenarios/legacy_scenarios.cpp
        legacy/scenarios/model_c7.cpp
        legacy/scenarios/model_flare.cpp
        legacy/scenarios/model_gentle.cpp
        legacy/scenarios/analytic_canopy.cpp
        legacy/scenarios/pfss_field_line.cpp
        legacy/two_temp/closure.cpp
        legacy/two_temp/hydro.cpp
        legacy/two_temp/conduction.cpp
        legacy/scenarios/model_column_2t.cpp
    )
    target_compile_definitions(chromosphere PUBLIC CHROMO_ENABLE_LEGACY=1)
endif()
target_include_directories(chromosphere PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_executable(chromo_main src/main.cpp)
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

if(NOT CHROMO_ENABLE_LEGACY)
    add_test(NAME legacy_scenario_disabled
             COMMAND ${CMAKE_COMMAND}
                     -DEXE=$<TARGET_FILE:chromo_main>
                     -DOUT=${CMAKE_CURRENT_BINARY_DIR}/legacy_disabled_unused.dat
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_legacy_disabled.cmake)
endif()

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
