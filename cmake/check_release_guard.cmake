# The release solver is a single-fluid common-temperature equilibrium mixture.
# Every legacy two-fluid / multi-temperature / finite-rate-ionization request
# must be rejected loudly rather than silently ignored.
if(GUARD STREQUAL "enable_te")
    set(ENV_ARGS ENABLE_TE=1)
    set(ION_ARG no-ionization)
    set(EXPECTED "ENABLE_TE is a legacy two-fluid setting")
elseif(GUARD STREQUAL "single_fluid")
    set(ENV_ARGS SINGLE_FLUID=0)
    set(ION_ARG no-ionization)
    set(EXPECTED "SINGLE_FLUID is a legacy two-fluid setting")
elseif(GUARD STREQUAL "two_fluid")
    set(ENV_ARGS ISO_TWO_FLUID=1)
    set(ION_ARG no-ionization)
    set(EXPECTED "ISO_TWO_FLUID is a legacy two-fluid setting")
elseif(GUARD STREQUAL "ionization")
    set(ENV_ARGS ISO_NS=4)
    set(ION_ARG ionization)
    set(EXPECTED "finite-rate ionization was requested")
else()
    message(FATAL_ERROR "unknown release guard: ${GUARD}")
endif()

execute_process(
    COMMAND ${CMAKE_COMMAND} -E env
            GAMMA_TABLE=${TABLE} ISO_NS=4 ${ENV_ARGS}
            ${EXE} ${OUT} full ${ION_ARG} model_column - 0 no-cooling
    RESULT_VARIABLE STATUS
    OUTPUT_VARIABLE STDOUT
    ERROR_VARIABLE STDERR)
set(COMBINED "${STDOUT}${STDERR}")
if(STATUS EQUAL 0)
    message(FATAL_ERROR "release guard ${GUARD} unexpectedly accepted the request")
endif()
string(FIND "${COMBINED}" "${EXPECTED}" MATCH_POSITION)
if(MATCH_POSITION EQUAL -1)
    message(FATAL_ERROR
        "release guard ${GUARD} failed for the wrong reason; expected '${EXPECTED}', got: ${COMBINED}")
endif()
message(STATUS "release guard ${GUARD} rejected with the expected diagnostic")
