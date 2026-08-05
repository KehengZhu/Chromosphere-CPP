if(GUARD STREQUAL "enable_te")
    set(ENV_ARGS ENABLE_TE=1 SINGLE_FLUID=1)
    set(ION_ARG no-ionization)
    set(EXPECTED "ENABLE_TE=1 is incompatible")
elseif(GUARD STREQUAL "single_fluid")
    set(ENV_ARGS ENABLE_TE=0 SINGLE_FLUID=0)
    set(ION_ARG no-ionization)
    set(EXPECTED "explicit SINGLE_FLUID=0 is incompatible")
elseif(GUARD STREQUAL "ionization")
    set(ENV_ARGS ENABLE_TE=0 SINGLE_FLUID=1)
    set(ION_ARG ionization)
    set(EXPECTED "finite-rate ionization was requested")
else()
    message(FATAL_ERROR "unknown gamma guard: ${GUARD}")
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
    message(FATAL_ERROR "gamma guard ${GUARD} unexpectedly accepted the request")
endif()
string(FIND "${COMBINED}" "${EXPECTED}" MATCH_POSITION)
if(MATCH_POSITION EQUAL -1)
    message(FATAL_ERROR
        "gamma guard ${GUARD} failed for the wrong reason; expected '${EXPECTED}', got: ${COMBINED}")
endif()
message(STATUS "gamma guard ${GUARD} rejected with expected diagnostic")
