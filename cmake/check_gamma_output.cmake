file(REMOVE "${OUT}" "${OUT}.gamma_diag")
execute_process(
    COMMAND ${CMAKE_COMMAND} -E env
            GAMMA_TABLE=${TABLE} SINGLE_FLUID=1 ENABLE_TE=0 ISO_NS=48
            ${EXE} ${OUT} full no-ionization model_column - 0 no-cooling
    RESULT_VARIABLE STATUS
    OUTPUT_VARIABLE STDOUT
    ERROR_VARIABLE STDERR)
if(NOT STATUS EQUAL 0)
    message(FATAL_ERROR "gamma output smoke run failed: ${STDOUT}${STDERR}")
endif()
if(NOT EXISTS "${OUT}.gamma_diag")
    message(FATAL_ERROR "gamma diagnostic sidecar was not created")
endif()
file(SHA256 "${TABLE}" EXPECTED_SHA)
file(READ "${OUT}.gamma_diag" CONTENT)
foreach(EXPECTED IN ITEMS
        "EOS_MODE=gamma_table"
        "GAMMA_TABLE_SHA256=${EXPECTED_SHA}"
        "potential=cell_center_mean_of_phi_g_imh_phi_g_iph"
        "columns=rho_total v_cm T x_eq n_e n_HI p_total Gamma1 kappa_physical kappa_solver")
    string(FIND "${CONTENT}" "${EXPECTED}" POSITION)
    if(POSITION EQUAL -1)
        message(FATAL_ERROR "gamma sidecar missing '${EXPECTED}'")
    endif()
endforeach()
string(REGEX MATCH "^([0-9]+) 10" HEADER_MATCH "${CONTENT}")
if(NOT HEADER_MATCH)
    message(FATAL_ERROR "gamma sidecar data-width header is not the expected 10 columns")
endif()
file(REMOVE "${OUT}" "${OUT}.gamma_diag")
message(STATUS "gamma diagnostic sidecar contains physical fields and EOS provenance")
