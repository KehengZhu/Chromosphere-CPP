# A build without -DCHROMO_ENABLE_LEGACY=ON must reject every legacy
# scenario immediately, with a message that names the option.
#   cmake -DEXE=<chromo_main> -DOUT=<unused output path> -P check_legacy_disabled.cmake
foreach(SCENARIO IN ITEMS model_c7 model_flare model_gentle analytic_canopy pfss_field_line model_column_2t)
    execute_process(
        COMMAND ${EXE} ${OUT} full ionization ${SCENARIO} - 0 cooling
        RESULT_VARIABLE STATUS
        OUTPUT_VARIABLE STDOUT
        ERROR_VARIABLE STDERR)
    set(COMBINED "${STDOUT}${STDERR}")
    if(STATUS EQUAL 0)
        message(FATAL_ERROR "legacy scenario ${SCENARIO} unexpectedly ran in a non-legacy build")
    endif()
    string(FIND "${COMBINED}" "-DCHROMO_ENABLE_LEGACY=ON" MATCH_POSITION)
    if(MATCH_POSITION EQUAL -1)
        message(FATAL_ERROR
            "legacy scenario ${SCENARIO} failed for the wrong reason; expected a "
            "message naming -DCHROMO_ENABLE_LEGACY=ON, got: ${COMBINED}")
    endif()
endforeach()
file(REMOVE "${OUT}")
message(STATUS "every legacy scenario is rejected with a message naming CHROMO_ENABLE_LEGACY")
