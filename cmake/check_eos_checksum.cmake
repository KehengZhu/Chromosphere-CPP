if(NOT DEFINED ROOT)
    message(FATAL_ERROR "ROOT must point to the Chromosphere2026 source tree")
endif()

set(TABLE "${ROOT}/data/eos/gamma1_hydrogen_v1.dat")
set(SIDECAR "${TABLE}.sha256")
if(NOT EXISTS "${TABLE}" OR NOT EXISTS "${SIDECAR}")
    message(FATAL_ERROR "EOS production table or SHA-256 sidecar is missing")
endif()

file(READ "${SIDECAR}" SIDECAR_TEXT)
string(REGEX MATCH "^[0-9A-Fa-f]+" EXPECTED_HASH "${SIDECAR_TEXT}")
string(LENGTH "${EXPECTED_HASH}" EXPECTED_LENGTH)
if(NOT EXPECTED_LENGTH EQUAL 64)
    message(FATAL_ERROR "Invalid SHA-256 sidecar: ${SIDECAR}")
endif()

file(SHA256 "${TABLE}" ACTUAL_HASH)
string(TOLOWER "${EXPECTED_HASH}" EXPECTED_HASH)
string(TOLOWER "${ACTUAL_HASH}" ACTUAL_HASH)
if(NOT ACTUAL_HASH STREQUAL EXPECTED_HASH)
    message(FATAL_ERROR
        "EOS table checksum mismatch: expected ${EXPECTED_HASH}, got ${ACTUAL_HASH}")
endif()

message(STATUS "EOS production table SHA-256 OK: ${ACTUAL_HASH}")
