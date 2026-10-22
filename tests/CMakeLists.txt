include(CTest)

set(UNIT_TESTS
    test_logical_sequence
    test_integer
    test_big_unsigned
    test_dyadic
    test_rational
    test_decimal
    test_extended
    test_real
    test_number
    test_operators
    test_rounding_mode
    test_float
    test_float_rounding
    test_rounded_float
    test_float_approximation
    test_float_error
    test_float_lower_upper_bound
    test_float_bounds
    test_float_ball
    test_complex
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne-numeric numeric)
    add_test(NAME ${TEST} COMMAND $<TARGET_FILE:${TEST}>)
    if(WIN32)
        set_tests_properties(${TEST} PROPERTIES
            ENVIRONMENT_MODIFICATION
                "PATH=path_list_prepend:$<TARGET_FILE_DIR:ariadne-numeric>"
        )
    endif()
endforeach()

# Absolute source path: out-of-source and Windows CTest need no copied fixtures.
target_compile_definitions(test_float_rounding PRIVATE
    "ARIADNE_FLOATAP_ROUNDING_ORACLE_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/data/floatap_rounding.tsv\""
    "ARIADNE_FLOATAP_ADD_SUB_ORACLE_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/data/floatap_add_sub.tsv\""
    "ARIADNE_FLOATAP_EXP_ORACLE_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/data/floatap_exp.tsv\""
    "ARIADNE_FLOATAP_POW_ORACLE_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/data/floatap_pow.tsv\""
    "ARIADNE_FLOATAP_SQRT_ORACLE_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/data/floatap_sqrt.tsv\""
    "ARIADNE_FLOATAP_LOG_ORACLE_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/data/floatap_log.tsv\""
    "ARIADNE_FLOATAP_CENTERED_EXP_ORACLE_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/data/floatap_centered_exp.tsv\""
    "ARIADNE_FLOATAP_ASIN_ORACLE_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/data/floatap_asin.tsv\""
    "ARIADNE_FLOATAP_ACOS_ORACLE_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/data/floatap_acos.tsv\""
    "ARIADNE_FLOATAP_ATAN_ORACLE_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/data/floatap_atan.tsv\""
)

add_library(check_numeric OBJECT check_numeric.cpp)
target_link_libraries(check_numeric PRIVATE ariadne-numeric numeric)

add_custom_target(tests)
add_dependencies(tests ${UNIT_TESTS} check_numeric)
