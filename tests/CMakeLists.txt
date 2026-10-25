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

# Absolute fixture base directory for the oracle-using test targets.
# The TSV file names live beside their consumers in the C++ test files.
# Absolute path preserves out-of-source builds and Windows CTest.
foreach(TEST IN ITEMS test_big_unsigned test_integer test_rational test_float_rounding test_float test_rounded_float test_real)
    target_compile_definitions(${TEST} PRIVATE
        "ARIADNE_TEST_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/data\"")
endforeach()

add_library(check_numeric OBJECT check_numeric.cpp)
target_link_libraries(check_numeric PRIVATE ariadne-numeric numeric)

add_custom_target(tests)
add_dependencies(tests ${UNIT_TESTS} check_numeric)

# Standalone C++20 generators are independent of unit-test builds.
set(ORACLE_GENERATORS
    generator_big_unsigned_gcd
    generator_big_unsigned_b76
    generator_big_unsigned_decimal
    generator_big_unsigned_views
    generator_exact_small
    generator_big_unsigned_power_div
    generator_big_unsigned_division_remaining
    generator_rational_exact
    generator_exact_b72_r17
    generator_integer_structural
    generator_big_unsigned_mul
    generator_big_unsigned_dispatch
    generator_integer_random
    generator_rational_binary64
    generator_rational_shared_stream
    generator_rational_r31
    generator_rational_random
    generator_big_unsigned_random_mul
    generator_floatmp_directed_bounds
    generator_floatmp_pi
    generator_floatmp_get_d
    generator_floatmp_power_exact
    generator_floatmp_sqrt
    generator_floatmp_pow_short
    generator_rounded_float_directed
    generator_floatmp_exp
    generator_floatmp_log
    generator_floatmp_sin_cos
    generator_floatmp_centered_exp
    generator_floatmp_atan
    generator_floatmp_tan
    generator_floatmp_asin_acos
    generator_floatmp_tanh
    generator_floatmp_rounding
    generator_big_unsigned_sqrt
)
set(ORACLE_GENERATOR_COMMANDS)
set(ORACLE_VERIFY_COMMANDS)
foreach(GENERATOR IN LISTS ORACLE_GENERATORS)
    add_executable(${GENERATOR} EXCLUDE_FROM_ALL "data/${GENERATOR}.cpp")
    target_compile_features(${GENERATOR} PRIVATE cxx_std_20)
    list(APPEND ORACLE_GENERATOR_COMMANDS
        COMMAND "$<TARGET_FILE:${GENERATOR}>" --write --data-dir "${CMAKE_CURRENT_SOURCE_DIR}/data")
    list(APPEND ORACLE_VERIFY_COMMANDS
        COMMAND "$<TARGET_FILE:${GENERATOR}>" --verify --data-dir "${CMAKE_CURRENT_SOURCE_DIR}/data")
endforeach()

# Build only the independent generator executables.
add_custom_target(generators)
add_dependencies(generators ${ORACLE_GENERATORS})

# Explicitly overwrite only the TSVs covered by versioned generators.
# Commands run sequentially; the first failure stops further generation.
add_custom_target(generate-data
    ${ORACLE_GENERATOR_COMMANDS}
    COMMENT "Regenerating source-tree test data using standalone generators"
    VERBATIM
)
add_dependencies(generate-data generators)

# Non-destructive byte-for-byte checks of every versioned generator fixture.
add_custom_target(verify-data
    ${ORACLE_VERIFY_COMMANDS}
    COMMENT "Verifying independent source-tree test data"
    VERBATIM
)
add_dependencies(verify-data generators)
