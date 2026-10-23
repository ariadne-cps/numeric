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

# Independent standalone C++20 TSV generators: compiled with tests, never auto-run.
set(ORACLE_GENERATORS
    generator_big_unsigned_gcd
    generator_big_unsigned_b76
    generator_exact_small
)
foreach(GENERATOR IN LISTS ORACLE_GENERATORS)
    add_executable(${GENERATOR} "data/${GENERATOR}.cpp")
    target_compile_features(${GENERATOR} PRIVATE cxx_std_20)
endforeach()
add_dependencies(tests ${ORACLE_GENERATORS})
