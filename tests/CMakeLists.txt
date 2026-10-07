include(CTest)

set(UNIT_TESTS
    test_logical_sequence
    test_integer
    test_dyadic
    test_rational
    test_extended
    test_real
    test_number
    test_operators
    test_rounding_mode
    test_float
    test_float_rounding
    test_rounded_float
    test_float_approximation
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

add_library(check_numeric OBJECT check_numeric.cpp)
target_link_libraries(check_numeric PRIVATE ariadne-numeric numeric)

add_custom_target(tests)
add_dependencies(tests ${UNIT_TESTS} check_numeric)
