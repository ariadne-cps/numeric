set(BENCHMARKS
    benchmark_arithmetic
    benchmark_big_unsigned_multiply
    benchmark_big_unsigned_arithmetic
    benchmark_big_unsigned_division_profile
    benchmark_compare_multiple_precision_backends
    benchmark_compare_double_precision_backends
    benchmark_compare_integer_backends
    benchmark_integer_alias_dispatch
    benchmark_integer_decimal_profile
    benchmark_integer_small_power_profile
    benchmark_integer_small_operations_profile
    benchmark_compare_rational_backends
    benchmark_integer
    benchmark_rational
)

if(NOT TARGET benchmarks)
    add_custom_target(benchmarks)
endif()

foreach(BENCHMARK ${BENCHMARKS})
    add_executable(${BENCHMARK} ${BENCHMARK}.cpp)
    target_link_libraries(${BENCHMARK} PRIVATE ariadne-numeric numeric GMP::GMP logging LOGGING_SRC)
endforeach()

add_dependencies(benchmarks ${BENCHMARKS})
