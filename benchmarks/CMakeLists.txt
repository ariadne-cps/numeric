set(BENCHMARKS
    benchmark_arithmetic
)

foreach(BENCHMARK ${BENCHMARKS})
    add_executable(${BENCHMARK} ${BENCHMARK}.cpp)
    target_link_libraries(${BENCHMARK} PRIVATE ariadne-numeric numeric GMP::GMP ${GMPXX_LIBRARY})
endforeach()

add_custom_target(benchmarks)
add_dependencies(benchmarks ${BENCHMARKS})
