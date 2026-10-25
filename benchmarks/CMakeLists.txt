set(BENCHMARKS
    benchmark_floatmp
    benchmark_compare_double_precision_backends
    benchmark_dyadic
)

if(NOT TARGET benchmarks)
    add_custom_target(benchmarks)
endif()

foreach(BENCHMARK ${BENCHMARKS})
    add_executable(${BENCHMARK} ${BENCHMARK}.cpp)
    target_link_libraries(${BENCHMARK} PRIVATE ariadne-numeric numeric logging LOGGING_SRC)
endforeach()

add_dependencies(benchmarks ${BENCHMARKS})
