set(BENCHMARKS
    benchmark_arithmetic
    benchmark_floatmp_backends
)

if(NOT TARGET benchmarks)
    add_custom_target(benchmarks)
endif()

foreach(BENCHMARK ${BENCHMARKS})
    add_executable(${BENCHMARK} ${BENCHMARK}.cpp)
    target_link_libraries(${BENCHMARK} PRIVATE ariadne-numeric numeric GMP::GMP)
endforeach()

if(BUILD_FLOATNP_COMPARISON)
    target_compile_definitions(benchmark_floatmp_backends PRIVATE ARIADNE_NUMERIC_ENABLE_FLOATNP)
endif()
add_dependencies(benchmarks ${BENCHMARKS})
