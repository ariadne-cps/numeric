set(BENCHMARKS
    benchmark_floatmp
    benchmark_compare_double_precision_backends
    benchmark_dyadic
    benchmark_big_unsigned
    benchmark_rational
)

if(NOT TARGET benchmarks)
    add_custom_target(benchmarks)
endif()

# Numeric's standalone library bundles the Numeric, Foundation and Utility
# object targets. Benchmarking is an INTERFACE target in consumer builds;
# its Threading/Logging implementations are not otherwise linked into the
# executable. Include only their object targets not already bundled by
# ariadne-numeric to avoid duplicate definitions, especially on Windows.
get_target_property(_benchmark_thread_objects
    threading ARIADNE_PROJECT_OBJECT_TARGETS)
get_target_property(_numeric_objects
    numeric ARIADNE_PROJECT_OBJECT_TARGETS)
if(NOT _benchmark_thread_objects OR NOT _numeric_objects)
    message(FATAL_ERROR "Benchmarking requires registered Numeric/Threading objects")
endif()
foreach(_object IN LISTS _numeric_objects)
    list(REMOVE_ITEM _benchmark_thread_objects "${_object}")
endforeach()

foreach(BENCHMARK IN LISTS BENCHMARKS)
    add_executable(${BENCHMARK} ${BENCHMARK}.cpp)
    target_link_libraries(${BENCHMARK} PRIVATE
        ariadne-numeric numeric ariadne-benchmarking)
    foreach(_object IN LISTS _benchmark_thread_objects)
        target_sources(${BENCHMARK} PRIVATE
            $<TARGET_OBJECTS:${_object}>)
    endforeach()
endforeach()

add_dependencies(benchmarks ${BENCHMARKS})
