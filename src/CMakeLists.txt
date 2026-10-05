add_library(NUMERIC_SRC OBJECT
    operators.cpp
    int.cpp
    double.cpp
    integer.cpp
    dyadic.cpp
    decimal.cpp
    rational.cpp
    real.cpp
    complex.cpp
    number.cpp
    floatdp.cpp
    floatmp.cpp
    float_approximation.cpp
    float_lower_bound.cpp
    float_upper_bound.cpp
    float_bounds.cpp
    float_ball.cpp
    float_error.cpp
    float_literals.cpp
    concepts.cpp
    logical_sequence.cpp
)

target_include_directories(NUMERIC_SRC PRIVATE
    ${PROJECT_SOURCE_DIR}/include/numeric
)

target_compile_definitions(NUMERIC_SRC PRIVATE ARIADNE_NUMERIC_BUILD)
