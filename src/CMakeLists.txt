add_library(NUMERIC_SRC OBJECT
    operators.cpp
    int.cpp
    double.cpp
    integer.cpp
    big_unsigned.cpp
    dyadic.cpp
    decimal.cpp
    rational.cpp
    real.cpp
    complex.cpp
    number.cpp
    floatdp.cpp
    floatmp.cpp
    floatap.cpp
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

# The optional x86-64 kernel uses baseline 64x64 multiplication only.
# Other architectures and unsupported toolchains use the C++20 reference.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$"
        AND CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang|MSVC)$")
    target_sources(NUMERIC_SRC PRIVATE big_unsigned_x86-64.cpp)
    target_compile_definitions(NUMERIC_SRC PRIVATE ARIADNE_BIG_UNSIGNED_X86_64=1)
endif()
