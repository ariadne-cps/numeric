add_library(NUMERIC_SRC OBJECT
    operators.cpp
    int.cpp
    double.cpp
    integer.cpp
    big_unsigned.cpp
    big_unsigned_karatsuba.cpp
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

# Architecture kernels use baseline 64x64 multiplication only.
# Unsupported architecture/compiler combinations keep the C++20 reference.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$"
        AND CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang|MSVC)$")
    target_sources(NUMERIC_SRC PRIVATE big_unsigned_x86-64.cpp)
    target_compile_definitions(NUMERIC_SRC PRIVATE ARIADNE_BIG_UNSIGNED_X86_64=1)
    # The fused 128-bit accumulator is compiler-sensitive: the same
    # algorithm was faster with GCC 14 and not consistently with Clang.
    # Keep the B3 schoolbook selector on Clang/AppleClang/MSVC and allow
    # a same-host GCC control build by setting this option to OFF.
    option(NUMERIC_BIG_UNSIGNED_GCC_FUSED_MULTIPLY
           "Select fused 128-bit BigUnsigned multiplication on x86-64 GCC" ON)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
            AND NUMERIC_BIG_UNSIGNED_GCC_FUSED_MULTIPLY)
        target_compile_definitions(NUMERIC_SRC PRIVATE
            ARIADNE_BIG_UNSIGNED_GCC_FUSED=1)
    endif()
endif()

# Experimental Apple Silicon/Linux AArch64 fused multiply. AppleClang/Clang
# and GNU support unsigned __int128, lowered to baseline AArch64 MUL/UMULH.
# Keep the portable production selector until same-host benchmarks establish
# a speedup; the forced mul_into_fused entry point is always available here.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64|ARM64)$"
        AND CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_sources(NUMERIC_SRC PRIVATE big_unsigned_arm64.cpp)
    target_compile_definitions(NUMERIC_SRC PRIVATE ARIADNE_BIG_UNSIGNED_ARM64=1)
    option(NUMERIC_BIG_UNSIGNED_ARM64_FUSED_MULTIPLY
           "Select fused BigUnsigned multiplication on ARM64 (experimental)" OFF)
    if(NUMERIC_BIG_UNSIGNED_ARM64_FUSED_MULTIPLY)
        target_compile_definitions(NUMERIC_SRC PRIVATE
            ARIADNE_BIG_UNSIGNED_ARM64_FUSED=1)
    endif()
endif()
