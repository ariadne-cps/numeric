# FloatNP migration: static coexistence and MPFR performance reference

The reference remains `FloatMP = Float<MultiplePrecision>` (MPFR).
The new implementation is named `FloatNP = Float<NativePrecision>`.
Both are concrete types, distinguished by their precision tags.
Neither introduces a virtual dispatch, a heap-allocated interface,
nor a function-pointer switch in scalar arithmetic.

`include/numeric/floatnp.hpp` reserves the new names and declares the
precision tag. **It does not yet implement Float<NativePrecision>.**
The dual-backend build is intentionally rejected until the native
implementation source `src/floatnp.cpp` is added. The default tests
and benchmark executable exercise the MPFR reference baseline.

## Baseline

```sh
cmake -S . -B build-mpfr -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_PYTHON_BINDINGS=OFF -DBUILD_NUMERIC_MIGRATION_BENCHMARKS=ON
cmake --build build-mpfr --target test_floatmp_backends benchmark_floatmp_backends
ctest --test-dir build-mpfr -R '^test_floatmp_backends$' --output-on-failure
./build-mpfr/benchmarks/benchmark_floatmp_backends 300 7 > mpfr_baseline.csv
```

The benchmark runs 53, 128, 256 and 1024-bit cases. It reports
`backend,operation,bits,iterations,samples,median_ns`. Deterministic
preconstructed nontrivial operands, a warm-up, repeated measurements,
a median, and a compiler escape prevent the most common measurement
pitfalls. Cases include constructors, copy, comparison, conversions,
directed arithmetic, FMA, sqrt and transcendental functions.

The test verifies exact rational equality of rounded outputs against
MPFR under nearest, down, up and toward-zero rounding; it does not
compare via `double` (which would hide errors above 53 bits).

## Dual backend (once FloatNP exists)

Add a definition of `Float<NativePrecision>`, with the same named
operations as the MPFR specialization, plus `src/floatnp.cpp`.
Enable dual compilation in an isolated Release build:

```sh
cmake -S . -B build-dual -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_PYTHON_BINDINGS=OFF -DBUILD_NUMERIC_MIGRATION_BENCHMARKS=ON \
  -DBUILD_FLOATNP_COMPARISON=ON
cmake --build build-dual --target tests benchmark_floatmp_backends
ctest --test-dir build-dual --output-on-failure
./build-dual/benchmarks/benchmark_floatmp_backends 300 7 > dual.csv
python3 benchmarks/compare_floatmp_backends.py dual.csv
```

The existing parameterized `test_float` and `test_float_rounding`
also instantiate the `NativePrecision` type in a dual build, while
retaining their MPFR cases. For new tests, prefer one templated test
body instantiated for both `FloatMP` and `FloatNP`; use
explicitly native-only assertions only for implementation details.

The comparison script prints FloatNP/MPFR ratios; values above 1
indicate a slowdown. Benchmark on the same CPU, compiler, flags,
load and power/frequency conditions; repeat full runs and consider
allocation counts. No numerical performance claims are made before
both backends have been measured.

## Correctness requirements still outstanding

The deterministic differential cases are a starting point, not a
proof of correctly rounded transcendental operations. Before removing
MPFR, add randomized differential cases, exact-halfway ties,
exponent-range extremes, cancellation, signed zero, NaN, infinity,
and hard-to-round transcendentals at all target precisions.
Retain MPFR as a test oracle until correctness parity has been shown.

When parity and performance are acceptable, change the public
`FloatMP` mapping separately. This is an ABI-affecting decision,
not a rename to perform while testing both implementations.
