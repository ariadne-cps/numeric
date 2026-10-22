# FloatAP migration: static coexistence and MPFR performance reference

The reference remains `FloatMP = Float<MultiplePrecision>` (MPFR).
The new implementation is named `FloatAP = Float<ArbitraryPrecision>`.
Both are concrete types, distinguished by their precision tags.
Neither introduces a virtual dispatch, a heap-allocated interface,
nor a function-pointer switch in scalar arithmetic.

`include/numeric/floatap.hpp` reserves the new names and declares the
precision tag. **It does not yet implement Float<ArbitraryPrecision>.**
The dual-backend build is intentionally rejected until the native
implementation source `src/floatap.cpp` is added. The default tests
and benchmark executable exercise the MPFR reference baseline.

## Baseline

```sh
cmake -S . -B build-mpfr -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_PYTHON_BINDINGS=OFF -DBUILD_NUMERIC_MIGRATION_BENCHMARKS=ON
cmake --build build-mpfr --target test_float_multipleprecision_backends benchmark_compare_multiple_precision_backends
ctest --test-dir build-mpfr -R '^test_float_multipleprecision_backends$' --output-on-failure
./build-mpfr/benchmarks/benchmark_compare_multiple_precision_backends 300 7 > mpfr_baseline.csv
```

The benchmark defaults to 53, 128, 256 and 1024-bit cases; explicitly
pass `--bits 53,128,256,1024,2048,4096` to include 2048 and 4096 bits. It reports
`backend,operation,bits,iterations,samples,median_ns`. Deterministic
preconstructed nontrivial operands, a warm-up, repeated measurements,
a median, and a compiler escape prevent the most common measurement
pitfalls. Cases include constructors, copy, comparison, conversions,
directed arithmetic, FMA, sqrt and transcendental functions.

The test verifies exact rational equality of rounded outputs against
MPFR under nearest, down, up and toward-zero rounding; it does not
compare via `double` (which would hide errors above 53 bits).

## Dual backend (once FloatAP exists)

Add a definition of `Float<ArbitraryPrecision>`, with the same named
operations as the MPFR specialization, plus `src/floatap.cpp`.
Enable dual compilation in an isolated Release build:

```sh
cmake -S . -B build-dual -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_PYTHON_BINDINGS=OFF -DBUILD_NUMERIC_MIGRATION_BENCHMARKS=ON \
  -DBUILD_FLOATAP_COMPARISON=ON
cmake --build build-dual --target tests benchmark_compare_multiple_precision_backends
ctest --test-dir build-dual --output-on-failure
./build-dual/benchmarks/benchmark_compare_multiple_precision_backends 300 7 > dual.csv
```

The existing parameterized `test_float` and `test_float_rounding`
also instantiate the `ArbitraryPrecision` type in a dual build, while
retaining their MPFR cases. For new tests, prefer one templated test
body instantiated for both `FloatMP` and `FloatAP`; use
explicitly native-only assertions only for implementation details.

## Slowdown summary entirely in C++

The benchmark calculates its own FloatAP/FloatMP ratios. No Python
postprocessing step or additional runtime dependency is necessary.
It writes raw nanoseconds per operation to CSV on **stdout** and
the per-precision and overall slowdown summary to **stderr**.

To see **just the summary**, without creating a CSV file:

```sh
./build-dual/benchmarks/benchmark_compare_multiple_precision_backends 300 7 --summary-only
```

To print the individual FloatAP/FloatMP ratios as well:

```sh
./build-dual/benchmarks/benchmark_compare_multiple_precision_backends 300 7 --bits 53,128,256,1024,2048,4096 --summary-only --details
```

To fail the process if **any** (operation, precision) runs slower
than, say, 1.50 times the MPFR reference:

```sh
./build-dual/benchmarks/benchmark_compare_multiple_precision_backends 300 7 --summary-only --max-slowdown 1.50
```

The comparison uses two distinct aggregate metrics:
- `geomean`: geometric mean of FloatAP/FloatMP time ratios across
  distinct (operation, precision) pairs, with **equal weight** each.
- `total_time`: sum of FloatAP median ns/op divided by the sum of
  FloatMP median ns/op. This represents **one execution of each
  benchmark case** and weights expensive operations more heavily.
- `worst`: largest observed slowdown, with its operation and precision.

Ratios above 1 indicate slowdown; ratios below 1 indicate a speedup.
The 'ALL' summary is a synthetic equal-frequency workload, **not** a
prediction of application-wide performance. Real Ariadne workloads
may execute different operations at different frequencies.

A reference-only build still produces MPFR timing CSV, but explicitly
reports that no FloatAP slowdown can be computed yet; supplying
`--max-slowdown` in such a build is an error.

Benchmark on the same CPU, compiler, flags, load and power/frequency
conditions; repeat full runs and consider allocation counts. No
numerical performance claims are made before both backends have been
measured.

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
