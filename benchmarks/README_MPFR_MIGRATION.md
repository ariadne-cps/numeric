# FloatVP migration: static coexistence and MPFR performance reference

The reference remains `FloatMP = Float<MultiplePrecision>` (MPFR).
The new implementation is named `FloatVP = Float<VariablePrecision>`.
Both are concrete types, distinguished by their precision tags.
Neither introduces a virtual dispatch, a heap-allocated interface,
nor a function-pointer switch in scalar arithmetic.

`include/numeric/floatvp.hpp` reserves the new names and declares the
precision tag. **It does not yet implement Float<VariablePrecision>.**
The dual-backend build is intentionally rejected until the native
implementation source `src/floatvp.cpp` is added. The default tests
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

## Dual backend (once FloatVP exists)

Add a definition of `Float<VariablePrecision>`, with the same named
operations as the MPFR specialization, plus `src/floatvp.cpp`.
Enable dual compilation in an isolated Release build:

```sh
cmake -S . -B build-dual -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_PYTHON_BINDINGS=OFF -DBUILD_NUMERIC_MIGRATION_BENCHMARKS=ON \
  -DBUILD_FLOATVP_COMPARISON=ON
cmake --build build-dual --target tests benchmark_floatmp_backends
ctest --test-dir build-dual --output-on-failure
./build-dual/benchmarks/benchmark_floatmp_backends 300 7 > dual.csv
python3 benchmarks/compare_floatmp_backends.py dual.csv
```

The existing parameterized `test_float` and `test_float_rounding`
also instantiate the `VariablePrecision` type in a dual build, while
retaining their MPFR cases. For new tests, prefer one templated test
body instantiated for both `FloatMP` and `FloatVP`; use
explicitly native-only assertions only for implementation details.

The dual-backend executable prints all per-operation measurements as
a machine-readable CSV on **stdout**, and an immediate slowdown summary
on **stderr** (also shown when stdout is redirected to `dual.csv`).
A baseline-only executable explicitly says that no slowdown can be
computed without FloatVP, rather than presenting a misleading number.

The postprocessing script prints a detailed per-operation comparison
plus the same two summary metrics for every precision and all operations:

```sh
python3 benchmarks/compare_floatmp_backends.py dual.csv --summary-only
```

- `geomean_x`: geometric mean of the FloatVP/FloatMP time ratios across
  distinct (operation, precision) pairs, giving **equal weight** to each.
- `total_time_x`: sum of FloatVP's median ns/op divided by the sum of
  FloatMP's median ns/op; this assumes one execution of each benchmark
  case, and **weights expensive operations more heavily**.
- `worst_x`: largest observed slowdown, with its operation and precision.

Ratios above 1 indicate slowdown; ratios below 1 indicate a speedup.
The 'ALL' summary is a synthetic equal-frequency mix of the listed
operations **not** a prediction of application-wide performance: real
Ariadne workloads may exercise a different operation mix.

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
