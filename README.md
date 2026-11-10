# Ariadne Numeric

[![macOS Status](https://github.com/ariadne-cps/numeric/actions/workflows/macos.yml/badge.svg)](https://github.com/ariadne-cps/numeric/actions/workflows/macos.yml)
[![Ubuntu Status](https://github.com/ariadne-cps/numeric/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/ariadne-cps/numeric/actions/workflows/ubuntu.yml)
[![Windows Status](https://github.com/ariadne-cps/numeric/actions/workflows/win.yml/badge.svg)](https://github.com/ariadne-cps/numeric/actions/workflows/win.yml)
[![Sanitizers Status](https://github.com/ariadne-cps/numeric/actions/workflows/sanitizers.yml/badge.svg)](https://github.com/ariadne-cps/numeric/actions/workflows/sanitizers.yml)
[![Coverage Status](https://github.com/ariadne-cps/numeric/actions/workflows/coverage.yml/badge.svg)](https://github.com/ariadne-cps/numeric/actions/workflows/coverage.yml)
[![codecov](https://codecov.io/gh/ariadne-cps/numeric/branch/main/graph/badge.svg)](https://codecov.io/gh/ariadne-cps/numeric)

Ariadne Numeric is the standalone C++20 numeric layer used by Ariadne. It provides exact, floating-point, approximate and validated number types together with the arithmetic and logical operations used by the higher-level Ariadne libraries. Although independently maintained and usable, this repository is also part of the [Ariadne repository structure](https://github.com/ariadne-cps/ariadne/blob/main/doc/structure.md).

## Features

- Exact arithmetic with integers, naturals, dyadics, decimals and rationals.
- Double- and multiple-precision floating-point arithmetic.
- Explicit rounding modes and directed rounding.
- Validated lower and upper bounds, intervals, balls and error bounds.
- Approximate, exact, validated and effective generic number abstractions.
- Constructive real-number computations controlled by computational effort.
- Complex arithmetic over supported numeric types.
- Generic numeric traits, factories and arithmetic concepts used throughout Ariadne.

## Dependencies

Numeric depends directly on:

- [ariadne-cps/foundation](https://github.com/ariadne-cps/foundation), included as a Git submodule.

Foundation provides [ariadne-cps/utility](https://github.com/ariadne-cps/utility) transitively. Build configuration is shared through [ariadne-cps/build](https://github.com/ariadne-cps/build).

Standalone development additionally pins [Ariadne Testing](https://github.com/ariadne-cps/testing) and [Ariadne Benchmarking](https://github.com/ariadne-cps/benchmarking). Both are independent, opt-in development dependencies: their gitlinks are deliberately excluded from production dependency-commit consistency checks. Logging is supplied by optional Benchmarking, not by Numeric's production library. Neither Testing nor Benchmarking is linked to Numeric production targets or installed with Numeric.

## Build

Clone the repository together with its Git submodules:

```bash
git clone --recurse-submodules https://github.com/ariadne-cps/numeric.git
cd numeric
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --parallel --target ariadne-numeric bindings
cmake --install build-release --prefix install
```

A C++20 compiler and CMake are required. Python 3.11+ with development
headers is needed for the Python bindings; disable them with
`-DBUILD_PYTHON_BINDINGS=OFF` when not needed.

The CI builds the standalone C++ tutorial against the installed library,
runs both C++ and Python tutorials, and verifies uninstallation independently
of the test suite.

## Tests

Configure Debug separately and run both the C++ and Python tests:

```bash
python -m pip install pytest
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug -DPYTHON_ENFORCE_TESTS=ON
cmake --build build-debug --parallel --target tests
ctest --test-dir build-debug -LE python --output-on-failure
ctest --test-dir build-debug -L python --output-on-failure
```

## Benchmarks

The optional Benchmarking module registers the benchmark executables in the standalone `benchmarks` target. Build a Release configuration before drawing performance conclusions:

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_BINDINGS=OFF
cmake --build build-release --parallel --target benchmarks
./build-release/benchmarks/benchmark_big_unsigned --bits 53,128 --iterations 1000 --samples 7 --warmups 1 --summary min,max,median,mad
./build-release/benchmarks/benchmark_rational --bits 53,128 --iterations 1000 --samples 7 --warmups 1 --summary min,max,median,mad
./build-release/benchmarks/benchmark_dyadic --bits 53,128 --iterations 1000 --samples 7 --warmups 1 --summary min,max,median,mad
./build-release/benchmarks/benchmark_floatmp --bits 53,128 --iterations 1000 --samples 7 --warmups 1 --summary min,max,median,mad
./build-release/benchmarks/benchmark_compare_double_precision_backends --bits 53 --iterations 1000 --samples 7 --warmups 1 --details
```

BigUnsigned and Rational are **absolute, noncomparative** benchmarks. They
measure the native implementations without GMP, reference backends, speedup
ratios or slowdown thresholds. BigUnsigned includes the historical unsigned
arithmetic and square-root operations plus native division, GCD, multiplication,
power and decimal-conversion timings. Rational preserves all 15 original
operation names, including construction, arithmetic, extraction, rounding and
binary64 conversion. Immutable 32-element operand pools are constructed before
timing, and algebraic correctness invariants are checked before samples.
These checks do not assert performance equivalence with the former GMP runs;
their timings must not be compared across changed workloads.

Measurements use Benchmarking's `Runner` (preparation, timed operation, verification, explicit warm-ups and normalized nanosecond samples), `Statistics` and reporting facilities. FloatDP/FloatMP comparisons use a typed `ComparisonBattery` and paired/aggregate reports. Existing operand pools, operation names and sample-count validation remain in the owning Numeric repository. The Dyadic pool retains its exact-arithmetic verification before timing; the original FloatMP backend drivers did not have per-round numerical verification callbacks, so the migration does not misrepresent the empty verification callback as an additional correctness oracle. The --summary option now selects the displayed statistics from retained samples; A/B ratios still use median by default. The new shared CLI defaults to zero warm-ups; specify `--warmups 1` to reproduce the earlier benchmark drivers' implicit warm-up count. The comparison requires 53-bit precision because FloatDP has fixed double precision.

## Sanitizers

The Sanitizers workflow runs two independent Ubuntu/Clang Debug jobs:
ThreadSanitizer for data races, and AddressSanitizer with
UndefinedBehaviorSanitizer for memory errors and undefined behavior.
Both jobs disable Python bindings; they are separate from Coverage.

## Coverage

Configure a separate Debug build with coverage enabled:

```bash
mkdir build-coverage
cd build-coverage
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCOVERAGE=ON
cmake --build . --parallel --target coverage
```

On Ubuntu coverage is generated with GCC/lcov. On macOS it is generated with AppleClang/LLVM coverage tools. CI uploads the macOS coverage report to Codecov.

## Contribution guidelines

If you would like to contribute to Numeric, please contact the developer:

- Luca Geretti <luca.geretti@univr.it>

## Testing framework

The `main` branch uses [Ariadne Testing](https://github.com/ariadne-cps/testing)
**only in standalone tests**, pinned in `ss/tst`. All existing CTest
executables and the compile-only `check_numeric` target remain registered;
production targets, public headers, Python bindings and install rules do not
link to the new testing library.

All 20 CTest executables have been migrated to Ariadne Testing and link
`ariadne-testing` directly, with no dependency on the retired `macro`
submodule. The test suites preserve their numeric oracles, precision
variants and mathematical evidence checks; `--concise` selects compact
reporting. The compile-only `check_numeric` target remains unchanged.

```bash
git clone --recurse-submodules https://github.com/ariadne-cps/numeric.git
cd numeric
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug -DBUILD_PYTHON_BINDINGS=OFF
cmake --build build-debug --parallel --target tests
ctest --test-dir build-debug --output-on-failure
```

For subsequent changes, validate the original checks and outcomes,
including indeterminate evidence and exceptional paths, on the full
GCC/Clang/MSVC CI matrix. Frozen TSV-oracle replay,
precision-parametrized cases and state restoration remain covered;
testing migration does not regenerate reference datasets.

The design journals now reside directly in `doc/`: `big_unsigned_design.md`,
`integer_design.md`, `rational_design.md`, `dyadic_design.md`,
`floatdp_design.md` and `floatmp_design.md`.

### Float-rounding suite migration

The full `test_float_rounding` executable is migrated to the canonical
Ariadne Testing API. It retains the legacy three-precision order, all 16
phases per precision, and their original numerical loops/oracles. Each
phase is now a separate case, allowing later phases to execute after an
exception. The `floatmp_sin/cos/tan` oracle checks gain nested precision,
original-input and rounding-mode diagnostic scopes; the frozen TSV files
and independent generators remain untouched. Legacy print/notification
macros become non-asserting structured diagnostics. Compare all baseline
assertions and failures in GCC/Clang/MSVC CI before accepting the migration.

### Additional migrated suites

The `test_extended` phases become three independent cases. Decimal's static
concept assertion remains compile-time, its eight expected exceptions remain
checks, and every existing frozen-oracle row still performs three comparisons;
scopes identify the oracle and row ID in diagnostics. `test_float_error`
retains its original two-precision sequence with one case per backend.
Neither numerical expectation nor reference TSV is changed.

### Float-approximation migration

All eight phase/precision executions are registered as independent Testing
cases. Both compile-time concept assertions and all mathematical identity,
equality, and predicate checks are retained. Fixture-local checks use the
implicit-context shorthand macros; free-function tests retain `_CTX`.
The shared `test_floats.hpp` header contains no legacy testing
dependency. Test fixtures use the new framework consistently.
