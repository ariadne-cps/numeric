# Ariadne Numeric

[![macOS Status](https://github.com/ariadne-cps/numeric/actions/workflows/macos.yml/badge.svg)](https://github.com/ariadne-cps/numeric/actions/workflows/macos.yml)
[![Ubuntu Status](https://github.com/ariadne-cps/numeric/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/ariadne-cps/numeric/actions/workflows/ubuntu.yml)
[![Windows Status](https://github.com/ariadne-cps/numeric/actions/workflows/win.yml/badge.svg)](https://github.com/ariadne-cps/numeric/actions/workflows/win.yml)
[![Sanitizers Status](https://github.com/ariadne-cps/numeric/actions/workflows/sanitizers.yml/badge.svg)](https://github.com/ariadne-cps/numeric/actions/workflows/sanitizers.yml)
[![Coverage Status](https://github.com/ariadne-cps/numeric/actions/workflows/coverage.yml/badge.svg)](https://github.com/ariadne-cps/numeric/actions/workflows/coverage.yml)
[![codecov](https://codecov.io/gh/ariadne-cps/numeric/branch/main/graph/badge.svg)](https://codecov.io/gh/ariadne-cps/numeric)

Ariadne Numeric is the standalone C++20 numeric layer used by Ariadne. It provides exact, floating-point, approximate and validated number types together with the arithmetic and logical operations used by the higher-level Ariadne libraries.

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

Foundation provides [ariadne-cps/utility](https://github.com/ariadne-cps/utility) transitively. Build configuration is shared through [ariadne-cps/configuration](https://github.com/ariadne-cps/configuration).

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

## License

Numeric is released under the GNU General Public License v3.0.

## Testing framework experiment (`testing` branch)

The `testing` branch introduces [Ariadne Testing](https://github.com/ariadne-cps/testing)
**only in standalone tests**, pinned in `ss/tst`. All existing CTest
executables and the compile-only `check_numeric` target remain registered;
production targets, public headers, Python bindings and install rules do not
link to the new testing library.

The first two migrated consumers are `test_logical_sequence` (explicit
mathematical evidence expectations, including possibly/not-possibly) and
`test_rounding_mode` (real backend rounding modes and equality checks).
All other executables still use the legacy testing macros. This is a
controlled migration experiment, **not** evidence that the entire Numeric
suite has been migrated.

```bash
git clone --branch testing --recurse-submodules https://github.com/ariadne-cps/numeric.git
cd numeric
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug -DBUILD_PYTHON_BINDINGS=OFF
cmake --build build-debug --parallel --target tests
ctest --test-dir build-debug --output-on-failure
```

Before merging any migration, compare the old and new checks and outcomes,
including indeterminate evidence and exceptional paths, with a green
GCC/Clang/MSVC CI matrix. Subsequent Numeric experiments should cover
frozen TSV-oracle replay, precision-parametrized cases and state restoration;
no reference datasets are regenerated as part of testing migration.
