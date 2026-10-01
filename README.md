# Ariadne Numeric

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Unix Status](https://github.com/ariadne-cps/numeric/actions/workflows/unix.yml/badge.svg)](https://github.com/ariadne-cps/numeric/actions/workflows/unix.yml)
[![Windows Status](https://github.com/ariadne-cps/numeric/actions/workflows/win.yml/badge.svg)](https://github.com/ariadne-cps/numeric/actions/workflows/win.yml)
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

- [ariadne-cps/paradigm](https://github.com/ariadne-cps/paradigm), included as a Git submodule.
- [GMP](https://gmplib.org/).
- [MPFR](https://www.mpfr.org/).

Paradigm provides [ariadne-cps/utility](https://github.com/ariadne-cps/utility) transitively. Build configuration is shared through [ariadne-cps/configuration](https://github.com/ariadne-cps/configuration).

## Build

Clone the repository together with its Git submodules:

```bash
git clone --recurse-submodules https://github.com/ariadne-cps/numeric.git
cd numeric
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --parallel
ctest --output-on-failure
```

A C++20 compiler, CMake, GMP and MPFR are required.

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
