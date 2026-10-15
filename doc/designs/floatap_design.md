# FloatAP — arbitrary-precision floating point

This is the living design and user manual for the **experimental** `FloatAP`
backend. Keep this document synchronized with every implementation change.
The migration is tracked on branch `remove-mpfr`.

## Names and coexistence

- **`FloatMP = Float<MultiplePrecision>`**: existing MPFR reference.
- **`FloatAP = Float<ArbitraryPrecision>`**: new native arithmetic, independent
  of MPFR once complete.
- `ArbitraryPrecision` is a bit count, **not** a multiple of IEEE double
  precision. Examples: 7, 53, 113, 127, 256, 1000 bits.
- Both are distinct concrete types, selected at compile time. No virtual
  scalar interface, backend pointer, or per-operation dynamic dispatch.

The name 'arbitrary precision' describes user-selectable significance,
not an infinite exponent range or necessarily unbounded resource usage.

## Planned representation

A finite nonzero number is represented as

```
value = (-1)^negative * significand * 2^exponent
```

with an integer significand and explicit precision in bits. A normalized
nonzero value has exactly the requested number of significant bits.
The exponent is a signed 64-bit integer in the initial kernel; checked
overflow raises an exception. Separate flags represent positive/negative
zero, positive/negative infinity, and NaN.

For the first implementation, GMP `mpz_t` supplies integer storage.
**MPFR is not used by the native kernel.** This is an intentionally
temporary dependency: the long-term objective remains replacing GMP
and MPFR. It also enables direct, exact rational comparison with MPFR
in the tests.

## Rounding contract

Supported initial modes: nearest-even, toward negative infinity,
toward positive infinity, and toward zero. A rational numerator /
denominator is converted by:

1. determining the binary magnitude/exponent;
2. scaling to a `p`-bit significand;
3. integer quotient and remainder;
4. comparing twice the remainder with the divisor for nearest-even;
5. applying directed rounding according to sign;
6. normalizing on significand carry and adjusting the exponent.

Exact values should remain exact. Special-value rules must follow
the floating-point contract (including signed zeros, NaN, infinities,
and cancellation). The implementation must not round twice; `fma`
requires one rounding after exact multiply-and-add.

## Performance roadmap

The exact-integer baseline is deliberately simple and may allocate
large temporaries. Optimization must preserve bit-identical directed
results:

- Fast path for single- and double-limb significands;
- preallocated/scratch temporaries to reduce allocations;
- guard/round/sticky bits for bounded-work alignment when exponents
  differ greatly (rather than shifting by the entire exponent gap);
- multiplication/division specialized for fixed requested precision;
- avoid rational canonicalization in ordinary arithmetic;
- profile separately at 53, 128, 256 and 1024 bits.

The first implementation should prefer correctness and transparent
measurements over speculative optimizations.

## Verification and benchmarks

`benchmark_compare_multiple_precision_backends` is the sole C++
benchmark executable. It measures both concrete backends and produces
`FloatAP / FloatMP` ratios when `BUILD_FLOATAP_COMPARISON=ON` and
the full `FloatAP` class is implemented:

```sh
./build-dual/benchmarks/benchmark_compare_multiple_precision_backends \
  300 7 --summary-only
```

`--details` reports operation-level ratios; `--max-slowdown 1.5`
fails if **any** operation exceeds that ratio. Geometric-mean slowdown
weights each (operation, precision) pair equally; time-sum slowdown
assumes one call to each operation; neither is an application workload
profile. All reports and comparisons are in **C++**, not Python.

The shared `test_float`, `test_float_rounding`, and differential
tests must instantiate for both backends once `FloatAP` meets their
API requirements. Correct rounding of transcendental operations
must be validated against MPFR, including hard-to-round cases.
The MPFR library remains a *test oracle*, not an implementation
dependency of native arithmetic.

## Status and limitations

- The `FloatAP` name/precision tag and dual-backend testing hooks
  are scaffolding, not a production floating-point implementation.
- The GMP-only kernel lives in `src/floatap.cpp`, with the internal
  interface in `include/numeric/floatap.hpp` and its C++
  test in `tests/test_floatap_kernel.cpp`. The test covers randomized
  rational-reference add/sub/mul/div/fma under four rounding modes,
  plus special values, with 7/53/113/128/256/1024-bit precision.
  **The public `Float<ArbitraryPrecision>` specialization is NOT complete
  or benchmarkable.**
- The kernel is initially exact-integer based; exponent-difference shifts can
  allocate excessively. This requires optimized guard/sticky-bit arithmetic.
- The complete public `FloatAP` API, range/error semantics,
  transcendental functions and performance optimization remain work.
- Until both versions build and run on the same machine, **no actual
  slowdown number exists**.

## Changing the architecture

Update this manual whenever representation, rounding, exponent limits,
allocation strategy, tests, benchmark methodology, or public names
change. Do not silently replace `FloatMP`; migration of the public
alias is a separate, ABI-relevant change after correctness and
performance review.

## Coverage-driven regression tests

The macOS LLVM coverage report dated 2026-10-15 identified untested
integer conversion, exact halfway rounding, negative denominators,
invalid rounding modes, exceptional rational conversion, assignment,
swap, signed zeros, and NaN/infinity arithmetic. These are now exercised
by the dedicated `test_floatap_edge_cases` C++ executable, alongside
`test_floatap_kernel` (randomized exact arithmetic).

The implementation is declared in the sole `numeric/floatap.hpp` header;
there is no `detail` directory and no separate kernel header. The
prototype `FloatAPKernel` class remains distinguishable from the
not-yet-implemented public specialization `Float<ArbitraryPrecision>`.

In an isolated Clang/LLVM instrumented build with both C++ test
executables, `src/floatap.cpp` achieved 100% line coverage and about
98.7% branch coverage. These figures do not describe the entire Ariadne
project, and some infeasible platform-dependent overflow branches remain.

## Coverage integration

The internal implementation is compiled into `NUMERIC_SRC` unconditionally;
`test_floatap_kernel` and `test_floatap_edge_cases` are in the normal CTest test list, so the CI `coverage`
target exercises and instruments `src/floatap.cpp`. This does **not**
mean the public `FloatAP` API is complete: CI still uses MPFR for `FloatMP`.

Build and run the kernel test:

```sh
cmake --build build --target test_floatap_kernel
ctest --test-dir build -R '^test_floatap_kernel$' --output-on-failure
```

The kernel is a prototype using GMP limbs. Comparison with MPFR and
hard-to-round transcendental validation are still pending.

## Coverage report of 15 October 2026, 21:32

The attached LLVM report `Archivio 74(1).zip` confirms:

| Source | Functions | Lines | Branches |
| --- | --- | --- | --- |
| `src/floatap.cpp` | 16/16 (100%) | 118/118 (100%) | 152/154 (98.70%) |
| `include/numeric/floatap.hpp` | 23/27 (85.19%) | 63/75 (84%) | 19/26 (73.08%) |

The four uncovered header functions are the `ArbitraryPrecision`
constructors, `bits()` accessor and equality operator. The
`precision_tag()` checks in `test_floatap_edge_cases` now exercise both
constructors, equal and unequal instances, and non-power-of-two bit counts.
The new tests require a fresh coverage run to confirm their effect.

Remaining uncovered kernel paths concern checked exponent arithmetic,
oversized shifts and integer size limits. Forcing extraordinary memory
allocation just to exercise platform-specific overflow is not a
correctness test. Coverage is not evidence of fully correct arithmetic:
additional independent boundary and MPFR differential checks remain
necessary.
