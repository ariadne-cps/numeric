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
  test in `tests/test_floatap.cpp`. The test covers randomized
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
`test_floatap` (randomized exact arithmetic).

The implementation is declared in the sole `numeric/floatap.hpp` header;
there is no `detail` directory and no separate kernel header. The
public `FloatAP` class remains distinguishable from the
not-yet-implemented public specialization `Float<ArbitraryPrecision>`.

In an isolated Clang/LLVM instrumented build with both C++ test
executables, `src/floatap.cpp` achieved 100% line coverage and about
98.7% branch coverage. These figures do not describe the entire Ariadne
project, and some infeasible platform-dependent overflow branches remain.

## Coverage integration

The internal implementation is compiled into `NUMERIC_SRC` unconditionally;
`test_floatap` and `test_floatap_edge_cases` are in the normal CTest test list, so the CI `coverage`
target exercises and instruments `src/floatap.cpp`. This does **not**
mean the public `FloatAP` API is complete: CI still uses MPFR for `FloatMP`.

Build and run the kernel test:

```sh
cmake --build build --target test_floatap
ctest --test-dir build -R '^test_floatap$' --output-on-failure
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

## Coverage report of 15 October 2026, 21:37

The next uploaded coverage report, `Archivio 75(1).zip`, confirms
that the new precision-tag tests cover all 27/27 functions in
`include/numeric/floatap.hpp`. Header lines improved from 63/75
(84%) to 69/75 (92%), while branch coverage stayed at 19/26
(73.08%). The remaining six uncovered header lines are checks for
exceptional 64-bit exponent/shift and huge-operand bounds that are
unreachable or resource-intensive via the ordinary public kernel API.

`src/floatap.cpp` remains 118/118 covered lines, 16/16 functions,
and 152/154 branches. The uncovered branch in `fma` involving two
infinities with the same sign can be tested without artificial inputs;
one other reported branch is related to the zero-remainder rounding
path. Coverage alone does not establish correct rounding or stability.

## Follow-up coverage tests (15 October 2026)

The updated `test_floatap_edge_cases` exercises same-sign infinite FMA,
classification of finite nonzero/zero/NaN, precision validation and
positive/negative exponent overflow. These represent genuine arithmetic
contracts, not fabricated coverage hooks.

An isolated Clang coverage run of the equivalent test cases together with
the randomized tests gave 100% line coverage and 99.35% branch coverage
for `src/floatap.cpp`, versus 98.70% previously. Header branch coverage
improved from 73.08% to about 80.77%; the remaining paths include
`checked_sub` with the minimum signed exponent, shifts outside the
supported type range and bit lengths exceeding 2^63-1. Some of those
conditions cannot arise with feasible inputs on 64-bit builds; they must
be distinguished from missing tests, not bypassed merely for coverage.
The CI report on this commit has not yet been run.


## Coverage report of 15 October 2026, 21:51

The attached `Archivio 76(1).zip` reports 100% line/function coverage for
`src/floatap.cpp` and 153/154 covered branches (99.35%).
For `include/numeric/floatap.hpp`, coverage is 27/27 functions,
69/75 lines (92%) and 21/26 branches (80.77%).

A missing branch in `checked_sub` is reachable through a finite,
normalized `FloatAP` value with exponent `INT64_MIN`.
`test_floatap_edge_cases` now tests construction of that value, correct
finite division yielding an exponent near `INT64_MAX`, and exceptional
cases when exact rational conversion or another division would overflow.
These cases require only tiny GMP significands, not enormous allocations.

The uncovered source branch inside `if (mpz_sgn(r))` is shown twice
by LLVM's instrumentation at the same source location. The reported
zero-count arm may be a duplicate branch region; ordinary exact and
inexact quotients are already exercised extensively. We do not mutate
correct arithmetic solely to increase the reported percentage.

The shift-bound and bit-length guards cannot be triggered with feasible
operands on the current 64-bit target: `mp_bitcnt_t` spans at least the
entire nonnegative int64 exponent range, and an actual GMP integer with
more than `INT64_MAX` bits would require unattainable storage.


## Coverage report of 15 October 2026, 21:59 (Archivio 77)

Source-based LLVM coverage from the provided archive:

| Source | Functions | Lines | Regions | Branches |
| --- | --- | --- | --- | --- |
| `src/floatap.cpp` | 16/16 (100%) | 118/118 (100%) | 254/255 (99.61%) | 153/154 (99.35%) |
| `include/numeric/floatap.hpp` | 27/27 (100%) | 71/75 (94.67%) | 74/77 (96.10%) | 22/26 (84.62%) |

This confirms the `INT64_MIN` tests added previously increased
header coverage from 69/75 to 71/75 lines, and from 21/26 to
22/26 branches. The uncovered locations have been inspected in the
LLVM-generated HTML and `branches.txt`:

- `floatap.hpp:64-66` — `shift_amount()` throws for negative values
  or values exceeding `mp_bitcnt_t`. Its callers pass nonnegative
  checked exponent differences. On the 64-bit target `mp_bitcnt_t`
  can represent the full nonnegative signed-int64 range, so neither
  exception arm has a feasible caller through the tested API.
- `floatap.hpp:71-73` — `bit_length()` would throw only if a GMP
  integer had more than `INT64_MAX` binary digits. Such an operand
  is not practicable in memory. These are the **only four uncovered
  executable header lines**.
- `floatap.hpp:114` — two coverage records for the same `is_zero()`
  expression show the zero case absent in one instrumented instance
  but present **40 times** in another. This is not an untested
  arithmetic case: zero classification is exercised.
- `floatap.cpp:23` — two coverage records exist for the same
  `if (mpz_sgn(r))` expression. One reports 0 true hits; the other
  reports approximately **134,000 true** and **21,700 false** hits.
  The true/false cases of the implemented rounding condition are
  demonstrably covered; the missing record appears to be a
  multi-instance instrumentation artifact.

**Decision:** Keep these safety checks and do not introduce impossible
huge operands or weaken invariants to satisfy instrumentation counters.
The realistic kernel correctness tests and coverage pass are now ready
for the next development stage. This does **not** certify the unimplemented
public `Float<ArbitraryPrecision>` interface or correctly rounded
transcendental functions.

## Direct class integration (October 2026)

The temporary kernel class was removed. `Float<ArbitraryPrecision>` now
owns GMP significand storage, precision, exponent, sign and special-value
state directly; `FloatAP` aliases the specialization. The existing tests
exercise its elementary arithmetic methods and special values without a
wrapper. Public constructors for integers/rationals, exact rational
conversion and comparisons have been introduced. They are not yet a
complete replacement for the `FloatMP` interface: transcendental
functions, advanced rounding session controls and optimized extreme
exponent comparisons are outstanding. The old coverage numbers describe
the pre-refactor kernel, so a fresh report is required. Final goal remains
100% line, function and branch coverage without weakening safety checks.


## Coverage report 78 (15 October 2026, 22:27)

After replacing `FloatAPKernel` with `Float<ArbitraryPrecision>`, the
uploaded `Archivio 78(1).zip` reports:

| File | Functions | Lines | Branches |
| --- | --- | --- | --- |
| `src/floatap.cpp` | 23/25 (92.00%) | 150/175 (85.71%) | 165/194 (85.05%) |
| `include/numeric/floatap.hpp` | 29/46 (63.04%) | 73/113 (64.60%) | 22/46 (47.83%) |

This is an expected coverage drop from adding API entry points, not
evidence that earlier randomized arithmetic tests ceased running. The
new `public_comparisons()` and `public_rounding()` C++ tests exercise
ordered and unordered `cmp`, all six comparison operators, ADL-visible
value classifiers, precision-tag NaN/infinity constructors, and explicit
rounded add/sub/mul/div/fma with both operand precision orderings.
The numbers above describe the **incoming** report, not coverage after
this change. The next report must confirm actual changes.

The remaining unrealistically-large shift/bit-length branches need
careful review, and cannot be made reachable through ordinary API inputs
without allocating impossible memory. A final 100% functions/lines/
branches target is recorded, but no current result is being called 100%.


## Coverage report 79 and integer-constructor cleanup (15 October 2026)

The received `Archivio 79.zip` measures:
- `src/floatap.cpp`: functions 25/25, lines 169/175, branches 192/194.
- `include/numeric/floatap.hpp`: functions 46/46, lines 109/113, branches 42/46.

The six uncovered source lines belong to the `int64_t` constructor's
impossible decimal parse failure and its manual cleanup exception path.
The constructor now imports the binary magnitude using GMP, with RAII
cleanup and sign treatment that safely handles `INT64_MIN`. This
eliminates unnecessary string conversion and redundant error handling,
rather than manufacturing test failures just to boost coverage.
Tests exercise zero, `INT64_MIN` and `INT64_MAX` as exact rationals.

The remaining header guards reject impractically large shifts and
significands on 64-bit targets. They are preserved until an appropriate
testability/design solution is established. The requested **100% function,
line and branch** coverage is still the final acceptance criterion,
not a currently achieved measurement. Existing numbers require a fresh
coverage run after the change.

## Windows MSVC and GMP predicate portability

On Windows x64, C++ `long` is 32 bits, whereas `std::int64_t` is 64 bits.
The GMP function `mpq_cmp_si` takes a `long`, so passing `INT64_MIN`
or `INT64_MAX` loses information and triggers MSVC C4244 under `/WX`.
The boundary tests now compare exact rationals against independently
parsed decimal values via `mpq_cmp`, without narrowing. The source also
uses `mpz_cmp_ui(..., 0)` to test nonnegative GMP remainder and
significand for zero. This is equivalent under the current representation
invariant, avoiding incidental sign checks from the `mpz_sgn` macro.
The change requires a new Windows build and coverage run to confirm.


## Coverage report 81 (16 October 2026, 09:43)

The uploaded `Archivio 81(1).zip` contains the following LLVM coverage
figures for the FloatAP files:

| Source | Functions | Regions | Lines | Branches |
| --- | --- | --- | --- | --- |
| `src/floatap.cpp` | 27/27 (100%) | 265/265 (100%) | 337/340 (99.12%) | 196/197 (99.49%) |
| `include/numeric/floatap.hpp` | 46/46 (100%) | 162/166 (97.59%) | 119/124 (95.97%) | 43/47 (91.49%) |

**Important:** The table columns follow LLVM's output order: functions,
regions, lines, branches. Both files have 100% function coverage.
The source has 100% region coverage but still **three uncovered lines**
within a macro expansion; the exact meaning of the counters is determined
from the LLVM HTML, not the previous summary table.

The last missing source branch is an internal implementation path in
`mpz_cmp_ui(remainder, 0)`. Likewise, `is_zero()` in the header
has a duplicate macro expansion branch. Both the zero and nonzero
conditions are already exercised by the tests. They are changed to
`mpz_size(value)` tests, which are equally exact for GMP zero
classification and have no internal signed-comparison branch. The
operation is O(1); it reads the number of limbs. A small Clang
coverage instrumented C++ program confirms both zero/nonzero outcomes
and compilation with warnings as errors. CI coverage on the updated
repository is pending.

Remaining uncovered header lines correspond to the exception guards in
`shift_amount` and `bit_length`. Their safety checks are preserved:
negative/oversized shifts and integers with more than `INT64_MAX` bits
cannot be produced by normal tested paths on 64-bit macOS. Achieving
100% functions/lines/branches must never require weakening these guards
or allocating impossible operands. A separate design review for testing
the checked conversions will be needed.


## Review of coverage 82 and meaningful bounds checks (16 October 2026)

The attached `Archivio 82(1).zip` reports:

| File | Functions | Lines | Branches |
| --- | --- | --- | --- |
| `src/floatap.cpp` | 27/27 (100%) | 333/333 (100%) | 194/194 (100%) |
| `include/numeric/floatap.hpp` | 46/46 (100%) | 115/117 (98.29%) | 41/44 (93.18%) |

The three remaining header branches all concern `shift_amount()` or
`bit_length()`. Their safety arguments are **not** interchangeable:

- The shift helper is called with nonnegative values as an invariant of
  the arithmetic and checked exponent differences. The negative-input
  check was redundant. The check that the nonnegative shift fits
  `mp_bitcnt_t` matters only when that type has fewer value bits than
  the signed exponent. It is now conditional at compile time.
- GMP's `mpz_sizeinbase` yields an unsigned `size_t` length. On 64-bit
  targets this value can exceed `INT64_MAX` in type semantics even
  though allocating an actual integer that large is unrealistic.
  The narrowing protection remains and is extracted into
  `floatap_exponent_from_bit_length(size_t)`, exercised directly by C++
  boundary tests without allocating a giant GMP integer.

This is a correctness-driven simplification, not a blanket removal of
checks to satisfy coverage. The next full CI report must verify the
updated totals on macOS and Windows. The goal remains 100% for
functions, lines, and branches across the completed FloatAP implementation.


## First comparable benchmark slice (16 October 2026)

`BUILD_FLOATAP_COMPARISON=ON` now enables direct, monomorphic
`FloatAP`/MPFR comparison for **only the implemented common API**:
rational construction, copying, comparison, directed add/sub/mul/div,
fused multiply-add, conversion to double and exact rational conversion.
These run at 53, 128, 256 and 1024 bits and use the same operand pool
for each backend. The C++ benchmark emits per-operation timings and
geometric-mean/weighted slowdown. MPFR's transcendental, square-root
and power benchmarks are retained in the MPFR-only configuration.

The differential conformance test similarly compares the two backends
on exact rational results for elementary arithmetic and four rounding
modes. The generic full-`FloatMP` test suites are not enabled for
`FloatAP` until the advanced operations exist. No shim calls MPFR
from `FloatAP`, and there is no dynamic backend dispatch.

The initial native implementation still uses exact GMP integers,
may allocate large temporaries and may be slower than MPFR. This
benchmark is a baseline for optimization, **not a claim of parity
for transcendentals or the entire floating-point interface**.

Build:

```sh
cmake -S . -B build-dual -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_PYTHON_BINDINGS=OFF \
    -DBUILD_NUMERIC_MIGRATION_BENCHMARKS=ON \
    -DBUILD_FLOATAP_COMPARISON=ON
cmake --build build-dual --target benchmark_compare_multiple_precision_backends test_float_multipleprecision_backends
ctest --test-dir build-dual -R '^test_float_multipleprecision_backends$' --output-on-failure
./build-dual/benchmarks/benchmark_compare_multiple_precision_backends 300 7 --summary-only
```

The dual configuration is not yet verified against the full Ariadne
CMake environment by this development session; CMake must not compile
`src/floatap.cpp` twice. After building, recheck all coverage metrics
because every added source branch must eventually be covered.


## Native square root and integer powers (16 October 2026)

After the preceding all-green 100% test and coverage report, the next
implementation block introduces native `sqrt(rounding_mode, FloatAP)`
and `pow(rounding_mode, FloatAP, Nat)`:

- `sqrt` constructs exact dyadic lower/upper bounds using GMP integer
  division and `mpz_sqrt`. It increases working precision until both
  certified endpoints round identically to the destination format.
  It never uses MPFR or an unverified `double` result. A bounded
  resource guard throws rather than returning an uncertain result.
- `pow` computes the exact numerator and denominator powers with GMP
  and rounds the resulting rational once; zero/infinity/NaN and signed
  zero have explicit cases.

These two operations are included in the dual-backend differential
test and benchmark slice alongside elementary arithmetic. More
advanced transcendentals (`exp`, `log`, trigonometric functions, π)
are **not implemented or benchmark-comparable yet**. The
`BUILD_FLOATAP_COMPARISON` option remains temporary until those
native implementations and tests exist. Consequently benchmark
measurements cannot yet be described as the complete set requested.

Coverage must be re-established for the newly added functions; the
previous 100% result covered the earlier implementation only.


## Coverage 83 and exact square-root algorithm (16 October 2026)

The uploaded `Archivio 83(1).zip` shows `floatap.hpp` at 100%
functions (47/47), lines (113/113) and branches (40/40). For
`src/floatap.cpp`, functions are 35/35 (100%), lines 402/406
(99.01%), and branches 220/228 (96.49%). The uncovered cases arise
principally from the first iterative implementation of `sqrt` and
additional special-value combinations in `pow`.

The iterative rational enclosure square root was replaced with direct
`mpz_sqrtrem` on a suitably shifted positive normalized significand.
Its remainder gives exact information about whether the root is integral;
nearest-even is decided by comparing `4*N` with `(2*q+1)^2`. The
midpoint equality cannot occur for integral `N`, since the right side
is odd. Upward increments for nonzero remainder, while downward and
toward-zero truncate; any significand carry is normalized. Exponent
halving uses mathematical floor and checked integer arithmetic. Signed
zero, infinity, NaN and negative input cases remain explicit.

An independent standalone C++ GMP property harness executed 237,000
normalization and directed-enclosure cases across 2..80-bit precision.
It is not a substitute for the complete Ariadne build or full MPFR
differential regression. The previous arbitrary 2048-bit iteration
limit, convergence failure and excessive precision guard are no
longer part of the algorithm. Edge tests now include multiple rounding
modes and signed special-value cases. The full coverage result must
be regenerated for this updated version.


## Coverage report 84 (16 October 2026, 10:58)

The uploaded `Archivio 84(1).zip` reports full function coverage for
`FloatAP`. `include/numeric/floatap.hpp` reaches 47/47 functions,
113/113 lines and 40/40 branches. The new `src/floatap.cpp` has
32/32 functions, 356/361 lines and 241/246 branches covered.

The missing executable lines are the invalid square-root rounding-mode
exception and the significand carry normalization. The new edge-case
tests exercise both, together with inexact rounding toward zero. In
particular, sqrt(255) at eight bits rounds upward to 16, necessitating
the carry and exponent normalization; the precise boundary was verified
with a standalone C++/GMP integer-square-root check. Three additional
uncovered branch records in the Rational constructor were duplicated
LLVM counters with zero executions while the corresponding other
instances already demonstrate both sides (NaN, infinities and signs).
No production code is altered merely to inflate such duplicate counters.
A fresh full CI coverage run must confirm the new totals.


## Native exp/log with exact rational enclosures (16 October 2026)

The first transcendental block adds native `exp` and `log` in
`src/floatap_exp_log.cpp`, registered in `NUMERIC_SRC`.
They do not invoke MPFR and do not use floating-point approximations
to determine the rounding bit. Taylor sums and their rigorous remainder
bounds are represented by GMP rationals.

For `exp(x)`, repeatedly halve the argument until `|x| <= 1/2`;
the Taylor remainder is bounded by twice the first omitted term.
Positive interval squaring then reverses the argument reduction.
For `log(x)`, scale the positive argument by powers of two into
`[1/2,2]` and use the exact atanh series
`log((1+z)/(1-z)) = 2 sum(z^(2n+1)/(2n+1))` for `|z| <= 1/3`,
combining an enclosure for `log(2)` when necessary.

Working precision increases until both interval endpoints round to the
same `FloatAP` at the requested precision and rounding mode. This
makes the result certifiably correctly rounded whenever the procedure
terminates. Certain extreme inputs are intentionally rejected by
resource guards instead of returning an uncertified value. This first
correctness-oriented implementation allocates large GMP rationals and
is **not yet performance optimized**.

The differential MPFR oracle test now includes `exp` and `log` for
four rounding modes and the benchmark includes both operations.
`sin`, `cos`, `tan`, `atan`, `asin`, `acos` and `pi` remain
missing, so the dual-backend opt-in must not yet be removed. Newly added
branches will need a fresh coverage report.


## Permanent dual-backend benchmark configuration (16 October 2026)

The benchmark infrastructure is now **unconditional**: both
`FloatMP` and `FloatAP` are compiled, timed and subjected to
differential checks in every standalone build. No CMake variables
`BUILD_FLOATAP_COMPARISON` or
`BUILD_NUMERIC_MIGRATION_BENCHMARKS`, and no
`ARIADNE_NUMERIC_ENABLE_FLOATAP` preprocessor macro, are used by
the benchmark or the dedicated differential test. The configuration
submodule registers `benchmarks/` once.

The C++20 compile-time `HasTrigonometricOperations<FloatAP>`
capability check is a **temporary completeness gate**, not a
build-configuration flag. Until native `sin`, `cos`, `tan` and
`atan` exist, the benchmark measures both implementations on
the same common operation set (including `sqrt`, `exp`,
`log`, and integer `pow`), and emits a diagnostic explaining that
the complete trigonometric comparison is not yet available.
When those native signatures exist, the benchmark automatically
includes the operations for **both** backends. The differential
test uses the analogous `HasAdvancedFunctions` check for
`sin`, `cos`, `tan`, `asin`, `acos`, `atan`, and π.

This change does not claim that missing transcendental functions
have been implemented. No dummy implementations or MPFR forwarding
have been introduced into `FloatAP`.

Typical build, without any benchmark/comparison option:

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --target benchmark_compare_multiple_precision_backends test_float_multipleprecision_backends
ctest --test-dir build-release -R '^test_float_multipleprecision_backends$' --output-on-failure
./build-release/benchmarks/benchmark_compare_multiple_precision_backends 300 7 --summary-only
```

Historical references to build flags in earlier sections describe
prior states of this development branch, not the current build contract.


## Complete native trigonometric interface and unconditional comparison (16 October 2026)

The native `Float<ArbitraryPrecision>` specialization now implements the
entire set of operations requested by the comparison benchmark:
`sin`, `cos`, `tan`, `atan`, `asin`, `acos`, and `pi`, in
addition to the existing `sqrt`, `exp`, `log`, `pow`, elementary
arithmetic, conversions and comparisons. Their definitions live in
`src/floatap_trigonometric.cpp`, not in an MPFR wrapper. The benchmarks
and differential tests no longer have temporary capability gates:
they instantiate **both** backends and exercise the full operation set.

Native trigonometry uses GMP rational enclosures and certified remainders:
- Taylor series for sine and cosine on `[-1, 1]`, with error bounded by
  the first omitted alternating term. Larger arguments are halved and
  reconstructed by double-angle identities with interval arithmetic.
- Arctangent is reduced to `[-1/2, 1/2]` with the exact identities
  `atan(x) = pi/4 + atan((x-1)/(x+1))` and
  `atan(x) = pi/2 - atan(1/x)`; each reduced series has a rigorously
  bounded alternating remainder.
- Pi uses Machin's identity
  `pi = 16 atan(1/5) - 4 atan(1/239)`.
- Inverse sine uses `asin(x) = 2 atan(x/(1+sqrt(1-x^2)))`, where
  the rational square-root interval is computed by integer division
  and `mpz_sqrt`. Inverse cosine uses `pi/2 - asin(x)`.
- Tangent divides certified sine/cosine intervals, explicitly rejecting
  an unresolved denominator enclosure containing zero.

For every finite nontrivial input, the implementation increases its
internal precision until **both** endpoints round to the same FloatAP
under the requested mode. No approximate final result is accepted if
the enclosure straddles a rounding boundary. Domain errors and special
values are treated explicitly; `atan(±inf)`, `asin(±1)`, and
`acos(-1)` use the same certified pi/2 or pi enclosures, not large
finite surrogates. Resource guards may throw for extremely large inputs
rather than return uncertified values.

This initial algorithm is correctness-oriented, with expensive rational
temporaries and an adaptive precision loop. It is not optimized for
benchmark speed, and does not promise termination for every theoretical
hard-to-round input. Full Ariadne CI compilation and MPFR differential
testing are required **before performance numbers can be trusted**.
The previously achieved 100% coverage applies only to the earlier
implementation; all newly introduced functions, lines and branches must
be covered in subsequent reports.


## Independent FloatMP and FloatAP test targets (16 October 2026)

The former dual-backend `test_float_multipleprecision_backends` target
is replaced with `test_floatmp_conformance`, which imports only the
FloatMP API. The independent `test_floatap` and
`test_floatap_edge_cases` targets use the native FloatAP and exact
GMP/Rational references, never `floatmp.hpp` or MPFR. This structure
allows the FloatMP tests and library backend to be deleted later without
modifying FloatAP tests. The explicit dual-backend *performance*
comparison remains the benchmark's responsibility.

Clang compilation failures in `test_floatap_edge_cases.cpp` caused by
unqualified `16_bits` literals were fixed by explicitly constructing
`ArbitraryPrecision(std::size_t{16})`. Missing default handling for
new `half_pi` and `negative_half_pi` switch cases in
`floatap_trigonometric.cpp` was fixed as well.


## Benchmark-driven optimization: direct binary comparison (16 October 2026)

The initial benchmark at 53 bits measures `cmp` at roughly 194 ns
for FloatAP versus 12 ns for FloatMP (about 16x), and at 256 bits
roughly 307 ns versus 12 ns (about 26x). The previous implementation
converted both operands to GMP rationals and compared their canonical
numerators and denominators. This is unnecessary because FloatAP
stores signed binary significands and exponents directly.

The optimized comparison now classifies NaN/infinities/signed zeros,
compares signs, then compares the most significant binary positions
from the significand bit lengths and exponents. Matching positions
are resolved by truncating the lower-exponent significand and
examining any discarded nonzero low bits; otherwise no GMP temporary
is required. The exponent-gap calculation uses unsigned subtraction,
avoiding signed overflow across the full int64 range, and large
exponent gaps never trigger enormous shifts. The result is inverted
for negative operands. This is independent of MPFR and does not
convert to `Rational`.

An independent GMP rational oracle test compares thousands of
randomized mixed-sign, mixed-scale, mixed-precision operands against
the new direct comparator, including signed zeros and equality.

This addresses one measurable overhead only. The remaining large
slowdowns in `exp`, `log`, `atan`, `asin` and `acos` come
primarily from full-precision `mpq_t` operations within each Taylor
series iteration. The next optimization stage should use bounded
fixed-point `mpz_t` arithmetic with certified interval error,
binary-splitting / faster convergence, and caching reusable constants
such as pi (without sacrificing rounding correctness).
No post-change benchmark improvement or full coverage is claimed
until verified by a fresh build and run.


## Broad precision-regression suite and GCC namespace portability (16 October 2026)

The uploaded coverage report 86 shows that the earlier 100% goal is not
yet achieved after the transcendental implementation. In LLVM's report,
the fields are functions, regions, lines, and branches respectively:

| File | Functions | Lines | Branches |
| --- | --- | --- | --- |
| `src/floatap.cpp` | 31/33 | 446/480 | 267/290 |
| `src/floatap_exp_log.cpp` | 27/27 | 113/136 | 39/58 |
| `src/floatap_trigonometric.cpp` | 43/43 | 189/219 | 94/122 |
| `include/numeric/floatap.hpp` | 47/47 | 113/113 | 40/40 |

A dedicated `test_floatap_precision` executable covers precisions
2, 7, 24, 53, 113, 256, and 1024 bits for elementary operations under
all four rounding modes, with independently constructed exact-rational
answers. It also checks selected transcendental results at 7, 24, 53,
and 113 bits for ordered directed rounding and exact identities. These
tests are native FloatAP tests, not FloatMP dependency tests. Directed
rounding ordering alone is *not* a proof of correctly rounded
transcendentals; more independent reference vectors remain necessary.

On GCC 13 with `-Wshadow -Werror`, nested enumerators `downward`,
`upward`, `toward_zero` and `infinity` conflict with similarly
named Ariadne namespace variables when the FloatMP header is present.
The enum class values are renamed to `ToNearest`, `Downward`,
`Upward`, `TowardZero` and `Finite`, `Infinity`, `NaN`.
The public `ROUND_*` constants retain their names and behavior.
The change is a source API rename for clients who used the nested
enumerators directly; benchmark calls through `ROUND_*` are unaffected.

The full CI build and a fresh coverage report are required to
validate this expanded test suite and renaming. The final 100%
coverage and full-precision correctness targets remain outstanding.


## Reusing the existing Ariadne Float test suites (16 October 2026)

The authoritative behavioral contract for the MPFR replacement is the
unchanged, macro-based `tests/test_float.cpp` and
`tests/test_float_rounding.cpp` (`ARIADNE_TEST_*`). New, smaller
FloatAP-only tests supplement these suites but are **not substitutes**
for instantiating the historical `TestFloat<ArbitraryPrecision>` and
`TestFloatRounding<ArbitraryPrecision>` classes. The existing
`ARIADNE_NUMERIC_ENABLE_FLOATAP` guards in those legacy files cannot
simply be enabled yet: the native specialization still lacks some
constructors, global rounding-mode overloads, mixed arithmetic,
`RawFloat` registration, bounds/balls, and conversion/streaming APIs.
Do not remove or skip assertions to obtain a green test run.

The first API-compatibility block introduces exact conversions from
`Dyadic`, `ExactDouble`, and `TwoExp`, mutable precision with
round-to-nearest semantics, raw access, a thread-local current rounding
mode, stream insertion, exponent shifts and elementary unary helpers.
These methods exist to satisfy contracts already expected by Ariadne,
not as special test shims. They require compilation and regression
verification on macOS, Ubuntu, and Windows. The historical tests have
**not yet been enabled**; the complete API parity task remains open.

The shared `test_float.cpp` now also contains direct
`ARIADNE_TEST_*` checks for this first compatibility block. All
existing FloatDP/FloatMP assertions remain in place. These checks are
**not** presented as full test parity: the generic
`TestFloat<ArbitraryPrecision>` and
`TestFloatRounding<ArbitraryPrecision>` instantiations still require
more production APIs. No test is marked passed merely because its
instantiation is not yet enabled.


## Correcting class-scope lookup in rounded helper methods (16 October 2026)

The first API-parity block introduced inline hidden friends
`sqr(mode,x)`, `rec(mode,x)`, `med(mode,a,b)`, and
`rad(mode,a,b)`. The original unqualified calls to
`mul(mode,...)`, `div(mode,...)`, `add(mode,...)`, and
`sub(mode,...)` mistakenly resolved to class-static implementation
functions, which require four parameters, rather than to the desired
namespace-scope hidden friends. The error prevented compilation of
every translation unit including `floatap.hpp`.

These helpers now explicitly invoke the class-static four-parameter
kernels with the appropriate destination precision and rounding mode.
For mixed-precision operands, `med` and `rad` retain the minimum
operand precision used by the public three-argument arithmetic friends.
All original tests remain unchanged; the Ariadne-macro compatibility
checks in `test_float.cpp` are extended for the four helpers.

In the same new test block, the invalid two-integer `Dyadic`
constructor was replaced with a `Dyadic` constructed from
`ExactDouble(1.25)`. Full build and the historical templated
`TestFloat<ArbitraryPrecision>` and
`TestFloatRounding<ArbitraryPrecision>` parity checks are still
outstanding; the new helpers alone do not establish API parity.


## Two incorrect regression expectations exposed by CTest (16 October 2026)

The first build with the new compatibility checks completed compilation,
but `test_float` failed because `rec(5/4)` was compared for *exact*
equality to `4/5`. This rational has an infinite binary expansion, so
no finite-precision FloatAP result can equal it exactly. The original
test now checks all four rounding modes against exact Rational
inequalities: downward < 4/5 < upward, toward-zero equals downward for
this positive result, and nearest equals upward since the upper
representable neighbor is closer to 4/5.

The new precision test also expected `sqrt(9) = 3` at just two
significand bits. At this precision, the input **9 first rounds** to
8 in nearest, downward and toward-zero modes, and to 12 in upward
mode. Thus `sqrt(8)` rounds to 3 in nearest and 2 in downward or
toward-zero mode, whereas `sqrt(12)` rounds to 4 in upward mode.
The test checks each of these exact mode-specific expectations,
retains `sqrt(9) = 3` only when 9 is exactly representable (four or
more bits), and verifies `sqrt(4) = 2` across every tested precision.
Additional exact Rational square inequalities independently enclose
`sqrt(2)`.

Neither correction changes production numerical algorithms or deletes
a meaningful behavioral test. Both distinguish a test expectation
error from an implementation bug. Full Ariadne CTest and updated coverage
are pending after this commit.

The mode-specific two-bit expectations and nearest rounding for
`rec(5/4)` at 128 bits were checked independently with a small C++/GMP
arithmetic program. That check is not a full Ariadne compile or CTest
run. The regression is pending CI validation.


## Removing the noncanonical compatibility test (16 October 2026)

The ad-hoc `test_floatap_compatibility()` block has been removed
from `tests/test_float.cpp`: it was not the established Ariadne
behavioral contract and introduced expectations inconsistent with
finite binary precision. No historical `TestFloat<...>` assertions
were modified. The existing templated suites
`TestFloat<ArbitraryPrecision>` and
`TestFloatRounding<ArbitraryPrecision>` are the only acceptable
end-to-end compatibility criterion; they are still gated off because
current `FloatAP` does not yet meet their constraints.

In particular, `IsARawFloat<FloatAP>` is not yet specialized in
`numeric/float.decl.hpp`, and the three-way rounded arithmetic
and bounds/ball integration are not complete. Activating the whole
suite before those implementation changes would introduce a known
build failure, not meaningful test coverage. Implement the missing
APIs in production sources, then instantiate the unchanged legacy
suites unconditionally and confirm them by actual CI.


## Historical Float test activation (16 October 2026)

The existing `tests/test_float.cpp` and
`tests/test_float_rounding.cpp` instantiate their **original**
`TestFloat<ArbitraryPrecision>` and
`TestFloatRounding<ArbitraryPrecision>` templates unconditionally,
with all `ARIADNE_TEST_*` assertions untouched. The former
`ARIADNE_NUMERIC_ENABLE_FLOATAP` test guards are removed, including
the previously guarded cross-precision conversion specialization.
`FloatAP` is registered as an `ARawFloat` in `float.decl.hpp`
so the existing Ariadne generic floating-point abstractions can
participate in overload resolution.

This change **enables the historical suites but does not establish
that they compile or pass**. Notably, the native backend still lacks
some of the existing FloatDP/FloatMP contract: conversions and
assignments, round-tag overloads, bounds/ball operations, parsing,
minimum/maximum exponent semantics, and several transcendental
helpers. Each missing API must be implemented rather than excluding
or weakening a test. The full build and CI results are the criterion
for completing this integration.


## First compilation failures with the historical suites (16 October 2026)

The user-supplied Clang/macOS build log shows that enabling the unmodified
`TestFloat<ArbitraryPrecision>` and
`TestFloatRounding<ArbitraryPrecision>` template instantiations reveals
genuine missing production interfaces. Among the first errors are:
rounding tags (`down`, `up`, `near`, `approx`, `rounded`) not
convertible to FloatAP's private rounding abstraction; missing
constructors from `FloatDP` and from `FloatAP` at a different precision;
missing exact comparisons with `Dyadic` and `Rational`, integer
casting, assignments, and `operator Dyadic`. Compiler diagnostics are
truncated to 20 errors per translation unit, so the full remaining
interface gap is not yet known.

The first integration pass introduces a small FloatAP-specific rounding
value type that accepts Ariadne's existing tag classes while retaining
the native GMP kernel's four rounding choices. Current-rounding tokens
read the thread-local rounding state; explicit rounding-tag operations
do not change it. It also adds exact cross-precision constructors,
`FloatDP` conversion, assignments, exact cross-type comparisons and
integer/Dyadic casts. Dyadic conversion uses a sufficiently precise
GMP `mpf_t`, not `double`, so finite large or high-precision values
remain exact.

The existing tests and CMake selection remain untouched. This pass
is deliberately limited to the compile errors actually identified;
it does **not** claim that the full historical tests compile or pass.
Later diagnostics will likely involve the still-missing validated
bounds/ball operations and float exponent-range API. No test gates
or weakened assertions were reintroduced.


## Consolidation into the original Float test suites (16 October 2026)

`tests/test_float.cpp` now shares its multiple-precision conversion,
arithmetic and comparison paths with both `MultiplePrecision` and
`ArbitraryPrecision`: generic operands and expected values use
`Float<PR>` and `PR`, while actual MPFR-only internals such as
`get_mpfr()`, `MPFR_RNDN`, raw pointers, built-in MPFR precision
conversions and MPFR default-precision management remain in explicit
`Same<PR,MultiplePrecision>` subbranches. The test assertions for
generic behavior have not been disabled.

`tests/test_float_rounding.cpp` shares the cross-precision conversion
test using `Float<PR>` for both multi-precision backends. Its original
templated classes, and every original rounding assertion, remain
enabled. A `test_precision_regressions()` method using Ariadne's
`ARIADNE_TEST_*` macros integrates the key tests formerly in
`test_floatap.cpp`, `test_floatap_edge_cases.cpp`,
`test_floatap_precision.cpp` and `test_floatmp_conformance.cpp`:
seven precisions between 2 and 1024 bits, all four rounding modes,
exact Rational oracles for basic arithmetic and FMA, low-precision
square-root edge cases, transcendental enclosures and identities,
varied-scale comparisons and nonfinite values. These tests execute
for FloatMP **and** FloatAP; they are not fallback FloatAP-specific suites.
Coverage-only tests of internal implementation details are not
substitutes for this shared mathematical contract.

The four dedicated test targets and their source files are removed
from CMake and the branch. The only primary float tests are the
original `test_float` and `test_float_rounding`, which keep the
historical DoublePrecision and MultiplePrecision instantiations as well
as the ArbitraryPrecision instantiations.

This changes the test arrangement but does not assert a passing build:
`FloatAP` still lacks parts of the historic API, including formatting,
input parsing, certain mixed operations, and round-tag helpers. The
original tests should continue to expose these as compile failures
until the production API is implemented. Do not restore gated or
weakened assertions for an artificially green regression.
