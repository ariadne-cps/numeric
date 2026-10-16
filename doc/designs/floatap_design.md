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
