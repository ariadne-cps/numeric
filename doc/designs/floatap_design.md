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

## Representation

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
large temporaries. Optimization must preserve **correctly rounded**
results under all four supported rounding modes, including difficult
cases; no optimization may silently weaken the enclosure/certification
contract.

Priority based on the 16 October 2026 same-host benchmark:

1. **Transcendentals first:** replace repeated exact `mpq_t` Taylor
   operations with fixed-point `mpz_t` series evaluation and explicitly
   propagated outward errors; accelerate argument reduction and the
   `atan`/`asin`/`acos` dependency chain.
2. Evaluate binary/rectangular splitting, reusable certified constants
   (particularly `pi` and `log(2)`), and efficient argument reduction
   *after* obtaining reliable per-operation profiles. Reuse must not
   depend on mutable global rounding state.
3. Optimize `add`, `mul`, `div`, `fma`, conversion, and comparison
   using limb-oriented fast paths, bounded exponent alignment,
   preallocated scratch integers, and fewer temporary allocations.
4. Profile separately at 53, 128, 256, and 1024 bits, including
   exceptional, exact, boundary, and hard-to-round cases.

The first fixed-point step for `log` is implemented, but an observed
speedup is **not** yet established; see the decision record below.
The long-term intent remains to avoid GMP as well as MPFR in the native
backend. MPFR is permitted in oracle tests and benchmarks, **not**
in native `FloatAP` arithmetic.

## Verification and benchmarks

`benchmark_compare_multiple_precision_backends` is the C++
dual-backend comparison executable. It exercises both concrete
backends unconditionally in current builds and reports
`FloatAP / FloatMP` timings for individual operations:

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

## Status and limitations (updated 16 October 2026)

- `FloatAP = Float<ArbitraryPrecision>` is a direct, usable public
  specialization with native basic arithmetic, conversions, square root,
  exponential, logarithm, trigonometry and `pi`, along with substantial
  compatibility with Ariadne's shared Float suites. This is **not**
  yet a drop-in replacement for every FloatMP application.
- Implementation now resides in **one translation unit**:
  `src/floatap.cpp`. The former `floatap_exp_log.cpp` and
  `floatap_trigonometric.cpp` were merged and removed from CMake.
  The transcendental code shares one exact `Fraction` and `Interval`
  representation while retaining function-specific internal namespaces.
- The original kernel is based on GMP integers; most transcendental
  series are still based on expensive GMP rational intermediates.
  The first exception is the new fixed-point `log_series`.
  Extreme exponents, range conventions, correctness of directed
  conversions and hard-to-round behavior need continued verification.
- The dual-backend comparison is already running. Timings supplied on
  16 October 2026 show approximately **5x** slower basic additions
  at 53 bits and **352x / 1,488x / 8,987x** slower `log` at
  53 / 128 / 256 bits **before** the fixed-point change.
  These are baseline measurements, not a measured improvement.
- CI for commit `8eb5c97` has succeeded on the listed Linux and
  macOS configurations, but Windows Release and Codecov checks are
  failing; full cross-platform validation is **not** complete.
- Historical sections below describe decisions *at the date made*.
  A historical statement that an operation was absent or lived in
  a separate file must not be mistaken for current implementation status.

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


## Compiler diagnostics after consolidation (16 October 2026)

The 3595-line macOS Clang build log revealed a common blocking
regression in the native `Rounding` wrapper: in-class constant
initialization of the wrapper's four `ROUND_*` objects was rejected by
Clang and GCC before the nested wrapper class was complete. The
wrapper also exposed implicit conversion to its underlying enum and
an equality operator taking two wrappers; these created ambiguous
equality expressions in `floatap.cpp`.

The first targeted fix declares the four class constants inside
`FloatAP` and defines them as inline `constexpr` after the class
definition. It replaces the implicit enum conversion with explicit
`value()` access, used by the native kernels' switches; equality
against enum values is then unambiguous. Forward declarations of
Ariadne's rounding tags use the matching `struct` keyword, avoiding
cross-platform mismatched-tag diagnostics.

The historical rounding suite additionally had an accidental use
of `Float<PR>` inside a member method already declaring `using
Float=RawFloat<PR>`; in this name lookup context the spelling is
invalid. The conversion tests now use the established alias `Float`
without changing what they assert.

The attached log also reports missing `FloatAP` production APIs
(string construction/parsing, default precision, exponent limits,
in-place operations, generic float concepts, and others). These
remaining diagnostics are not resolved by this targeted compilation
fix. Both original float test suites remain enabled without
suppressing assertions. A complete project rebuild is still required.


## Performance redesign: certified fixed-point logarithm (16 October 2026)

**Change:** `8eb5c97722fa41532174470022cc9cddeae454af`,
`src/floatap.cpp`. This is the first response to the benchmark
finding that the GMP-rational transcendental prototype is far slower
than MPFR. Earlier `log_series` iterations used `mpq_t` for every
Taylor term and sum; each rational multiplication/addition could
grow numerator and denominator and require canonicalization.

**Decision:** keep the externally visible `FloatAP` representation,
all four rounding modes, the certified-enclosure refinement loop,
and the reduction `x = 2^k y` with `y in [1/2, 2]`. Replace only the
internal evaluation of

```text
log((1+z)/(1-z)) = 2 * sum(n=0..infinity) z^(2n+1)/(2n+1)
z = (y-1)/(y+1), |z| <= 1/3
```

with **bounded fixed-point GMP integer arithmetic**. Work at
`B = work + 32 + bit_width(N)` fractional bits with
`N = floor(work/3) + 4` terms. Enclose the nonnegative magnitude
`|z|` between directed floor/ceiling scaled integers; maintain
lower/upper fixed-point enclosures for its odd powers using directed
squaring, multiplication, and integer division. Sum the nonnegative
terms with directed division, apply the factor of two, and restore
the sign by swapping/negating endpoints for negative `z`.

**Remainder argument:** for `|z| <= 1/3`, the omitted part of the
doubled atanh series after `N` terms is bounded by

```text
2 * sum(n=N..infinity) |z|^(2n+1)/(2n+1)
    <= (9/4) * 3^(-2N-1)
    < 2^(-work), when N = floor(work/3) + 4.
```

The implementation adds `2^(-work)` to its fixed-point upper
endpoint; the directed integer operations enclose the included
terms. A zero argument returns an exact zero enclosure. The outer
`certified` loop still refines until **both** endpoints round to the
same requested result, rather than accepting an approximate answer.

**Trade-offs and current limits:** the series is now linear in a
bounded count of fixed-point steps, rather than manipulating exact
rational denominators at each step. This does **not** yet add fast
argument reduction, binary splitting, reuse of `log(2)`, or a
scratch allocator. The `log_interval` reconstruction still handles
powers of two and `log(2)` using rational enclosures; other
transcendentals still use the exact rational prototype. No
post-change speed ratio is available yet, and the bound/protocol
must still be stress-tested against MPFR at all rounding modes and
near rounding boundaries, including cases that trigger refinement.

**Performance baseline supplied before this change** (median ns,
300 iterations, seven samples, identical benchmark operations):

| Operation | Bits | MPFR (ns) | FloatAP (ns) | FloatAP / MPFR |
|---|---:|---:|---:|---:|
| add, nearest | 53 | 69.03 | 335.28 | 4.9x |
| log, nearest | 53 | 660.28 | 232180.14 | 352x |
| log, nearest | 128 | 1551.39 | 2309173.75 | 1488x |
| log, nearest | 256 | 2671.25 | 24007044.58 | 8987x |
| atan, nearest | 53 | 667.50 | 232492.08 | 348x |
| atan, nearest | 128 | 2683.61 | 2245674.17 | 837x |

These timings justify prioritizing the transcendental algorithms;
they are **not** measurements of `8eb5c97`. Retain the original
same-backend, same-input and same-rounding benchmark methodology,
and validate output equivalence before claiming an optimization.

**First post-change benchmark (16 October 2026).** These are user-supplied
same-workload median timings with 300 iterations and seven samples.
The log speedup is computed against the preceding supplied baseline;
the MPFR comparisons use the contemporaneous new run:

| Bits | Old FloatAP log (ns) | New FloatAP log (ns) | Speedup | New MPFR log (ns) | Remaining slowdown |
|---|---:|---:|---:|---:|---:|
| 53 | 232180.14 | 3870.00 | 60.0x | 657.50 | 5.9x |
| 128 | 2309173.75 | 6245.14 | 369.8x | 1520.56 | 4.1x |
| 256 | 24007044.58 | 10135.28 | 2368.7x | 2609.31 | 3.9x |

The baseline FloatAP `atan` times remain approximately 237687 ns
at 53 bits and 2262583 ns at 128 bits; `asin`/`acos` remain even
more expensive. These were **not** modified by the logarithm change.
The subsequent fixed-point arctangent decision is recorded below.
No speedup should be inferred for other functions from this table.

**Cross-platform correction:** commit `aa31f47` guards the series
index in `exp_series` by explicitly comparing unsigned indices to
a safely converted signed maximum minus one; this avoids MSVC
C4389 under `/WX` and ensures that `n+1` can be converted to
`long`. The compile and regression result of this commit must
still be confirmed on Windows.

**Next investigations:** replace `small_atan` and the
`asin`/`acos` dependency chain with rigorously bounded fixed-point
methods; examine Machin/pi caching, fast angle reduction and
rectangular/binary splitting; then optimize the ordinary arithmetic
kernel's GMP allocations and repeated `from_ratio` path.
Every future algorithmic change must update this decision record
with its proof obligations, measurements, and regressions.

## Performance redesign: certified fixed-point arctangent (16 October 2026)

**Change:** `src/floatap.cpp` function `small_atan`. Following the
successful fixed-point logarithm redesign, replace the previous
`mpq_t` Taylor sum, which repeatedly constructed and reduced exact
rational terms, with fixed-point GMP integer lower and upper bounds.

**Contract and scope:** for `|x| <= 1/2`, including the Machin
arguments `1/5` and `1/239`, evaluate

```text
atan(x) = x - x^3/3 + x^5/5 - ...
```

using the magnitude `|x|` with `B = work + 32 + bit_width(N)`
fractional bits, where `N = floor(work/2) + 4` is an iteration
bound. Form `floor(|x| 2^B)` and `ceil(|x| 2^B)` once, then
propagate both bounds through the squared argument and odd powers
using GMP fixed-point multiplication with outward truncation.
For each positive term, take the floor and ceiling of the scaled
power divided by its odd denominator. Alternate subtraction and
addition in **opposite directions** so that cancellation never
narrows the enclosure spuriously.

**Certified alternating remainder:** after accumulating terms
`n = 0, ..., k`, the exact tail has the sign of term `k+1`,
and its magnitude is bounded by
`|x|^(2k+3)/(2k+3)`. The algorithm obtains a conservative scaled
integer upper bound from the upper power and uses it to widen
**only the appropriate endpoint**. It stops early when the bound
is at most `2^(-work)`, or after `N` terms with the actual
certified tail (even if this still exceeds the target width).
The sign of a negative argument is restored by negating and
reversing the endpoints. Integer outward rounding guarantees
inclusion independently of the number of guard bits; the guard
precision serves efficiency, not correctness.

The existing `certified` refinement loop remains unchanged: it
returns only after both enclosure endpoints round identically in
the requested rounding mode. Thus the change is compatible with
the four rounding modes and does not forward numerical work to MPFR.

**Expected propagation:** `pi_interval` uses two calls to
`small_atan` in Machin's identity; `atan_interval` uses it
after argument reduction; `asin_interval` calls `atan_interval`
twice on the endpoints of its rational square-root enclosure;
`acos` depends on `asin` and `pi`. All therefore benefit
without changing their identities. This update does **not**
accelerate the rational Taylor evaluation of `exp`, `sin`,
or `cos`, and does not cache `pi`.

**Pre-commit verification:** a standalone arithmetic check compared
the enclosing interval against high-precision arctangent values
for 208 rational arguments, including `+/-1/2`, `1/5`,
`1/239`, and very small values, with work precision 32, 53,
85, 128, 256, and 512 bits. This is an algorithm-level
cross-check, **not** a complete Ariadne build or MPFR differential
regression. Numerical speedup and the final cross-platform
regression remain to be measured.

**Literature and next steps:** Arb documents fixed-point Taylor
arctangent evaluation, table-based argument reduction, and
rectangular splitting at low/medium precision; see
https://arblib.org/arb.html and the MPFR algorithm notes at
https://www.mpfr.org/algo.html. Investigate reducing transformed
arguments to below `1/16`, caching certified Machin constants
at appropriate working precision, and rectangular/binary splitting
when iteration cost dominates. Retain hard-to-round and all-mode
differential comparisons as acceptance criteria.


## Performance redesign: certified fixed-point exponential (16 October 2026)

**Decision:** replace the per-iteration GMP `mpq_t` Taylor sum in
`FloatAPExpLogDetail::exp_series` by bounded fixed-point `mpz_t`
arithmetic. Leave the outer argument-halving and interval-squaring
strategy, special-value handling, all four rounding modes, and the
`certified` endpoint agreement loop intact.

For `|x| <= 1/2`, compute `exp(|x|)` from its **positive-term**
series `sum(a^n/n!, n >= 0)` at fixed-point scale
`Q = 2^B`, with `B = work + 32 + bit_width(N)` and a conservative
iteration cap `N = floor(work/2)+64`. Enclose `a` once using
`floor(aQ)` and `ceil(aQ)`, and propagate floor/ceiling bounds
through each multiplication and division by `n`. The lower and
upper accumulators enclose the included partial sum. For
`0 <= a <= 1/2`, the infinite positive tail after term `n`
is bounded by **twice the first omitted term**
`2 a^(n+1)/(n+1)!`. The implementation computes an outward-rounded
upper bound of this next term and adds twice that bound to its
upper accumulator. It can stop when this certified remainder is
at most `2^(-work)`; it also stops at the conservative cap using
the *actual* certified remainder rather than assuming it is small.

For negative input, form the exact reciprocal interval
`[1/upper, 1/lower]` of the positive result. Positivity makes
this order reversal sound, and repeated interval squaring after
argument reduction remains monotone. The final endpoints are
rounded only by the preexisting certification loop. Consequently,
fixed-point truncation does **not** reduce the correctness contract;
the guard bits affect efficiency, not the validity of the enclosure.

**Independent algorithm check:** 209 rational arguments in
`[-1/2,1/2]`, including `+/-1/2`, tiny arguments and randomized
fractions, were checked against high-precision reference evaluations
for `work = 8, 32, 53, 85, 128, 256, 512, 1024`. All
evaluations were enclosed. This is not a substitute for the C++
MPFR differential suites or a full cross-platform regression.

**Pre-change benchmark provided by the user:** `exp` takes
25.08 us at 53 bits, 70.03 us at 128 bits, 250.45 us at
256 bits, and **11.56 ms at 1024 bits**; MPFR takes
0.79, 1.12, 1.62 and 6.61 us respectively. These are baseline
figures measured *before this change*, not achieved speedups.
A new same-machine C++ benchmark and all-mode MPFR comparisons
are required before accepting performance claims.

**Follow-up:** apply a jointly evaluated certified fixed-point
sine/cosine series; then optimize argument reduction and the
rational interval squaring operations. Arb/FLINT documents
fixed-point Taylor/rectangular-splitting exponential methods
(`_arb_exp_taylor_naive`, `_arb_exp_taylor_rs`) and more
aggressive reduction of small exponential arguments:
https://flintlib.org/doc/arb.html . No MPFR function is called
by the new native code.


## Performance redesign: certified simultaneous sine and cosine (16 October 2026)

**Implementation:** `src/floatap.cpp`, `small_sine_cosine`.
Replace the rational GMP Taylor iteration by a **single pair of
outward-rounded fixed-point integer series** for both functions. The
previous `mpq_t` approach created and normalized exact fractions
at every step, causing a steep increase in cost with precision.

**Series and argument domain:** the current argument-halving path
reduces to `|x| <= 1` before invoking the kernel. Evaluate the
positive magnitudes of both Taylor series simultaneously:

```text
sin(a) = sum(k >= 0) (-1)^k a^(2k+1)/(2k+1)!
cos(a) = sum(k >= 0) (-1)^k a^(2k)/(2k)!,   0 <= a <= 1
```

Use a binary fixed-point scale `Q = 2^B`, with
`N = floor(work/2)+48` and `B = work+32+bit_width(N)`.
Enclose `a = |x|` by floor/ceiling integers and its square
by directed integer multiplication and division by `Q`.
Maintain **both** lower/upper bounds for the nonnegative sine
and cosine term magnitudes; generate successive terms by the
factor `a^2 / ((2k+2)(2k+3))` and
`a^2 / ((2k+1)(2k+2))`, respectively, always rounding outward.
The loop accumulates signed terms using the proper opposite
endpoints for subtraction. The only rational conversions happen
at the end, not in the Taylor loop.

**Error certification:** for `0 <= a <= 1`, each of the
alternating sine and cosine series has monotonically decreasing
positive term magnitudes. Following term `k`, the exact
remaining tail has the sign of term `k+1` and its absolute
value does not exceed the **first omitted term**. The upper
fixed-point bound of that term therefore widens only the
appropriate endpoint of each interval. Stop once both first
omitted upper bounds are at most `2^(-work)` (or at
the conservative maximum, still including the actual
enclosure of the remainder). For negative arguments, negate and
swap the sine endpoints; cosine is even. A zero argument
returns the exact pairs `(0,1)`.

The original argument-halving, double-angle interval propagation,
`tan = sin/cos` enclosure, all four rounding modes, and
endpoint-agreement refinement remain unchanged. The algorithm
is still native GMP, **without MPFR calls**. It handles
`|x| > 1` via the existing reduction, although large
argument reduction remains an important independent limitation.

**Validation:** a standalone exact-rational/fixed-point
implementation was checked on **219 rational arguments**, at
work precisions 8, 32, 53, 85, 128, 256, 512, 1024 bits,
against high-precision sine and cosine references; all
intervals contained their reference results. This is a useful
check of the interval recurrence, not a substitute for the
Ariadne regression suite. New C++ regression cases compare
`FloatAP` with `FloatMP` on dyadic inputs, including
negative and reduced arguments, in **all four rounding modes**
at 53, 128, 256 and 1024 bits.

**Performance baseline (user benchmark before this change):**

| Bits | FloatAP sin (ns) | MPFR sin (ns) | FloatAP cos (ns) | MPFR cos (ns) |
|---|---:|---:|---:|---:|
| 53 | 19531.67 | 704.44 | 19196.11 | 521.67 |
| 128 | 66202.36 | 1061.67 | 64755.83 | 691.11 |
| 256 | 252700.97 | 1541.53 | 251945.56 | 1150.14 |
| 1024 | 11626473.33 | 6159.58 | 11562120.00 | 4843.89 |

These are **pre-change measurements**, not demonstrated
speedups. The outer `tan` operation was 28.94 ms at
1024 bits versus 6.77 us for MPFR, motivating the
joint kernel. Rebenchmark identical inputs and verify
bit-exact results under every rounding mode before
claiming the optimization is accepted.

**State of the art and next steps:** FLINT/Arb provides
a pair of fixed-point Taylor kernels and a rectangular-splitting
variant, as well as higher-performance argument reduction:
https://flintlib.org/doc/arb.html .
Once these fixed-point series are validated and benchmarked,
the next priorities are certified constant caching (especially
`pi`), improved reduction for large arguments, fewer
rational interval temporaries during double-angle formulas,
and limb-level arithmetic fast paths.


## Performance redesign: cached certified pi (16 October 2026)

**Benchmark motivation:** after fixed-point sine/cosine, `pi_nearest`
became the largest *relative* discrepancy: repeated calls to
`FloatAP::pi` took approximately 5.9, 8.2, 12.7, and **84.3 us**
at 53, 128, 256, and 1024 bits, versus 49, 49, 50, and
51 ns for MPFR, respectively. The benchmark calls `pi`
300 times per sample after warm-up, so it measures a **hot
constant**, not one-off computation. MPFR caches its pi
constant; see https://www.mpfr.org/mpfr-current/mpfr.html
(section on memory handling) and `src/const_pi.c`.
FLINT's Arb `arb_const_pi` also caches constants:
https://flintlib.org/doc/arb.html .

**Decision:** implement two independent, bounded **thread-local**
caches, without delegating any numeric work to MPFR:

1. `FloatAPTrigonometricDetail::pi_interval(work)` retains the
   most recently constructed certified Machin interval at the
   **largest requested work precision** on each thread. It may
   safely return a tighter interval when a caller subsequently
   requests fewer bits. `atan_interval` and `acos` also use
   this shared interval, avoiding redundant evaluations.
2. `FloatAP::pi(mode, precision)` retains the already correctly
   rounded `FloatAP` value keyed by the **exact requested
   precision and rounding mode** (four cache slots for the
   most recently requested precision). Returning a copy
   preserves value semantics. A different precision invalidates
   the four slots; each rounding mode has its own entry.

**Safety and correctness:** no cached approximation is treated
as exact pi. The interval cache contains lower/upper rationals
that are certified by the existing Machin series and interval
arithmetic; the rounded-value cache is populated **only after**
the existing endpoint-agreement refinement certifies the requested
rounding mode. A higher-precision interval is a valid enclosure
for a lower-precision request. Thread-local state avoids races
and cross-thread mutations, with memory bounded by one interval
and four rounded values **per thread** (largest attained interval
precision may grow). Both caches are independent of the mutable
default rounding mode. Invalid rounding modes are rejected
before indexing the slots.

**Testing:** the shared float rounding suite now interleaves
53, 256, 7, 1024, 128, and 53 bits; it requests pi in all
four modes and compares the exact rational value with the
corresponding `FloatMP` result. A repeated same-key call
tests the hot cache. The existing `atan`/`asin`/`acos`
tests also exercise pi-dependent paths.

**Benchmark interpretation:** the next repeated-pi benchmark
should speed up markedly, but that result measures **cache
hits**. Measure first-call cost separately if relevant, as
well as mixed precision, concurrency, and allocation overhead.
The current benchmark timings are pre-cache baseline only.
All-mode numerical equivalence and CI must pass before
the cache is considered validated. Future work should speed
cold pi generation (e.g. Brent-Salamin/AGM) when profiling
demonstrates a need.

## Performance redesign: direct FloatAP-to-double conversion (16 October 2026)

**Measured motivation:** after the certified pi cache, `to_double`
was the worst relative benchmark: 136 ns at 53 bits and 324 ns
at 1024 bits, versus 18 ns for MPFR at both precisions
(approximately 7.5x and 17.8x slower). The earlier
`FloatAP::get_d()` constructed a full `Rational` and called
GMP `mpq_get_d`: for negative binary exponents this allocates
a power-of-two denominator and canonicalizes an exact rational
that is not otherwise needed.

**Decision:** convert the native positive GMP significand and
signed binary exponent directly to a hardware `double`. GMP
`mpz_get_d_2exp` extracts a normalized, **truncated**
double significand with at most 53 significant bits without
building a quotient or denominator. The original stored
significand remains available for exact guard, tie, parity,
and sticky-bit tests via `mpz_tstbit`/`mpz_scan1`. The
converter supports nearest-even, upward, downward, and
toward-zero using the **current FloatAP rounding mode**.
This corrects the earlier semantic discrepancy: the former
`Rational::get_d()`/GMP `mpq_get_d` always truncated
toward zero, while `FloatMP::get_d()` calls MPFR with
the current rounding mode.

**Boundary correctness:** determine the binary exponent from
the exact stored significant-bit count, guarding signed overflow.
For normal doubles, compute the 53-bit integer significand
and scale by an exact power of two. For subnormal doubles,
discard the additional bits below the `2^-1074` quantum;
ties are rounded to even, including the half-min-subnormal
case. A carry out of the significand adjusts the exponent or
causes overflow; directed modes return the finite maximum
where appropriate rather than indiscriminately returning
infinity. A nonzero underflow rounds to signed zero or the
minimum subnormal according to the direction. NaN,
signed infinity, and signed zero are handled without GMP
conversions. All finite results are exactly representable
scalings, avoiding a second rounding in `ldexp`.

**Independent verification:** the integer-extraction and
rounding recurrence was compiled with GMP and compared
against an exact GMP rational oracle on **54,880 cases**,
spanning 14 significand widths from 2 to 1024 bits,
14 exponent categories (including subnormals and normal
boundaries), positive/negative values, and all four modes.
Each compared the full 64-bit double representation. This
is a standalone algorithm check, not a complete C++/MPFR
regression. The Ariadne rounding suite now additionally
compares against actual `FloatMP::get_d()`, including
large finite values and cases that overflow `double`.
The benchmark speedup remains to be measured.

**Contract and follow-up:** the new implementation uses only GMP
for significand inspection; no MPFR arithmetic is introduced.
The mode-aware numerical behavior is an intentional
compatibility correction and must be validated against the
shared regression suite under all rounding modes. Benchmark
normal, subnormal, overflow, and exact-conversion workloads
separately; do not infer cold or worst-case latency from the
normal-input microbenchmark alone.

Reference: GNU GMP manual, integer conversion functions
(`mpz_get_d_2exp`), https://gmplib.org/manual/Converting-Integers;
GNU MPFR manual, `mpfr_get_d`, https://mpfr.org/mpfr-current/mpfr.html .

## Performance redesign: adaptive log series and log2 cache (16 October 2026)

**Evidence:** in the benchmark supplied after the direct `get_d()`
optimization, `FloatAP::log` at 1024 bits took 81.06 us versus
8.53 us for MPFR (9.51x), making it the worst relative case.
The benchmark's rational inputs are near `1/2`, so the old
reduction interval `[1/2,2]` left them near `1/2` and the
series argument `z=(x-1)/(x+1)` near `-1/3`. The previous
fixed-point implementation always evaluated `floor(work/3)+4`
terms, even for inputs whose transformed argument was tiny.

**Decisions:**

1. Normalize `x=2^k y` into **`3/4 <= y <= 3/2`**.
   Consequently `z=(y-1)/(y+1)` lies in `[-1/7,1/5]`.
   Inputs close to powers of two now have a very small `z`,
   reducing the work needed for the Taylor series.
2. At the end of each iteration, `power_hi` is an outward-rounded
   fixed-point upper bound for `|z|^(2N+1)`, where `N`
   terms have been included. Since `|z| <= 1/3` both for
   the reduced main input and for `log(2) = 2 atanh(1/3)`,
   the positive doubled tail obeys
   `tail <= (9/4)|z|^(2N+1) <= 3*power_hi/2^B`.
   Stop as soon as `3*power_hi <= 2^(B-work)`, or at
   the original conservative maximum count. Add the
   **actual certified tail bound** to the upper endpoint;
   no assumption of negligible error is required.
   The signed result uses interval negation and reversal.
3. Keep the highest-work **certified interval** for
   `log(2)` in a bounded, thread-local cache. It is valid
   for subsequent lower-work requests without changing its
   interval endpoints. Cache entries are independent of
   rounding mode; the existing `certified` loop still
   performs the final correctly rounded conversion.
   This is not a cache of approximate scalar values.

**Correctness invariants:** the fixed-point term and power
bounds are always outward-rounded. Early termination uses
an analytic bound on the complete omitted tail, including
the factor of two in `log((1+z)/(1-z))`. Neither the
normalization nor constant caching changes the exact
identity `log(x)=log(y)+k*log(2)`. The constant is never
treated as exact. The same monotone endpoint arithmetic
continues to support all four rounding modes.

**Regression added:** a differential oracle test compares
all four rounded modes of `FloatAP` against `FloatMP`
at 53, 128, 256, and 1024 bits, using exact binary inputs
near 1/2, 3/4, 1, 3/2 and 2, and small/large powers of two.
This is a test addition, not a claim that the CI has passed.

**Performance status:** the stated times are the *pre-change*
baseline. Future measurements must distinguish cold versus
warm `log(2)` and mixed-precision calls; expected hot
speedups must not be presented as measured results. This
targets both series complexity and redundant evaluation
of the same transcendental constant.

**State-of-the-art context:** FLINT/Arb uses table-based log
argument reduction and fixed-point Taylor evaluation via
rectangular splitting for low/medium precision:
https://flintlib.org/doc/arb.html . The present reduction
and cache are simpler, conservative improvements. More
aggressive log argument tables and rectangular splitting
remain possible future steps, after this regression.

## Arithmetic optimization: integer-denominator rounding (16 October 2026)

**Motivation:** the newest benchmark shows `FloatAP` still approximately
5x slower than MPFR on 53-bit addition, and 4.7x slower on 53-bit
multiplication, despite the fixed-point transcendental improvements.
Both arithmetic kernels assemble an **exact binary integer** and
invoke the general `from_ratio(numerator,denominator,scale,precision,mode)`
routine with `denominator = 1`. Previously, that routine allocated
five additional temporary GMP integers, normalized the rational and
performed a general quotient/remainder division.

**Decision and invariant:** detect the **positive unit denominator**
in `FloatAP::from_ratio`, and round the signed integer significand
directly. Copy its magnitude into the output `mpz_t`; use the exact
integer bit length to calculate the normalized binary exponent and
required shift. If the integer has fewer than the requested precision
bits, left-shift it exactly. Otherwise, inspect the discarded bit
pattern with `mpz_tstbit` and `mpz_scan1`, truncate directly with
`mpz_tdiv_q_2exp`, and increment for the requested mode when
necessary. Nearest uses **ties-to-even**: increment for a guard bit
and either nonzero lower discarded bits or an odd retained LSB.
Directed modes increment only for a nonzero discarded portion and
the appropriate result sign. If increment creates a carry beyond
the target precision, shift the significand right one bit and
increment the stored exponent.

The exact stored value remains
`(-1)^negative * significand * 2^exponent`. Preserve signed-zero
handling, overflow-checked exponent arithmetic, the original
general-rational division path (for genuine nonunit denominators),
all four rounding modes and all special-value handling. This is a
**single fast path**, not a change to the public `FloatAP` API.

**Verification before publication:** a standalone native GMP C++
harness compared the new integer-specialized path with the old
general denominator-one routine in **318,400 directed-rounding cases**,
covering precision from 2 to 2048 bits, significand lengths below,
at and above the target, signed operands, rounding carries, and all
four modes. Both normalized significand and exponent matched.
The Ariadne rounding suite additionally receives exact FloatAP vs
FloatMP comparisons at 7/53/128/256/1024 bits, checking integer
construction, add/sub/mul/fma under all four rounding modes. Neither
the harness nor the additional tests imply that full CI has passed.

**Performance status:** the reference benchmark provided on
16 October 2026 measured 53-bit `add_nearest` at 346.39 ns for
FloatAP versus 74.45 ns for MPFR, and 53-bit `mul_nearest` at
294.03 ns versus 67.22 ns. At 1024 bits `add_nearest` took
392.36 ns versus 70.55 ns, and `mul_nearest` 472.50 ns versus
177.22 ns. These measurements are **before** this change.
A same-host post-change benchmark is needed to quantify the
benefit and detect any new bottlenecks.

**Next options:** avoid constructing a temporary unit denominator
even at the call sites, and specialize ordinary direct limb addition
and multiplication where profiling justifies the complexity. This
fast path does not change rational-based transcendental enclosures
or borrow rounding logic from MPFR.

## Interval quotient optimization: sign-aware extrema (16 October 2026)

The benchmark after integer-denominator rounding showed 1024-bit
`tan_nearest` at 61.61 us versus 6.72 us for MPFR (9.17x).
The original certified trigonometric interval division built two
rational reciprocals and four rational cross-products, then selected
their minimum and maximum, even though endpoint signs determine the
extrema in advance.

**Decision:** preserve exact rational enclosures, but choose the
extremal endpoint pair using the signs of both intervals. For a
positive denominator `[c,d]` and numerator `[a,b]`, return
`[a/d,b/c]` if `a>=0`, `[a/c,b/d]` if `b<=0`,
and `[a/c,b/c]` if the numerator straddles zero. For a
negative denominator, return `[b/d,a/c]`, `[b/c,a/d]`,
or `[b/d,a/d]` for the same three numerator cases.
The computation now requires only **two exact rational divisions**,
rather than two reciprocals and four cross-products.

The old domain check is unchanged: a denominator touching or
straddling zero is rejected. This matters for tangent poles.
Because the two selected quotients are exactly the minimum and
maximum of all four endpoint combinations, no rounding or
enclosure quality is sacrificed. The improvement also applies
to other certified trigonometric paths using interval division,
including arcsine.

**Validation:** an independent exact-arithmetic enumeration checked
240,426 valid interval configurations over rational endpoints,
including negative denominators and numerators crossing zero;
the optimized pair matched the old four-candidate extrema in
every case. The shared rounding tests now additionally compare
`FloatAP::tan` with `FloatMP::tan` at 53, 128, 256, and
1024 bits in all four rounding modes. Negative arguments and
arguments with negative cosine are explicitly included.

**Performance status:** the benchmark values above precede this
change; a new run is needed to measure actual speedup and detect
any trade-offs. FLINT/Arb's interval arithmetic design remains
a useful reference for rigorous error propagation:
https://flintlib.org/doc/arb.html . Further performance work
should profile the exact-rational endpoint construction and
conversion overhead before undertaking more intrusive changes.

## Selective fixed-point trigonometric evaluation (16 October 2026)

**Observed bottleneck:** the latest benchmark reports 1024-bit
`cos_nearest` at 41.82 us versus 4.84 us for MPFR (8.63x),
and `sin_nearest` at 41.92 us versus 5.90 us (7.11x).
`tan_nearest` is 49.41 us versus 6.53 us (7.56x).
The fixed-point Taylor kernel previously evaluated both the sine
and cosine recurrence on every call, even when the public request
only required one result. Each unused recurrence performed two
GMP multiplications, shifts, directed divisions and accumulator
updates per iteration.

**Decision:** add a private `TrigOutputs` selector to the
fixed-point kernel and propagate it from the trigonometric
evaluation dispatcher. When the reduced argument is already
in `[-1,1]`, a request for `sin` skips the entire cosine
Taylor recurrence, and vice versa. This avoids building an
unused series and converting its interval endpoints to
rationals. For `tan`, always evaluate both functions.
If the original argument requires one or more halvings,
**both** series remain necessary: the double-angle identities
`sin(2x) = 2 sin(x) cos(x)` and
`cos(2x) = cos(x)^2 - sin(x)^2` must be propagated
using certified intervals. No change is made to the
existing argument-reduction threshold or identities.

**Certification:** the active recurrence retains its original
outward-rounded GMP fixed-point endpoints and the alternating
Taylor remainder bound. The stopping condition only consults
the bounds of the requested series (or both when necessary).
Thus a single-output evaluation may finish sooner without
changing the mathematical enclosure guarantee.
The unrequested element of the private return pair is a
placeholder and **must not** be used by callers; the public
dispatcher accesses only the requested result. The
double-angle and tangent paths always request both intervals.
No rounding mode or public interface is changed.

**Regression:** extend the differential sine/cosine/tangent
suite comparing FloatAP with MPFR at 53, 128, 256, and
1024 bits under all four rounding modes. New arguments
include exactly +1/-1 and +7/8/-7/8, in addition to values
above 1 and negative-cosine cases. These tests exercise both
the selective and joint kernel paths; correctness of the
rounded results is the acceptance criterion.

**Measurements:** the reported 41.82 us cosine and 41.92 us
sine at 1024 bits are the **pre-change** reference. The
latest benchmark was collected under visibly different
timing conditions from the preceding run (for example,
the MPFR comparison microbenchmark changed substantially).
Do not attribute inter-run changes exclusively to
algorithmic modifications without same-environment controls.
Rebenchmark identical inputs, including `tan`, and
measure cold and warm cases before recording any speedup.

**Relation to the state of the art:** FLINT/Arb explicitly
supports evaluating sine alone or both sine and cosine in
its rectangular-splitting Taylor kernel (the `sinonly`
option). See https://flintlib.org/doc/arb.html .
The present change removes clearly unnecessary computations;
it does not yet implement rectangular splitting, advanced
argument reduction, or limb-level intermediate reuse.


## Exponential reduction with bounded fixed-point reconstruction (16 October 2026)

**Motivation:** after the selective sine/cosine implementation, the
user's same-run benchmark reported 1024-bit `exp_nearest` at
54.97 us for FloatAP and 6.42 us for MPFR (8.57x), the largest
remaining relative slowdown. At 53 bits, it took 4.66 us versus
1.66 us. The existing `exp_series` already uses rigorous GMP
fixed-point Taylor intervals but reduces its argument only to
`|x| <= 1/2`. More aggressive halving reduces the number of
series terms, but the old *exact-rational* interval squaring
doubled the denominator bit length after each squaring. Thus
simply reducing to `1/2^16` would make reconstruction slower,
not faster.

**Decision:** select a target `|x| <= 2^(-d)` with
`d = min(16, max(6, floor(work/32)))`. Compute
`exp(x/2^k)` using the existing outward-rounded fixed-point
series with `work + k + 8` working bits. Convert the resulting
strictly positive rational enclosure once to **B-bit scaled
integers**, using floor for the lower endpoint and ceiling for
the upper, with `B = work + k + 64`. Recover `exp(x)` by
`k` interval squarings in the integer domain:

```text
L' = floor(L*L / 2^B)
U' = ceil (U*U / 2^B)
```

As squaring is monotone for nonnegative operands, the new
`[L', U']` encloses the exact squared interval. There are no
growing rational denominators in the reconstruction loop;
only the two final endpoints become `Fraction` objects.
The reduction is a general algorithmic choice based on
requested precision, **not** a branch for benchmark inputs.

**Rounding guarantee:** the reduced exponential series already
provides an outward certified enclosure. Converting to integers
rounds outward; every reconstruction step also rounds outward.
These are closed supersets of the mathematical result regardless
of whether `B` proves sufficient for a particular hard-to-round
case. The outer `certified` loop continues increasing precision
until both endpoints give the same rounded FloatAP value for the
requested rounding mode. Existing special cases, positive and
negative inputs, and all four rounding modes retain their
contracts. Overflow checks guard the working bit count and
GMP shift conversion.

**Prototype verification:** a standalone GMP C++ program compared
bounded fixed-point reconstruction with exact rational squaring
of **the same** certified reduced-series enclosure. On **1,040
cases** across 53, 128, 256, and 1024 bits, reduction limits
`1/2`, `1/4`, `1/8`, `1/16`, random rational inputs and
both signs, the new lower/upper endpoints enclosed the exact
reference endpoints in every case. This tests reconstruction
independently, not the full Ariadne build. Differential C++
tests have also been added for `exp` on positive/negative
dyadics, in all four modes at those four precisions.

**Prototype performance, not an Ariadne benchmark:** on a Linux
GMP C++ harness at 1024 working bits and an argument near 1/2,
the original half-reduced series and exact-rational squaring took
approximately 78 us; bounded reconstruction with a `1/2^16`
reduction took approximately 35 us. These figures are
environment-specific, useful to choose the algorithm, and are
**not** asserted to be the speedup of the committed code.
At 53 bits, aggressive reduction has a weaker benefit, hence
the precision-dependent reduction schedule.

**Further improvements:** optimize `exp` with rectangular
splitting or bit-burst at sufficiently high precision, and
profile fixed-point series iterations separately from argument
reduction and endpoint conversion. FLINT/Arb implements
`_arb_exp_taylor_naive`, `_arb_exp_taylor_rs`, and more
aggressive argument reduction:
https://flintlib.org/doc/arb.html .

## Trigonometric performance: derive cosine from a certified sine interval (16 October 2026)

**Benchmark motivation:** the latest same-run 1024-bit median for
`tan_nearest` is 49.82 us versus 6.68 us for MPFR (7.46x).
Although the preceding selective-Taylor change roughly halved
the time of separate `sin` and `cos`, tangent still requested
**both** complete fixed-point Taylor recurrences.

**State of the art:** FLINT/Arb's `sin_cos.c` deliberately computes
only sine and derives cosine by an integer square root when the
Taylor term count is sufficiently large; see the code path described
as "Compute cos(a) from sin(a) using a square root" in
https://fossies.org/linux/flint/src/arb/sin_cos.c and the
`_arb_sin_cos_taylor_rs` `sinonly` option documented at
https://flintlib.org/doc/arb.html . The technique is especially
attractive for rigorous fixed-point arithmetic.

**Decision:** for the **private** `TrigOutputs::both` mode,
evaluate only the certified sine recurrence; obtain a certified
cosine enclosure using `cos(a)=sqrt(1-sin(a)^2)`.
This is valid for the already reduced argument `|a|<=1`,
where sine is nonnegative and cosine is strictly positive.
Individual sine-only and cosine-only calls keep their current
selective Taylor implementations. All reduction/double-angle
and tangent callers retain their original interfaces and
continue requesting both results where needed.

**Exactly directed integer-square-root enclosure:** let the
fixed-point scale be `Q=2^B`, and let the original outward
sine enclosure be `[Slo/Q, Shi/Q]` for `sin(|a|)`.
Intersect the bounds with the proven `[0,1]` range,
giving integer endpoints `0<=L<=U<=Q`. Since
`sqrt(1-s^2)` decreases on `[0,1]`, the exact cosine
obeys

```text
floor(sqrt(Q^2 - U^2))/Q <= cos(a) <=
ceil(sqrt(Q^2 - L^2))/Q.
```

Compute both radicands **exactly in GMP integers**.
Use `mpz_sqrt` for the lower root and `mpz_sqrtrem`
for the upper root, incrementing the upper result only
when the integer remainder is nonzero. Thus the new
interval uses only exact integer operations, introduces
no additional unquantified error, and cannot silently
truncate cosine. Negative arguments negate only sine;
cosine remains even. The existing correctly-rounded
`certified` refinement is unchanged in every rounding
mode. The same procedure also applies before existing
double-angle reconstruction for larger arguments.

**Independent verification:** 7,344 cosine enclosures were
checked at eight fixed-point precisions (7, 24, 53, 113,
128, 256, 512, 1024 bits) against high-precision
`mpmath.cos` values, including input zero, one, 7/8,
1023/1024, tiny arguments, 300 random arguments in
`[0,1]` per precision and three choices of endpoint
slack. All intervals contained the reference value.
This validates the cosine-reconstruction inequality,
not the full C++ numerical kernel. The shared Ariadne
differential tests now include positive and negative
inputs near 1, very small arguments, arguments requiring
double-angle reconstruction, all four rounding modes,
and 53, 128, 256, 1024-bit target precisions.

**Acceptance:** rerun the complete C++ and MPFR differential
suites and measure the benchmark on the same hardware.
The expected improvement concerns `tan` and trigonometric
calculations requiring both outputs. The separate
`cos_nearest` kernel is not changed by this optimization,
and no improvement is claimed until measured. The
mathematical guarantee does not depend on performance.


## Inverse trigonometric optimization: single arctangent enclosure (16 October 2026)

**Benchmark motivation:** after the certified cosine-from-sine
optimization, the largest remaining relative slowdown was
`acos_nearest` at 1024 bits: **156.70 us** for FloatAP
versus **23.42 us** for MPFR (6.69x). At the same
precision, `asin_nearest` took 155.04 us versus 23.40 us.
The common `asin_interval` routine obtained an interval
for `z=x/(1+sqrt(1-x^2))` and computed `atan(z.lower)`
**and** `atan(z.upper)` separately with the certified
fixed-point arctangent series. The endpoints differ by only
the previous outward square-root uncertainty, so the second
arctangent is largely redundant.

**Decision:** replace the two independent arctangent calls with
**one** certified arctangent interval and an exact rational
width correction. For the certified transformed interval
`a <= z <= b`, with `w=b-a >= 0`, the arctangent is
monotone and globally **1-Lipschitz** because
`0 < atan'(t) = 1/(1+t^2) <= 1`. Therefore, given
`atan(a) in [L,U]`, all exact values of `atan(z)`
lie in `[L,U+w]`. Alternatively, from `atan(b) in [L,U]`,
all values lie in `[L-w,U]`. Select the endpoint with the
smaller absolute value, which is cheaper for the alternating
arctangent Taylor series, then double both interval endpoints
to obtain the rigorous enclosure for `asin(x)`.

**Mathematical and numerical invariants:** no approximation is
substituted for an exact input endpoint. All widths and bound
enlargements use GMP rational arithmetic with no rounding.
The derivative bound holds for **all real inputs**, so it is
also valid for negative `x` and near both ends of the
arcsine domain. The prior domain checks, signed behavior
and correctly-rounded final certification remain unchanged.
`acos(x)=pi/2-asin(x)` automatically benefits from the same
change; its `pi` interval remains certified and cached.

**Verification:** an independent 425-decimal-digit mpmath
check exercised **7,464** rational intervals across widths
corresponding to 7, 24, 53, 128, 256, and 1024 bits, with
both endpoint choices and positive/negative inputs. The
new intervals contained arctangent at each endpoint and at
the midpoint in all cases. This is numerical corroboration;
the rigorous guarantee follows from the analytic derivative
bound, not from finite testing. The shared C++ rounding suite
also now compares FloatAP `asin` and `acos` with MPFR at
53, 128, 256, and 1024 bits, using inputs near +/-1, zero,
both signs, and all four rounding modes.

**Expected performance:** this eliminates one arctangent
series evaluation from each ordinary `asin`/`acos` call.
The width correction introduces a small number of exact
rational operations and may very occasionally require more
working precision in the existing certified endpoint-agreement
loop; correctness is independent of the extra iterations.
The measurements above are **before** the optimization.
A same-environment post-change benchmark and a full CI pass
are required before any speedup is claimed. More extensive
changes to `atan` argument reduction remain a separate,
higher-risk opportunity.

## Exponential optimization: single certified fixed-point Taylor recurrence (16 October 2026)

**Bottleneck:** the 16 October benchmark after the one-atan optimization
showed FloatAP `exp_nearest` taking 5.37 us at 128 bits versus 1.04 us
for MPFR (5.15x), and 28.34 us at 1024 bits versus 6.42 us
(4.42x). The previous `exp_series` evaluated independent lower
and upper Taylor-term recurrences with two GMP big-integer products,
directed shifts, and divisions per term, as well as another upper
multiplication/division to bound the first omitted term.

**Decision:** compute only the downward-rounded positive Taylor
recurrence at fixed-point scale `Q=2^B`, with
`B=work+32+bit_width(floor(work/2)+64)`.
For `0 <= x <= 1/2`, define
`a=floor(Q*x)`, `t_0=Q`, and
`t_n=floor(t_(n-1)*a/(Q*n))`.
Every `t_n` is a lower bound for the exact scaled Taylor term
`T_n=Q*x^n/n!`; the two successive GMP floor operations are
equivalent to the single floor in this definition.

**Rigorous bound for *all* intermediate rounding errors:**
Write `e_n=T_n-t_n`. Because `0 <= Q*x-a < 1` and
`t_(n-1) <= Q*exp(x) < 2Q`,
the recurrence satisfies

```text
0 <= e_n < e_(n-1)*x/n + t_(n-1)/(Q*n) + 1
          < e_(n-1)/(2*n) + 2/n + 1.
```

Since `e_0=0`, induction gives `e_n<6` for **every**
`n>=1`, independently of target precision or Taylor length.
After including `n` terms, their aggregate downward error
is less than `6*n`. Furthermore, for `x<=1/2`,
the complete omitted tail satisfies
`sum_{k>n} T_k <= 2*T_(n+1) <= T_n <= t_n+6`.
Consequently, with `S=sum_{k=0}^n t_k`, the entire exact
exponential satisfies

```text
S / Q <= exp(x) <= (S + t_n + 6*(n+1)) / Q.
```

The upper bound is constructed in exact GMP integers; computing
`6*(n+1)` with GMP also avoids native unsigned overflow.
Termination uses `t_n <= 2^(B-work)` or the existing
conservative maximum-terms limit. All certified endpoint
rounding and extra-precision refinement stay unchanged.
Negative inputs use `exp(-x)=1/exp(x)` with swapped rational
interval endpoints. The subsequent bounded fixed-point
squaring reconstruction, introduced earlier, is unaffected.

**Independent verification:** a standalone compiled GMP C++ harness
tested 360 positive and negative rational arguments with precisions
53, 128, 256, and 1024 bits, comparing the resulting rational
intervals against lower/upper exact rational Taylor bounds. All
360 certificates enclosed the reference. The proof above provides
the general guarantee, not the finite test count. Additional
Ariadne differential tests compare FloatAP and FloatMP `exp`
for small arguments, plus inputs close to unity, at each of the
four target precisions and all four rounding modes.

**Prototype benchmark (not the Ariadne benchmark):** on the local
Linux GMP harness with a benchmark-like input near 1/2 and 32
extra working bits, the single-recurrence exponential completed
in approximately 2.57 us at 128 bits and 17.1 us at 1024 bits.
The earlier two-recurrence benchmark measured approximately
4.1 us at 128 bits and over 30 us at 1024 bits. Cross-run
differences are not controlled; evaluate the final outcome
using the same-host FloatAP/MPFR benchmark and full CI.

**Scope:** this is an algebraic interval bound, not a heuristic
or an MPFR-dependent shortcut. The additional interval width is
deliberately conservative; the existing certified loop retries
at higher work precision if both endpoints do not round to
the same target value. No special-value semantics, public API,
or functional test expectations are relaxed.

## Trigonometric optimization: one certified fixed-point recurrence (16 October 2026)

**Measured motivation (pre-change):** the 16 October benchmark after
the exponential single-recurrence optimization reported `sin_nearest`
at 1024 bits as 22.49 us (MPFR 5.90 us), `cos_nearest` as
22.16 us (MPFR 4.79 us), and `tan_nearest` as 30.70 us
(MPFR 6.47 us). FloatAP's sin/cos Taylor kernel still calculated
both lower and upper outward-rounded sequences of positive term
magnitudes, performing two big-integer products and two divisions
for every Taylor term. The existing tangent path already avoids
evaluating cosine's Taylor series, reconstructing positive cosine
from a certified sine interval using integer square roots.

**Decision:** compute only a **downward-rounded fixed-point
recurrence** for the positive magnitudes of the alternating
sine or cosine Taylor terms at the reduced argument `|x|<=1`.
Let `Q=2^B`, `a=floor(Q*|x|)`,
`b=floor(a*a/Q)`, and

```text
t_0 = a                 for sine
t_0 = Q                 for cosine
t_(n+1) = floor(t_n*b/(Q*d_n)),
d_n = (2n+2)(2n+3)      for sine,
d_n = (2n+1)(2n+2)      for cosine.
```

The two sequential GMP divisions by `Q` and `d_n` implement
the same floor as the combined rational quotient because all
intermediate quantities are nonnegative. The remaining
alternating sum uses the standard `+,-,+,-` signs of sine
and cosine. For a tangent request, compute only the sine
recurrence; reconstruct cosine with the existing exact integer
square-root bounds. Negative inputs negate only the resulting
sine interval. The existing halving/double-angle mechanism is
unchanged.

**Certification proof:** let `T_n` be the exact positive Taylor
term scaled by `Q`, and `e_n=T_n-t_n`. The initial error
is less than one for sine and zero for cosine. Since
`|x|<=1` and `a=floor(Q*|x|)`, we have

```text
0 <= Q*x^2 - b < 3,
0 <= t_n <= T_n <= Q,
0 <= e_(n+1) < e_n/d_n + 3/d_n + 1.
```

Every `d_n>=2`, so induction yields `0<=e_n<8`
for all `n`. After accumulating `n+1` terms with
alternating signs, the difference between the exact partial
sum and the integer sum `S_n` is bounded in magnitude by
`8*(n+1)` fixed-point units. The alternating-series
remainder has the sign of the **first omitted term** and
magnitude at most `T_(n+1)<t_(n+1)+8`. Consequently,
with `E=8*(n+1)` and `R=t_(n+1)+8`,

```text
n even: [S_n - E - R, S_n + E] / Q
n odd:  [S_n - E,     S_n + E + R] / Q
```

contains the exact sine or cosine of the positive reduced
argument, regardless of when the adaptive series terminates.
All error budgets use GMP integer arithmetic and avoid
native integer overflow. For a tangent request, intersect
the positive-sine enclosure with the mathematically known
`[0,1]` range and compute the cosine bounds as

```text
floor(sqrt(Q^2-U^2))/Q <= cos(x) <=
ceil(sqrt(Q^2-L^2))/Q.
```

This retains correct rounding in all modes via the existing
endpoint-agreement refinement. Certification does not rely
on testing or on a heuristic estimate of discarded bits.

**Independent numerical checks:** a Python exact-integer and
`mpmath` harness tested 160 nonzero arguments at each of
8, 53, 128, 256, and 1024 working bits, for a total of
**800 argument/precision combinations**. Both sine and
cosine enclosures contained the high-precision reference
values in all cases. The cosine-from-sine integer-square-root
reconstruction also enclosed the reference for all cases.
Arguments included +/-1, +/-1023/1024, +/-1/2, random
rationals across [-1,1], and values as tiny as +/-2^-1050.
These finite tests corroborate the mathematical argument;
they are not a replacement for the proof.

**Ariadne regression:** the shared differential suite
compares sine, cosine, and tangent from FloatAP and MPFR
at 53, 128, 256, and 1024 bits in all four rounding modes.
New arguments include +/-31/32, +/-17/16 and +/-2^-128.
The test also exercises the separately rounded cosine
series and the sine-derived cosine path. No previous
expectations were relaxed.

**Performance:** all timings above are a **pre-change**
baseline. The new kernel removes half of each active
fixed-point Taylor term recurrence and associated GMP
division work; however, actual impact must be measured
in the same-host benchmark after CI passes. If the
remaining bottleneck is factorial division or GMP
allocation overhead, additional improvements should
be driven by profiling rather than assumptions.

## Shared GMP kernel for addition and subtraction (16 October 2026)

**Motivation:** the preceding same-run benchmark measured FloatAP
`sub_nearest` at 170.14 ns versus MPFR 40.00 ns at 256 bits.
The old implementation copied and negated the second operand
before calling add; add then aligned two independent signed GMP
temporaries. This caused unnecessary allocations, copying and
sign changes.

**Design decision:** represent subtraction as a logical inversion
of the sign bit of the second operand, sharing the `add_sub`
kernel. Select the operand with larger binary exponent and
shift only its positive significand down to the smaller exponent.
Then perform a single GMP `mpz_add` or `mpz_sub` on unsigned
magnitudes and apply the effective sign of the larger-exponent
operand. The resulting integer numerator and binary scale are
algebraically identical to the former algorithm. The existing
exact normalization, nearest-even and directed rounding in
`from_ratio` are unchanged, including signed cancellation zero.
The exponent-gap fast path uses effective (possibly inverted)
operand signs. NaN, infinity and zero handling remains explicit.

**Standalone validation:** 200,000 pseudorandom signed GMP integer
alignments with variable exponent differences and significand
sizes up to 1024 bits exactly matched the previous kernel.
On the same Linux executable, the isolated GMP alignment
microbenchmark decreased from 41-45 ns to 11-12 ns per call.
This is **not** an Ariadne benchmark and excludes final rounding.

**Regression:** expand FloatAP-versus-MPFR exact-result testing
for add and sub across 7, 53, 128, 256, and 1024 bits,
all four rounding modes, both signs, cancellation, exponent
gaps and signed zero. CI and the same-environment benchmark
must verify the integrated implementation before attributing
an end-to-end speedup.

## Logging integration and benchmark reporting (16 October 2026)

**Dependency:** add `ariadne-cps/logging` under `ss/logging` as a
Git submodule. The pinned logging revision is
`3df461f6aabe6e87eb7cabfbabc5e93e0f5be87c`.
Reuse the existing `utility` dependency, with identical
upstream revisions of `configuration` and `utility` across
the numeric/foundation/logging submodule graph. The CMake
benchmark targets link to `ariadne-logging` and `logging`;
the numeric library's public API does **not** acquire a
logging dependency.

**Output conventions:** human-readable benchmark output uses
`LOGGING_PRINTLN` with aligned columns, wide separators,
and whitespace between precision groups and benchmark
phases, in the style of the integration project's benchmarks.
For the FloatAP/MPFR comparison, table output is the default.
Machine-readable CSV is retained behind `--csv` or
`--format=csv`, with the same column names as before and
without logger prefixes or mixed summary data.
The median-of-samples measurement procedure and existing
FloatAP/MPFR comparison arithmetic are unchanged.

**Command-line configuration:** comparison benchmark accepts
the old positional `iterations odd_samples` arguments as
well as `--iterations N`, `--samples N`,
`--bits 53,128,256,1024`, `--summary-only`,
`--details`, and `--max-slowdown RATIO`.
The arithmetic benchmark accepts the existing positional
exponents and runs plus `--vector-exponent N`,
`--repetition-exponent N`, and `--runs N`, including
`--help`. Defaults remain backward compatible.
The CSV option is intended for scripted performance
regressions and spreadsheet comparison.

**Independent functional failure exposed by the supplied CI
logs:** the previous large-exponent-gap shortcut in the new
shared `add_sub` kernel produced wrong directed subtraction
results versus MPFR (in `test_precision_regressions`).
Remove the shortcut and always compute the exact signed
integer before performing the **single final rounding**.
The differential tests remain intact, so this fixes the
production implementation, not the test expectations.
A conservative exact path is preferable to a faster but
numerically incorrect shortcut; later optimizations must
prove directed-rounding behavior over large exponent gaps.

**Verification requirements:** the uploaded Release logs
compiled successfully on Linux GCC, Linux Clang and AppleClang
but failed the FloatAP/MPFR subtraction differential checks.
The ThreadSanitizer build completed, but a full TSAN test
outcome was not established. Neither post-fix correctness nor
logging-integrated CMake builds are claimed until CI completes.

## Benchmark reporting via the Ariadne logger (16 October 2026)

Following `ariadne-cps/integration`, both numeric benchmarks explicitly
set `Logger::instance().configuration().set_verbosity(1u)` on entry:
logging defaults to verbosity zero, which suppresses every
`LOGGING_PRINTLN` call. All comparison measurements, headings,
separators, progress notifications, slowdown summaries, help and error
messages now pass through the logger. **There is no CSV mode**: the
previous `--csv` and `--format=csv` options and the direct
`std::cout` output have been removed. This requirement supersedes the
earlier CSV-output design recorded above. Progress is reported before
and after each precision group so long benchmark runs are not silent;
the formatted result table and aggregate slowdown are emitted after
all requested groups complete. The timed operations and sample
aggregation are unchanged.

## Benchmark command-line parsing via Utility (16 October 2026)

Numeric benchmarks use the shared Ariadne `utility/command_line_interface.hpp`
API, matching `ariadne-cps/integration`: `CommandLineInterface`,
`CommandLineContribution`, `make_valued_argument_parser` and
`make_unvalued_argument_parser`. The former ad hoc per-benchmark argument
scanning was removed. Numeric's existing foundation -> utility dependency
provides both the interface and its compiled implementation; no new
dependency or private copy of the command-line parser is introduced.

The common `benchmarks/benchmark_command_line.hpp` adapter translates
the previous positional invocations (comparison: iterations, samples;
arithmetic: vector exponent, repetition exponent, runs) into named arguments
before calling the same parser. The adapter does not parse values and cannot
bypass normal CLI validation or duplicate-option detection. `--help` and
`-h` use the built-in CLI help, as in integration, including its built-in
diagnostics. Explicit options are `--iterations/-i`, `--samples/-s`,
`--bits/-b`, `--summary-only`, `--details`, `--max-slowdown`,
`--vector-exponent/-v`, `--repetition-exponent/-e`, and `--runs/-r`
(as relevant to each executable). Invalid values are rejected before any
measurements run. Measurement reports continue to use the Ariadne logger,
not CSV or standard output.

## Built-in logging CLI options (16 October 2026)

Both benchmarks now add the upstream `logging_command_line()` contribution
from `logging/command_line.hpp` to the common Ariadne
`CommandLineInterface`, as the kernel CLI does. Consequently `--scheduler/-s`
(immediate, blocking, nonblocking), `--theme/-t` (none, light, dark), and
`--verbosity/-v` are available and appear in standard `--help` output.
The benchmark-specific short names have been made non-conflicting:
comparison samples use `-n`, arithmetic vector exponent uses `-x`.
The long option names are unchanged. The argument normalization adapter
preserves logging option values, including `--scheduler immediate`.

The logging scheduler requires an attached thread registry when changed;
standalone numeric benchmarks attach a long-lived registry that reports
no additional application threads, then set the default report verbosity
of 1 before CLI acquisition. The logging contribution can override this
with `--verbosity 0` to mute reports. The numeric library itself has no
logging dependency; this concerns the benchmark executables only.

## Signed-exponent-gap arithmetic at the native range limits (16 October 2026)

**Observed regression:** the user-supplied CI logs for PR #28 at
`3ffefc11` report `test_float_rounding` failing in both
`TestFloatRounding<ArbitraryPrecision>` instantiations (64 and 192 bits)
on Linux GCC 13, Clang 16, AppleClang 21, and Linux ThreadSanitizer.
The existing assertion `add(up,max,min)==inf_` throws
`FloatAP exponent overflow` instead. `FloatAP::max` and
`FloatAP::min` use `INT64_MAX` and `INT64_MIN` respectively;
subtracting their exponents in the exact-addition alignment exceeds
the signed exponent type before numeric rounding occurs.

**Design:** keep exact GMP alignment for ordinary exponent gaps.
When the unsigned, overflow-safe difference of two exponents exceeds
the sum of both significand widths, the destination precision, and two
bits, the smaller operand
cannot cross any noncoincident rounding boundary at the public result
precision (which may be the destination precision for in-place addition).
Its sign alone
selects the appropriate side of an exactly representable value or a
tie. Instead of making a GMP integer with billions of bits, substitute
a signed one-bit perturbation `2^(-high_precision-result_precision-2)` relative to the
larger significand and perform **one** correctly directed rounding via
`from_ratio` with a power-of-two denominator. Unlike the former
large-gap shortcut, this preserves both directed results and nearest-even
tie breaking for opposite-sign operands. A rounded significand carry
beyond `INT64_MAX` now yields the correctly signed infinity; other
checked-exponent arithmetic still throws on nonrepresentable exponents.

**Regression:** extend the historical `test_limits()` for FloatAP at
both 64 and 192 bits, covering positive/negative maximal finite values,
both signs of the minimum positive finite operand, overflow, adjacent
finite results, higher-precision in-place addition, and
round-to-nearest/downward/upward behavior. The
existing shared assertion remains unchanged. Independent exact
rational models checked 116,000 comparisons between direct and substituted
rounding, including halfway and endpoint significands and cases where the
requested result precision exceeds the dominant operand precision. This is not a
full Ariadne build, and no new CI pass or coverage percentage is claimed
until the published commit is tested on the requested branch.

## Native-precision FloatAP rescaling at extreme exponents (16 October 2026)

**CI follow-up:** the attached `logs_105616396908.zip` for PR #28
(checkout merge `11c9323` incorporating `f1c78c18`) shows that the
original `add(up,max,min)==inf_` regression and all other new
far-exponent arithmetic checks now pass. Nevertheless,
`test_float_rounding` fails on AppleClang 21, GCC 13, Clang 16 and
ThreadSanitizer. The TSAN stack identifies GMP's
`__gmp_overflow_in_mpz -> __gmpz_mul_2exp` via
`FloatAP::to_rational()`, `operator Rational()`, and
`FloatAP(FloatAP const&, RoundingModeType, ArbitraryPrecision)`.
The triggering expression is `FloatAP wide_high(max,wider)` in the
new historical limits test. This is a production bug, not an invalid
test expectation: exponent `INT64_MAX` must not be expanded into a
`Rational` merely to change FloatAP precision.

**Design fix:** when the precision is unchanged, copy the representation;
otherwise handle nonfinite/signed-zero values directly and rescale a
finite significand using `from_ratio(signed_significand, 1,
original_exponent, target_precision, rounding_mode)`. GMP now works
only with the existing bounded significand and its destination-precision
rounding bits rather than allocating `2^INT64_MAX` limbs. Using a
signed significand retains directed-rounding behavior for negative
operands. Keep the existing exponent representability checks.

**Regression:** extend `test_limits()` to verify exact widening and
narrowing at `INT64_MAX`, preservation of sign, and the already present
in-place wider-precision addition. Tests compare FloatAP values directly;
they do not convert extreme values to Rational for their assertions.
The incoming logs are not a successful test report. A fresh build,
CI execution and coverage report for the published correction are needed.

## Coverage 91: Rational-left FloatAP comparisons (16 October 2026)

**Incoming baseline:** The user-supplied `Archivio 91.zip` is the
source-based LLVM coverage report for the `remove-mpfr` work.
It reports for `src/floatap.cpp`: 181/199 functions, 1476/1753 lines,
695/978 branches; for `include/numeric/floatap.hpp`: 84/102 functions,
197/252 lines and 56/110 branches. The related Coverage workflow #358
at `d2a0e0d0` is green, but does not assert complete coverage.
Regions are secondary and not the acceptance criterion.

**Observed, per-source gap:** In `include/numeric/floatap.hpp`, all
six `Rational`-left comparison overloads (`==`, `!=`, `<`, `>`,
`<=`, `>=`, lines 417-436 of this baseline) have zero execution
counts. Both outcomes of both branch records at lines 431 and 435
also have zero counts. By contrast, `FloatAP`-left comparisons already
run. In `src/floatap.cpp`, the unordered return path of both
`cmp(FloatAP,Rational)` and `cmp(Rational,FloatAP)` is also unexecuted
(the conditions have false counts but no true counts).

**Classification and hypothesis:** These are reachable public operations,
not impossible implementation branches or compiler-only artifacts.
Testing finite rational operands strictly below, equal to and strictly
above the FloatAP operand exercises the six left-hand overloads and both
short-circuit outcomes of `<=` and `>=`. Testing a NaN on either
side of each `cmp` verifies the unordered outcome without converting
nonfinite values to exact rationals.

**Test-only intervention:** Extend the existing
`TestFloatRounding<ArbitraryPrecision>::test_comparison()` with
explicit observable assertions for each ordered outcome and for the
four operand-order/NaN combinations. Do not change production arithmetic,
weaken expectations, or add arbitrary template instantiations.
The next user-supplied coverage report must establish actual changes
in official per-file counters and confirm the report's base commit.
No new coverage improvement or new CI success is claimed for this change
until measured. This is the first coverage hypothesis for Archive 91.

## Coverage report 92 and FloatAP/FloatDP mixed comparisons (17 October 2026)

**Baseline and measured delta:** The user-supplied `Archivio 92.zip`
contains the LLVM coverage report following `02fca0f4`.
The published `remove-mpfr` branch remains at that commit before this
new intervention. The official per-file counters are:

| File | Archive 91 functions / lines / branches | Archive 92 functions / lines / branches |
| --- | --- | --- |
| `src/floatap.cpp` | 181/199, 1476/1753, 695/978 | 181/199, 1478/1753, 699/978 |
| `include/numeric/floatap.hpp` | 84/102, 197/252, 56/110 | 90/102, 207/252, 64/110 |

Thus the preceding Rational-left comparison intervention added six
covered functions, twelve covered lines and twelve covered branch
outcomes across the two files. Archive 92's per-line records show
`floatap.hpp:417-436` executing six times per operator; at lines
431 and 435 both true and false outcomes are now covered on both
short-circuit comparisons. The unordered paths in the two
FloatAP/Rational `cmp` functions account for the source delta.
A green CI is not a 100% coverage claim.

**Next per-symbol gap:** Archive 92's `src/floatap.cpp` records show
the dedicated `cmp(FloatDP, FloatAP)` function unexecuted at
lines 1013-1015; of the twelve mixed-comparison operators at
1017-1039, only `operator==(FloatAP,FloatDP)` has executions (four).
The four compound `<=`/`>=` expressions at 1022, 1026, 1034,
1038 therefore contain eight branch records with zero true/false
counts. The corresponding `FloatAP`-left `cmp` function does run,
but that does not exercise the other public operators.

**Classification and targeted test:** These are valid public,
statically distinguished overloads, not impossible invariants.
In the existing `test_comparison()` function, only for
`ArbitraryPrecision`, compare a FloatAP holding -1.25 against
FloatDP values -2.25, -1.25 and +1.25. Exercise all six operator
families in both operand orders and assert their mathematical truth
values in each of the three ordered outcomes. The operands are exact
in binary, so the test has no rounding-dependent assumptions.
The predicted effect is execution of eleven previously unused
mixed-operator functions and `cmp(FloatDP,FloatAP)`, plus true and
false outcomes for all eight short-circuit branch records.
The actual aggregate effect must be verified in the next LLVM report;
no coverage improvement or new CI success is claimed yet.

**Delivery rule:** The only source modification is to the test;
no FloatAP production code or shared assertions have been weakened.
Keep the journal and test together on `remove-mpfr`.

## New-machine performance baseline and direct GMP result rounding (17 October 2026)

**Reference:** the user changed benchmark machines; therefore the
provided `Markdown incollato(20261017-185756).md` is the ONLY
performance baseline for this optimization cycle. It contains
300 iterations and 7 timed samples per operation at 53, 128, 256,
and 1024 bits, with 27 operation classes per precision. The total
108-point geometric-mean slowdown is **2.160x**, time-weighted ratio
**2.771x**, and worst case **7.367x** (`atan_nearest`, 53 bits).
70 of 108 ratios are >=2x. The required aspirational target is a
global worst-case FloatAP/FloatMP ratio **strictly below 2x** with
geometric mean as close as possible to 1x. All before/after timings
must be collected on the same host, Release configuration, inputs,
iterations and sample count. Older machine timings are not comparable.

**Evidence and bottleneck:** FloatAP's add/sub kernel previously
constructed two independent GMP temporaries (`exact` and `one`),
aligned and combined them, then invoked `from_ratio(exact,1,...)`,
which allocated its own FloatAP and copied the exact integer a
second time into the result significand. The multiplication kernel
followed the same pattern. On the new machine, `sub_nearest` is
4.11x, 4.73x, 4.40x, 3.85x slower at the four precisions.
Nearest addition is 2.68x, 3.68x, 3.21x, 3.42x; nearest
multiplication is 2.79x, 2.15x, 1.64x, 1.41x.
The old general-rational path is not needed for these exact binary
integer intermediates.

**Optimization hypothesis:** split the existing integer-denominator
rounding into `round_scaled_integer(FloatAP&,Exponent,Rounding)`.
This helper works directly on a nonzero, already-owned positive GMP
significand, with its sign stored as a separate bit; it retains
the same guard/sticky/tie-to-even checks and checked exponent arithmetic.
`from_ratio(n,1,...)` calls this same helper after copying the
numerator exactly once, preserving all other call sites.
`add_sub` and `mul` instead construct the output first, do the
GMP accumulation directly inside its significand and round in-place,
eliminating the extra `mpz_t` temporaries and one integer copy.
The very-large-exponent signed-sticky branch remains unchanged.
Exact cancellation preserves the pre-existing signed-zero rule
(downward yields negative zero). For opposite-sign additions, the
sign is adjusted if the subtraction of magnitudes underflows.
NaN/infinity/zero cases, exception checks and precision propagation
are untouched.

**Validation requirements:** existing differential FloatAP-versus-MPFR
tests at 7, 53, 128, 256 and 1024 bits cover all four rounding modes,
including exact cancellation and far exponent gaps.
Rebuild the full project and run all affected tests/CI before accepting
this production refactor. Use a fresh same-machine 300x7 benchmark to
record per-operation and aggregate deltas, including possible regressions
in the already-fast operations. The hypothesis predicts improvement
primarily in add/sub/mul; no result or attainment of <2x is claimed
until measured. Treat transcendental, division, and square root as
separate follow-up algorithmic work, not as effects of this change.

## Warmed-up benchmark and non-unit quotient storage (17 October 2026)

**Benchmark replacement:** The user provided a further same-machine
measurement `Markdown incollato(20261017-190359).md` after the prior
warm-up run `Markdown incollato(20261017-190253).md`. This latest
report is the performance reference for subsequent improvements.
It measures 27 operations x 4 precisions, 300 iterations and 7
samples. Previous tables remain historical, not calibration targets.

| Bits | New geometric mean | New time-weighted ratio | Worst case |
| --- | ---: | ---: | --- |
| 53 | 1.202x | 2.088x | 3.851x (`atan_nearest`) |
| 128 | 1.938x | 2.857x | 4.290x (`exp_nearest`) |
| 256 | 1.600x | 2.013x | 3.674x (`exp_nearest`) |
| 1024 | 1.471x | 2.948x | 3.663x (`acos_nearest`) |
| ALL | **1.530x** | **2.698x** | **4.290x** |

Across 108 operation/precision pairs, **45** remain >=2x and
**18** are >=3x (computed from reported median times). Compared to
the 19:02 warmed-up measurement, the geometric mean decreased from
1.816x to 1.530x and 15 additional cases fell below 2x; however,
the worst-case ratio remained above 4x. Notable FloatAP absolute
changes: 53-bit `add_nearest` 127.223->62.780 ns, `sub_nearest`
125.557->51.667 ns, 128-bit `mul_nearest` 76.250->37.363 ns,
256-bit `sub_nearest` 118.887->53.333 ns, and 1024-bit
`add_nearest` 100.557->63.057 ns. These results are consistent with
the direct-GMP-result refactor `5bbfd92e` but neither report records
an exact measured commit, and CPU frequency/warm-up effects also
altered MPFR timings. Do not treat the differences as a controlled
performance attribution. The current success criterion remains all
108 ratios **strictly below 2x** with mean as close to 1x as feasible.

**Next targeted bottleneck:** Even after the arithmetic refactor,
`div_down` and `div_nearest` at 128 bits remain approximately
4.006x and 3.858x slower than FloatMP. The general `from_ratio`
path (nonunit denominator) previously made separate GMP
absolute-value copies of numerator and denominator and a temporary
quotient, then copied that quotient into the return value. Unlike
the denominator=1 path, division still used five GMP temporaries.
The new nonunit path reuses positive operands directly, builds an
absolute-value scratch only for a negative operand, computes the
quotient directly in `result._significand`, and compares the exact
remainder against the active divisor. Both quotient normalization,
nearest-even tie handling, directed rounding and exponent carry
checks remain algebraically identical. The same behavior applies
to Rational construction and signed division; negative denominators
remain supported at the internal helper level.

**Differential regression:** Extend the historical shared
FloatAP-vs-FloatMP precision regression with both operand signs,
denominators 3/7/16 (inexact and exact), precisions 7,53,128,256,1024,
and all four rounding modes, testing both division and Rational
construction against MPFR. This does not relax any old assertion.
CI and the next same-host benchmark must verify both correctness
and any actual speedup. No performance delta is claimed for this
nonunit-quotient change before measurement. Follow-up needs separate
work on exp, acos/asin, and general trigonometric certification.

## Benchmark 19:08 and dyadic certified-arithmetic fast path (17 October 2026)

**New reference:** the warmed-up user report `Markdown incollato(20261017-190811).md`
replaces the 19:03 measurement for subsequent experiments. It has
27 operations at 53/128/256/1024 bits, each measured at 300 iterations
and seven samples. Ratios (FloatAP / FloatMP) are:

| Bits | Geomean | Time-weighted | Worst |
| --- | ---: | ---: | --- |
| 53 | 1.134x | 2.001x | 3.799x (`atan_nearest`) |
| 128 | 1.852x | 2.717x | 4.019x (`exp_nearest`) |
| 256 | 1.532x | 1.967x | 3.601x (`exp_nearest`) |
| 1024 | 1.438x | 2.936x | 3.653x (`acos_nearest`) |
| ALL | **1.467x** | **2.661x** | **4.019x** |

Across all 108 pairs, 42 ratios are >=2x, 13 >=3x and 1 >=4x.
The 19:03 report had 1.530x geomean, 2.698x total-time ratio,
4.290x worst, 45 ratios >=2x and 18 >=3x. Selected FloatAP division
improvements: nearest at 53 bits 231.250 -> 176.803 ns; at 128 bits
229.860 -> 188.057 ns; at 1024 bits 365.000 -> 312.917 ns.
53-bit construction from rational 314.167 -> 226.390 ns.
These observations are consistent with the nonunit-rational change
in `c2a0db3d`, but neither report independently records the
measured commit; CPU warm-up/frequency and MPFR variation also matter.

**Cross-cutting optimization:** Certified endpoints returned by
exp/log/trigonometric algorithms are often exact dyadic rationals.
Previously `FloatAP::from_rational` treated their denominator as
general and performed GMP quotient/remainder. For a canonical
dyadic `n/2^k`, the same exact value is `n * 2^-k`, so
detect the single-set-bit denominator and reuse the
already-tested `round_scaled_integer` path for the signed
integer numerator and scale `-k`. A shift larger than the
representable Exponent falls back to the existing path.
The existing nearest-even, directed rounding, carry-overflow and
special-value contracts are unchanged; non-dyadic denominators
still use the general division implementation.

**Argument reduction:** Repeated `x = x / Fraction(2)` in
`exp_interval`, `log_interval`, and `sine_cosine_intervals`
creates temporary rationals and invokes general division.
Replace only these operations with a new in-place
`Fraction::halve()` that uses `mpq_div_2exp`. It is exact,
including negative inputs, and leaves all argument thresholds,
Taylor recurrences, certified error bounds and reconstruction
formulas unchanged. A standalone GMP harness compared 961 canonical
signed rational cases divided by powers of two to denominator scaling,
confirming identical exact fractions. Standalone exploration of
different exponential reduction limits yielded workload-dependent
results, so the reduction thresholds remain unchanged.

**Semantic regression:** New FloatAP/FloatMP directed-rounding
checks exercise dyadics with numerator immediately above or below
a precision boundary, both signs, denominator shifts 1/7/31/130,
precision 7/53/128/256/1024, and all four modes. The existing
shared differential exp/log/trigonometric suite must pass.
Full CI and a fresh 300x7 benchmark on the warmed-up user's machine
are required before accepting a measured speedup. None is claimed.
Global target: worst <2x, geomean approaching 1x.

## Benchmark 19:16 and one-recurrence arctangent (17 October 2026)

**New benchmark:** `Markdown incollato(20261017-191618).md`, same new
machine, 300 iterations, seven samples and four precision groups.
Latest result: geomean **1.463x**, time-weighted **2.617x**, worst
**4.121x** (`atan_nearest`, 53 bits). By precision: 53 bits
1.204x / 2.143x / worst 4.121x; 128 bits 1.753x / 2.631x /
worst 3.924x (`tan_nearest`); 256 bits 1.506x / 1.880x /
worst 3.574x (`tan_nearest`); 1024 bits 1.441x / 2.871x /
worst 3.612x (`acos_nearest`), where each row lists geomean,
time-weighted, and worst. Out of 108 points **49** are >=2x,
**11** >=3x, and **1** >=4x; the prior 19:08 report had 42,
13, and 1 respectively. Aggregate geomean and weighted ratio
improved very slightly from 1.467x and 2.661x, but the global
worst worsened from 4.019x. The dyadic/halving change in
`f93538c3` is consistent with a reduction in measured
`exp_nearest` FloatAP time by 25%, 17%, 20% and 14%
at 53/128/256/1024 bits respectively; benchmark metadata does
not include measured commit identity. Unix and Coverage CI for
`f93538c3` passed at this inspection; Windows was still running.

**Next measured bottleneck:** `atan_nearest` at 53 bits costs
3555.277 ns against 862.780 ns for MPFR (4.121x). The
previous `small_atan` maintained upper and lower fixed-point
power recurrences separately, requiring two independent GMP
multiplications and directed shifts for each Taylor term.
The same kernel is used by atan, pi interval computation,
and inverse-trigonometric functions.

**Certified single-recurrence replacement:** for `|x|<=1/2`,
set `Q=2^B`, `a=floor(Q|x|)`, `b=floor(a*a/Q)`,
`p_0=a`, and `p_(n+1)=floor(p_n*b/Q)`.
Let `P_n=Q|x|^(2n+1)`, `e_n=P_n-p_n`. Since
`0<=Q*x*x-b<2` and `p_n/Q<=|x|<=1/2`,
`0<=e_(n+1)<e_n/4+2` with `e_0<1`, and thus
`0<=e_n<4` for all n. The exact scaled Taylor term
`P_n/(2n+1)` exceeds the truncated integer
`floor(p_n/(2n+1))` by less than 5 units. Therefore an
alternating partial sum of n+1 such integer terms differs
from the exact corresponding sum by at most `5*(n+1)`.
The complete alternating tail has the sign of the next
term and magnitude at most `P_(n+1)/(2n+3)`,
which is below `p_(n+1)+4`. Subtract/add the
term-error bound to both endpoints and the tail to the
appropriate endpoint depending on parity. Divide endpoints
by Q. For negative x, negate and reverse the endpoints.
These are rigorous inequalities, not heuristic truncation
error estimates; the original guard bits and stopping rule
are retained.

**Independent check:** 595 pseudo-random and boundary rational
arguments across work=48,85,160,288,1056 were compared
against high-precision mpmath evaluations; every generated
interval enclosed the numerical reference. This is numerical
support, NOT a formal proof of containment; the inequalities
above establish mathematical correctness. A standalone GMP
microbenchmark has not been used to claim a speedup.

**Ariadne regression:** compare the fully rounded FloatAP
and FloatMP values for atan at 7/53/128/256/1024 bits in
all four rounding modes, both signs, very small arguments
and values around the 1/2 and 2 argument-reduction thresholds,
without altering existing transcendental tests.
The new single-recursion algorithm and test expectations
must pass Unix, Windows and Coverage CI. The performance
effect awaits a fresh same-host benchmark: neither a
<2x worst ratio nor a new speedup is claimed yet.
