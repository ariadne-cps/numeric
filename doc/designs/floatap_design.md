# FloatAP — arbitrary-precision floating point

This is the living design and user manual for the **experimental** `FloatAP`
backend. Keep this document synchronized with every implementation change.
The historical MPFR-removal work is tracked on `remove-mpfr`. The
GMP-free continuation is tracked on `remove-gmp`; see
[`integer_rational_design.md`](integer_rational_design.md) for the shared
unsigned limb core and the `Integer`/`Rational` migration.

## Active implementation policy — architecture-specific substrate (18 October 2026)

**CURRENT, superseding earlier architecture-neutrality restrictions
ONLY for the internal native multiprecision substrate:**
independently developed, explicitly documented ARM64/AArch64 and
x86-64 fast paths are now **permitted**. Keep a fully correct
portable C++20 implementation and test it as a first-class
reference, not merely as a rarely compiled fallback.
This is the user's deliberate exception to the original
"no hardware-architecture-specific implementation" rule
recorded earlier in this journal. Earlier historical
commit notes remain historical, not current prohibitions.

No third-party multiprecision engine may be integrated or copied.
**The low-level unsigned limb substrate is SHARED** across
the future GMP-free `Integer` (replacing `mpz_t`),
`Rational` (replacing `mpq_t`) and `FloatAP`:
it belongs to the common internal numeric arithmetic layer
described in `integer_rational_design.md`, NOT solely
inside `FloatAP`. Shared generic primitives may have
their own properly documented internal source files.
**Only FloatAP-specific** representation, rounding,
fixed-point certification and floating-point algorithms
must stay inside the existing `FloatAP` class and its
`floatap.hpp`/`floatap.cpp` implementation.
No independent third-party engine, public backend abstraction,
or duplicated FloatAP-specific source files.
The optimization boundary is **shared limb primitives**,
not architecture-dependent floating-point semantics:
one exact mathematical and rounding contract for all
implementations.
Architecture-specific code must pass the detailed documentation,
dispatch, correctness, and performance requirements in
"Architecture-specific limb kernels: mandatory engineering
contract" later in this document.

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

## Benchmark 19:22 and one-root certified cosine for tangent (17 October 2026)

**Latest warmed-up same-machine measurement:** `Markdown incollato(20261017-192212).md`,
300 iterations, 7 samples, 27 operations at 53/128/256/1024 bits.
This report replaces 19:16 as the performance reference for subsequent
work. The 108-pair geometric mean is **1.398x** (previous 1.463x);
time-weighted ratio is **1.805x** (previous 2.617x); global worst is
**3.719x**, `tan_nearest` at 128 bits (previous 4.121x).
42/108 ratios remain >=2x and 10/108 >=3x, versus 49 and 11
previously. By precision (geomean, time-weighted, worst):
53 (1.145x,1.963x,3.292x atan), 128 (1.765x,2.293x,3.719x tan),
256 (1.440x,1.496x,3.561x tan), 1024
(1.314x,1.806x,3.138x tan). Measurements do not embed an exact
code revision, so do not claim an isolated causal attribution.
Comparing the 19:16 and 19:22 FloatAP times suggests the
single-recurrence atan change greatly helped atan (19-52%),
asin (13-47%) and acos (14-46%). Unix and Coverage on
`bd412a26` succeeded; Windows was in progress on inspection.

**Discarded optimization:** A fast tangent ratio computed by
outward-rounded integer division of dyadic rational endpoints
was investigated first. An independent GMP harness verified its
enclosure across 120,000 random dyadic interval divisions.
However, native GMP microbenchmarks showed that it was slower
than the existing exact `mpq_div`, particularly at 1024 bits
(roughly 462ns versus 141ns for one representative ratio).
No such dyadic-quotient shortcut was merged.

**Next targeted kernel:** `small_sine_cosine(...,both)` is used
by `tan` and by double-angle reconstruction. It currently computes
both cosine bounds from two large integer square roots, even when
the sine enclosure is very tight.

**Certified one-root construction:** Let `Q=2^fixed_bits` and
`[L,U]` be the already certified scaled sine magnitude interval,
clamped to `[0,Q]`. For `U<=8Q/9`, define
`f(s)=sqrt(Q^2-s^2)`. Its derivative magnitude satisfies
`|f'(s)|=s/sqrt(Q^2-s^2) <= 8/sqrt(17) < 2`
on `[0,U]`. Therefore
`f(U) <= Q*cos(|x|) <= f(L) <= f(U)+2*(U-L)`.
With `c=floor(f(U))`, the integer interval
`[c, c+1+2*(U-L)]` rigorously encloses `Q*cos(|x|)`,
including exact-square cases. It needs only one `mpz_sqrt`
rather than two square roots and a remainder. When
`U>8Q/9`, keep the original exact two-root construction.
The guard is an exact GMP integer comparison (9U<=8Q),
not a floating-point heuristic. Sine, tangent pole handling,
rounding decisions, Taylor error bounds and argument reduction
are unchanged. This is an *outward* enclosure, never an
uncertified cosine approximation.

**Independent check:** A standalone native GMP harness generated
80,000 integer sine intervals at scales corresponding to
53,128,256,1024 bits, including near-sin(1) values, and verified
the new cosine interval always contains the previous exact
two-root interval; the lower endpoints agreed. In that
isolated harness, computing the cosine bounds was roughly
1.5-2x faster with one root. This is not an end-to-end
Ariadne speedup claim.

**Regression and acceptance:** Extend FloatAP-vs-FloatMP tangent
differential rounding tests for both signs, small arguments,
values at or near |x|=1 and beyond (double-angle reconstruction),
precisions 7/53/128/256/1024, and all four rounding modes.
No historical test or exception contract was weakened.
Run Unix, Windows and Coverage CI, then compare a fresh
warmed-up 300x7 benchmark on the user's machine to this
19:22 report. The numerical and performance outcome of this
new commit is not yet established in Ariadne. Global objective:
worst strictly <2x, with geometric mean near 1x.

## Controlled 1000x11 tan benchmark, CI failure and rollback (17 October 2026)

**Measurement protocol:** The user explicitly identified
`Markdown incollato(20261017-192803).md` as the benchmark of the code
**before** `368348c0`, namely `bd412a26` if the checkout was clean,
and `Markdown incollato(20261017-192855).md` as the benchmark
**after** the one-root cosine change. Both reports use the same
machine and `iterations=1000, samples=11`. This new 1000x11
pre-change report supersedes older 300x7 runs as the benchmark
baseline; do not compare these protocols to attribute performance
changes. Neither report embeds its own commit SHA, so the precise
code identities are provided by the user rather than machine-verified.

| Metric (all 108 pairs) | Before 368348c0 | After 368348c0 |
| --- | ---: | ---: |
| Geometric mean FloatAP / FloatMP | 1.389544x | 1.388666x |
| Total-time FloatAP / FloatMP | 1.811418x | 1.817842x |
| Worst ratio | 4.666x (atan, 53 bits) | 4.742x (atan, 53 bits) |
| Ratios >= 2x | 36/108 | 36/108 |
| Ratios >= 3x | 13/108 | 13/108 |
| Ratios >= 4x | 1/108 | 1/108 |

**Measured `tan_nearest` per precision:**

| Bits | FloatAP ns before -> after | FloatAP reduction | Slowdown before -> after |
| --- | --- | ---: | --- |
| 53 | 2592.750 -> 2537.250 | 2.14% | 3.300x -> 3.236x |
| 128 | 3524.166 -> 3420.083 | 2.95% | 3.737x -> 3.385x |
| 256 | 4785.541 -> 4671.333 | 2.39% | 3.339x -> 3.299x |
| 1024 | 18061.542 -> 17363.250 | 3.87% | 3.175x -> 3.029x |

The 128-bit ratio improvement is partly due to FloatMP itself
slowing from 943.167 to 1010.375 ns (+7.13%). The absolute
FloatAP gain is ~3%, not ~9%. The 53-bit atan worst ratio
worsened 4.666x -> 4.742x; this is outside the modified
kernel and is not evidence that the one-root code caused it.
These observations fall far short of the global worst <2x goal.

**CI blocker:** GitHub Actions for `368348c0`
(Unix workflow run `39062535889`, Windows
`39062535867`, Coverage `39062535783`) all finished
**failure** on 17 October 2026. The failing steps were
`Test C++`, `Test`, and `Build, Test and Generate Code Coverage`,
respectively, across macOS, Linux GCC/Clang/TSAN, and Windows.
GitHub Actions for its parent `bd412a26` had succeeded
on Unix, Windows and Coverage. This establishes a regression
in the new change set, although the available GitHub connector
does not expose the detailed failing assertions or full job logs.
Do not misdiagnose the exact test or assert a confirmed root cause
without that evidence.

**Engineering decision:** Revert **both** implementation and tests
introduced by `368348c0`, restoring the last CI-green
`bd412a26` file contents for `src/floatap.cpp` and
`tests/test_float_rounding.cpp`. Retain this journal,
the recorded 1000x11 before/after performance experiment,
and the one-root mathematical exploration as *rejected pending
debugging*, not as accepted code. The restored tree has
the same numerical implementation as the pre-change benchmark.
A new CI run of the rollback commit is still required before
claiming it passes.

**Next priority:** a measurable algorithmic speedup rather than
a marginal hot-path micro-optimization: `atan` at 53 bits
remains worst at 4.666x, followed by `tan` 128 bits
at 3.737x (pre-change). Investigate with targeted profiling,
preserve certified rounding and the full old tests, and
repeat the same warmed-up 1000x11 protocol for comparisons.

## Postmortem: one-root tangent CI failure is a test fixture error (17 October 2026)

**New evidence:** The user provided `logs_105860848323.zip` with the
full failing GitHub Actions logs for `368348c0`. Unix GCC 13,
Clang 16, AppleClang 21 and ThreadSanitizer all report one failing
C++ executable, `test_float_rounding` (18 of 19 CTest entries pass).
They identify the same exception inside `test_precision_regressions()`:

```
FloatAP: Dyadic is not exactly representable at the requested precision
tests/test_float_rounding.cpp:143: calling test_precision_regressions()
```

**Root cause (test fixture, not proven RTL/production error):**
The tangent regression introduced by `368348c0` iterates precisions
`{7,53,128,256,1024}` and includes
`Dyadic(1023,10u)`, `Dyadic(1025,10u)` as well as values near
the `|x|=1` reconstruction boundary. Its construction expression
`FloatAP a(argument,ap)` requests **exact** representation of
each dyadic input. Values with 10 or more significant bits cannot
be represented exactly with a 7-bit FloatAP significand; the
constructor correctly throws. This is not evidence that the
one-root cosine enclosure violated its mathematics, and the
reported historical `atan` slowdown was unrelated to this path.
The earlier rollback to `9ed3a462` was conservative given the
absence of logs, but its generic description of a production
regression is superseded by this diagnosed test-fixture failure.

**Corrected test:** Restore the *same* one-root production
implementation from `368348c0`, restore its full tangent test
matrix (no removal of values or precisions), but explicitly
round the dyadic input using
`FloatAP a(Rational(argument), FloatAP::ROUND_TO_NEAREST, ap)`.
Construct the MPFR oracle from `static_cast<Rational>(a)`
at the identical precision. All four rounding modes of the
`tan` operation remain tested, independent of this input
pre-rounding. Do not change the exact `Dyadic` constructor
contract or weaken any original differential assertions.

**Performance evidence and decision:** Controlled same-host
1000-iteration, 11-sample benchmark before/after `368348c0`
found 2.1%, 3.0%, 2.4%, 3.9% faster FloatAP tangent at
53/128/256/1024 bits, respectively. That is a modest real
optimization, not a solution to the worst-case threshold;
geometric mean 1.390x -> 1.389x, aggregate total ratio
1.811x -> 1.818x, worst 4.666x -> 4.742x due to `atan`
53-bit variation. Reintroducing the one-root optimization is
justified **only if the corrected tests pass CI**. Keep all
historical tests and the high-precision rounding contracts
intact. On failure distinguish implementation error, test
fixture and CI/environment before another rollback.
This corrective commit does not establish CI success yet.

## CI confirmation for corrected tangent regression (17 October 2026)

The corrective commit `b3a5522da833ab7db4bd37347f7cf1d3e49e5ba7`
restored the one-root cosine kernel with the original full tangent
differential regression matrix, changing only the test fixture to
explicitly round dyadic inputs before constructing FloatAP.
**All three GitHub Actions workflows for this exact commit passed**
on 17 October 2026:
- Unix: run 39063273221 — success.
- Windows: run 39063273235 — success.
- Coverage: run 39063273251 — success.

This validates the correction against the supported CI matrix.
The prior failure at `368348c0` was due to requesting exact
construction of dyadics that exceed 7-bit precision; the
one-root implementation is retained. No production code was
changed between `368348c0` and `b3a5522d`, so the user-supplied
post-one-root `1000 iterations / 11 samples` benchmark at 19:28:55
remains the latest performance reference: geomean **1.389x**,
time-weighted **1.818x**, worst **4.742x** (`atan_nearest`,
53 bits), 36/108 ratios >=2x. The earlier 19:28:03 report is
the corresponding pre-one-root baseline. Do not infer any
additional performance benefit from the test-only correction.
The next performance work must target much larger gains on
`atan` 53-bit and `tan`/`cos` at 128 bits while preserving
MPFR differential tests and all four rounding modes.

## Certified atan(1/2)-centered argument reduction (17 October 2026)

**Starting point:** CI-green source at `b1602707` (same numerical
implementation as `b3a5522d`). The latest user-supplied benchmark
is `Markdown incollato(20261017-192855).md`, using 1000 iterations,
11 samples on the warmed-up new machine. It reports 108 paired
observations, geomean **1.389x**, time-weighted **1.818x**, worst
**4.742x** (`atan_nearest`, 53 bits) and 36 ratios >=2x.
This is the pre-intervention baseline; no new timing is yet available.

**Why the current atan algorithm is expensive:** The benchmark's
32 rational operands at each precision are
`(2^bits + 13 + 2*i) / (2*2^bits + 17 + 3*i)`,
rounded to the target FloatAP/FloatMP format. They cluster just
*above* 1/2. The previous `atan_interval` switches at 1/2 to
`atan(x) = pi/4 + atan((x-1)/(x+1))`; the reduced argument
is near -1/3. An independent GMP implementation of the existing
fixed-point recurrence found 27, 50, 91, 333 Taylor steps
for the 32 exact benchmark operands at 53,128,256,1024 bits,
respectively (work=bits+32), even after pi is cached.
Centering instead at 1/2 reduces these to **one step each**.
These are recurrence iteration counts, not end-to-end timings.

**General exact transformation:** For positive
`x in [1/4,5/7]`,
`atan(x) = atan(1/2) + atan((2*x-1)/(x+2))`.
The reduced argument is monotone in x, taking values in
`[-2/9,3/19]`, hence safely within the domain `|r|<=1/2`
of the existing certified single-recurrence Taylor kernel.
The atan addition identity has no branch ambiguity for
positive x in the chosen interval. For negative x use oddness,
negating and reversing the certified interval. Values with
`|x|<1/4` use the prior direct Taylor series; positive
`x>5/7` uses the unchanged pi/4 or pi/2 reduction.
All interval operations remain exact.

**Certified reusable constant:** Add `atan_half_interval(work)`
with a `thread_local` cache of the *interval*, not a rounded value,
recomputing only if the requested work precision exceeds the cached
precision. Its endpoints are obtained from the already validated
`small_atan(1/2,work)` kernel; reusing a higher-precision enclosure
at lower precisions is always safe. Summing both certified intervals
preserves outward bounds. The outer `certified` loop retains
exactly the previous test: rounded lower and upper endpoints must
agree in the requested mode, otherwise increase work precision.
This preserves correctly rounded results and all four modes,
including ties, directed rounding and NaN/infinity handling.

**Independent checks:** A GMP C++ harness verified 50,000 exact
rational reduction identities (including the boundaries), and
the bound on the reduced magnitude. A separate GMP C++ probe
reproduced the 32 benchmark arguments and old-versus-new Taylor
iteration counts at all four advertised precisions.
Neither check replaces an end-to-end Ariadne build or a controlled
performance benchmark.

**Ariadne regression:** Extend the historical shared
`test_precision_regressions()` for `ArbitraryPrecision`
with operands just below, equal to and just above 1/4,
1/2 and 5/7, plus low magnitudes and the existing reduction
region; both signs; six requested precision steps
(53,7,128,1024,256,53) to exercise a cache upgrade and
reuse; and all four IEEE rounding modes. Construct FloatAP
inputs via explicit rounding to avoid violating the exact Dyadic
constructor, then give the MPFR oracle the identical binary input.
The existing sine/cosine/tan/asin/acos and pi tests remain unchanged.

**Acceptance:** Validate Unix, Windows and Coverage CI, then
collect a same-machine, same-build, fully warmed-up benchmark with
1000 iterations and 11 samples; compare all 108 ratios, especially
`atan_nearest` 53/128/256/1024 and derived
`asin`/`acos`. The aspirational global target is worst
strictly below 2x and geomean near 1x. No improvement is claimed
until the new benchmark is measured; if regressions occur, classify
implementation, test and environment causes explicitly.

## Benchmark 19:48: atan acceleration and inverse-trigonometric regression (17 October 2026)

**Incoming measurement:** `Markdown incollato(20261017-194856).md`
is the first warmed-up 1000-iteration/11-sample benchmark after
`53f00b55`. Its pre-change comparator is
`Markdown incollato(20261017-192855).md` (same protocol and host).

| Metric (108 operation/precision pairs) | Before | After |
| --- | ---: | ---: |
| Geometric mean FloatAP/FloatMP | 1.389x | **1.351x** |
| Time-weighted FloatAP/FloatMP | 1.818x | **1.580x** |
| Worst | 4.742x atan@53 | **4.043x acos@53** |
| Ratios >=2x | 36/108 | **36/108** |
| Ratios >=3x | 13/108 | **13/108** |

The new half-centered `atan` reduction is highly effective for
standalone `atan_nearest`: FloatAP times at 53/128/256/1024 bits
changed from 2601.000/3792.208/5639.750/36879.416 ns to
1718.750/2080.542/2223.834/2572.000 ns, reductions of
33.9%/45.1%/60.6%/93.0%. FloatAP atan is now faster than MPFR
at 128, 256 and 1024 bits for these benchmark operands.

**Unacceptable unintended costs:** `asin` became slower by
27.8%/24.0%/16.0% at 53/128/256 bits, and `acos` by
23.4%/22.1%/14.2%; both improve slightly at 1024 bits.
The latest global worst is therefore `acos` at 53 bits.
Unix, Windows and Coverage CI on `53f00b55` all passed, so
this is a measured performance regression, not a known
numerical-correctness failure.

**Root-cause hypothesis:** `asin(x)=2*atan(x/(1+sqrt(1-x*x)))`
passes an argument close to `2-sqrt(3)=0.267949...` for
`x` near 1/2, the current 32-input benchmark pool. The
centered atan identity activates for this argument because
`1/4 < 2-sqrt(3) < 1/3`. A standalone GMP fixed-point
iteration probe compared the direct kernel and the centered
reduced kernel at 53/128/256/1024 bits (using work=bits+32):
direct requires 22/42/76/278 iterations, whereas the
centered reduced kernel needs 19/35/63/231. The savings
at moderate work precision are insufficient to compensate
for interval-cache combination, extra fraction arithmetic
and certified-constant handling; at 1024 bits, the saved
iterations outweigh the overhead. This is a strong inference
from source and measurements, not an instrumented per-function
profile.

**Next hypothesis: precision-sensitive crossover.** Prefer the
existing direct certified Taylor kernel for `|x|<1/3`
when `work<=512` bits, retaining the current `|x|<1/4`
direct region for higher work precision. In the centered
region `[1/3,5/7]`, the transformed argument lies in
`[-1/7,3/19]`, narrower than before. The standalone atan
benchmark operands near 1/2 stay on the fast centered branch
at every precision. `asin`/`acos` operands near 1/2
should return to the direct branch at moderate precision,
without special-casing the benchmark or changing certified
interval semantics. The 512-bit crossover is a performance
hypothesis to be validated, not a numerical threshold.

**Regression:** Extend the existing FloatAP-versus-FloatMP
atan checks to either side of 1/3, retaining 1/4, 1/2 and
5/7 boundaries and all four rounding modes, signed inputs,
and multi-precision cache ordering. Add signed inverse-
trigonometric checks immediately around 1/2 to the existing
asin/acos suite. No previous assertions are removed.
All extrema are still computed as outward-certified intervals
and the final equality-of-rounded-endpoints criterion is
unchanged. The next CI and warmed-up 1000x11 benchmark must
validate the precision-sensitive dispatch and any performance
claim; do not attribute a speedup to this patch in advance.

## Benchmark 19:54: validated 1000x11 baseline and centered sine/cosine trial (17 October 2026)

**Latest usable baseline:** `Markdown incollato(20261017-195459).md`
is the user's rerun on the same warmed-up machine with the CPU no
longer occupied by Coverage. It is the valid **after-`e9eeb095`**
measurement (1000 iterations, 11 samples, 27 operations at each of
53/128/256/1024 bits). **Discard** the immediately preceding
`Markdown incollato(20261017-195337).md` for quantitative performance
decisions: Coverage was simultaneously using the CPU, as the user
explicitly confirmed. The valid pre-`e9eeb095` comparator is the
19:48:56 1000x11 measurement.

| Aggregate | 19:48 before | 19:54 after |
| --- | ---: | ---: |
| Geomean FloatAP/FloatMP | 1.351x | **1.375x** |
| Time-weighted ratio | 1.580x | **1.541x** |
| Worst slowdown | 4.043x (acos 53) | **3.538x (cos 128)** |
| Ratios >=2x | 36/108 | **44/108** |
| Ratios >=3x | 13/108 | **13/108** |

`asin` FloatAP median reductions versus 19:48 were
22.4%, 17.7%, 12.3% at 53/128/256 bits (2.6% slower at 1024);
`acos` decreased 19.1%, 16.7%, 10.9% (5.4% slower at 1024).
The intended precision-aware atan dispatch recovered the earlier
inverse-trig regressions at the low/medium precisions without losing
the massive 1024-bit atan improvement (~2.57 us vs 24.30 us for MPFR).
Different reference MPFR timings, notably at 53 bits, still move the
ratios: monitor absolute FloatAP time and both backends.
**Unix, Windows and Coverage on `e9eeb095` are all green.**

**Next bottleneck:** cos@128 remains 3.538x, tan@128 3.507x,
tan@256 3.357x and log@53 3.333x. All are much larger than
the 2x worst-case target. The current trigonometric kernel
evaluates the sine/cosine Taylor series around zero on the
benchmark's 32 arguments clustering near 1/2; for tangent,
it derives cosine from sine via one exact integer root.
The new approach exploits a *general mathematical neighborhood*
of +/-1/2, never matching preselected benchmark indices.

**Certified local sine/cosine transformation:** Store certified
intervals `S = sin(1/2)`, `C = cos(1/2)` in a thread-local cache,
recomputing them when greater precision is requested. For
`y=|x|=1/2+d`, Taylor's theorem and the global bounds
`|sin''(t)|<=1`, `|cos''(t)|<=1` give

    sin(y) in S + C*d + [-d*d/2,+d*d/2]
    cos(y) in C - S*d + [-d*d/2,+d*d/2].

Each interval operation is performed with exact GMP rational
arithmetic, including sign-aware multiplication by d.
For a requested work precision `work`, take the shortcut
only when `|d|<=2^-ceil(work/2)` and there were no previous
halvings. Then `d*d/2<=2^(-work-1)`; center enclosures are
computed at `work+16` and certified, so the local error is
smaller than the requested work unit. For other inputs, and
whenever certified refinement increases work sufficiently that
the condition no longer holds, use the old kernel unchanged.
For negative x, sine is odd and cosine even. Crucially, this
is an *outward* interval, not an approximation rounded to the
requested result. No correct-rounding criterion is relaxed.

**Independent check and limits:** A standalone rational-bound
test of the linearized intervals over 600 near-center positive
and negative deviations at 53/128/256/1024-bit scales enclosed
365-digit mpmath references for both sin and cos. The guarantee
comes from Taylor's theorem, not the numerical sample. No end-to-end
speedup is claimed for the shortcut. It may make low-precision
operations slower due to exact-rational overhead; the next 1000x11
benchmark must test all 108 cases and the per-operation absolute
timings before accepting it. A bad result warrants rollback.

**Regression matrix:** Add FloatAP/FloatMP differential sin/cos/tan
checks around +/-1/2 at step widths matching the precision, plus
farther neighbors and out-of-fast-path values, using
7/53/128/256/1024 bits and all four rounding modes. Revisit a
precision after updating the cache to check its reuse.
Keep the historical tests and one-root tangent implementation.
Run Unix, Windows and Coverage CI for this commit before relying
on it. The global target remains worst strictly <2x and geomean
approaching 1x.

## Benchmark 20:00: isolate the centered trigonometric fast path (17 October 2026)

**Measured report:** `Markdown incollato(20261017-200024).md`
is the user's new same-host 1000-iteration, 11-sample benchmark
after `4c75d3a7`, compared with the valid 19:54:59 report
before `4c75d3a7`. The prior 19:53 run remains invalid for
comparison because concurrent Coverage occupied the CPU.

| Aggregate metric | Before 4c75d3a7 | After 4c75d3a7 |
| --- | ---: | ---: |
| Geometric mean (108 ratios) | 1.375x | **1.300x** |
| Total-time ratio | 1.541x | **1.455x** |
| Worst-case ratio | 3.538x (cos@128) | **4.934x (tan@128)** |
| Ratios >=2x | 44/108 | **35/108** |
| Ratios >=3x | 13/108 | **13/108** |

**Direct FloatAP timing effects of 4c75d3a7:**

| Precision | sin_nearest | cos_nearest | tan_nearest |
| --- | ---: | ---: | ---: |
| 53 | +21.3% | +12.0% | **+36.1%** |
| 128 | +5.9% | -0.9% | **+34.0%** |
| 256 | **-14.3%** | **-18.0%** | **+30.4%** |
| 1024 | **-67.5%** | **-67.9%** | -1.8% |

A negative percentage means a faster absolute FloatAP time.
The local cached sin/cos enclosure therefore provides a major
high-precision benefit, but the `TrigOutputs::both` tangent
path is an unacceptable regression. The likely reason is that
combining cached dyadic enclosures with a small exact rational
`delta` and the exact `delta^2/2` remainder enlarges
intermediate denominators. `tan` divides two such endpoints
as general GMP rationals, whereas the previous full Taylor
kernel delivered compact fixed-point dyadics and used one
integer square root for cosine. This is a source-grounded
performance hypothesis, not yet a per-function profiler
trace. On the measured host, the effect is clearly negative
for tangent from 53 to 256 bits.

**Targeted follow-up:** Keep the certified centered
sin/cos fast path ONLY when `halvings==0`,
`work>=256` and a single result (`sine` or `cosine`)
is requested. For `TrigOutputs::both`, retain the
previous Taylor kernel and one-root cosine construction
at every precision; this both prevents expensive local
rational interval division and preserves the same behavior
for trigonometric double-angle reconstruction. The
`work>=256` threshold avoids small-precision cache/rational
overhead and corresponds approximately to precision
>=224 bits at the initial 32-guard-bit call. These
thresholds affect only computational cost, never output
semantics or certification. All existing certified bounds
and final rounded-endpoint equality rules are unchanged.

**Acceptance:** The performance hypothesis predicts
preserving most of the 256/1024-bit sin/cos benefit while
returning tangent to roughly the pre-`4c75d3a7`
timings (subject to machine variation). Never claim
the predicted measurements as facts until a new
warmed-up, unloaded, same-machine 1000x11 report is
available. In the differential rounding test matrix,
include widths 192 and 224 bits around the dispatch
crossover, retaining all existing precisions, signs and
four rounding modes. The previous tests are not weakened.
Confirm Unix, Windows and Coverage CI for the follow-up
commit. Later, optimize tangent independently, potentially
by outward rounding its interval endpoints to compact
dyadics before exact division, subject to a rigorous
error-width analysis and an isolated microbenchmark.

## Benchmark 20:06 and unresolved Windows compile failure (17 October 2026)

**Controlled performance measurement:** User-supplied
`Markdown incollato(20261017-200644).md` is the first
1000-iteration/11-sample benchmark after `5116d004`;
the comparator is `Markdown incollato(20261017-200024).md`
for its direct parent `4c75d3a7`, on the user's warmed-up
machine. The earlier concurrent-coverage 19:53 run is invalid
for performance comparisons. Values are FloatAP/FloatMP
slowdowns on the same 108 operation/precision pairs.

| Measure | 20:00 parent | 20:06 current |
| --- | ---: | ---: |
| Geometric mean | 1.299815x | **1.300768x** |
| Time-weighted total | 1.454818x | **1.419272x** |
| Worst observed ratio | 4.934x (tan@128) | **3.702x (tan@128)** |
| Slowdowns >=2x | 35/108 | **34/108** |
| Slowdowns >=3x | 13/108 | **8/108** |
| Slowdowns >=4x | 3/108 | **0/108** |

**Measured tangent recovery (FloatAP absolute median):**

| Precision | Before ns | After ns | AP runtime change |
| --- | ---: | ---: | ---: |
| 53 bits | 3431.209 | 2518.917 | **-26.6%** |
| 128 bits | 4714.917 | 3478.375 | **-26.2%** |
| 256 bits | 6289.458 | 4744.625 | **-24.6%** |
| 1024 bits | 17708.375 | 17037.167 | **-3.8%** |

As intended, limiting the centered sin/cos cache to *single*
outputs and `work>=256` removes nearly all of tangent's
earlier slowdown at 53/128/256 bits. The high-precision
single-output advantage remains: at 1024 bits, FloatAP
`sin_nearest` is 3362.792 ns versus 5271.625 ns for MPFR
(0.638x), and `cos_nearest` is 3358.250 ns versus 4260.667 ns
(0.788x). Most low/medium single-output timing changes are
small; exact percentages depend on per-backend timing drift.
The 20:06 performance report is the new reference for
subsequent numeric optimization, pending platform validation.

**CI blocking issue discovered in follow-up:** GitHub Actions
for `5116d004`: Unix run `39065259290` succeeded,
Coverage run `39065259252` succeeded, but Windows run
`39065259277` failed **during Build** (Windows 2025 Release,
job 117491463543, step 13). The *parent* `4c75d3a7`
also failed during Windows Build (run `39064807829`,
job 117490109874, step 13); the preceding `e9eeb095`
was green on Unix, Windows and Coverage. Therefore the
portability/compilation failure first appears with the
centered trigonometric code or tests introduced in
`4c75d3a7`, not with the filtering change in `5116d004`.
Detailed MSVC diagnostics are not available through the
current GitHub connector (run-level log endpoint empty;
job-level logs and check annotations are not exposed).
**Do not claim an identified Windows root cause or infer
a mathematical failure** from the Build result. The new
optimization remains pending Windows validation.

**Next action:** obtain the Windows Build diagnostic,
fix the underlying portable C++ issue without changing the
certified rounding contracts, verify all three workflows,
then resume profiling `tan` at 128/256 bits and
`cos` at 128 bits. Do not stack another speculative
numeric optimization on top of an unresolved failing
platform build. Distinguish GitHub Coverage workflow
success from separate Codecov PR patch thresholds.

## Windows MSVC C4459 shadowing fix (17 October 2026)

**Confirmed root cause:** User-supplied Windows 2025 [Release] build
diagnostic for `5116d004`, `src/floatap.cpp:2153`:

```
error C2220: the following warning is treated as an error
warning C4459: declaration of 'near' hides global declaration
include/numeric/rounding.hpp(286): note: see declaration of 'Ariadne::near'
```

**Classification:** C++ portability / compiler warning promoted to
an error by the existing `/W4 /WX /permissive-` flags, **not** an
arithmetic or rounding failure. The local variable in
`sufficiently_near_half` was introduced by the centered
sin/cos fast path in `4c75d3a7` and persisted through
`5116d004`; both commits failed Windows at the Build step.
Other GitHub Unix and Coverage workflows on `5116d004`
passed. Earlier inability to identify the root cause was due to
missing detailed Windows logs; this new diagnostic supersedes it.

**Minimal fix:** Rename local `const bool near` to
`const bool within_half_threshold` and update the corresponding
return. The boolean expression, cleanup, certified enclosures,
dispatch, tests, compiler warning settings and numeric behavior
remain **byte-for-byte identical** apart from the local name and
line wrapping. Do not suppress C4459, relax `/WX`, or change
production mathematics to satisfy the compiler.

**Verification:** The modified source contains no
`const bool near` at the affected site; the exact new
identifier and return were checked against the source before
publishing. A new Windows 2025 [Release] build on this commit,
as well as Unix and Coverage, are required before claiming
all-platform success. The last performance baseline remains the
valid `Markdown incollato(20261017-200644).md` (1000x11):
geomean 1.301x, time-weighted 1.419x, worst 3.702x
(`tan_nearest` 128 bits), 34/108 ratios >=2x. Renaming a
local variable does not establish any new speedup.

## Cached certified tan(1/2) with local tangent evaluation (17 October 2026)

**Starting state:** `5555f262`, which passed all three GitHub Actions
workflows (Unix run 39065677444, Windows run 39065677386, Coverage
run 39065677381). The user-confirmed valid 1000-iteration/11-sample
baseline is `Markdown incollato(20261017-200644).md`, taken after
`5116d004` and before this new change:

| All 108 operation/precision pairs | Current |
| --- | ---: |
| Geometric mean FloatAP/FloatMP | **1.301x** |
| Time-weighted ratio | **1.419x** |
| Worst slowdown | **3.702x** (`tan_nearest` at 128 bits) |
| Slowdowns >=2x | **34/108** |
| Slowdowns >=3x | **8/108** |

`tan_nearest` FloatAP times are 2518.917/3478.375/4744.625/
17037.167 ns at 53/128/256/1024 bits, compared with FloatMP
817.000/939.542/1436.209/5716.167 ns. Another significant
outlier is `cos_nearest` at 128 bits (3.660x).
Earlier attempts to compute tangent from local cached sine and
cosine independently caused ~25-35% slowdown at low/medium
precisions because their general rational endpoints were
divided exactly on every call. The verified `5116d004`
therefore restricts the local sine/cosine path to
single-output functions. This tangent improvement is
**independent** of that route.

**New, general identity for values very near 1/2:** Let
`y=|x|=1/2+delta` and `T=tan(1/2)`. Then by Taylor's
theorem,
`tan(y)=T+(1+T^2)*delta+R`. For
`|delta|<=2^-ceil(work/2)` with `work>=32`,
`y in [1/4,3/4]`. Over that interval
`|sin(y)|<=3/4` and
`cos(y)>=1-y*y/2>=23/32>1/2`.
Consequently `0<=tan(y)<=3/2`,
`sec(y)^2<=4`, and
`|tan''(y)|=|2 tan(y) sec(y)^2|<=12`.
Therefore the rigorous remainder bound is
`|R|<=6 delta^2`. Negative arguments use oddness,
with reversed endpoints. This is not a heuristic
approximation; it certifies an entire open neighborhood
around the center and is independent of benchmark inputs.

**Cached certified endpoint construction:** The existing
`half_trigonometric_intervals(work+48)` gives proven sine
and cosine enclosures at `1/2`, both with positive cosine.
Perform their exact interval division once per thread and
precision upgrade to enclose `T`. Quantize both endpoints
**outward** to dyadics with denominator `2^(work+32)`
using `mpz_fdiv_q` for the lower endpoint and
`mpz_cdiv_q` for the upper endpoint. Cache both T and the
outward-rounded derivative interval
`[1+T.lower^2, 1+T.upper^2]`; its monotonicity is valid
because `T>0`. The resulting per-call operations are
limited to cached exact dyadic additions/multiplications
by delta and the explicit rational remainder interval.
No division of separately reconstructed sine/cosine
intervals occurs on repeated calls. A precision upgrade
recomputes the cache; a tighter cached interval is valid
for lower work precisions. The existing `certified` loop
still compares BOTH outward endpoints rounded in the
requested mode and refines if they differ. Outside the
neighborhood, retain all existing trigonometric paths
and exception contracts, including pole detection.

**Independent numerical check:** An external Python
fraction-arithmetic + mpmath model tested 1,550
containments across 7/53/128/256/1024 bits and both
signs, checking interval enclosure at the requested
working precision. The proof is the derivative bound
above, not the finite numerical sample. The repository
C++ suite and Windows CI have not yet been run for this
new commit.

**Differential tests:** Extend the existing FloatAP/FloatMP
regression to check tangent at both signs, center, points
inside/outside the precision threshold and values
outside the local fast path for 7/53/128/256/1024 bits
and all four rounding modes; use identical rounded
FloatAP inputs for the MPFR oracle. Revisit 53 bits
after a 1024-bit cache upgrade to test reuse. No old
tests or rounding conditions are relaxed.

**Acceptance:** Require Unix, Windows and Coverage success.
Rebenchmark on the user's warmed-up unloaded host with
1000 iterations, 11 samples, same input pool and Release
flags. Expect a major reduction in `tan_nearest`, but
no actual speedup is claimed until measured. Quantify
all 108 ratios and inspect any shifts in the worst case;
the hard objective remains worst <2x with geomean near 1x.

## Benchmark 20:19 and dyadic fixed-point centered sin/cos endpoints (17 October 2026)

**User-measured result:** `Markdown incollato(20261017-201956).md`,
1000 iterations, 11 median samples, warmed-up and unloaded same host.
Compare solely with the immediately preceding
`Markdown incollato(20261017-200644).md` at `5555f262`.
The new report measures the certified centered tangent cache
introduced at `4ae35e8c`; all 108 operation/precision pairs are
present. The new latest baseline is:

| Statistic | Before 4ae35e8c | After 4ae35e8c |
| --- | ---: | ---: |
| FloatAP/FloatMP geometric mean | 1.300768x | **1.234574x** |
| Total-time ratio | 1.419272x | **1.294434x** |
| Worst-case | 3.702x (tan@128) | **3.613x (cos@128)** |
| Pairs >=2x | 34/108 | **31/108** |
| Pairs >=3x | 8/108 | **7/108** |

The new tangent kernel delivered dramatic reductions in
absolute FloatAP `tan_nearest` time:
53 bits 2518.917->1415.083 ns (-43.8%, ratio 1.732x);
128 bits 3478.375->1759.375 ns (-49.4%, ratio 1.805x);
256 bits 4744.625->1879.917 ns (-60.4%, ratio 1.347x);
1024 bits 17037.167->3169.375 ns (-81.4%, ratio **0.553x**).
All four tangent ratios are now <2. The benchmark shift
also exposes `cos_nearest` at 128 bits (3.613x) as the
new worst; `exp_nearest` at 128 bits is 3.272x,
`acos_nearest` at 53 bits 3.255x, `atan_nearest`
at 53 bits 3.165x, and `pow7_nearest` at 53 bits 3.103x.
The aggregate objective worst <2x is NOT attained.
GitHub Unix run 39066193945, Windows run 39066193956,
and Coverage run 39066193962 all ultimately succeeded
for commit `4ae35e8c` (status verified during this intervention).

**Performance source inspection:** The 128-bit cosine still
takes the direct fixed-point Taylor kernel: the old
centered sin/cos cache was enabled only for `work>=256`,
while the initial requested work is `precision+32=160`.
Earlier use of exact GMP `Fraction` arithmetic to center
cosine at this precision had only marginal benefit.
The expensive path performs exact dyadic-rational
multiplications, endpoint additions, repeated canonicalization
and construction. The latest tangent optimization demonstrated
that a cache of *outward-rounded dyadic integer* endpoints
can avoid this overhead.

**New proposed local kernel:** For `y=|x|=1/2+d` and
`|d|<=2^-ceil(work/2)`, use the unchanged certified Taylor
bounds
`sin(y) in S+C*d+[-d²/2,+d²/2]` and
`cos(y) in C-S*d+[-d²/2,+d²/2]`.
On a canonical dyadic `d=n/2^k`, quantize cached
certified half-angle intervals outward to a common
`Q=2^(work+16)` using floor/ceiling integer division.
Store four scaled center endpoints per thread and
upgrade them only when `work` increases. Evaluate
`Q*C*d` and `Q*S*d` with signed GMP integer
multiplication and directed floor/ceiling shifts by k,
choosing reversed coefficient endpoints when n<0.
Since `|d|<=2^-ceil(work/2)`,
`d²/2 <= 2^(-work-1)`; at this scale a
conservative certified remainder is the exact integer
`2^(Q_bits-work-1)`. Add/subtract it outward.
Construct each final interval endpoint as a dyadic
`Fraction(scaled_integer,Q_bits)`. Negative inputs
negate the sine interval, cosine remains even.
This is a rigorous entire-neighborhood enclosure:
there is NO weakening of the rounded-endpoint
agreement test. Keep the previous exact-rational
centered code for non-dyadic offsets and the old
general Taylor kernel for all other inputs.
Do not affect the `TrigOutputs::both` path
or the separately cached centered tangent kernel.

**Dispatch:** Lower the single-output centered-kernel
threshold from `work>=256` to `work>=128`,
thereby exercising it at 128-bit target precision
(initial `work=160`) while leaving 53-bit
(initial `work=85`) on the older, measured-faster path.
The 96/128/192/224/256/1024-bit regressions
cover the threshold and precision-cache upgrades.
The threshold is solely a performance hypothesis,
not a change of function semantics.

**Independent checks (not a substitute for Ariadne CI):**
A high-precision Python/mpmath dyadic fixed-point
model checked 1260 combinations across working
precisions 85,128,160,224,256,288,1056;
both sine and cosine outward intervals contained
the high-precision evaluation for signed dyadic
offsets. The mathematical Taylor bound above
provides the proof, not these finite samples.
A standalone GCC/GMP microbenchmark comparing
the local rational and integer endpoint routines,
including GMP allocations and final rational
canonicalization, measured 550 vs 250 ns at
work=160, 567 vs 272 ns at work=288, and
1014 vs 444 ns at work=1056, for ~2.1-2.3x
local-kernel speedups on the runner. These are NOT
end-to-end FloatAP results and were not measured
on the user's macOS benchmark host.

**Validation and review:** No prior test is relaxed.
Require Unix, Windows and Coverage CI, including
all four rounding modes and cache-upgrade order,
before claiming correctness on supported platforms.
Require a new 1000x11 benchmark on the same
unloaded host; compare all 108 ratios and
FloatAP absolute medians, especially cos/sin at
128,256,1024 and prior-fast tangent. If the
global maximum does not improve or 1024-bit
sin/cos regress, re-evaluate or revert this
isolated change. Do not mistake the standalone
microbenchmark for an official improvement.

## Benchmark 20:25 and signed division without GMP numerator copying (17 October 2026)

**New 1000x11 warmed-up baseline:** user-supplied
`Markdown incollato(20261017-202537).md` is the first full measurement
after the cached fixed-point centered sine/cosine optimization
`04e1b6f3`, compared against
`Markdown incollato(20261017-201956).md` before that commit.
All 108 operation/precision pairs are included. The results are:

| Statistic | Before | After |
| --- | ---: | ---: |
| Geometric mean FloatAP/FloatMP | 1.235x | **1.212x** |
| Total-time ratio | 1.294x | **1.248x** |
| Worst-case | 3.613x (cos@128) | **3.364x (div_up@128)** |
| Cases >=2x | 31/108 | **32/108** |

**Measured sin/cos FloatAP times (ns):**

| Bits | sin before -> after | cos before -> after |
| --- | --- | --- |
| 53 | 1321.167 -> 1384.834 | 1408.458 -> 1454.750 |
| 128 | 1852.541 -> 1275.292 | 1928.292 -> 1269.500 |
| 256 | 2098.375 -> 1350.292 | 2092.459 -> 1345.791 |
| 1024 | 3434.209 -> 1499.209 | 3349.917 -> 1500.417 |

The 128/256/1024-bit sine and cosine gains range from
approximately 31% to 56% in absolute FloatAP time, consistent
with the fixed-point centered interval path. The 53-bit
initial work precision is below its dispatch threshold,
so the 3-5% slower reported times at that precision
cannot be directly attributed to the new path. At 1024 bits,
the absolute FloatAP/FloatMP ratios for sine and cosine are
0.284x and 0.355x. Tangent remains below 2x at all four
precisions. The remaining slowdowns >=2x number 32 (53:8,
128:13, 256:9, 1024:2). New top outliers are
`div_up` 128 bits 146.333/43.500=3.364x,
`acos_nearest` 53 bits 4150.417/1263.292=3.285x,
`atan_nearest` 53 bits 1783.625/548.917=3.249x,
`exp_nearest` 128 bits 2812.208/881.459=3.190x
and `div_down` 128 bits 133.958/43.042=3.112x.
At inspection of `04e1b6f3`, Unix and Coverage
workflow runs `39066603650` and `39066603557`
had passed; Windows run `39066603660` was still running.
Check Windows before claiming three-platform success.

**Next cross-cutting optimization:** `FloatAP::div` has separate
sign metadata in both operands and strictly nonnegative stored
significands. Nonetheless it previously allocated a fresh
`mpz_t`, copied `a._significand`, and negated this copy
if the result was negative. It then passed the temporary to
`from_ratio`, which may allocate another absolute magnitude.
Eliminate the copy by letting the existing private
`from_ratio` helper accept a defaulted `invert_sign`
boolean (default false for EVERY unchanged call site).
Compute the result sign from numerator and denominator signs,
then XOR with `invert_sign`. Pass the original immutable
positive significands directly from `div` with
`invert_sign = a._negative != b._negative`. The entire
quotient, remainder, exact nearest-even and directed rounding,
exponent normalization, exceptional operand handling and
signed zero cases are left unchanged. This is internal
API plumbing, not a new division algorithm or an external
signature change. The bit-exact behavior follows from
moving the sign operation from integer magnitude to
metadata before the same existing rounding operation.

**Regression and measurement:** Existing FloatAP-vs-FloatMP
rational construction/division tests span precisions
7,53,128,256,1024 and all four rounding modes. Expand
the signed-divisor matrix from {3,7,16} to
{-16,-7,-3,3,7,16}, preserving every existing test.
The generated C++ CI and a new unloaded same-machine
1000x11 benchmark are needed before accepting the
change as correct and faster. Predict fewer GMP allocations
and copies, especially for signed division, but DO NOT
claim a measured division speedup yet. The user objective
remains every worst-case <2x and geometric mean near 1x.

## Benchmark 20:28: division mixed results and Windows Build blocker (17 October 2026)

**Benchmark identity:** `Markdown incollato(20261017-202844).md`
is the latest user-provided warmed-up, same-host 1000-iteration,
11-sample report, following `89d93735`. Compare it only to the
20:25:37 1000x11 report measured before that commit. Each contains
27 operations at 53,128,256,1024 bits (108 operation/precision
pairs).

| Metric | 20:25 prior | 20:28 current |
| --- | ---: | ---: |
| Geometric mean AP/MP | 1.212040x | **1.237244x** |
| Total-time ratio | 1.248318x | **1.281647x** |
| Worst slowdown | 3.364x (`div_up`, 128) | **3.446x (`acos`, 53)** |
| Cases >=2x | 32/108 | **34/108** |
| Cases >=3x | 9/108 | **6/108** |
| Cases >=4x | 0/108 | **0/108** |

**Direct division observations:** The sign-metadata refactor in
`89d93735` removed a GMP numerator copy. Its FloatAP
absolute timings improved in 10/12 division/precision pairs.
In nanoseconds (before -> after):
- 53 bits: nearest 122.208 -> 114.000 (-6.7%);
  downward 95.625 -> 79.625 (-16.7%);
  upward 97.500 -> 81.959 (-15.9%).
- 128 bits: nearest 181.541 -> 168.709 (-7.1%);
  downward 133.958 -> 137.625 (+2.7%);
  upward 146.333 -> 140.417 (-4.0%).
- 256 bits: nearest 213.041 -> 206.459 (-3.1%);
  downward 162.458 -> 167.208 (+2.9%);
  upward 177.209 -> 170.334 (-3.9%).
- 1024 bits: nearest 313.625 -> 292.250 (-6.8%);
  downward 275.042 -> 266.708 (-3.0%);
  upward 276.750 -> 268.541 (-3.0%).

The 128-bit ratios for nearest/downward/upward improve to
2.754x/2.799x/2.888x, but all remain above the 2x objective.
The modest regression in overall aggregates involves timing
variation across operations not changed by this commit, so
it does not establish a causal global slowdown. Retain the
sign-plumbing optimization pending a repeated controlled run.
High-priority 20:28 ratios also include
`atan_nearest`@53 3.426x, `exp_nearest`@128 3.297x,
`pow7_nearest`@53 3.279x, `log_nearest`@53 3.252x,
and `cos_nearest`@53 3.236x.

**CI blocker:** Both `04e1b6f3` and `89d93735`
pass Unix and Coverage but FAIL Windows 2025 Release at
GitHub Actions **Build** step 13:
`04e1b6f3` Windows run 39066603660, job 117495485037;
`89d93735` Windows run 39066824427, job 117496125176.
The last confirmed fully green numerical commit was
`4ae35e8c`. These observations locate the first failure
to the fixed-point centered sin/cos change or its tests in
`04e1b6f3`, not the later division sign refactor.
The available connector does not expose the MSVC diagnostic.
Do NOT assert the exact root cause, suppress warnings,
or weaken the tests without that diagnostic.

**Potential next algorithmic direction (NOT merged):**
The current `asin_interval` evaluates
`sqrt(1-x*x)`, a rational interval division and
`atan`, and `acos` subtracts the result from
`pi/2`. For inputs very near +/-1/2, exact center
identities can bypass this chain:
`asin(1/2)=pi/6`, `acos(1/2)=pi/3`,
`acos(-1/2)=2*pi/3`, and the slope magnitude is
`D=2/sqrt(3)`. For `y=|x|=1/2+delta`
with `|delta|<=2^-ceil(work/2)` and work>=39,
the entire segment lies inside [1/4,3/4], where
`|asin''(y)|=y/(1-y*y)^(3/2)<=24/7<4`.
Thus `|R|<=2*delta*delta`, giving certified
`asin(y) in pi/6+D*delta +/-2*delta*delta`,
`acos(y) in pi/3-D*delta +/-2*delta*delta`,
and `acos(-y) in 2*pi/3+D*delta +/-2*delta*delta`.
For negative asin use oddness. Cached outward dyadic
intervals at work+24, sourced from existing certified pi
and sqrt3 routines at work+40, could make per-input
cost low. The general evaluator would remain the
fallback for all other values. A 500-digit independent
mpmath probe checked 1500 signed evaluations with only
numerical-evaluation tolerance; the analytic derivative
bound, not that sample, supplies the proof. A candidate
implementation and differential-test changes were prepared
but deliberately **NOT committed** because the parent
Windows Build remains broken. This avoids stacking an
independent feature over an unvalidated platform baseline.

**Next acceptance gate:** Obtain the Windows Build diagnostic
for run 39066824427 (and if necessary 39066603660),
classify the actual MSVC error, and restore a passing
Unix/Windows/Coverage matrix first. Only then introduce
the centered asin/acos experiment, with precision/boundary/
sign/directed-rounding differential tests and a controlled
new 1000x11 benchmark. The numerical target is unchanged:
worst AP/MP strictly <2x over all 108 pairs with
geometric mean approaching 1x; never trade away exact
rounding or historical test coverage to meet that target.

## MSVC C4459 fix for dyadic centered sin/cos (17 October 2026)

User-provided Windows 2025 Release build log for `89d93735`
(also failing on `04e1b6f3`) conclusively identifies:
`src/floatap.cpp(2173): warning C4459: declaration of 'upward'
hides global declaration`, referring to
`Ariadne::upward` in `include/numeric/rounding.hpp:280`.
Because the build uses `/W4 /WX`, the warning causes
`error C2220` and stops the Build step. This is a
**name-shadowing compiler failure**, not evidence of
incorrect numerical rounding.

Rename only the Boolean lambda parameter in
`centered_trigonometric_dyadic` from `upward` to
`round_toward_positive_infinity`, and the corresponding
conditional. All exact GMP arithmetic, outward-directed
integer rounding, selected endpoints, tests, and other
numerical algorithms are unchanged. Keep MSVC `/WX`
enabled. The new Windows CI, alongside Unix and Coverage,
must pass before all-platform acceptance. The latest
unloaded 1000x11 performance baseline remains the user
20:28 report: geomean 1.237x, total-time 1.282x,
worst 3.446x (acos@53), 34/108 cases >=2x.
No performance speedup is claimed for this identifier
rename. The future certified asin/acos center expansion
remains a separate unmerged design idea.

## Certified half-centered inverse trigonometric fast path (17 October 2026)

**Verified green starting point.** Commit `86901ead` passed Unix
(run 39067423601), Windows (run 39067423650), and Coverage
(run 39067423654). The MSVC C4459 issue from the local
`upward` lambda parameter is resolved without relaxing /WX.
The latest user-measured warmed-up, unloaded same-host 1000x11
benchmark is `Markdown incollato(20261017-202844).md` (after
`89d93735`, before this change):
geomean **1.237x**, weighted **1.282x**, worst
**3.446x (acos_nearest@53)**, and **34/108 >=2x**.
At 53 bits, FloatAP asin/acos take
3727.209/4431.208 ns versus MPFR 1339.375/1285.875 ns.
At 1024 bits they take about 40.7/42.1 us, despite the
fast centered atan kernel. The existing general
`asin_interval` computes a certified
`sqrt(1-x*x)`, divides two rational intervals, and
evaluates atan of one quotient endpoint; `acos`
also subtracts from certified `pi/2`. These operations
remain appropriate outside a small well-defined neighborhood.

**New mathematical reduction, applicable to any real input
in the neighborhood (not benchmark-index matching).** For
`y=|x|=1/2+d`, and
`|d|<=2^-ceil(work/2)`, work>=4,
the entire segment from 1/2 to y lies in [1/4,3/4].
On this segment
`|asin''(t)|=t/(1-t*t)^(3/2) <24/7<4`:
`t<=3/4`, `1-t*t>=7/16`, and
`sqrt(1-t*t)>1/2`. Thus Taylor's theorem gives the
**rigorous outward remainder** `|R|<=2*d*d`.
The exact center identities and derivative are
`asin(1/2)=pi/6`, `acos(1/2)=pi/3`,
`acos(-1/2)=2*pi/3` and
`D=asin'(1/2)=2/sqrt(3)`. Accordingly,
`asin(y) in pi/6 + D*d +/-2*d*d`,
`asin(-y)=-asin(y)`,
`acos(y) in pi/3 - D*d +/-2*d*d`, and
`acos(-y) in 2*pi/3 + D*d +/-2*d*d`.
All interval operations preserve exact outward bounds.
Outside the precision-dependent neighborhood, the old
general inverse-trigonometric logic remains unchanged.

**Implementation:** Compute certified `pi` at work+40
and `sqrt(3)` at work+40 once per thread/cache-precision
upgrade. Form exact intervals for `pi/6`, `pi/3`,
`2*pi/3` and `2/sqrt(3)` (positive denominator).
Quantize these four endpoint pairs **outward** to dyadics
at work+24 using the existing `outward_dyadic_interval`
helper, not ordinary rounded floats. A cached interval
computed at higher precision remains valid for lower
work. Per call, multiply the certified slope interval
by the exact signed offset, add/subtract the certified
center, then widen both endpoints by exactly
`2*d*d`. Only accept a final FloatAP result when
*both* outward endpoints round to the same result
in the selected mode; otherwise refine work as before.
No MPFR dependency is added to production code.
The general arcsine/arccosine paths, NaN/domain/pole
handling, signed zero and existing external APIs
are left unchanged.

**Tests and acceptance:** Extend the existing FloatAP/
FloatMP differential suite (do not remove other cases)
at precisions 7,53,128,256,1024, and revisit 53
after a higher-precision cache update. For each
precision check both signs, exact +/-1/2 and several
ULP-scale offsets, plus out-of-neighborhood points
1/4,3/4,1/3,2/3. Check both asin and acos, nearest,
downward, upward and toward-zero using the *same rounded
FloatAP input* as the MPFR oracle. All 3 CI workflows
must succeed, including MSVC /WX on Windows.
A new warmed-up, unloaded same-host 1000x11 benchmark
is needed before reporting *any* new speedup. The
global target remains worst strict <2x for all
108 cases and geomean approaching 1x.

## Benchmark 20:43 and next general certified-center candidates (17 October 2026)

**User-supplied benchmark:** `Markdown incollato(20261017-204344).md`
is the 1000-iteration/11-sample warm benchmark on the same host,
following `ecc13dcc`. The matching pre-`ecc13dcc`
baseline is `Markdown incollato(20261017-202844).md`.
Each run measures 27 operations at 53,128,256,1024 bits,
108 FloatAP/FloatMP pairs.

| Statistic | 20:28 baseline | 20:43 with centered asin/acos |
| --- | ---: | ---: |
| Geomean AP/MP | 1.237x | **1.090x** |
| Time-weighted AP/MP | 1.282x | **0.591x** |
| Worst ratio | 3.446x (acos@53) | **3.335x (exp@128)** |
| Cases >=2x | 34/108 | **25/108** |
| Cases >=3x | 6/108 | **6/108** |
| Geomean at 1024 bits | 0.992x | **0.809x** |

**Measured local inverse-trig gains** in absolute FloatAP runtime:

| Bits | asin before -> after (ns) | acos before -> after (ns) |
| --- | --- | --- |
| 53 | 3727.209 -> 1658.833 | 4431.208 -> 1667.458 |
| 128 | 5080.542 -> 2081.750 | 5844.000 -> 1984.834 |
| 256 | 7800.833 -> 2071.458 | 8662.667 -> 2071.833 |
| 1024 | 40680.791 -> 3453.292 | 42136.875 -> 3451.417 |

The single largest relative improvements (approx. 12x) were
at 1024 bits. Crucially, this benchmark pools 32 values
near 1/2; it demonstrates strong speed on that input
distribution, **not** an all-domain speedup. All certified
local paths must maintain mathematically proven bounds,
exact rounding, and the previous general fallback.

**Outstanding AP/MP slowdowns >=3x (2026-10-17 20:43):**
`exp_nearest`@128 3.335x;
`pow7_nearest`@53 3.332x;
`cos_nearest`@53 3.229x;
`atan_nearest`@53 3.151x;
`log_nearest`@53 3.050x;
`div_up`@128 3.011x.
These set sensible engineering priorities beyond global
geomean, since the hard objective remains every ratio <2x.

**Next certified-center candidate 1: exp near +/-1/2.**
For `x=c+delta`, `c=+/-1/2`,
`|delta|<=2^-ceil(work/2)`, `work>=4`,
the entire segment lies inside [-3/4,3/4].
`exp(x) = E_c*(1+delta)+R`,
`E_c=exp(c)`. Since `exp''(t)=exp(t)<=exp(3/4)<3`,
Taylor gives the conservative rigorously outward remainder
`|R|<=2*delta^2`. Compute the certified center
`E_+=exp_interval(1/2,work+40)` using the unchanged
general algorithm, and `E_-=1/E_+` using monotone
endpoint reversal. Cache outward dyadic enclosures
at `work+24` per thread and precision upgrade.
For each call, `1+delta>=3/4>0` permits just two
interval-endpoint products plus the remainder. Keep the
current range-reduced Taylor/squaring path for all other
inputs. This targets the current worst at 128 bits and
exp outliers at 256 and 1024 bits without weakening
the rounded-endpoint agreement rule.

**Candidate 2: log near +1/2.** The existing `log2_interval`
already caches an outward interval for ln2.
For `x=1/2+delta` in [1/4,3/4],
`log(x)=-ln2+2*delta+R`. As
`|log''(t)|=1/t^2<=16`, the certified remainder is
`|R|<=8*delta^2`. Reuse the cached certified ln2
(or an outward dyadic endpoint cache) and return directly
without reconstructing `z=(x-1)/(x+1)` or evaluating
the general log series. Keep all the old cases outside
the local neighborhood and preserve the existing
NaN/negative/zero/one/infinity behavior. This directly
targets `log_nearest`@53 (3.05x).

**Candidate 3: atan near +/-1/2.** Reuse the already
certified cached `atan(1/2)` interval. The derivative
`atan'(1/2)=4/5` is EXACT rational, with no second
nontrivial constant to compute.
`atan(1/2+delta)=atan(1/2)+(4/5)*delta+R`.
On [1/4,3/4],
`|atan''(t)|=2t/(1+t*t)^2<=3/2`,
so `|R|<=delta^2` is a certified bound.
Apply oddness on negative arguments and retain the
existing Möbius argument reduction elsewhere.
This could remove fraction reduction and a tiny
Taylor kernel from `atan_nearest`@53 (3.15x).

**Other avenues:**
- The new integer/dyadic centered sin/cos path is only
  enabled at work>=128. A controlled small-precision
  crossover experiment may assess activating it at
  the initial 53-bit work=85, where cos@53 is 3.23x.
  The old *fraction-based* small-precision centered
  path regressed, so activation must not be assumed
  faster without measured isolated comparisons.
- `pow(mode,x,Nat)` still converts the finite dyadic
  operand to a generic Rational, exponentiates both its
  numerator and denominator, then divides and rounds.
  Instead directly compute `mpz_pow_ui` on the exact
  nonnegative significand, multiply the binary exponent
  by the integer degree with checked overflow, retain
  the odd-power sign bit, and call the existing
  `round_scaled_integer` ONCE. This removes unnecessary
  large denominator/exact-rational work without
  intermediate rounding. Carefully test exponent range,
  signed zero/infinity/NaN, degree 0/1 and sign, every
  rounding mode, and values near round-to-even ties.
  This is a complementary **architectural** optimization,
  not a Taylor-center trick. It targets pow7@53=3.332x.
- Division@128, compare@1024 and exp outside center
  require independent profiling and input-distribution
  testing; avoid replacing rigorous one-round behavior
  with successive rounded operations.

**Verification and release gate:** A standalone 4000-bit
mpmath check validated the local exp/log Taylor remainder
inequalities on 105 cases at work=39,85,160,288,1056,
both exp center signs, and positive log offsets.
The *proof* is the global second-derivative bound, not
the numerical probe. An implementation of exp/log
centered paths with a shared exact GMP offset predicate,
outward dyadic caches, and matching MPFR differential
tests for signs, local inclusion thresholds, 7/53/128/
256/1024 bits and four modes has been prepared as
UNMERGED candidate source/test blobs. Do not claim any
benchmark benefit or CI success for this candidate.

**Blocking issue:** `ecc13dcc` passes GitHub Unix
(run 39067889123) and Coverage (run 39067889148),
but its Windows 2025 [Release] run 39067889104
(job 117499230326) FAILED at **Build step 13**.
The GitHub connector cannot retrieve the detailed MSVC
compiler log, so no exact Windows root cause is established.
Do not stack the independent exp/log candidate on top of
an unresolved Windows build failure or suppress /WX.
Obtain the Windows compiler diagnostic, fix the actual
cause, verify all three workflows, then introduce
the candidates in separately measurable commits.
Separately, Codecov patch check on ecc13dcc reported
79.60% diff coverage versus a 99.20% threshold; this
is not the same as the successful Coverage workflow,
and tests must be extended semantically rather than
adjusted just to make a metric pass.

## Windows MSVC C4267 fix in inverse-trig regression tests (17 October 2026)

**Confirmed compiler diagnostic** supplied for GitHub Actions
Windows 2025 Release run 39067889104 (commit `ecc13dcc`):
`tests/test_float_rounding.cpp(603)`, error C2220,
warning C4267: initializing `const Ariadne::Nat`
from `size_t` may lose data. The build uses MSVC /W4 /WX,
so this warning stops compilation. This identifies the
concrete Windows failure that was previously documented
without accessible compiler logs. It is unrelated to
FloatAP arithmetic or the certified inverse-trig bounds.

**Minimal portable test fix:** Explicitly convert the
known small, fixed precision values with
`const Nat offset_power(static_cast<Nat>(ap.bits()));`.
The regression matrix is unchanged: 7,53,128,256,1024
bits, plus a repeat at 53 to exercise cache reuse.
The power is fed into the same exact Integer(2)^Nat
ULP construction; no existing test is removed or weakened.
A source search found only this `const Nat` initialization
directly from `ap.bits()` in the test source.

**Validation:** The previously merged `ecc13dcc` passed
Unix and Coverage, and failed Windows only during Build.
A fresh full Unix/Windows/Coverage matrix must be
verified for this commit before considering the
centered asin/acos optimization green on all platforms.
Codecov patch coverage is a separate status and must
be assessed via semantic test coverage, not metric hacks.
The performance baseline remains the user's 20:43
1000-iteration/11-sample benchmark: geomean 1.090x,
time-weighted 0.591x and worst 3.335x (exp@128).
Do not infer a new speedup from a test-only type correction.
The independently planned exp/log certified-center work
remains unmerged until the new CI is validated.

## Certified centered exponential near +/-1/2 (17 October 2026)

**Validated starting commit:** `f75a4548`, which passed
GitHub Unix run `39068524458`, Windows run `39068524529`,
and Coverage run `39068524358`. The MSVC C4267 test-only
narrowing issue is resolved without weakening `/WX`.
The most recent measured benchmark is
`Markdown incollato(20261017-204344).md`, 1000 iterations,
11 samples on the same unloaded host:
AP/MP geometric mean **1.090x**, time-weighted **0.591x**,
worst **3.335x** (`exp_nearest` at 128 bits), and
**25/108** operation/precision ratios >=2x. In that report
FloatAP `exp_nearest` absolute medians at
53/128/256/1024 bits are 1979.417, 2872.625, 3678.167,
12275.709 ns respectively; matching FloatMP medians are
1063.625, 861.416, 1387.875, 6127.000 ns.
This is the **pre-optimization** benchmark, not an outcome.

**Certification proof:** For `x=c+delta` with
`c` exactly +1/2 or -1/2, take the local route only when
`work>=4` and `|delta|<=2^-ceil(work/2)`.
Then `|delta|<=1/4`, so the entire segment from c
to x lies in [-3/4,3/4]. Here
`exp''(t)=exp(t)<=exp(3/4)<e<3`.
Taylor's theorem therefore gives
`exp(c+delta)=E_c*(1+delta)+R` with
`|R|<=2*delta*delta`; this uses the exact center
derivative `exp'(c)=E_c`.
Because `1+delta>=3/4>0`, multiplying a certified
interval of E_c by 1+delta preserves endpoint order.
Then subtract/add `2*delta^2` from/to the result.
The condition is checked with exact GMP integer
comparison, with overflow-safe `mp_bitcnt_t` handling.
Inputs outside this mathematically defined window
retain the exact prior `exp_interval` algorithm;
it does not key on benchmark array position.

**Certified caching:** One thread-local cache per work
precision maximum holds BOTH `E_+=exp(1/2)` and
`E_-=exp(-1/2)` intervals. At a cache miss/precision
increase, evaluate existing certified
`exp_interval(1/2,work+40)`; the negative center is
the exact outward reciprocal interval
`[1/E_+.upper,1/E_+.lower]`.
Convert both intervals to dyadic endpoints at
`work+24` using exact outward GMP floor/ceiling,
including negative-value handling in the helper.
Cache enclosures, **never caller-rounded constants**:
all four rounding modes share a thread-safe local
cache. If work decreases, use the tighter existing
bounds. The enclosing `certified` refinement loop
still requires both endpoints to round to the same
target result, otherwise it increases precision.
No approximation or different rounding criterion
is introduced.

**Regression matrix:** Add FloatAP/FloatMP differential
tests for exact shared input values around both
+1/2 and -1/2 at bits 7,53,128,256,1024, with
53 revisited after a high-precision cache upgrade.
Exercise +/- one ULP, +/- four ULPs, the exact
`2^-ceil(work/2)` window boundary and +/- twice
the boundary, and fallback values 1/4, 3/4, 1/3,
2/3, 1, 2. Preserve tests for all four rounding
modes and all earlier scenarios.
No previous test assertions are removed.

**Measurement and acceptance:** This commit isolates
exponential optimization only: `log`, `atan`,
`pow` and trigonometric kernels are unchanged.
Do not claim any speedup until Unix/Windows/Coverage
CI pass and a same-machine, CPU-idle, warmed-up
1000/11 benchmark reports all 108 cases.
Track both AP and MP absolute medians, overall
geomean, time-weighted ratio, worst pair and >=2x.
If newly expensive cache startup dominates
short-lived workloads, investigate that separately
without tuning solely to the 32-value benchmark pool.
The separately proven local log and atan reductions
and integer-power refactor remain future *independent*
experiments, not part of this commit.

## Benchmark 20:58 and exact integer-significand pow optimization (17 October 2026)

**Latest controlled baseline:** `Markdown incollato(20261017-205819).md`
was measured on the same host with the established 1000 iterations,
11 median samples, 27 operations and four precisions (53,128,256,
1024 bits), immediately after commit `1d5adf8c`.
Compare with the preceding valid benchmark
`Markdown incollato(20261017-204344).md` after `ecc13dcc`.
The 19:53 benchmark with simultaneous coverage remains excluded.

| Aggregate | 20:43 before exp | 20:58 after exp |
| --- | ---: | ---: |
| Geomean FloatAP/FloatMP | 1.090x | **1.063x** |
| Total-time ratio | 0.591x | **0.492x** |
| Worst slowdown | 3.335x (exp@128) | **3.302x (pow7@53)** |
| Geomean, 1024 bits | 0.809x | **0.765x** |

**Measured FloatAP exponential gains:**
53 bits: 1979.417->1101.333 ns (-44.4%),
128 bits: 2872.625->1387.375 ns (-51.7%),
256 bits: 3678.167->1478.625 ns (-59.8%),
1024 bits: 12275.709->2591.083 ns (-78.9%).
The centered exponential therefore achieved the predicted
strong improvement *on this near-1/2 workload*.
At 1024 bits, exp FloatAP/FloatMP=2591.083/5563.042
=0.466x. The worst-case now occurs in integer power:
pow7@53 262.916/79.625=3.302x; the next high ratios
include div_up@128 132.875/47.375=2.805x,
div_nearest@128 170.459/63.291=2.693x and
sqrt@53 144.916/49.125=2.950x.
The user hard objective all 108 ratios <2x is NOT met.
Keep measuring all operations, not only the improved exp.

**Starting CI:** GitHub Unix run `39068877587`, Windows run
`39068876991`, and Coverage run `39068877551` all
completed successfully on `1d5adf8c`. The 20:58 report is
the last measured baseline before the present power change.

**New architectural power path:** A FloatAP finite nonzero value
is represented exactly as
`(-1)^negative * significand * 2^exponent` with nonnegative
GMP integer significand. Thus an integer power n>=0 has exact
magnitude `significand^n * 2^(exponent*n)` and sign
`negative && (n odd)`. There is no need to build a generic
GMP `mpq_t`, canonicalize a huge dyadic denominator, compute
numerator and denominator powers separately, and pass their
quotient to `from_ratio`. Instead, take a single
`mpz_pow_ui` of the stored significand; use exponentiation
by doubling with `checked_add` to calculate exponent*n
without signed overflow; then call the existing
`round_scaled_integer` once. This has the same exact
mathematical input to the single rounding step as the old
rational division. For n=1 return the existing FloatAP
value unchanged. Keep the entire original rational path
as a fallback only if the intermediate exponent product
does not fit Exponent: its normalization can exceptionally
cancel an enormous negative intermediate exponent,
so changing that case to an overflow exception would be
a semantic regression. Preserve the previous ordering
of n=0, NaN, infinity and signed-zero checks.

**Tests:** Extend FloatAP/FloatMP differential checks
over 7/53/128/256/1024 bits; both signs; diverse exact
fractional and integer inputs; powers
0,1,2,3,7,8,13,31; nearest and all three directed
rounding modes; all MPFR oracle values derive from
the already rounded binary FloatAP operand.
Add separate explicit signed-zero, signed-infinity,
NaN and exponent-zero assertions. Do not modify any
previous test. A preexisting round-to-even test covers
integer powers indirectly; new cases increase semantic
coverage for exact vs inexact results and directed
rounding on both signs.

**Next validation:** The change is algorithmic and
has not yet been benchmarked or passed CI. Require Unix,
Windows and Coverage success (and review the distinct
Codecov patch status), then another unloaded same-machine
1000x11 report comparing all 108 ratios, absolute
FloatAP and MPFR medians and the new global worst.
The power optimization is general to any finite input,
not a shortcut keyed to the benchmark's 32 near-half
arguments. If a latency improvement comes with changed
rounding behavior, treat that as a bug rather than
removing or tailoring the failing test. Other planned
log/atan-centered optimizations remain independent.

## Benchmark 21:01 and certified half-centered atan (17 October 2026)

**New measured 1000x11 baseline:** User-supplied
`Markdown incollato(20261017-210150).md` is the unloaded warmed-up
same-host report after commit `aa2b61c8`, versus
`Markdown incollato(20261017-205819).md` before it.
All 108 FloatAP/FloatMP operation/precision pairs are present
(27 operations at 53,128,256,1024 bits).

| Statistic | 20:58 pre-pow | 21:01 post-pow |
| --- | ---: | ---: |
| Geometric mean AP/MP | 1.063x | **1.024x** |
| Time-weighted total ratio | 0.492x | **0.488x** |
| Worst slowdown | 3.302x (pow7@53) | **3.257x (atan@53)** |
| 53-bit geometric mean | 1.046x | **0.972x** |
| Cases >=2x | count not recorded here | **27/108** |

**Direct measured power latency (FloatAP ns):**

| Precision | 20:58 | 21:01 | Reduction | 21:01 AP/MP |
| --- | ---: | ---: | ---: | ---: |
| 53 | 262.916 | 60.250 | -77.1% | 0.735x |
| 128 | 392.042 | 97.292 | -75.2% | 0.421x |
| 256 | 476.750 | 152.917 | -67.9% | 0.341x |
| 1024 | 1912.417 | 1408.000 | -26.4% | 0.385x |

This verifies substantial speedups on the tested
`pow7` workload, without changing the single-round
contract. Do not generalize the timings to all exponents.
The measured worst cases are now `atan_nearest`@53
1772.750/544.333=3.257x, `log_nearest`@53
2088.208/643.833=3.243x, `cos_nearest`@53
1470.292/454.417=3.236x, `div_up`@128
146.291/46.041=3.177x, and `sqrt_nearest`@53
146.459/49.500=2.959x. The 27 ratios >=2x
split 4/12/8/3 across 53/128/256/1024 bits.
The global acceptance goal worst strictly <2x
is NOT yet met. Unix and Coverage workflows for
`aa2b61c8` passed; Windows was still in progress
when this follow-up was prepared.

**Certified local atan improvement (independent change).**
For `y=|x|=1/2+delta` with
`|delta|<=2^-ceil(work/2)` and `work>=4`,
the entire Taylor segment lies in [1/4,3/4],
where `atan'(1/2)=4/5` exactly and
`|atan''(t)|=2t/(1+t*t)^2<=3/2`.
Taylor's theorem gives the rigorous
`atan(y) in atan(1/2)+(4/5)*delta+[-delta^2,+delta^2]`.
For negative inputs, use oddness and reverse bounds.
This works for *any* input in a mathematical neighborhood,
not just the benchmark's sequence, and the old
argument reduction/Taylor and infinity paths remain
the fallbacks elsewhere.

**Interval cache and rounding:** Reuse the existing
`atan_half_interval` certified core at `work+40`,
then quantize its endpoints OUTWARD to a dyadic scale
`work+24`. Cache those enclosures thread-locally,
upgrading work as needed. The per-call linear term
`(4/5)*delta` and `delta^2` remainder are computed
as exact GMP rationals. The `certified` loop still
accepts a result only if BOTH bounds round to the
same value in the selected rounding mode; otherwise
it increases work. This is an independent first
implementation with a likely strong speedup on
near-half inputs but no measured benefit yet.
Check if exact rational operations dominate at
53 bits before attempting an integer-only variant.

**Regression:** Append MPFR differential comparisons
at 7,53,128,256,1024 bits and revisit 53 after
a cache precision upgrade. Cover both signs, exact
1/2, +/-1 or 4 input ulps, exact initial-work
threshold and +/-2x threshold (the old general
fallback), and wider inputs 1/4,3/4,5/7,1/3,
2/3,1,2. Compare the identical rounded FloatAP
input and all four modes, preserving every old test.
No failure scenarios are removed or altered to
make CI pass.

**Acceptance:** Require Unix/Windows/Coverage success
on the new commit plus an unloaded warmed-up
1000-iteration/11-sample benchmark on the same
machine. Report all 108 ratios, FloatAP absolute
timings, geomean, weighted ratio and worst ratio;
revert/adjust if the small-precision rational
kernel increases latency. Future isolated
opportunities include cached certified log around
1/2, smaller-work fixed-point sin/cos, and
128-bit division/compare/add paths. A green
Coverage workflow is NOT the same as the
Codecov patch coverage check, which should be
reviewed separately.

## Independent FloatAP / FloatDP binary64 benchmark (17 October 2026)

Add a separate benchmark executable
`benchmark_compare_double_precision_backends` comparing FloatDP
(the existing hardware-double backend) against FloatAP at exactly
53 significand bits. The existing MPFR/FloatAP benchmark remains
`benchmark_compare_multiple_precision_backends`.
Both executables use the SAME templated 27-operation pool:
construct_rational, copy, compare; nearest and directed add,
multiply and divide; nearest subtract, FMA, sqrt, exp, log,
sin, cos, tan, atan, asin, acos, pi, pow7 and both conversions.
The deterministic 32-input exact-Rational pool, same
preconstructed input strategy, median-per-operation code and
one unmeasured warm-up sample are now in
`benchmarks/benchmark_backend_workload.hpp`. Changes to
these workloads therefore affect both comparators together
rather than subtly diverging.

FloatDP has fixed binary64 precision: reject any
`--bits` value other than 53 in the new benchmark;
FloatMP/FloatAP retain the existing 53/128/256/1024
precision options. Both executables now DEFAULT to
`--iterations 1000 --samples 11` (previous comparison
defaults were 300/7); explicit flags and the legacy
positional counts continue to work. Output labels,
matched-operation slowdown summaries and
`--summary-only`, `--details`, `--max-slowdown`
apply to FloatAP/FloatDP in the new executable and to
FloatAP/FloatMP in the original.

Build:
```sh
cmake --build build --target benchmark_compare_double_precision_backends benchmark_compare_multiple_precision_backends
./build/benchmarks/benchmark_compare_double_precision_backends --summary-only --details
./build/benchmarks/benchmark_compare_multiple_precision_backends --summary-only
```
Use the real local build directory as configured. Run Release
on an unloaded machine; timing results are not generated by
adding source code, and no speculative slowdown ratios are
claimed. Compare identical binary64 arithmetic semantics with
care: FloatDP is backed by hardware double and library
transcendentals, and the purpose here is speed on the same
battery, NOT to infer correctly rounded transcendental
equivalence from performance alone. Preserve existing tests
and require Windows/Unix/Coverage CI. The preexisting FloatAP
design/performance journal remains the source of truth.

## FloatDP benchmark: omit redundant ALL row (17 October 2026)

The FloatAP/FloatDP benchmark is limited to binary64/53-bit precision, so its `53 bits` summary already represents every one of the 27 operation ratios. Remove the redundant `ALL` summary row **only** from `benchmark_compare_double_precision_backends.cpp`. Retain the aggregated slowdown accumulator for `--max-slowdown` validation and keep all 27 timing rows, CLI flags and the 1000-iteration/11-sample defaults unchanged. The separate FloatAP/FloatMP benchmark still measures 53/128/256/1024 bits and continues to print its useful cross-precision `ALL` row. This is an output-only change; no numerical backend, operation or measurement method is modified.

## Benchmark 21:11 and certified half-centered log optimization (17 October 2026)

**User report:** `Markdown incollato(20261017-211153).md`
is a 1000-iteration/11-sample full benchmark of the existing
27 FloatAP/FloatMP operations at 53/128/256/1024 bits.
The previous comparable 21:01 report is
`Markdown incollato(20261017-210150).md`.

| Metric | 21:01 before centered atan | 21:11 after centered atan |
| --- | ---: | ---: |
| Geometric mean, 108 FloatAP/FloatMP pairs | 1.024x | **1.014x** |
| Total-time AP/MP ratio | 0.488x | **0.483x** |
| Worst slowdown | 3.257x (atan@53) | **3.181x (log@53)** |
| Pairs >=2x | 27/108 | **20/108** |

The 20 pairs >=2x split across precisions as 4 at 53,
10 at 128, 5 at 256 and 1 at 1024 bits. The 53-bit
`atan_nearest` median decreased 1772.750->1353.792 ns
(~23.6% faster) and its slowdown fell to 2.505x.
The present worst is `log_nearest`@53:
FloatAP **2109.875 ns**, FloatMP **663.375 ns**,
ratio **3.181x**. Other 53-bit outliers are
`cos_nearest` 3.073x and `sqrt_nearest` 2.836x.
The worst at 128 is `div_up` 3.164x, at 256
`div_nearest` 2.220x, and at 1024
`compare` 2.190x. The time-weighted 0.483x
ratio is for the synthetic benchmark operation
distribution, not a general application-wide claim.

**Release status:** Commit `3a4993da` introducing the
centered atan algorithm was green on Unix (run
39069436882), Windows (39069436958), Coverage
(39069436892). Commit `1513659f` adding the separate
FloatAP/FloatDP benchmark was likewise green on all
three workflows (Unix 39069736274, Windows
39069736273, Coverage 39069736336).
At the start of this log work, the next output-only
commit `f6b19da7` was still executing all 3 CI
workflows. No output or correctness claim about it is
made beyond the source changes.

**New certified logarithm reduction:** For positive
`x=1/2+delta`, require `work>=4` and
`|delta|<=2^-ceil(work/2)` using the already proven
exact GMP predicate `exp_inside_half_window`.
Then x and the entire segment from 1/2 to x lie in
[1/4,3/4]. Taylor's theorem with
`log'(1/2)=2` (EXACT) and
`|log''(t)|=1/t^2<=16` proves
`log(x)=-log(2)+2*delta+R`,
`|R|<=8*delta*delta`.
The constant -log(2) comes from the already existing
certified thread-local cache `log2_interval(work+40)`,
which computes an interval using the proven
`log_series(1/3)` kernel. Negate and reverse its
endpoints and round them OUTWARD to dyadics at
`work+24` using `exp_outward_dyadic` (which
performs floor for the lower endpoint and ceiling
for the upper endpoint, including negatives).
Cache the result thread-locally, upgrading whenever
work increases. Computing the per-call linear term
and the exact symmetric Taylor remainder requires
no series, argument scaling, or rational
Mobius reduction. On arguments outside this exact
neighborhood, preserve `log_interval` unchanged.
No other operation (including exp, atan, pow, sin,
cos or the comparison benchmark) changes.

**Rounding and regression:** The main certified loop
still compares both outward bounds in the requested
rounding mode and refines work whenever needed.
No intermediate result is rounded to a FloatAP.
Keep the existing NaN, sign/domain, signed-zero,
infinity and exact-log(1)=0 branches unchanged.
Add FloatAP/FloatMP differential checks around the
exact center, +/-1/4 ULP patterns, both sides of the
work-dependent inclusion cutoff, and fallback
1/4,3/4,1/3,2/3,1,2 at precision
7/53/128/256/1024 bits. All nearest/down/up/
toward-zero modes are retained and 53-bit cache
reuse is checked after a 1024-bit upgrade.
No existing test is removed or weakened.

**Acceptance:** The mathematical remainder proof
guarantees enclosure; the differential tests and
new GitHub Unix/Windows/Coverage CI must additionally
validate implementation portability and historical
correctness. DO NOT claim any runtime benefit until
another same-host idle-CPU warmed-up 1000x11 report
compares all 108 pairs with the 21:11 baseline.
The global target remains all 108 ratios <2x and a
geomean approaching 1x. Re-evaluate the center
if its 53-bit cost does not beat the general
log path, especially outside the benchmark pool.
## Benchmark 21:31 and MPFR/GMP-independent optimization priorities (17 October 2026)

**Baseline and provenance.** The user-supplied
`Markdown incollato(20261017-213109).md` contains the subsequent
same-machine FloatAP/FloatMP 1000-iteration/11-sample medians for
27 operations at 53, 128, 256 and 1024 significand bits.
It follows `20187a89` (certified log around 1/2); the immediately
preceding 21:11 benchmark is recorded in the previous section.
The 21:31 times are external measured results, not a benchmark rerun
performed by this document-only change.

| Bits | Geomean AP/MP | Total time AP/MP | Worst AP/MP | Pairs >=2x |
| --- | ---: | ---: | --- | ---: |
| 53 | 0.961x | 1.481x | cos_nearest 3.163x | 3/27 |
| 128 | 1.332x | 1.120x | div_nearest 2.902x | 10/27 |
| 256 | 1.035x | 0.579x | div_nearest 2.275x | 6/27 |
| 1024 | 0.759x | 0.250x | compare 2.104x | 1/27 |
| ALL | **1.001x** | **0.451x** | **cos_nearest@53 3.163x** | **20/108** |

The measured sum of per-operation medians (AP vs MPFR, ns):
53: 12229.165 vs 8256.917; 128: 14633.414 vs 13068.625;
256: 15394.955 vs 26570.290; 1024: 25717.373 vs 102807.790.
The preceding 21:11 global geomean/total/worst were 1.014x,
0.483x and 3.181x (log@53); this run improves them to 1.001x,
0.451x and 3.163x (cos@53). Ratios >=2x remain 20/108,
but their distribution moves from 4/10/5/1 to 3/10/6/1.
In particular, log@53 drops from 2109.875 ns at 21:11 to
1102.541 ns (approximately 47.7% less AP time) and from
3.181x to 1.675x. This is consistent with the new
half-centered log dispatch, not isolated causal proof for
all operation changes between benchmark runs.

**New acceptance contract.** The objective is simultaneously
(1) every individual AP/MPFR median ratio strictly below 2x,
(2) geometric mean strictly below 1x **at each precision**,
and (3) ratio of the sum of the 27 AP medians to the sum of
the same MPFR medians strictly below 1x **at each precision**.
The global aggregates are additional diagnostics, not
substitutes for the per-precision constraints. This release
still fails (1) at all precisions, (2) at 128 and 256 bits,
and (3) at 53 and 128 bits. No performance improvement may
sacrifice correctly rounded results in nearest-even, downward,
upward or toward-zero mode.

**Performance budgets and impact of previous work.**
At 53 bits the net excess is 3972.248 ns; only 850.79 ns
would be saved by capping the current >=2x outliers at
exactly 2x, so reaching the total-time goal requires another
3121.46 ns of savings beyond that (slightly more for strict
inequalities). The six sin/cos/tan/atan/asin/acos cases alone
contribute a net +3899.3 ns; exp+log another +568.3 ns,
offset by faster nontranscendental operations. At 128 bits,
the net excess is 1564.789 ns; correcting only the >=2x
ratios to 2x saves about 279.91 ns, still leaving 1284.88 ns.
The 128-bit excess is distributed across trigonometry
(+570.2 ns), exp/log (+761.4 ns) and basic arithmetic
including sqrt (+578.7 ns), offset by faster conversions
and other operations. Thus improving only the worst-ratio
division will not fix the 128-bit total-time objective.

**Historically high-leverage commits, with attribution limits.**
The comparable new-host cycle started at a 300x7 baseline
of geomean 2.160x, total 2.771x, worst 7.367x and 70/108
pairs >=2x (documented at 18:57). This is evidence of the
direction of progress, NOT a controlled 300x7 to 1000x11
speedup measurement. The direct-result rounding series
(`5bbfd92e`, `c2a0db3d`, `aa2b61c8`) removes intermediate
copies/denominator work; the 19:03 benchmark showed reductions
in add/sub/mul, though that comparison was not perfectly
commit-isolated. The fixed-point elementary kernels and
cached constants (16 October series of commits) chiefly
address the former multi-microsecond transcendental bottleneck.
The certified centered asin/acos change (`ecc13dcc`) has
stronger direct before/after evidence at 20:28 versus 20:43:
asin@1024 40680.791 -> 3453.292 ns and acos@1024
42136.875 -> 3451.417 ns, on this near-half workload.
The direct integer power change (`aa2b61c8`) reduced
pow7@53 262.916 -> 60.250 ns and pow7@128
392.042 -> 97.292 ns in the recorded 20:58/21:01
measurements. The centered atan change (`3a4993da`)
reduced atan@53 1772.750 -> 1353.792 ns in the
21:01/21:11 reports. The subsequent centered log change
(`20187a89`) has the measured reduction above.
These are workload-local observations; do not generalize
center optimizations to arbitrary arguments without
an independent corpus.

**Proposed optimization sequence (NOT implemented here).**

1. **P0 — remove trivial comparison overhead.** In
   `cmp(FloatAP,FloatAP)` the input significands are
   normalized to exactly their declared bit precision.
   After special values, signed zeros and signs have been
   checked, if the precisions agree but the binary exponents
   differ, the exponent order alone determines positive
   magnitude order. This allows an allocation-free return
   before the current `mpz_sizeinbase`/shift machinery.
   Keep the generic differing-precision and equal-exponent
   paths. Test signs, equal/different precision, extreme
   exponent gaps and exact equal magnitudes. Current
   compare ratios are 2.140x@128, 2.245x@256 and
   2.104x@1024, albeit only ~2.7-3.0 ns excess each.
   This is an inexpensive worst-case-ratio improvement,
   not an important total-time fix. It makes less
   reliance on GMP, rather than more.

2. **P1 — remove exact-square threshold construction in sqrt.**
   For `mpz_sqrtrem(root,remainder,scaled)`, write
   `scaled=root^2+remainder` with integer
   `0<=remainder<=2*root`. In nearest-even mode, the
   existing condition
   `4*scaled > (2*root+1)^2` is exactly equivalent
   to `remainder > root` because the remainder and
   root are integers and equality at the half-integer
   threshold is impossible. The current code constructs
   `threshold`, doubles and increments the root,
   squares it, shifts `scaled`, and compares; replace
   just that branch with the integer remainder/root
   comparison, without touching directed rounding or
   normalization. Verify the identity and MPFR differential
   behavior at 7/53/128/256/1024 bits, rounding modes,
   exact squares, small/large exponents and carry.
   The four-limb scratch allocation and the integer
   square itself still need profiling. Present sqrt ratios:
   2.828x@53, 1.948x@128, 2.039x@256 and 1.738x@1024.
   This is a simple, generic speed opportunity before
   considering a replacement integer-root algorithm.

3. **P2 — division on one/two machine limbs.** The
   `FloatAP::from_ratio` non-unit path is still
   generic GMP: four initialized scratch integers,
   normalization shifts, `mpz_tdiv_qr`, remainder
   comparisons and result finalization. Measured division
   at 128 is 176.791/60.917=2.902x nearest,
   126.584/44.083=2.871x downward, and
   128.458/44.625=2.879x upward; all three at
   256 bits also exceed 2x. Design a portable fixed
   one-/two-limb quotient/remainder kernel with exact
   guard, remainder, ties-to-even, signed directed rounding
   and exponent carry. A 53-bit division may need up to a
   128-bit intermediate; 128-bit division may need a
   four-limb numerator. Do not require GCC/Clang
   `__int128` on MSVC: use portable limb algorithms or
   conditional accelerated primitives plus portable tests.
   Keep the large-precision fallback while designing the
   limb API to outlive GMP. Avoid another GMP-only
   micro-optimization as the primary architectural fix.

4. **P3 — shared small-limb storage for add/sub/FMA and
   construction.** At 128 bits, add nearest/down/up are
   2.209x/2.321x/2.475x, sub 2.239x, FMA 2.106x;
   at 256 bits FMA is 2.004x. `add_sub` still creates
   a general `mpz_t` result and performs full alignment,
   while `fma` produces a scratch `product`, copies it
   into a temporary exact FloatAP, then calls `add`.
   Replace per-call heap allocation for <=1/2/4 limbs
   with inline significand and stack scratch capacity.
   Construct a single exact product-plus-addend accumulator
   and round once; DO NOT compute separately rounded
   multiply followed by add. Refactor behind reusable
   unsigned-integer limb helpers (shifts, add/sub,
   mul, divrem, sticky/guard, bit length), not new
   GMP-specific APIs. Preserve mixed precision, cancellation,
   signed zero, huge exponent gaps and FMA semantics.
   An eventual full GMP removal should reuse this work.

5. **P4 — reduce rational/interval overhead in the
   certified transcendental hot paths.** The 32 inputs
   in `benchmark_backend_workload.hpp` are specifically
   near 1/2 (operand a) and 1/3 (operand b), with
   almost all transcendental tests operating on a.
   They naturally exercise the exact half-centered
   caches and do not represent arbitrary-domain costs.
   The 53-bit trig entry work is 85 bits
   (`bits+32`), whereas `sine_cosine_intervals`
   currently enables the dyadic center optimization
   only for `work>=128` and single-output calls.
   Therefore cos@53 (1441.791 vs 455.833 ns, 3.163x)
   and sin@53 (1339.541 vs 674.125 ns, 1.987x)
   still use the general Taylor kernel initially.
   Experiment with a measured crossover at 53 bits,
   but DO NOT reenable the earlier slower fractional
   centered path blindly. At 128 bits the center path
   already exists, yet cos is 2.184x, so dispatch
   threshold alone cannot be the full solution.
   The broader bottleneck is repeated conversion
   `FloatAP -> Rational -> Fraction`, creation of
   rational intervals, and two independent
   `FloatAP::from_rational` conversions in `certified`.
   Introduce a read-only dyadic significand/exponent
   view, outward fixed-point interval endpoints,
   and direct agreement-to-rounded-result from those
   endpoints. Use precomputed certified constants,
   preserve Ziv-style refinement and all four modes,
   and retain the general correctly rounded evaluator.
   The same design should accelerate cos, atan, tan,
   exp/log and near-half inverse trigonometry without
   baking 1/2 into the architecture.

**Performance methodology and correctness gates.**
Keep the canonical 27-operation/4-precision 1000x11
report and record absolute AP and MPFR times as well
as all 108 ratios. Benchmark operations on the same
unloaded Release host; pair repeated A/B runs or
alternate backend order to detect frequency, cache
and thermal drift. Instrument heap allocations,
temporary integer sizes, Ziv refinement counts,
cached-center hits, time spent in rational
conversion and exact division versus kernel work.
Add an independent domain-diverse corpus (signs,
large/small exponents, away from 1/2, hard-to-round
cases, both cold and warm certified constants);
never redefine the existing 27-operation corpus
just to obtain favorable ratios. A 53-bit FloatDP
comparison is a separate speed reference, not
a correct-rounding oracle for every libm function.

The implementation branch at `20187a89` had successful
Unix, Windows and Coverage workflow runs
39070083926, 39070083947 and 39070083939.
A passing Coverage workflow does not itself
certify the separate Codecov patch-coverage status.
This section is a **design/analysis-only** journal
addition: none of P0-P4 is merged, no source
regressions are altered, and no new CI or
benchmark speedups are claimed for proposed changes.

## P0 implemented: normalized equal-precision comparison (17 October 2026)

The implementation of `cmp(FloatAP,FloatAP)` now returns immediately
when both finite nonzero operands have matching precision but distinct
binary exponents. A normalized positive significand always occupies
exactly its declared precision in bits, so exponent order is magnitude
order for equal precision; an already-checked negative sign reverses it.
Special values, signed zeros, opposite-sign operands, equal exponents
and mixed precisions retain their existing paths.

Regression tests in the canonical `test_float_rounding` suite compare
FloatAP ordering with exact Rational ordering across all pairings of
positive/negative same- and differing-exponent values at
7/53/128/256/1024 bits, together with reversed mixed-precision
comparisons and comparisons at extreme representable exponents.
This does not require or introduce any additional GMP API.

**Validation state at authoring:** source/tests reviewed and committed;
a complete CI run and same-host benchmark must be checked before
claiming correctness on Unix/Windows or any measured speedup.
The 21:31 benchmark above remains the last measured baseline.

## P1 implemented: sqrt remainder-based nearest rounding (17 October 2026)

For the exact integer square root `mpz_sqrtrem(root,remainder,scaled)`
with `scaled=root*root+remainder`, the old nearest-rounding
criterion `4*scaled>(2*root+1)^2` simplifies algebraically to
`4*remainder>4*root+1`. Since both values are integers,
this is equivalent to `remainder>root` (the half-integer
tie is impossible). The nearest branch now uses a direct
`mpz_cmp(remainder,root)`: no temporary threshold, threshold
square, or quadrupled `scaled`. The zero-remainder case,
directed rounding, result carry, exponent checks, sign and
special values are unchanged. The scratch allocations shrink
from four to three initialized `mpz_t` objects.

A separate integer model checked the old/new predicates on
about 60,000 root/remainder boundary samples as a sanity check;
the algebra above, not the sample count, provides the proof.
The shared precision regression adds FloatAP/MPFR differential
sqrt tests at 7,53,128,256,1024 bits and in nearest-even,
downward, upward and toward-zero modes. Inputs include exact
squares, non-squares, low-scale fractions, both sides of
1 and 4 and large positive/negative exponent shifts. The
reference evaluates the exact same rounded FloatAP input.

**Validation state at authoring:** code and regression added
as an independent P1 commit after the P0 commit.
Unix/Windows/Coverage status must be checked and a new
unloaded same-host 1000x11 measurement performed before
attributing any speedup. The original 21:31 benchmark is
still the performance baseline. No MPFR or GMP dependency
has been added; the long-term independent-limb design
remains the intended route for P2/P3.

## P0/P1 benchmark evaluation and independent P1 re-run (17 October 2026)

**Measured reports and identification:** The original
`Markdown incollato(20261017-213109).md` precedes P0.
`Markdown incollato(20261017-214146).md` is the first full
1000-iteration, 11-sample report after P0 commit `d9e9b348`.
`Markdown incollato (2).md` is the initial analogous report
after P1 commit `bb38b4b1`.
The independent repeat after P1 is
`Markdown incollato(20261017-214445).md`. These are
user-supplied, not benchmarks run by this journal change.
All four contain 27 operations per precision at
53/128/256/1024 bits, paired FloatAP and MPFR timings.
Do not change the workload to target particular ratios.

| Bits | Metric | Before P0 | After P0 | First P1 | Repeat P1 |
| --- | --- | ---: | ---: | ---: | ---: |
| 53 | geomean AP/MP | 0.961x | 0.932x | 0.974x | 0.937x |
| 53 | total AP/MP | 1.481x | 1.460x | 1.486x | 1.462x |
| 53 | worst AP/MP | 3.163x cos | 3.020x cos | 3.179x cos | 3.343x div_nearest |
| 128 | geomean AP/MP | 1.332x | 1.279x | 1.336x | 1.299x |
| 128 | total AP/MP | 1.120x | 1.087x | 1.301x | 1.138x |
| 128 | worst AP/MP | 2.902x div_nearest | 2.895x div_up | 2.984x div_up | 2.930x div_up |
| 256 | geomean AP/MP | 1.035x | 1.015x | 1.011x | 1.012x |
| 256 | total AP/MP | 0.579x | 0.588x | 0.573x | 0.588x |
| 256 | worst AP/MP | 2.275x div_nearest | 2.228x div_up | 2.443x div_up | 2.192x div_nearest |
| 1024 | geomean AP/MP | 0.759x | 0.755x | 0.754x | 0.739x |
| 1024 | total AP/MP | 0.250x | 0.252x | 0.251x | 0.248x |
| 1024 | worst AP/MP | 2.104x compare | 1.975x add_nearest | 1.974x add_down | 2.024x add_nearest |
| ALL | geomean AP/MP | 1.001x | 0.977x | 0.998x | 0.977x |
| ALL | total AP/MP | 0.451x | 0.455x | 0.467x | 0.451x |

The AP/MP ratios >=2x per precision are, in this
same run order:
53: 3/3/4/3; 128: 10/9/10/8;
256: 6/4/6/3; 1024: 1/0/0/1;
ALL: 20/16/20/15. The second P1 run therefore has
15 of 108 pairs >=2x, but the stable <2x target is
not achieved. In particular, the 1024-bit worst
moves from just under to just over 2x, showing why
a small-margin pass should not be treated as robust.

**P0 and P1 intended kernels have independently
clear latency improvement, unlike aggregate variation.**
FloatAP `compare` before P0 -> first P0,
in ns: 53 4.958->2.750; 128 5.083->2.792;
256 5.333->2.750; 1024 5.083->2.959,
decreases of 41.8-48.4%. These gains persist
after P1: in the repeated P1 run, compare takes
2.750, 2.792, 2.750, 2.750 ns respectively.
The FloatAP `sqrt_nearest` latency after P0 ->
first P1 -> repeated P1, ns:
53 138.959->100.167->97.833,
128 203.083->133.292->127.125,
256 244.292->165.000->161.917,
1024 514.708->362.750->334.250.
The P1 repeat ratios for sqrt are
53 1.995x, 128 1.343x, 256 1.438x,
1024 1.210x. At 53, 1.995x is a
*borderline single-run pass*; the initial P1
run was 2.343x because the MPFR reference
time differed. Preserve both observations.

**The suspicious 128-bit `log`/`sin` slowdown is
not reproduced.** Their FloatAP medians (ns):

| Operation @128 | After P0 | First P1 | Repeat P1 |
| --- | ---: | ---: | ---: |
| log_nearest | 1430.042 | 3205.166 | 1439.041 |
| sin_nearest | 1259.709 | 1831.708 | 1236.875 |

The initial P1 run had a 2389.334 ns rise
in the sum of all 27 FloatAP @128 medians;
log (+1775.124 ns) and sin (+571.999 ns)
explained 98.2% of it despite neither kernel
being changed by P1. The repeat puts the
FloatAP sum back near the post-P0 baseline:
P0 14068.414 ns; first P1 16457.748 ns;
repeat P1 14084.247 ns. Thus the 128-bit
aggregate 1.301x to 1.138x recovery has
a concrete explanation in non-reproduced
log/sin outliers. It does not identify their
root cause; CPU load, clock and cache effects,
or measurement noise remain unverified hypotheses.
The P0 128-bit MPFR sum was 12945.249 ns
versus 12371.458 ns on P1 repeat, so even
with roughly equal absolute FloatAP total,
the relative total ratio changed from 1.087x
to 1.138x because the reference backend
timings changed. Never attribute this
difference automatically to P1.

**Additional single-run comparator anomaly @53.**
The P1 repeat global worst is
`div_nearest`@53 = 108.792 / 32.542
= 3.343x. The MPFR division median of
32.542 ns is unusually low compared to
the before-P0, P0, and initial P1
MPFR medians 60.167, 60.167 and
59.667 ns respectively. FloatAP division,
by contrast, was 109.875, 109.875,
97.709 and 108.792 ns across the
four runs. This mismatch may be
measurement noise or benchmark effects,
NOT proven speed regression. Under the
ordinary ~60 ns MPFR reference, the
observed 108.792 ns would be ~1.81x.
Do not simply discard the recorded 3.343x:
repeat the 53-bit MPFR/FloatAP division
measurement until the reference discrepancy
is understood, retaining raw absolute times.

**Verified CI status.** Both implementation
commits passed Unix, Windows and Coverage:
P0 `d9e9b348`: Unix 39071901197,
Windows 39071901178, Coverage 39071901183;
P1 `bb38b4b1`: Unix 39071981314,
Windows 39071981324, Coverage 39071981337.
The separate Codecov patch status has
not been verified here. These benchmark
comparisons are not application-wide
performance claims, and MPFR is still
used as differential oracle / benchmark
rather than as native FloatAP arithmetic.

**Decision:** Keep P0 and P1. Neither of
the anomalous 128-bit log/sin results
reappeared, and the target operation
latencies improve. Do not alter or roll
back correctly rounded source code
to chase unisolated benchmark variance.
P2 (portable, GMP-independent interface
for small-limb division and exact
quotient/remainder) is still the next
important architectural optimization,
especially for 128 and 256 bits. P4
(dyadic fixed-point certification to reduce
rational allocations) remains important
for the 53-bit time-total objective.
Before calling any <2x outcome stable,
repeat cases with narrow margins, including
sqrt@53 and add_nearest@1024, and isolate
MPFR div_nearest@53.

## P2 phase 1: portable small-limb exact division at <=64 bits (17 October 2026)

**Implementation.** New internal `src/floatap_small_division.hpp`
defines fixed-capacity unsigned integers in 32-bit little-endian
limbs. It provides bit length, compare, shifts, increment and a
Knuth-style normalized quotient/remainder implementation with
a 64/32 single-divisor-limb fast path. It uses only standard C++20
`uint32_t`/`uint64_t` and stack arrays; NO MPFR, GMP or mandatory
`__int128` is involved in the integer arithmetic. The only GMP
bridge, `import_small_magnitude`, extracts up to 64 bits of
magnitudes from existing `mpz_t` storage on both 32- and
64-bit GMP limb configurations. Results are imported into the
final `FloatAP::_significand` once. Every other configuration
retains the original exact GMP division. No FloatAP public API,
precision semantics or benchmark workload is changed.

The `from_ratio` dispatch is deliberately restricted to
`precision <= 64`, numerators and denominators with <=64-bit
magnitudes and supported GMP limb layouts. Exact binary magnitude
normalization is followed by aligned integer division, remainder
comparison (including nearest-even parity ties), sign-aware
directed rounding and a single carry adjustment. The original
checked exponent combination with the dyadic scale remains
in force. If the fast path is inapplicable or scratch capacity
is insufficient, the existing GMP path remains the fallback;
zero/one denominator special handling is unchanged.

**Why not enable 128-bit here?** The pure 32-bit-limb prototype
passed local quotient/remainder and rounding tests up to
128-bit inputs, but isolated Linux x86-64 timing with the
same workload shape indicated ~80 ns versus ~135-165 ns
for an existing-GMP-style 53-bit nearest/directed ratio,
whereas at 128 bits it was ~190-215 ns versus
~178-208 ns for the GMP-style path. These microbenchmarks
include harness-side allocation/rounding and are NOT
a full library performance result or the user's Apple-host
1000x11 benchmark. Activating 128 bits before proving a
faster algorithm would risk regressions in the worst area.
The portable primitive may be reused in a faster wide-limb
implementation (or on future GMP-free storage), but 128-bit
dispatch stays unchanged in this commit.

**Independent arithmetic validation:** Standalone locally
compiled C++20/GMP oracle test passed 300000 random
full multi-limb divrem cases, 100000 random shifts,
and 500000 positive/negative ratio roundings at
7/24/53/113/128 bits, across all four rounding modes,
including varied operand-length ratios. These exercise
the shared pure header and integration-equivalent
normalization/rounding algorithm. The committed canonical
`test_float_rounding` suite further compares the AP
result to MPFR at 7,53,64,128,256 bits with 32
nontrivial input pairs, both signs, all four modes,
both Rational construction and FloatAP division on
identical rounded inputs. Previous tests remain intact.

**Release/performance gate:** This is a P2 *partial*
implementation, not completion of 128/256-bit division
specialization. Require Unix/Windows/Coverage success
(including compiler warning-as-error checks) and a
same-host unloaded 1000x11 user report before claiming
a FloatAP timing benefit. Inspect 53-bit division and
Rational construction absolute times, 128-bit division
for stability, all 108 AP/MPFR ratios, per-precision
geomean/total and worst; compare against the 21:44
P1 repeat. The future P2 wide-limb kernel should
remain GMP-independent and must *earn* its 128-bit
dispatch with isolated measurements, not optimistic
extrapolation.

## Correctness audit of the performance comparison (17 October 2026)

**Critical measurement limitation identified in source review.**
The current `benchmarks/benchmark_backend_workload.hpp`
constructs 32 exact Rational inputs independently as FloatMP and FloatAP
at nominal precision 53/128/256/1024 bits. The `report` lambda
only calls `median_ns`; each operation's return value is kept
alive but never numerically compared with the other backend.
`benchmark_compare_multiple_precision_backends.cpp` pairs only
the **timing** records. Consequently all previous speed ratios,
including the 21:44 P1 repeat, are **timing-only evidence**,
NOT confirmation that every benchmark operation computes the
same correctly rounded numerical value. Nor do they prove
that the nominal target precision is achieved on all inputs.

**Distinguish representation from mathematical accuracy.** The
FloatAP invariant is a p-bit normalized integer significand for
finite nonzero values, along with a signed power-of-two exponent.
This representation alone does not guarantee p bits of accuracy:
a low-accuracy result padded with trailing zeros might also have
a p-bit representation. For transcendental functions, the
certified-enclosure/Ziv-refinement algorithm *is intended* to
guarantee correctly rounded p-bit results IF every computed
interval really encloses the mathematical value and the rational
endpoint rounding is correct. These premises require independent
validation. Existing differential FloatAP-vs-MPFR tests in
`tests/test_float_rounding.cpp` include multiple precisions,
rounding modes and boundary cases, but are not an exhaustive
oracle for every benchmark input, every operation or the whole
real domain.

**Mandatory next validation gate, before further optimization.**
Introduce a correctness preflight into
`benchmark_compare_multiple_precision_backends`, reusing
the exact same shared workload's operations and 32 inputs.
Outside ALL timed regions, compare at each requested precision:

- All 32 converted `a` and `b` input values across backends
  as exact dyadic `Rational` values.
- For every one of the 27 operation names and 32 input indices,
  require exact numerical equality between independently
  evaluated FloatAP and FloatMP results. Float outputs compare
  through exact `Rational` conversions (no conversion to
  `double` and no tolerance); comparison outcomes compare
  directly; `to_double` compares the exact binary64 bit
  patterns, including sign bits; `to_rational` compares exactly.
- Assert the requested precision in both floating backends and
  the p-bit significand normalization for each nonzero finite
  FloatAP output. This is a structural invariant check,
  *not* a substitute for numerical agreement.
- Fail explicitly with precision, operation, input index,
  rounding mode and exact differing values. No benchmark
  slowdown aggregates are valid on a failed preflight.
- Run the preflight *before* measurements, preferably in
  a separate phase for BOTH backends, and ensure all expensive
  conversions/checks lie outside `median_ns`. Preserve the
  1000x11/27-operation timing workload, one warm-up sample and
  the independent FloatDP-vs-FloatAP speed benchmark.
- Keep the complete four-mode differential test suite and extend
  domain-diverse/hard-rounding tests independently of the
  32-element benchmark pool; a successful benchmark preflight
  establishes equality only on that finite pool.

There are 27 * 32 * 4 = 3456 operation/precision/index
comparisons, plus 2 * 32 * 4 = 256 input comparisons,
before any special-value or adversarial expansion.
This numerical preflight is REQUIRED before treating
<1x or <2x performance objectives as meaningful
correct-rounding comparisons.

**Independent CI blocker on P2.** As checked after commit
`4e23c5db`, Unix run 39072784107 (GCC/Clang/macOS
and sanitizer jobs), Windows run 39072784090, and
Coverage run 39072784093 all have conclusion FAILURE,
with the named Build step failing. The GitHub connector
does not expose the detailed build diagnostic through
its allowed job-log endpoints, so the exact compiler
error is NOT known. Do not guess a root cause or claim
P2 functional correctness, portability or speedup.
Resolve the actual compiler diagnostics first, rerun
Unix/Windows/Coverage and separately review Codecov
patch coverage. P0/P1 precede this failed P2 commit
and were green. The validator design above is not yet
implemented: this journal change makes no new source,
test or benchmark performance claims.

## P2 compile repair: countl_zero signedness (17 October 2026)

The user supplied the actual macOS AppleClang -Werror diagnostic from P2:
`std::countl_zero(uint32_t)` returns `int`, producing two
`-Wsign-conversion` failures in `src/floatap_small_division.hpp`:
`bit_length` and the division normalizing shift. Explicit, range-safe
casts of the [0,32] count to `std::size_t` and `unsigned` fix those
specific errors without changing any integer operation or rounding.
Verify Unix, macOS, Windows and Coverage before claiming overall CI green.

## Independent correctly-rounded literal oracles (17 October 2026)

The precision regression now includes **fixed exact expected dyadic
literals**, not values computed by FloatMP or by FloatAP. Twelve
rows cover `sqrt(2)`, `exp(1/2)`, `log(2)` and `pi` at
128, 256 and 1024 bits. Each row contains the exact nearest-even
integer significand and binary exponent plus independently certified
integer offsets for downward and upward rounding (toward zero equals
downward for these positive values). The test checks exact Rational
equality for all four modes as well as reported result precision.

The source constants were generated offline using Python's standard
library `fractions.Fraction` and integer arithmetic ONLY:
`sqrt(2)` by floor integer square root of `2*2^(2*1300)`;
`exp(1/2)` by a positive Taylor series with an explicit geometric
tail bound after decreasing terms; `log(2)` by the positive atanh
series `2*atanh(1/3)` with a bounded remainder; `pi` from the
Machin identity `16*atan(1/5)-4*atan(1/239)`, using rational
alternating-series bounds. Each lower and upper rational endpoint,
rounded separately by an independently implemented integer
quotient/remainder nearest-even and directed rounding routine,
yielded the SAME target dyadic literal at 128, 256 and 1024 bits.
Thus the constants assert the expected mathematically correctly
rounded result for the chosen arguments; they are NOT merely
comparisons with MPFR, p-bit significand width, or assertions
that an answer lies between directed endpoints.

This is strong but finite coverage, not a proof of universal correct
rounding for all transcendental inputs. The previously existing
MPFR differential tests remain valuable auxiliary regression checks,
without constituting the source of truth for these fixed oracles.
The benchmark still lacks a numerical preflight and no benchmark
performance number is affected by this test-only addition.

## Compile repair of independent expected-literal tests (17 October 2026)

The 128/256/1024-bit hard-coded exact dyadic literal regression was
initially written with `Bits(row.bits)`, but `Bits(unsigned long
long)` is deliberately private: only the `_bits` literal factory
may instantiate this tag. The public `ArbitraryPrecision(std::size_t)`
constructor explicitly supports dynamic precision counts. Construct
`ArbitraryPrecision p{static_cast<std::size_t>(row.bits)}` instead.
No expected mantissa, rounding offset, case or comparison changes.
Full CI pass must be confirmed before calling the test valid.

## P2 compiler portability follow-up: GMP signed limb indices (17 October 2026)

GCC 13 with `-Wsign-conversion -Werror` reported that
`mpz_getlimbn(input,i/2u)` passes `std::size_t` to the signed
GMP index type `mp_size_t`. The same latent diagnostic exists in
the 32-bit GMP-limb branch `mpz_getlimbn(input,i)`. Both call sites
now use explicit `static_cast<mp_size_t>` on the index. In this
adapter `bits <= 64`, so `words <= 2`, `i < words` and the
indices are always 0 or 1, making the conversion lossless on all
supported GMP platforms. No division algorithm, rounding branch,
correctness oracle, expected literal or timed workload is changed.
This is a diagnostic-specific source fix: verify that the post-fix
Unix/Windows/Coverage CI passes before claiming portability.

## P2 Windows compiler repair: shadowed Ariadne rounding name (17 October 2026)

MSVC 14.51 with `/W4 /WX` reported C4459 in `src/floatap.cpp`: the `small_rounded_ratio` output argument `rounded` shadows the existing namespace-level `Ariadne::rounded` declaration (`include/numeric/rounding.hpp`). The output parameter has been renamed to `rounded_significand`, including its assignment, without touching calling syntax, algorithm, expected results, or the fallback. This is a source-level diagnostic correction, **not evidence that Windows or the full CI is green** until the next workflow completes.

## P2 first full green CI and new-host 1000x11 performance (17 October 2026, 22:14)

**Source and validity.** The user-provided
`Markdown incollato(20261017-221453).md` is the first
1000-iteration, 11-sample report after the P2 small-limb division
change (`4e23c5db`) and the subsequent compiler/portability fixes,
with final implementation HEAD `272b6f59`. Its timing measurements
are from the user's benchmark host, not an execution here.
The immediately preceding baseline is the *repeated P1* report
`Markdown incollato(20261017-214445).md`, with HEAD `bb38b4b1`
(no P2 fast path). Same four precisions, 27 workloads per precision,
both backends and original exact-rational input generator.

**Build/test gate now green:** Verified directly from GitHub Actions:
Unix 39074145329 SUCCESS, Windows 39074145320 SUCCESS,
Coverage 39074145321 SUCCESS, all at commit `272b6f59`.
These runs include the committed, fixed literal-oracle tests
for `sqrt(2)`, `exp(1/2)`, `log(2)` and `pi` at
128/256/1024 bits; a green test run is evidence for those
particular cases and the other test-suite cases, not a
universal mathematical proof. Separate Codecov patch status
is not asserted. The benchmark still compares **timings** only:
it does not compare output values on all 32 inputs per operation.

| Bits | P1 repeat geomean | P2 geomean | P1 repeat total | P2 total | P1 repeat worst | P2 worst | P2 ratios >=2x |
| --- | ---: | ---: | ---: | ---: | --- | --- | ---: |
| 53 | 0.937x | **0.860x** | 1.462x | **1.436x** | 3.343x div_nearest | 3.714x cos | 3 |
| 128 | 1.299x | **1.268x** | 1.138x | **1.117x** | 2.930x div_up | 2.902x div_down | 7 |
| 256 | 1.012x | **0.999x** | 0.588x | **0.579x** | 2.192x div_nearest | 2.185x div_nearest | 4 |
| 1024 | 0.739x | **0.723x** | **0.248x** | 0.251x | 2.024x add_nearest | **1.940x add_nearest** | 0 |
| ALL | 0.977x | **0.942x** | 0.451x | 0.454x | 3.343x div_nearest | 3.714x cos | 14 |

The first post-P2 report satisfies all three constraints at
1024 bits. At 256 bits, the geomean 0.999x has an
approximately 0.1% margin: this is a *borderline*
single-run pass, not yet a stable achievement. At
128 bits neither mean nor total target is achieved;
at 53 bits the total is still substantially too slow.
At 53/128/256 bits individual operations remain >=2x.

**P2 intended-path timings (FloatAP, nanoseconds):**

| 53-bit operation | P1 repeat ns | P2 ns | Change | P2 AP/MPFR ratio |
| --- | ---: | ---: | ---: | ---: |
| construct_rational | 157.625 | 123.667 | -21.5% | 0.488x |
| div_nearest | 108.792 | 78.709 | -27.7% | 1.114x |
| div_down | 74.958 | 74.750 | -0.3% | 1.133x |
| div_up | 82.625 | 74.083 | -10.3% | 1.110x |

The elapsed-time reductions in construction and nearest
division are evidence of the P2 small-limb dispatch
helping those measured paths. It does **not** provide
a consistent ~30% improvement to all division modes.
At 53 bits the earlier anomalously low MPFR nearest
division timing of 32.542 ns recovers to 70.625 ns,
removing the spurious 3.343x worst; the actual FloatAP
nearest division also improves independently.
There is no new 128/256/1024 fast-path dispatch:
the changes to ratios at those precisions should
not be assigned to P2. In particular the 128-bit
division ratios after P2 remain
2.794x nearest, 2.902x downward, 2.869x upward,
and at 256 bits
2.185x nearest, 2.043x downward, 2.059x upward.

**Do not confuse ratio changes with AP regressions.**
At 53 bits, `cos_nearest` changes only 1427.084 ->
1442.250 ns (+1.1%) for FloatAP, but the MPFR time
drops 463.417 -> 388.333 ns (-16.2%). Hence its
ratio rises 3.079x -> 3.714x largely due to the
reference variation, not a comparably large AP slowdown.
The new 53-bit `sin_nearest` is 2.011x and
`atan_nearest` 2.462x. The newly observed
AP total summed across all operations at 53 bits
*increases* 11939.789 -> 12184.251 ns; MPFR
increases even more, 8164.170 -> 8485.290 ns.
Accordingly the improved **relative** total ratio
1.462x -> 1.436x does not mean absolute AP total
fell. The 128-bit log/sin outlier of the first P1
report does not recur in the new run (log 1438.625,
sin 1280.458 ns). Retain raw absolute measurements
and rerun near-threshold cases.

**Priority by absolute target budget.** At 53 bits
the sum of individual AP medians exceeds MPFR
by 3698.961 ns. `cos`, `atan`, `sin`, `tan`
alone account for excesses of approximately
1053.917, 808.750, 681.959 and 662.875 ns,
respectively (total +3207.501 ns). Reducing
these *certified transcendental* costs is essential
to meeting the <1x **total** objective at 53 bits.
At 128 bits the excess is +1520.209 ns, led by
`tan` (+793.166), `cos` (+705.458),
`sin` (+531.083), `exp` (+508.458)
and `log` (+222.834), although faster atan,
asin and acos partly offset these. Meanwhile
division remains the clear short-latency
worst-ratio issue at both 128 and 256 bits.

**Next work, contingent on the correctness gate.**
Keep the original P2 <=64-bit fast path and
portable GMP-free unsigned-limb kernel.
Do not raise the 128-bit fast-path threshold:
the earlier isolated 32-bit-limb prototype was
not faster than GMP there. Design a better
portable wide-limb algorithm with fallback and
measure it before dispatching. In parallel,
the P4 fixed-point interval/certified-rounding
work for 53 and 128 bits is now more important
for *total-time* acceptance than additional
53-bit division micro-optimizations. Preserve
the independent exact expected-literal test
approach for correctly rounded results, extend
hard-to-round/transcendental domains and
add a correctness preflight to the benchmark
outside all timed regions. No kernel or test
source is modified by this journal-only entry.

## P4 phase 1: remove duplicate exact-rational argument materialization (18 October 2026)

**Change.** The shared transcendental input adapter previously performed
`FloatAP -> Rational -> Fraction`: `FloatAP::operator Rational()`
allocates/canonicalizes an `mpq_t`, after which `Fraction(mpq_srcptr)`
allocates another pair of GMP integers and copies the numerator and
denominator. A new `Fraction(FloatAP const&)` instead calls the existing
`FloatAP::to_rational(mpq_ptr)` once into the owned `mpq_t`.
`rational_argument()` uses this constructor, eliminating one rational
instance and its deep copy for **exp, log, sin, cos, tan, atan, asin,
acos**, without changing the exact dyadic value, the interval
algorithms or the certified Ziv-refinement criterion.
This remains an incremental allocation improvement, NOT the complete
GMP-independent interval backend; the new code adds no GMP calls or
persistent dependencies beyond those already used by FloatAP.

**Independent expected-literal extension.** The shared
`test_float_rounding` suite also checks explicit dyadic integer
significands for `sin(1/2)` and `cos(1/2)`, at 53, 128, 256, and
1024 bits in all four rounding modes. Expected values were generated
offline, without MPFR or FloatAP: exact `fractions.Fraction` terms
of the alternating sine/cosine Taylor series at x=1/2, with the
first omitted signed term as the rigorous interval bound.
For each result and rounding mode, **both rational interval endpoints
independently round to the same dyadic literal** using exact
quotient/remainder and nearest-even parity. All values are compared
exactly as rationals; tests do not accept loose tolerances,
rounding-direction containment alone, or arbitrary zero padding.
The direct-argument change also remains subject to the unchanged
existing MPFR differential and domain-diverse suites.

**Acceptance gate.** The source/expectations have not been timed
on the user's benchmark host yet. Check Unix, Windows, Coverage and
all literal oracles first, then repeat the canonical unchanged
1000x11 benchmark at 53/128/256/1024 bits. Examine absolute
times (particularly 53-bit cos/sin/atan/tan and 128-bit cases)
and report all geomean/total/worst ratios. A CI pass proves only
the exercised cases, not correctness over every transcendental input.
Do not claim a speedup or prematurely dispatch centered algorithms
at lower work precision without an A/B measurement.

## P4 phase 1, first same-host 1000x11 measurement (17 October 2026, 22:20)

**Source and reproducibility.** User-supplied
`Markdown incollato(20261017-222054).md` reports
27 operations, 53/128/256/1024 bits, 1000 iterations and
11 median samples per backend. It follows P4 commit
`72659ca7` (one exact FloatAP-to-Fraction materialization,
without intermediate Rational deep copy), and is compared
to the immediately preceding P2 measurement
`Markdown incollato(20261017-221453).md`
(implementation `272b6f59`). The report is a measured
timing comparison, **not a check of benchmark outputs**.
The independent literal-oracle tests are a separate
functional correctness gate.

| Bits | P2 -> P4 geomean AP/MPFR | P2 -> P4 total AP/MPFR | P2 -> P4 worst | P2 -> P4 count >=2x |
| --- | --- | --- | --- | --- |
| 53 | 0.860x -> 0.941x | 1.436x -> 1.431x | cos 3.714x -> cos 2.966x | 3 -> 2 |
| 128 | 1.268x -> 1.271x | 1.117x -> 1.068x | div_down 2.902x -> div_nearest 2.907x | 7 -> 9 |
| 256 | 0.999x -> 0.984x | 0.579x -> 0.551x | div_nearest 2.185x -> div_down 2.262x | 4 -> 4 |
| 1024 | 0.723x -> 0.733x | 0.251x -> 0.244x | add_nearest 1.940x -> add_up 1.964x | 0 -> 0 |
| ALL | 0.942x -> 0.964x | 0.454x -> 0.436x | cos 3.714x -> cos 2.966x | 14 -> 15 |

**All eight affected transcendental kernels become faster in
absolute FloatAP time at all four precisions** (32 of 32
operation/precision comparisons in this first report):

| Function | 53 bits change | 128 bits change | 256 bits change | 1024 bits change |
| --- | ---: | ---: | ---: | ---: |
| exp | -6.1% | -6.9% | -6.2% | -4.3% |
| log | -9.1% | -5.7% | -5.7% | -4.8% |
| sin | -5.4% | -5.5% | -8.0% | -6.2% |
| cos | -5.5% | -5.3% | -6.2% | -4.6% |
| tan | -6.7% | -8.1% | -4.7% | -4.6% |
| atan | -5.0% | -9.0% | -4.6% | -3.0% |
| asin | -5.6% | -10.0% | -4.5% | -2.9% |
| acos | -4.0% | -6.4% | -4.3% | -4.1% |

Sum of the *eight affected* FloatAP medians, P2 -> P4:
53: 11244.584 -> 10593.750 ns (-5.79%);
128: 12951.665 -> 12003.666 ns (-7.32%);
256: 13760.958 -> 13027.045 ns (-5.33%);
1024: 21012.999 -> 20148.415 ns (-4.11%).
The remaining 19 AP medians in aggregate are approximately
+1.34%, -2.47%, -1.78%, +0.14% respectively.
This affected/unaffected contrast is consistent with the
removed rational deep copy reducing fixed overhead; it
is *evidence*, not a controlled statistical proof of
causality from a single run.

The sums of all 27 AP medians P2 -> P4 are:
53 12184.251 -> 11546.002 ns (-5.24%);
128 14489.835 -> 13503.831 ns (-6.80%);
256 15557.417 -> 14791.462 ns (-4.92%);
1024 25162.207 -> 24303.540 ns (-3.41%).
The corresponding MPFR sums move:
53 8485.290 -> 8070.458 ns (-4.89%);
128 12969.626 -> 12638.169 ns (-2.56%);
256 26888.127 -> 26834.337 ns (-0.20%);
1024 100051.712 -> 99486.874 ns (-0.56%).
Thus worsening 53/128/1024 geometric ratios does NOT
mean the transcendental AP kernels grew slower.

**Worst-ratio interpretation:** At 53 bits
cos falls 1442.250 -> 1362.916 ns in AP (-5.5%),
while MPFR rises 388.333 -> 459.458 ns (+18.3%).
Consequently the large ratio drop
3.714x -> 2.966x is not solely attributable to
P4. The current >2x 53-bit cases are cos 2.966x
and atan 2.385x. At 128 bits the >2x cases
include all three divisions (2.863-2.907x),
add nearest/down/up (2.309-2.392x), sub nearest
2.267x, cos 2.293x and construct_rational
1.625x is **not** >2x; 9 measured operations
total are >=2x. At 256 bits the ratios >=2x
include all three divisions (2.160-2.262x)
and add_nearest 2.021x. At 1024 bits all
27 ratios are <2x, though add_up 1.964x is
near the boundary.

**Remaining absolute deficits:** FloatAP total
minus MPFR total after P4 is +3475.544 ns at
53 bits and +865.662 ns at 128 bits; at 256
and 1024 bits the totals favor FloatAP.
The direct Rational-copy removal helps but
cannot by itself achieve the 53-bit <1x total
objective. The next P4 changes should target
the certified interval representation and
decision-to-round paths, measuring the
cost of `Fraction` endpoint construction,
exact-rational arithmetic, cached half-center
operations, and two independent endpoint
`FloatAP::from_rational` calls per refinement.
Keep outward enclosure proofs and independent
exact literal expectations: do not replace them
with less accurate fixed-precision approximations.
For 128-bit worst ratios, a separate portable
wide-limb P2 division algorithm and arithmetic
add/sub optimization are still necessary.

**CI status when recorded:** For `72659ca7`,
Unix run 39074633880 SUCCESS and Coverage
39074633895 SUCCESS; Windows 39074633860 was
still IN PROGRESS on last check. No claim of
all-three-green or separate Codecov patch
coverage is made. Re-run the same-host 1000x11
workload before accepting small geometric/worst
margin changes. This entry changes **only**
the design journal; it does not alter test
workloads or numerical kernels.

## Final post-timing correctness check in the FloatMP/FloatAP benchmark (18 October 2026)

**Motivation and contract.** The timing benchmark previously
did NOT compare calculated numbers. As requested, it now
reruns the precise same workload **after completing all timing
loops** and asserts exact result agreement, *before* printing
aggregate slowdown ratios or accepting `--max-slowdown`.
There is no validation code inside `median_ns`, so the timing
samples themselves have no additional comparisons, allocations
or exact-rational conversions.

**No workload drift:** `benchmarks/benchmark_backend_workload.hpp`
defines one `BackendPool<Float>` constructor (same 32 exact
Rational pairs and same nearest-rounded backend input
construction) and one `for_each_backend_operation` visitor
(same 27 operation names and lambda definitions). Both
the original `run_backend` timing and the new verifier call
these helpers. The unchanged median warmup and 1000x11
default iteration/sample settings remain in effect.
FloatDP/FloatAP shares the same visitor and input pool;
it must continue to compile and benchmark unchanged.

**Exact comparisons:** On each selected 53/128/256/1024-bit
precision the verifier compares all 32 pairs of constructed
inputs and 27*32 operation results. The default four-precision
run therefore makes **256 paired input value checks and 3456
paired operation-result checks**. Finite float results,
including `construct_rational`, copy, arithmetic,
transcendentals and `pi`, are compared using exact
`Rational` values; `to_rational` also compares exactly.
`compare` uses the returned comparison outcome directly.
`to_double` compares all 64 bits of IEEE binary64 via
`std::bit_cast<uint64_t>`, not an epsilon tolerance.
The precision of every returned FloatMP/FloatAP object
is checked against the requested precision.
First mismatch throws with precision, operation,
input index and, where applicable, both exact values
or double bit patterns. Successful execution logs a
clear numerical PASS and retains the original report.

**Scope/limitations:** This agreement checks the exact
measured *input pool* and the specific nearest/down/up
operation modes present in the 27-op battery. It is
additional differential evidence, NOT an independent
oracle for mathematical correct rounding. The hard-coded
independently certified expected dyadic literals at
53/128/256/1024 bits, plus boundary/property tests,
remain necessary. Other domains, extreme exponents,
special values and unexercised rounding modes require
separate semantic tests. The verifier currently runs
once at the end of a timing invocation, not within
the timed loop, to avoid changing performance results.
CI/build and the user's benchmark run must validate
the refactor; no timing advantage is claimed from
this instrumentation.

## Fix benchmark refactor linkage and switch diagnostics (17 October 2026)

**Actual diagnostic:** In the shared `BackendPool<FloatDP>` constructor,
`std::vector<FloatDP>::emplace_back(exact,FloatDP::ROUND_TO_NEAREST,precision)`
passes the rounding constant through a forwarding reference.
This ODR-uses `FloatDP::ROUND_TO_NEAREST`, an in-class initialized
static const integral member without an out-of-class definition, and
caused an undefined symbol at link time on Apple arm64.
The constructor now materializes the same value in a local
`Float::RoundingModeType` and passes that variable to both
`emplace_back` calls. This preserves all 32 exact input rationals,
input rounding, precision, and the 27-operation timed workload,
while removing the unintended static-member reference.
The independent MP/DP backends still share the same template.

**Compiler warning:** The new `CheckedValue::Kind` switch now has
an explicit throwing default to satisfy `-Wswitch-default` and to
reject any invalid kind. It does not relax the exact Rational,
comparison-outcome, binary64-bit, or precision checks.

**Validation gate:** Verify that both benchmark executables link
on AppleClang/MSVC/GCC, and that the Unix, Windows and Coverage
workflows succeed. Neither fix changes numerical kernels or
the expected dyadic literals. The end-of-benchmark numerical
comparison must be exercised by actually running the benchmark:
successful compilation alone is not proof it passes the 3456
operation comparisons.

## P4 phase 2: certify and enable centered sin/cos at 53 bits (18 October 2026)

**Motivation:** P4 phase 1 removed unnecessary Rational materialization
from all eight transcendental functions. The 22:20 report still has
53-bit cosine about 2.97x MPFR, with an AP absolute median of ~1363 ns.
The existing centered half-angle dyadic algorithm was only enabled
when `work >= 128`; the first 53-bit refinement requests `work=85`
(`bits + 32`) and therefore *always* fell back to a full sine/cosine
Taylor series, even for dyadic arguments within 2^(-43) of 1/2.

**Change:** Lower the guard for the already-existing standalone
sine/cosine centered optimization from `work >= 128` to
`work >= 64`. No formula, truncation bound, rounding decision,
certified interval or fallback is weakened. The rigorous
`|delta| <= 2^(-ceil(work/2))` test is unchanged, so only arguments
very near +/-1/2 enter the cached outward-dyadic Taylor
linearization. Values outside this window, inputs requiring
argument halving, and tangent (which needs both results and its
separate certified path) still execute their original algorithms.
Cache entries are always outward intervals, never a rounded
approximation. This is a *precision-independent mathematical*
special case, not a shortcut returning a cached low-accuracy
floating number.

**Independent oracle coverage:** A new static set of **32**
mathematically certified dyadic literal entries checks sin/cos
on both sides of +/-1/2 at 53/128/256/1024 bits, all
four rounding modes, including sign-dependent directed
rounding and exactly representable arguments. These
expectations were produced offline using rational-arithmetic
alternating Taylor series with the exact first omitted
term as a rigorous bound. At every precision and mode,
rounding both interval endpoints independently produced
the same integer mantissa/exponent. Every test checks exact
Rational equality to that fixed literal and verifies the
requested precision. Both the existing half-center literal
suite and the unchanged MPFR differential test suite remain.

**Release gate and caveats:** CI must pass Unix, Windows and
Coverage, and the 1000x11 benchmark's new 3456-result
post-timing equality verification must pass. Then compare
actual absolute sine and cosine latency at 53 bits with
the 22:20 report and assess all four per-precision
geomean/total/worst values. The improvement is expected
but not yet measured. Because the benchmark inputs
concentrate near 1/2, retain separate tests for
non-centered inputs and do not infer general-purpose
trigonometric gains outside the certified window.
The most recent independently attached report in this
conversation is still `Markdown incollato(20261017-222054).md`;
the user's message "Ultimo risultato" did not provide
a new accessible benchmark attachment at commit time.

## Post-comparison benchmark: exact value agreement and measured ratios (17 October 2026, 22:29)

The subsequently accessible user-pasted `Markdown incollato(20261017-222901).md` contains an actual **numerical verification PASS**: 256 paired input values and 3456 paired operation results, including exact dyadic Rational equality, comparison outcomes, and 64-bit double conversion bit patterns. This is the first explicitly logged success of the new checker on the canonical 4-precision, 27-operation, 32-input 1000x11 workload. It confirms agreement **on those inputs**; it does not replace the fixed independently certified mathematical literal oracles or establish correctness on all inputs. Its measurements precede P4 phase 2 commit `a8290e20`, whose focused 53-bit centered trig optimization was committed separately after this report was produced.

| Precision | 22:20 P4 geomean -> 22:29 | 22:20 total -> 22:29 | 22:29 worst |
| --- | --- | --- | --- |
| 53 | 0.941x -> 0.883x | 1.431x -> 1.442x | cos_nearest 3.166x |
| 128 | 1.271x -> 1.260x | 1.068x -> 1.043x | div_up 3.098x |
| 256 | 0.984x -> 0.980x | 0.551x -> 0.558x | div_up 2.245x |
| 1024 | 0.733x -> 0.724x | 0.244x -> 0.244x | add_down 1.966x |
| ALL | 0.964x -> 0.943x | 0.436x -> 0.438x | cos_nearest 3.166x |

At 53 bits FloatAP sine changes 1283.709 -> 1299.417 ns and cosine 1362.916 -> 1407.333 ns between these two timing-only runs, despite no further mathematical kernel edit in that interval. MPFR cosine varies 459.458 -> 444.500 ns. These numbers **must not** be interpreted as the performance of the next centered sine/cosine threshold adjustment, which has no user benchmark yet. At 128 bits overall weighted ratio 1.043x is a further marginal improvement but still fails <1x. 1024-bit meets all three constraints on the supplied measurement; 256-bit worst still exceeds 2x.

As separately checked, commit `94c3cce1` (benchmark linking repair) passed Unix and Coverage; Windows workflow 39075150674 is recorded as FAILURE. We do not have the Windows failure diagnostic and must not claim it was green. A later commit `a8290e20` has Unix/Windows/Coverage workflows started but was still pending at time of this journal entry. The Windows error, if still present, requires diagnosis before considering P4 phase 2 portable. Further benchmark comparisons must preserve correct rounding and all exact expected literals.

## P4 phase 2: first measured 53-bit centered sin/cos improvement (17 October 2026, 22:34)

**Measured source:** User-supplied `Markdown incollato(20261017-223420).md`,
1000 iterations, 11 samples, 27 operations at each of 53, 128,
256 and 1024 bits. Compared to `Markdown incollato(20261017-222901).md`
taken before centered-trigonometric dispatch threshold was lowered
from 128 to 64 work bits (implementation commit `a8290e20`).
The latest source journal commit `104f9f02` does not alter arithmetic.
The new report explicitly prints:

```
Numerical verification PASS: 256 paired inputs and 3456 paired
operation results (exact values, comparison outcomes and binary64 bits).
```

This confirms exact FloatAP/FloatMP agreement on the complete
*benchmarked input pool*, not universally correct rounding.
Independently certified dyadic expected-literal tests for sin/cos
near +1/2 and -1/2, at 53/128/256/1024 bits in all four
rounding modes, remain part of the regression suite.

**Targeted 53-bit gain in absolute FloatAP latency:**

| Operation | 22:29 ns (pre-P4 phase 2) | 22:34 ns (post-P4 phase 2) | Change | AP/MPFR at 22:34 |
| --- | ---: | ---: | ---: | ---: |
| sin_nearest | 1299.417 | 1009.458 | -22.31% | 1.450x |
| cos_nearest | 1407.333 | 1010.875 | -28.17% | 2.203x |
| tan_nearest | 1440.541 | 1437.167 | -0.23% | 1.801x |
| atan_nearest | 1331.541 | 1329.750 | -0.13% | 2.419x |

The sum of sine and cosine medians improves by
686.417 ns; the full 53-bit sum of 27 FloatAP
medians improves by 666.544 ns, from 11843.751
to 11177.207 ns (-5.63%). Thus the intended
kernels explain essentially all the absolute AP
saving, with the other operations net +19.873 ns.
The corresponding 53-bit MPFR sum changes from
8214.874 to 8403.711 ns (+188.837 ns).
The improvement in the *relative* total ratio
therefore reflects both real FloatAP savings and
between-run MPFR drift. The speedup is measured
on this concentrated near-1/2 workload and should
not be extrapolated to arbitrary trig arguments.

| Bits | P4 phase 1 -> phase 2 geomean AP/MPFR | P4 phase 1 -> phase 2 total | P4 phase 2 worst | Phase 2 >=2x count |
| --- | --- | --- | --- | ---: |
| 53 | 0.883x -> **0.835x** | 1.442x -> **1.330x** | 2.419x atan_nearest | 2 |
| 128 | 1.260x -> **1.239x** | 1.043x -> 1.068x | 2.969x div_up | 8 |
| 256 | 0.980x -> **0.946x** | 0.558x -> **0.547x** | 2.101x div_nearest | 2 |
| 1024 | 0.724x -> 0.725x | 0.244x -> 0.245x | 1.940x add_nearest | 0 |
| ALL | 0.943x -> **0.918x** | 0.438x -> **0.434x** | 2.969x div_up | 12 |

53-bit cosine slowdown is substantially reduced
from 3.166x to 2.203x but remains above the
strict 2x ceiling; atan is now worst at 2.419x.
The 53-bit absolute excess over MPFR summed over
all 27 operations is still +2773.496 ns
(11177.207 - 8403.711). The most important
remaining positive excesses include atan
(+779.958 ns), tangent (+639.042 ns),
cosine (+552.084 ns), logarithm (+415.250 ns),
acos (+373.583 ns), asin (+307.499 ns)
and sine (+313.458 ns). Thus reaching <1x
weighted total at 53 bits needs broader
certified-kernel improvements, not another
tiny arithmetic primitive optimization.

128-bit total becomes 1.068x in this run,
versus 1.043x immediately before; absolute
FloatAP total rises from 13595.208 to
13860.167 ns (+264.959 ns), chiefly from
asin/acos/tan timing changes, NOT the changed
53-bit centered trig dispatch. No persistent
128-bit regression is established from one run.
256-bit worst remains above 2x, though only
two division modes exceed the threshold.
1024 bits passes all three current performance
goals for this run (worst 1.940x).

**CI:** GitHub Actions for source commit
`a8290e20` (Unix 39075450493,
Windows 39075450467, Coverage 39075450477)
all finished SUCCESS. The immediately
following journal-only commit `104f9f02`
also completed Unix 39075496983, Windows
39075496973 and Coverage 39075496977
SUCCESS. As requested by the user,
the unrelated proposed Windows Ninja
link-pool/CMake workaround was **not applied**.
This entry changes only design documentation,
and no new performance numbers are claimed
without an additional host benchmark.

**Next algorithmic target:** 53-bit
`atan` (largest remaining worst ratio and
large absolute excess); then investigate the
certified `tan` path and fixed-point rounding
from certified endpoints. For 128-bit worst
ratio, retain separate work on portable
division and exact add/sub. Expand
independent hard-to-round/domain-diverse
literals before extending an optimized
shortcut. Keep the exact benchmark
post-timing equality verifier enabled.

## P4 phase 3: fixed-point certified centered atan near +/-1/2 (18 October 2026)

**Motivation.** The 22:34 benchmark after P4 phase 2 showed
`atan_nearest`@53 = 1329.750 ns, versus MPFR 549.792 ns
(2.419x), the worst remaining 53-bit ratio; it also contributed
about 780 ns to the 53-bit weighted-time excess.
The existing `centered_atan_half_interval` already uses
`atan(1/2)` and a certified first-order Taylor enclosure,
but it computes the linear term `4*delta/5`, remainder
`delta*delta` and endpoint additions in full rational arithmetic
on every invocation.

**Change.** Keep the same exact neighborhood gate
`|delta| <= 2^-ceil(work/2)`, work >=4, and the same
mathematically certified half-center interval. Cache
outward integer bounds of `atan(1/2)` on scale
`Q=2^(work+24)` (upgrading the thread-local cache only
when greater work precision is requested). For dyadic
`delta=n/2^k` with k <= scale_bits+2, obtain the
exact *outward* fixed-point linear-term bounds using
`floor(n*2^(scale_bits+2-k)/5)` and
`ceil(n*2^(scale_bits+2-k)/5)`, including for negative n.
Bound the quadratic Taylor remainder by the integer
`2^(scale_bits-work)`, as
`|R| <= delta^2 <= 2^-work` follows rigorously from
the preexisting neighborhood test and the second-derivative
bound on [1/4,3/4]. Add/subtract these certified
integer bounds, then convert the two final integer
endpoints to dyadic Fraction objects. Preserve the
original exact rational implementation as fallback
when the denominator is not a power of two or its
scale exceeds the cached fixed-point capacity.
Return the negated enclosure for negative x, exactly
as before; do not approximate the center with a
rounded floating scalar.

The new arithmetic uses GMP's existing integer operations,
but no new MPFR dependency or new public API. For eventual
GMP removal, the arithmetic is designed as quotient by
small unsigned 5, exact shifts and directed integer
rounding: a portable limb implementation can replace
the temporary GMP bridge later.

**Independent mathematical oracle tests.** 16 fixed exact
dyadic result literals now cover atan(+/-[1/2 +/- 2^-(p-3)])
for p=53,128,256,1024 and nearest-even, downward, upward,
toward-zero modes, including signed directed rounding.
Expected values were computed offline using only Python
`fractions.Fraction` and integer operations, NOT FloatAP
or MPFR. The certified rational enclosure for atan(1/2)
comes from an alternating Taylor sum (terms through
n=599 give the lower bound; including n=600 gives
the upper bound, interval width < 2^-1200).
The exact angle-addition identity
`atan(y)=atan(1/2)+atan((2*y-1)/(y+2))` reduces
each offset to a tiny rational r. The alternating
bound between `r` and `r-r^3/3` provides exact
lower/upper rational bounds for the small correction.
Each pair of rational endpoints was independently
rounded to all four target modes using integer quotient,
remainder and parity; the endpoints agreed on every
literal. Test comparisons are exact Rational equalities,
not comparisons to MPFR or coarse tolerances. Existing
domain-diverse tests and benchmark 3456-result post-timing
verifier are unchanged.

**Performance/acceptance gate.** No user-host timing
has yet been obtained for this P4 phase 3 change.
Unix, Windows and Coverage must pass first.
Then rerun the exact 1000x11 four-precision
benchmark and confirm 256 input / 3456 result
equality, plus absolute atan time and all per-precision
worst/geomean/total ratios. No CMake/Windows build
settings or benchmark workload are changed.

## P4 phase 3 build repair: avoid rounding-name shadowing (17 October 2026)

The user supplied an AppleClang `-Wshadow -Werror` build failure in
`tests/test_float_rounding.cpp` (lines 632-634 after commit
`c5e77c19`). The 16 newly added independently certified atan literal
cases declared local variables `near`, `down`, `up`, which shadow
namespace-scope `Ariadne::near`, `Ariadne::down`,
`Ariadne::up` in `numeric/rounding.hpp`. Rename the local
results to `nearest_result`, `downward_result`,
`upward_result`, and `toward_zero_result`. The exact expected
literals, mathematical comparisons, number of cases and modes, and
numerical implementation are unchanged. This addresses the supplied
compiler diagnostic only; full Unix/Windows/Coverage success and the
post-timing benchmark result verification still need confirmation.
The newly created future `remove-gmp` branch is not modified.

## Companion GMP-free core design (18 October 2026)

This journal is retained as the authority for FloatAP precision,
correct rounding, independently certified literals and per-precision
performance. On branch `remove-gmp`, it is maintained **alongside**
[`integer_rational_design.md`](integer_rational_design.md), which
owns the shared `BigUnsigned`, signed `Integer`, finite/extended
`Rational`, the transition away from `mpz_t`/`mpq_t`, and planned
Integer-vs-`mpz_t` and Rational-vs-`mpq_t` benchmarks. The branch
was forked from the all-green `remove-mpfr` commit `9823a3b1`.

Shared-limb changes that alter FloatAP behavior, runtime cost,
rounding certification, or its benchmark must update **both journals
in the same commit**. GMP and MPFR may remain in testing/benchmark
backends during the migration; no GMP-independent FloatAP status is
claimed until its production representation and algorithms cease to
call GMP. Preserve the 1000x11 FloatMP/FloatAP workload, exact
post-timing equality verification, and all existing independent
mathematical literal tests.

## G1 FloatAP-first: initial native unsigned significand kernel (18 October 2026)

**User priority change:** Complete the GMP removal from **FloatAP first**;
only afterwards examine whether the implementation can be refactored
to serve `Integer` and `Rational`. The earlier shared
`BigUnsigned`-first migration roadmap is **superseded**.
Do not impose a general-purpose numeric facade or shared allocator
on FloatAP's critical paths to manufacture reuse.

**Initial isolated implementation:**
`src/floatap_native_unsigned.hpp` introduces
`Ariadne::FloatAPNative::Unsigned`, a private, GMP-free
little-endian unsigned integer with 32-bit limbs, four inline limbs,
dynamic growth, canonical zero, bit manipulation, bit length,
comparison, add/sub with carry and borrow, shifts, multiplication,
and exact base-2^32 long division returning quotient and remainder.
All arithmetic uses standard `uint32_t`/`uint64_t`, not GMP,
MPFR or compiler-specific `__int128`. The type has arbitrary
size rather than the fixed 18-limb ceiling of the earlier
`FloatAPSmallDivision::Unsigned` optimization.
This initial version is **not** yet used by FloatAP itself;
its numerical representation, arithmetic and benchmark performance
remain unchanged. In particular, FloatAP still stores its
significand in `mpz_t` and still uses `mpq_t` for certified
intervals and exact conversions. No GMP-free FloatAP claim.

**Correctness tests:** add the functional CTest target
`test_floatap_native_unsigned` with fixed independent literal
oracles (32-bit carry, 64-bit boundary product, shifts, quotient/
remainder, zero, underflow and >1024-bit growth), plus 120
deterministic unsigned operand pairs up to 2048 bits cross-checked
against GMP `mpz_*` for add/sub/mul, multiple shifts and divrem.
These use a GMP backend **only as a test oracle**, outside
the header implementation. The test explicitly checks
`a=q*b+r` and `0<=r<b`. The benchmark targets remain separate
from CTest. Respect Ariadne GPL headers, test macros and formatting.

**Next stage:** after Unix, Windows and Coverage are green,
introduce a small, explicit GMP import/export adapter and
migrate a carefully selected FloatAP hot path to the native
significand, measuring absolute timing, the 53/128/256/1024-bit
worst/geomean/sum ratios and exact 3456-result benchmark check.
For the full removal, FloatAP construction from legacy
`Integer`/`Rational`, output `Rational` conversions and
certified transcendental `Fraction` intervals also require
GMP-free replacements or appropriate boundary adapters.
The old GMP-backed `Integer` and `Rational` may remain
temporarily in production outside FloatAP while the FloatAP
kernel is migrated. Do not conflate source independence
with removing GMP from every linked object.
## G1 FloatAP-native unsigned MSVC identifier correction (18 October 2026)

The first native unsigned kernel built and passed the Unix and Coverage
jobs at `dcab9609`, but its Windows/MSVC build failed in
`src/floatap_native_unsigned.hpp`, starting at the declaration
`std::array<Limb,inline_capacity> _inline{}` (line 57). Microsoft's
compiler supports `_inline` as a nonstandard synonym of the
`inline` keyword; use of that spelling as a field name conflicts
with the MSVC parser and cascades into errors in the member functions.
Rename the private array field consistently from `_inline` to
`_small_limbs`. No arithmetic, layout capacity, algorithm, test
cases, benchmark workload, public interface or rounding behavior
is changed. Recheck Unix, Windows and Coverage before claiming
portable build success.

## G1 FloatAP-native exact ratio path above 1024 bits (18 October 2026)

The newly introduced `floatap_native_ratio.hpp` implements GMP-free
quotient/remainder rounding for strictly positive unsigned magnitudes,
with correct nearest-even and directed modes, explicit binary exponent,
and a separate carry flag to preserve FloatAP's max-exponent infinity
behavior. `floatap.cpp` now imports legacy GMP values to
`FloatAPNative::Unsigned` and exports the rounded significand only
at its boundary, and **executes the native exact ratio calculation
when precision > 1024 bits**. The standard 53/128/256/1024-bit
benchmark paths and all other GMP-backed FloatAP operations remain
unchanged. This is a first integration milestone, **not yet removal
of GMP from FloatAP**; tests and same-host performance must be
verified before broadening the replacement. Exact unsigned ratio
oracles cover 1/3 at 2,53,128,1024,1100,1536 bits in four modes
and both signs, plus half-way nearest-even carry cases. The core
has no GMP references; conversion adapters do, temporarily.
## FloatAP-first consolidation: implementation in floatap.cpp (18 October 2026)

At the user's direction, the experimental unsigned arithmetic and
correctly-rounded rational division are **implementation details of
`floatap.cpp`**, not separate headers or an additional public numeric
type. Consolidated the old `floatap_small_division.hpp`,
`floatap_native_unsigned.hpp` and `floatap_native_ratio.hpp`
into private translation-unit namespaces, and deleted the three
standalone headers. `floatap.hpp` remains the sole FloatAP API/header.
This is a structural correction preparing for a direct future
replacement of `_significand` and the remaining GMP calls.

Removed `test_floatap_native_unsigned.cpp`, which depended on the
private representation, from the test tree and CTest. Added
mathematically expected results for ±1/3 at 1100/1536 bits and all
four rounding modes, plus FloatMP comparisons for wide non-dyadic
operands, to `test_float_rounding.cpp`. These exercise FloatAP's
*public* API and avoid coupling the test suite to the internals.
The 53/128/256/1024-bit benchmark remains unchanged, and all
four Integer/Rational benchmarks remain separate from CTest.
The earlier 645b2ed9 changes failed compilation on all CI
platforms, so CI validation is required for this consolidation
before extending the native implementation. No claim of total
GMP independence is made.

## G1 FloatAP native-ratio AppleClang Wshadow fix (18 October 2026)

At commit `a1ece3e3`, macOS AppleClang 21 with `-Wshadow -Werror`
reports that the local result variable `rounded` in
`src/floatap.cpp` (inside the precision > 1024 native ratio path)
shadows namespace-scope `Ariadne::rounded` declared by
`numeric/rounding.hpp`. Rename that local result to
`rounded_ratio` and update its three accesses. No changes to
the ratio algorithm, rounding modes, exponent carry, expected test
values, precision thresholds, or benchmark workloads. The fix needs
Unix, Windows and Coverage validation before further migration.

## G1 FloatAP native division at all precisions (18 October 2026)

The verified baseline for this step is `ae165e0b`, green on Unix,
Windows and Coverage. Replaced the generic GMP
`mpz_tdiv_qr` division and rational rounding fallback in
`FloatAP::from_ratio` with the GMP-free
`FloatAPNative::round_ratio` unsigned algorithm for **all** precisions.
Keep the existing allocation-free `FloatAPSmallDivision` fast path
for inputs that fit its <=64-bit optimization; its capacity failures
now fall through to the unrestricted native unsigned algorithm,
not GMP division. Keep the denominator=1 scaled-integer fast
path unchanged. GMP `mpz_import`/`mpz_export` still bridge the
legacy significand and input Integer/Rational representation;
`FloatAP` is not yet GMP-independent as an object or as a whole.

Extended public-API regression of exactly rounded ±1/3 from
only 1100/1536 bits to precisions
2,7,31,32,33,53,64,65,127,128,129,255,256,257,1023,1024,
1100 and 1536 bits, with four rounding modes and both signs.
The expected dyadic mantissas come from exact Integer division,
independently of FloatAP and FloatMP, and results are checked
against FloatMP. Added explicitly expected 5/4, 7/4,
-5/4 and -7/4 tie-to-even and directed-rounding cases at
2 bits, including significand carry. Keep the existing
independent certification tests, benchmark input workload,
the 256 input/3456 operation post-timing equality checks
and all other tests unchanged. Benchmarks stay outside CTest.

**Performance requirement:** collect the unchanged MacBook
1000x11 FloatMP/FloatAP benchmark and compare each
operation/precision ratio to the user-provided
*beginning-of-remove-gmp* MacBook baseline (overall
geomean 0.951x, total-time 0.404x and worst 2.971x;
per-precision 53: 0.966/1.266/2.378,
128: 1.242/1.020/2.971,
256: 0.960/0.535/2.329,
1024: 0.711/0.224/2.143, in
geomean/total/worst order). This change may affect
native division performance at 53/128/256/1024 bits;
do not claim an improvement without same-host measurements,
and do not change the benchmark workload to hide a regression.
CI status and MacBook timing are not yet established
for this commit.

## G1 wide-limb native division specialization (18 October 2026)

The user's first MacBook 1000x11 report after the all-native
`from_ratio` conversion (`96d571e4`) verifies the full
256-input / 3456-result post-timing comparisons but reveals
material regressions in division. In absolute FloatAP times,
`div_up` at 128 bits increased from 129 to 289 ns, at
256 bits from 158 to 435 ns, and at 1024 bits from 244 to
1654 ns. `construct_rational` similarly slowed because it
uses the same exact quotient/remainder routine. The ALL
slowdown worsened from 0.951x geomean / 0.404x weighted
total / 2.971x worst at the beginning-of-branch MacBook
baseline to 1.019x / 0.459x / 9.739x. MPFR timing on the
53-bit group varied materially between runs, so track
absolute times and per-operation data as well as ratios.

**Follow-up implementation:** keep all arithmetic in
`floatap.cpp`. On GCC/Clang with compiler-supported 128-bit
intermediates, dispatch large denominators to a native
base-2^64 Knuth quotient/remainder loop to halve the
number of quotient and divisor limbs. It stores exact
multiplication/subtraction in `unsigned __int128` and
corrects overestimates before returning an exact quotient
and remainder. GCC/Clang pedantic warnings are suppressed
only on the extension typedef; the implementation remains
portable by preserving the existing base-2^32 fallback
for MSVC and small denominators. No general-purpose
`BigUnsigned` class or new header is added.

**Independent validation before integration:** the
isolated wide division algorithm passed 16,640
deterministic GMP cross-checks spanning 1–4096-bit
numerators and 1–2048-bit denominators with exact
quotient/remainder equality, under strict Clang C++20
warning flags. A local x86 Linux microprobe gave around
730 ns for 1024-bit wide divrem; this cannot be
compared quantitatively with the ARM MacBook benchmark,
so no MacBook speedup is claimed yet. Add public-API
FloatAP/FloatMP differential cases for non-dyadic
53/128/256/1024/1536-bit operands, both signs and all
four rounding modes; retain independent fixed 1/3 and
nearest-even tests. Benchmark workload and CTest
separation are unchanged.

**Next validation:** cross-platform Unix/Windows/Coverage,
then the identical MacBook 1000x11 comparison. Do not
change the numerical workload, precision groups, or
rounding certification to improve timings. The production
FloatAP class still stores `mpz_t`; this is only a
division performance milestone, not complete GMP removal.

## G1 wide-division scratch optimization: one GMP-free path at a time (18 October 2026)

**Scope and checkpoint.** Continue optimizing the already-introduced
GMP-free `FloatAPNative::round_ratio` / `divrem_wide` path before
replacing the next FloatAP arithmetic operation. Do not treat removal
of GMP from the division **kernel** as removal from `FloatAP` as a
whole: `mpz_t` remains the owning significand, and import/export
still bridges legacy `Integer` / `Rational`.

**Frozen benchmark comparison supplied by the user (same 1000 iterations,
11 samples, 27 operations x 4 precisions).** The *better* supplied
report is the comparison baseline; its SHA was not supplied. The
*worse* supplied report is attributed to parent commit
`8f488f593f7e8f5d701cf85c00c9695a5172b457`.
Do not confuse this baseline with the earlier beginning-of-branch
baseline 0.951/0.404/2.971 documented above.

| FloatAP measurement | Better baseline | Worse parent commit |
| --- | ---: | ---: |
| ALL geomean FloatAP/FloatMP | 0.912x | 1.010x |
| ALL sum-weighted slowdown | 0.407x | 0.448x |
| ALL worst slowdown | 2.994x | 7.471x |
| 128-bit div_up (ns) | 129.375 | 339.625 |
| 256-bit div_up (ns) | 154.083 | 467.250 |
| 1024-bit div_up (ns) | 256.875 | 1188.334 |
| 128-bit construct_rational (ns) | 266.959 | 535.916 |
| 256-bit construct_rational (ns) | 291.583 | 598.375 |
| 1024-bit construct_rational (ns) | 484.375 | 1203.500 |

Both reports state PASS on 256 paired inputs and 3456 paired results.
They do **not** establish that this change restores baseline performance.

**Implementation.** The 64-bit Knuth `divrem_wide` previously built
packed `a` and `b`, normalized `u` and `v`, and separate
`q` and remainder dynamic word arrays (six vector objects with
dynamic buffers for these operand sizes). The optimized implementation
packs directly into normalized `u` and `v`, uses a single quotient
buffer, and denormalizes the remainder in place in `u`.
Each of the three buffers uses 64 inline 64-bit words with vector
fallback for larger operands. The exact quotient estimation,
correction, remainder and rounding semantics are unchanged;
the portable base-2^32 kernel on MSVC is unchanged.
No new public type or header, global scratch, hidden precision cap
or GMP dependency is introduced.

**Completed isolated validation (NOT integrated Ariadne validation).**
Using Clang 17/C++20 `-O3 -Wall -Wextra -Wpedantic` and GMP as
an independent oracle, a standalone kernel probe passed 35,000
deterministic randomized exact quotient/remainder comparisons among
the old kernel, new kernel and `mpz_tdiv_qr`. Input divisors ranged
from 4 to 259 32-bit words; dividends had 0 to 256 additional
32-bit words, exercising both inline and dynamic fallback buffers.
A local x86_64 Linux kernel-only timing run (not the user's Mac)
observed the following **non-portable and noisy** examples:

| Input family | Old wide kernel (ns) | Reduced-scratch kernel (ns) |
| --- | ---: | ---: |
| 128-bit divisor | 120-138 | 65-69 |
| 256-bit divisor | 147-150 | 93-116 |
| 1024-bit divisor | 505-516 | 439-473 |
| 2048-bit divisor | 1294-1363 | 1306-1414 |

These samples are only a reason to try the integration, not a
FloatAP speedup claim: they exclude the public API, rounding,
`mpz_import`, `mpz_export`, and real workload distribution.
No complete CTest/CI or same-host MacBook 1000x11 report has been
run at this checkpoint. The source change therefore **remains
unvalidated as a baseline recovery**.

**Rejected alternatives at this checkpoint.** An exploratory
two-64-bit-division quotient estimator and a normalized single-limb
reciprocal quotient estimator were each tested against 35,000 GMP
divisions. Neither was a consistent gain at 1024-4096-bit divisors
on x86_64, so neither is integrated. Such architecture-sensitive
tuning needs an Apple arm64 measurement rather than assumptions
about the compiler's 128/64 division code generation. Relevant
prior art: GNU MP division algorithms
(https://gmplib.org/manual/Division-Algorithms.html) and
the correctly-rounded MPFR benchmark contract
(https://www.mpfr.org/algo.html).

**Next gate.** Run regular public FloatAP/FloatMP regression, edge,
rounding and CI tests; then repeat exactly the frozen 1000x11
MacBook benchmark for the same 53/128/256/1024-bit inputs and
compare **absolute medians per operation** and geomean/sum/worst
ratios. If 1024-bit division remains materially below the baseline,
profile quotient-digit estimation vs. `Unsigned` import/export
and scratch traffic on Apple arm64; investigate preinversion and
native 64-bit owning significands on measured evidence.
Do not proceed to a different GMP-free substitution until this
division gate has been assessed, or a quantified exception is
documented.

## G1 capacity tuning after same-Mac scratch optimization benchmark (18 October 2026)

**Measured checkpoint supplied by the user.** Commit
`a944c9814aba106b40e7f9bc55dc1571cd7c0d91` (three local
wide-division scratch buffers instead of six heap buffers) has now
been measured with the unchanged `1000 11` benchmark. It passes
verification of 256 paired inputs and 3456 operation results.
The following are **FloatAP absolute medians in nanoseconds**;
the first column is the previously supplied better baseline,
the second is its worse pre-scratch successor `8f488f5`,
and the third is the measured `a944c981` result.

| Operation | Bits | Baseline | Pre-scratch | Scratch optimized |
| --- | ---: | ---: | ---: | ---: |
| div_up | 128 | 129.375 | 339.625 | 246.417 |
| div_up | 256 | 154.083 | 467.250 | 384.042 |
| div_up | 1024 | 256.875 | 1188.334 | 1094.542 |
| div_nearest | 1024 | 289.375 | 1226.667 | 1121.958 |
| construct_rational | 128 | 266.959 | 535.916 | 444.000 |
| construct_rational | 256 | 291.583 | 598.375 | 511.042 |
| construct_rational | 1024 | 484.375 | 1203.500 | 1086.209 |
| ALL geomean slowdown | -- | 0.912x | 1.010x | 0.985x |
| ALL cost-weighted slowdown | -- | 0.407x | 0.448x | 0.437x |
| ALL worst slowdown | -- | 2.994x | 7.471x | 6.580x |

These measurements demonstrate a real recovery, **not** restoration
of baseline. Large division and rational construction are still the
dominant *regressed* operations at 128/256/1024 bits; elementary
functions mostly remain close to their earlier medians.
Checks for `a944c981`: Windows release, AppleClang/macOS release
and coverage, GCC/Clang/Linux release and coverage, and Linux
ThreadSanitizer reported success; Codecov project and patch checks
reported failure (do not silently describe the entire CI as green).

**Next isolated change: native `Unsigned` inline capacity 4 -> 16.**
The old four 32-bit inline limbs cover only up to 128-bit magnitudes.
Normalizing and shifting the numerator in `round_ratio` frequently
increases the number of words, triggering allocation even when the
original inputs are small. Raise only the private inline capacity
to sixteen 32-bit limbs (512 bits); preserve `_large` for all
other sizes. All arithmetic, rounding decisions, exponent handling,
portable 32-bit Knuth division and 64-bit wide division are unchanged.

**Exploratory same-host kernel-adapter microprobe, not Ariadne.**
A locally compiled synthetic `Unsigned`-plus-GMP-import/export
probe ran repeatable 120,000-iteration cases and indicated
(with substantial platform/workload limitations):

| Input bits | Four inline limbs | Sixteen inline limbs | Sixty-four inline limbs |
| --- | ---: | ---: | ---: |
| 53 | 138.79 ns | 111.09 ns | 127.98 ns |
| 128 | 165.59 ns | 134.88 ns | 158.45 ns |
| 256 | 242.71 ns | 216.25 ns | 201.07 ns |
| 1024 | 619.96 ns | 645.10 ns | 602.62 ns |

This probe is **not a complete `round_ratio` test or performance
claim**; especially at 1024 bits, the sixteen-limb choice may
regress performance. Sixteen is a deliberately smaller and more
conservative first experiment for 53/128/256, avoiding wide-object
copy overhead. Re-run the identical MacBook 1000x11 public benchmark
before keeping or adjusting the capacity.

**Additional algorithm investigation not integrated:** a
normalized-leading-limb reciprocal computed once per wide division
can replace many 128/64 estimates by products and small corrections,
as in Granlund-Montgomery preinversion. An independent C++/GMP
prototype matched 80,000 full quotient/remainder cases and
10,000,000 isolated quotient digit cases. Isolated x86 timings
favored preinversion for large quotients only, and penalized
128-bit divisors. This remains a **candidate for a separately
measured follow-up**, not part of this change. The Apple arm64
behavior may differ because compiler-rt implements 128/64
division in software; reference:
https://github.com/llvm/llvm-project/blob/main/compiler-rt/lib/builtins/udivmodti4.c
and Granlund & Montgomery (PLDI 1994):
https://doi.org/10.1145/178243.178249

**Gate:** obtain standard CTest/CI results, then per-operation
same-Mac timing. Accept only improvements in reproducible aggregate
and critical division/rational cases, with explicit monitoring of
53-bit latency and 1024-bit regressions. If sixteen limbs is not
beneficial, revert rather than carry an unmeasured capacity choice.
Continue optimizing this GMP-free ratio path before moving to a
different GMP removal milestone.

## G1 portability decision: defer architecture-specific division optimization (18 October 2026)

**Decision requested by project owner.** Do **not** introduce an
ARM64-only or x86-only performance path while restoring the GMP-free
FloatAP ratio/division implementation to the frozen baseline.
Specifically, the proposed Apple ARM64-only preinverted
quotient-digit dispatch has **not been merged**, and must not be
added as a separate algorithm merely to improve Mac timing.
Keep the same high-level division/rounding algorithm across
architectures. Compiler-feature fallbacks already required for
portability (e.g. the portable 32-bit limb divider when a native
128-bit intermediate type is unavailable) are a different concern
from microarchitecture-targeted tuning; document their impact.

**Why.** A universal source algorithm alone cannot prove parity
of speed between ARM64 and x86_64. Benchmark both architectures
with the *same deterministic workload* and identical precision
groups before making performance claims. The supplied historical
baseline and subsequent measurements are from the same Mac;
there is **no comparable historical x86_64 baseline supplied**.
Do not compare absolute nanosecond figures across hosts.
Collect an x86_64 FloatMP/FloatAP baseline and record machine,
compiler, optimization flags, GMP/MPFR versions and commit SHA;
then evaluate subsequent portable changes against the frozen
same-host reference on each architecture.

**Acceptance gate before another GMP removal milestone:**
1. At 53, 128, 256 and 1024 bits, every operation affected by
   the current GMP-free ratio/division transition (at minimum
   `construct_rational`, `div_nearest`, `div_down`,
   `div_up`) must reach or beat the frozen same-Mac FloatAP
   absolute median, allowing only explicitly measured noise.
   Overall geomean or time-weighted slowdown alone does not
   satisfy the gate; track FloatMP ratios and unaffected operations
   as safeguards.
2. Preserve all rounding modes, exact quotient/remainder,
   special-value behavior, 256-input/3456-result benchmark
   verification, public regression tests and cross-platform CI.
3. The same portable implementation should be benchmarked on
   x86_64 as well as ARM64, with a separately recorded local
   x86_64 baseline and explicit analysis of any architecture gap.
   Missing x86 measurements must be reported, not inferred from
   Mac results.

**Next algorithm priorities (proposals, not implementations).**
Profile ratio scaling, `Unsigned` normalization/import/export,
quotient digit estimation, multiply-subtract, and allocation
overheads. Explore architecture-neutral normalization and rounding
that avoid unnecessary full-size quotient/remainder materialization;
consider reciprocal estimation only if applied to *both* x86_64
and ARM64 and proven beneficial in the targeted precision bands.
Prefer size/precision-based algorithm selection over CPU-target
selection. Avoid extra code paths without clear, reproducible
cross-platform evidence.

**Status:** The 16-limb experiment at `9343c309` remains the branch
head before this documentation-only decision. The supplied Mac
benchmark has not reached the gate, especially at 1024 bits.
No new division performance optimization is introduced here.

## G1 portable preinverted quotient digits, first complete integration candidate (18 October 2026)

**Target and policy.** The project owner requires *at least* reaching
the frozen better FloatAP baseline for all affected operations at
53/128/256/1024 bits before moving to any other GMP-free substitution.
Use ARM64 Mac as the primary benchmarking host until this gate is met;
only then assess x86-64-specific performance issues. **Do not**
introduce ARM64-only optimizations: the algorithm and its
operand-size-based threshold must be the same across CPU families.

**Current measured checkpoint, immediately before this experiment.**
The unchanged Mac 1000 iterations x 11 samples, 27 operations in four
precision groups, on `9343c309` reported PASS for 256 inputs and
3456 operation outputs. The all-operations geomean slowdown was
0.963x versus a frozen baseline 0.912x; time-weighted slowdown was
0.435x versus 0.407x; worst slowdown 6.542x versus 2.994x.
At 1024 bits: `div_nearest` was 1120.750 ns (baseline 289.375),
`div_up` 1092.833 ns (baseline 256.875), and
`construct_rational` 1088.583 ns (baseline 484.375).
The latest benchmark is **not** evidence of baseline recovery.

**Algorithm change, with no architecture-specific dispatch.**
The existing base-2^64 long divider repeatedly estimates each
quotient digit by dividing a 128-bit prefix by the normalized
leading divisor limb. When the quotient has at least four 64-bit
digits, compute the exact reciprocal
`floor((2^128-1)/top_divisor_limb)` once, and estimate each digit
using a high-half multiplication, additions, and two corrections.
For shorter quotients, retain exact built-in wide division.
The unusual case where the high dividend limb is at least the
top divisor limb uses an exactly saturated digit estimate as
before. Knuth quotient correction, full subtraction, exact
remainder and all rounding modes are unchanged.
Both x86-64 and ARM64 take the same branch for identical operands
when their toolchains support `__int128`, while the longstanding
base-2^32 fallback remains the portability path elsewhere.

**Local independent verification (Linux x86-64, not a Mac performance
measurement).** A standalone algorithm probe ran 80,000 deterministic
complete division cases: both previous and preinverted long
division returned identical exact quotient/remainder, independently
checked against GMP `mpz_tdiv_qr` (divisor 128-3840 bits, numerator
up to 4480 bits). A separate probe tested 3,000,000 normalized
128/64 digit cases plus 150 hand-picked boundary combinations,
with exact quotient/remainder matching built-in `__int128`.
The full-division probe also passed with Clang/ASan+UBSan.
Standalone (GCC 14 `-O3`) approximate timing with equally sized
dividend and divisor-plus-quotient bit widths:
512-bit divisor: old 291 ns, preinverse 259 ns;
1024-bit divisor: old 800 ns, preinverse 730 ns.
At 128/256-bit divisors the probes were approximately tied.
These numbers exclude FloatAP adapters, public rounding, and
Apple ARM64, so they do **not** demonstrate the requested recovery.

**Pending gates.** Validate integration with cross-platform CI
and exhaustive FloatAP public tests. Then benchmark
`1000 11` on the *same ARM64 Mac* with unchanged operands
and all modes. Compare every affected FloatAP operation against
the frozen absolute baseline medians, not merely the aggregate.
If insufficient (especially 1024 bits), profile quotient loop,
remainder construction, normalization, GMP boundary
import/export, and investigate a single-pass rounded-quotient
kernel that avoids materializing a full remainder, while
retaining correct nearest-even, directed rounding and ties.
Do not initiate a different GMP-free substitution while this
baseline gate remains open.

## G1 no-allocation ratio exponent comparison (18 October 2026)

**Frozen performance gate.** Continue optimizing only the existing
GMP-free native ratio/division kernel until affected absolute FloatAP
operation medians on ARM64 meet or beat the user's *better* original
baseline at each precision, within reproducible measurement noise.
Do not introduce architecture-targeted branches; cross-architecture
benchmarking is postponed until the ARM64 performance gate is met.

**New user-supplied same-host benchmark for preinversion commit
`4b517d3b25952cfb32f3cab697d7a2c7ee141f1f`.**
1000 iterations, 11 samples, 27 operations at 53/128/256/1024 bits;
numerical verification PASS on 256 input pairs and 3456 results.
All-geomean FloatAP/FloatMP 0.975x (previous 0.963x,
better frozen baseline 0.912x); weighted slowdown 0.436x
(previous 0.435x, baseline 0.407x); worst 5.729x
(previous 6.542x, baseline 2.994x).
The preinversion reduced the 1024-bit FloatAP division costs
by 16.7-18.9% relative to its parent, but the larger benchmark
still falls substantially short of the original baseline:

| FloatAP operation / bits | Original baseline ns | Previous ns | New report ns |
| --- | ---: | ---: | ---: |
| construct_rational / 128 | 266.959 | 294.834 | 300.709 |
| div_up / 128 | 129.375 | 189.167 | 188.750 |
| construct_rational / 256 | 291.583 | 396.791 | 407.542 |
| div_up / 256 | 154.083 | 297.250 | 271.750 |
| construct_rational / 1024 | 484.375 | 1088.583 | 1178.917 |
| div_nearest / 1024 | 289.375 | 1120.750 | 933.041 |
| div_down / 1024 | 260.167 | 1049.959 | 853.250 |
| div_up / 1024 | 256.875 | 1092.833 | 886.833 |

Do not interpret the 1024-bit `construct_rational` worsening as
proven solely caused by reciprocal estimation: the supplied
samples have visible cross-run variability and the FloatMP
timings drifted too. The change is worth independently measuring.

**This experiment.** In `FloatAPNative::round_ratio`, determining
`floor(log2(n/d))` compared one input with a newly allocated,
shifted full-size temporary operand, even though only the
comparison was needed. Introduce `compare_shifted`, which
computes the shifted words of the positive operand on demand
and compares directly against the other operand without
allocation, copy, or mutation. The logic handles arbitrary
positive shift lengths, 32-bit word boundaries and equal/
unequal bit lengths; `shift=0` forwards to `compare`.
After determining the exponent, quotient/remainder division,
rounding and the final native/GMP adapter are exactly unchanged.
No architecture-specific dispatch or public API change.

**Independent local correctness evidence, not integrated tests.**
A C++20/GMP deterministic randomized probe exercised 140,000
shifted comparisons with arbitrary shifts and 140,000 complete
quotient / remainder-rounding status cases, comparing against
GMP. Both passed. The fused remainder-rounding alternative was
**not included**: isolated timing was mixed and sometimes
slower at 1024 bits even though it was correct.
A vector-limb microprobe suggested shifted comparison is much
cheaper than materializing a shifted vector (usually ~2 ns
versus ~15-55 ns for the synthetic cases); this is an *isolated*
probe, not a FloatAP or ARM64 speed claim.

**Pending gate.** Run full cross-platform compile/test CI.
Then run unchanged `1000 11` benchmark on the same Apple
ARM64 host. Keep the experiment only if affected operations
progress towards their **absolute** original baselines
without a material slowdown elsewhere. In particular, the
1024-bit division and ratio construction require substantially
more than this memory allocation improvement. Do not
proceed to another GMP removal substitution yet.

## G1 architectural correction: native arithmetic belongs to FloatAP (18 October 2026)

**User architecture requirement.** All production FloatAP arithmetic
must be owned by the existing `FloatAP` class
(`Float<ArbitraryPrecision>`), using only its existing header
`include/numeric/floatap.hpp` and source `src/floatap.cpp`.
No parallel `FloatAPNative` or `FloatAPSmallDivision` namespace,
no extra header or CPP. Private implementation details are allowed,
but they are private nested types and static member functions,
not an independently maintained backend.

**Structure-only change.** Transform the two namespace-level
arithmetic groups into `FloatAP` private implementation members:
`FloatAP::SmallUnsigned`, `FloatAP::Unsigned`,
`FloatAP::WideScratch`, `FloatAP::RoundedRatio`,
`FloatAP::RatioRounding`, and their import/export, bit and
limb arithmetic, quotient/remainder, and rounding utilities.
The private types are only forward-declared in
`include/numeric/floatap.hpp`, with all definitions and
arithmetic implemented in the existing `src/floatap.cpp`.
Algorithms, public API, existing benchmark, precision thresholds,
and rounding semantics are intended to remain unchanged.

**Crucial distinction.** This refactoring does **not** make
FloatAP wholly GMP-free: the public class still owns
`mpz_t _significand`, `Integer`/`Rational` conversions
still import/export GMP magnitudes, and other basic operations
still execute GMP calls. Do not present namespace removal
as performance progress or completion of the original migration.

**Validation gate.** Compile on macOS/AppleClang, GCC, Clang,
MSVC and execute existing public FloatAP tests, rounding
oracles, sanitizer and benchmark equality checks. This
is a structural candidate subject to that validation;
no baseline performance benefit is asserted.
After structural regression checks, resume the same
ARM64-first, architecture-neutral speed recovery gate,
particularly for `construct_rational` and `div_*`
at 128, 256 and 1024 bits. No new GMP-free path before
those medians have recovered the frozen original baseline.

## G1 fused numerator scaling into wide division scratch (18 October 2026)

**Benchmark checkpoint received after `bfe14262ffc5438acbd73127c4725c430c3225f5`.**
The same user-supplied 1000-iteration, 11-sample, 27-operation,
53/128/256/1024-bit ARM64 FloatMP/FloatAP benchmark reports
PASS on 256 paired inputs and 3456 operation results. Its
ALL geomean is 1.001x (original better baseline 0.912x),
ALL weighted slowdown is 0.434x (baseline 0.407x),
and ALL worst is 5.657x (baseline 2.994x).
Absolute FloatAP medians at 1024 bits: `div_nearest` 977.625 ns
(baseline 289.375 ns), `div_up` 935.000 ns (baseline
256.875 ns), `construct_rational` 1288.291 ns (baseline
484.375 ns). The checkpoint therefore **fails the recovery gate**.
FloatMP timing drift shows that small run-to-run absolute
differences must not be confused with an algorithmic effect.
At this checkpoint macOS/Linux build-test jobs on `bfe1426`
succeeded; Windows was in progress at the last check;
Codecov project and patch checks failed.

**Bottleneck targeted.** Prior `round_ratio` built an
arbitrarily large `shift_left(n,shift)` temporary, then packed
that temporary into `divrem_wide`'s normalized `u`
scratch buffer. For precision 1024, scaling the numerator
often grows a 1024-bit input to around 2048 bits, so this
duplicates the significant digits across two objects,
incurring allocation, initialization, and copy overhead.

**Implementation.** Keep `FloatAP`'s private wide division
kernel in the existing `src/floatap.cpp` and extend its
private signature with `numerator_shift`. For nonnegative
ratio shifts and denominators with at least four 32-bit
limbs, calculate the normalized dividend directly into
the existing zero-initialized base-2^64 `WideScratch u` by
combining the public ratio's numerator shift and the
divisor's normalization shift. No shifted native numerator
is allocated or imported twice. The exact base-2^64
Knuth digit estimate, correction, borrow, quotient,
denormalized remainder and rounding decisions remain
unchanged. On small divisors, retain the old portable
base-2^32 divider with a shifted temporary; on negative
ratio shifts, shift the denominator as before. Compilers
without `__SIZEOF_INT128__` retain the complete prior
division code path. Dispatch depends on data size /
compiler integer capabilities, never ARM64 versus x86-64.

**Independent isolated validation.** A synthetic C++20
wide-limb (64-bit) probe implemented this combined
normalization and compared exact `q,r` against
`mpz_tdiv_qr` on 100,000 deterministic random input
pairs, with arbitrary divisor bit-length and shifts
including shifts not divisible by 64. All passed.
The same probe passed Clang ASan+UBSan. On the
available x86-64 host, its unsanitized GCC -O3
synthetic division results for a 1024-bit divisor
and approximately 1024-bit quotient were
old 1109.5 ns vs combined 665.8 ns
(128 bits 146.1 vs 107.5;
256 bits 183.7 vs 152.7).
**These are NOT benchmarked FloatAP or ARM64 timings**:
the probe uses `std::vector<uint64_t>`, omits
GMP import/export, construction, precision-aware rounding
and the real Ariadne benchmark's distributions.
A large speedup here does not prove the integrated
FloatAP gate is satisfied.

**Acceptance gate unchanged.** Compile and test on
macOS, Linux, Windows; verify mathematical edge
cases (especially shifted division, nearest-even,
inexact directed rounding, huge precisions) and
3456 benchmark results. Rerun the unchanged
1000x11 workload on the same ARM64 machine.
Retain only if all affected absolute FloatAP medians
move toward the better original frozen baseline.
Do not advance to any other GMP-free substitution
until `construct_rational` and all `div_*` reach or
beat the baseline at each precision; an aggregate
geomean alone cannot meet this criterion.

## Stable benchmark defaults: 20000 iterations and 31 samples (18 October 2026)

**User-requested benchmark-only configuration change.** Increase the
default measured iterations per sample from `1000` to `20000`
and the number of measured samples from `11` to `31`, in all
six backend comparison / exact arithmetic benchmark executables:
`benchmark_compare_multiple_precision_backends`,
`benchmark_compare_double_precision_backends`,
`benchmark_compare_integer_backends`,
`benchmark_compare_rational_backends`,
`benchmark_integer`, and `benchmark_rational`.
Their defaults are declared in the respective two floating
comparison sources and the shared exact arithmetic report header.
Explicit `--iterations/-i` and `--samples/-n` arguments,
including historical positional arguments, continue to override
these defaults. This change does **not** alter workloads,
operands, operation lists, arithmetic algorithms, reference
backends, verification, or reported median computation.
The prior `1000 11` results remain historical records and
should not be directly used as the performance baseline for
new default runs.

**Rebased Apple ARM64 baseline.** The user has measured three
commits on the same host with the new `20000 31` settings:
the beginning-of-branch baseline, previous `bfe1426`,
and current `ba2ef1a`. The reports passed numerical
verification (256 paired inputs, 3456 paired results).
The 1024-bit FloatAP `div_nearest` absolute median
was 292.404 ns (baseline), 945.806 ns (previous), and
893.652 ns (current); `construct_rational` was
508.450 ns, 1198.633 ns, and 1157.854 ns respectively.
The corresponding all-operation FloatAP/FloatMP slowdown
geomeans were 1.060x, 1.109x, and 1.100x.
The latest optimization improved the affected cases
but **did not reach** the required absolute baseline.
Before acting on unusual reports, compare independent
FloatMP control timings and unaffected FloatAP operations;
significant widespread control regressions may indicate
host contention and should trigger a user check **before
further code changes**.

**Acceptance gate unchanged:** match or beat the measured
`20000 31` baseline at every affected operation/precision
before proceeding to another GMP-free substitution. Measure
on idle Apple ARM64 first; investigate x86-specific
performance afterward. Do not loosen numerical oracles.

## G2: sparse divisor products in FloatAP wide division (18 October 2026)

**Problem and scope.** The 20000-iteration, 31-sample Mac ARM64
baseline still outperforms native FloatAP quotient/remainder at
128/256/1024 bits. Many large rationals have sparse divisors
(e.g. 2^(p+1)+small) after normalization. The existing Knuth
64-bit quotient-digit loop multiplies every divisor limb even
when the limb and pending product carry are zero.

**Change.** Inside the existing private FloatAP::divrem_wide
method, select the sparse algorithm when the normalized
divisor contains at least eight 64-bit limbs and at most
one third are nonzero. In that case, when the divisor
limb and multiplication carry both equal zero, skip the
128-bit product and propagate any subtraction borrow
by decrementing the dividend limb; its borrow continues
iff its old value was zero. The general loop, quotient
estimate, corrections, add-back and rounding remain
unchanged. This is data-size and sparseness based,
not architecture-specific. No other GMP paths changed.

**Validation before integration.** A standalone C++20
64-bit Knuth division implementation with the identical
sparse check and borrow handling matched GMP
mpz_tdiv_qr on 100000 deterministic randomized cases,
with both dense and sparse divisors. Clang ASan+UBSan
on the same random cases passed. Synthetic, x86-64
GCC -O3 sparse 1024-bit divisor kernel-only timing
decreased from 1177.04 to 544.21 ns, and at
2048 bits from 4318.08 to 1649.24 ns. These are
not Ariadne or ARM64 measurements; dense cases
showed no systematic improvement. Added public
FloatAP/FloatMP regression tests for sparse ratios,
positive/negative inputs, four rounding modes,
and 128/256/1024/1536-bit precision.

**Acceptance criterion unchanged.** Run CI on
macOS/Linux/Windows and the unchanged benchmark
with 20000 iterations and 31 samples on idle
Apple ARM64. Check FloatMP and unaffected FloatAP
controls for CPU contention before acting.
Require absolute baseline recovery for every
construct_rational and div_nearest/down/up precision
before attempting any additional GMP substitution.

## G3: separate sparse/dense Knuth subtract loops (18 October 2026)

**Motivation.** User-supplied, comparable Apple ARM64 20000x31
benchmarks for the G2 sparse-divisor fast path (`9e5ff09`)
versus `ba2ef1a` showed that sparse
`construct_rational` at 1024 bits improved from
1157.854 ns to 997.198 ns (13.9% faster), whereas
the dense FloatAP divisions regressed:
`div_nearest` 893.652 to 939.519 ns,
`div_down` 829.712 to 871.923 ns and
`div_up` 870.975 to 914.373 ns (approximately
5% slower). FloatMP and unaffected FloatAP timing
controls were stable; all 3456 numerical checks passed.
The G2 implementation tested the runtime `sparse_divisor`
flag on every inner-loop limb, even for dense divisors.

**Change.** In the existing private
`FloatAP::divrem_wide` method, move the runtime
sparse/dense decision to the quotient-digit level.
Use one local C++20 templated lambda instantiated
as two independent inner loops. `if constexpr (Sparse)`
compile-time-elides all zero-limb condition checks
from the dense specialization; the sparse version
retains the exact G2 zero-product/borrow optimization.
A shared body preserves exact product carry, subtraction
borrow, quotient correction, addback and remainder
behavior. No changes to thresholds, divisor
normalization, GMP adapters, rounding modes, work
inputs, public types, or any other arithmetic paths.
This is fully architecture-independent and stays
within the existing `FloatAP` class source.

**Verification.** The existing FloatAP regression
cases for both sparse and dense ratios across limb
boundaries, signed inputs, four rounding modes,
and independent exact mathematical oracles are
unchanged. Repository CI and an independent
Apple ARM64 `20000 31` benchmark are required
before accepting any performance claim for G3.
The goal is to keep G2's `construct_rational`
benefit and recover its roughly 5% division
regression, not to claim recovery of the original
baseline.

**Frozen acceptance gate.** The original `20000 31`
Apple ARM64 baseline still must be matched or
beaten in **absolute median ns** for all of
`construct_rational`, `div_nearest`,
`div_down`, and `div_up` at each precision
53, 128, 256, 1024 bits before any additional
GMP-free substitution. Check FloatMP and
unchanged FloatAP controls for anomalous
CPU contention and stop for user feedback if
those show widespread regressions.

## Architectural division research and revised performance strategy (18 October 2026)

### Why the preceding work is insufficient

On the comparable idle-host Apple ARM64 workload (20000 iterations,
31 median samples), the beginning-of-branch 1024-bit reference
FloatAP medians were 508.450 ns for construct_rational,
292.404 ns for div_nearest, 255.475 ns for div_down and
260.046 ns for div_up. At G2 (commit 9e5ff09),
the respective measurements were 997.198 ns,
939.519 ns, 871.923 ns and 914.373 ns.
G3 (5a66ee2) currently has no validated ARM64 performance
report; do not extrapolate its gain.
The native division remains roughly 3.2 to 3.5 times
slower than the initial baseline at 1024 bits, regardless
of small loop-level improvements. Stop optimizing the
same quadratic kernel incrementally without evidence
of a viable path to baseline parity.

### Source-informed alternatives

1. **GMP-style short/medium-length pre-inverted schoolbook kernel.**
   GMP's basecase long division uses a quotient estimate based
   on normalized high limbs and precomputed inverse, plus a
   dedicated fused multiply-and-subtract primitive. GMP
   mpn_tdiv_qr dispatches according to limb-count thresholds:
   schoolbook below DC_DIV_QR_THRESHOLD, recursive division
   above it, and blockwise inverse/Barrett only at much
   larger sizes. Reference: GNU MP 6.3.0 manual,
   https://gmplib.org/manual/Basecase-Division ,
   https://gmplib.org/manual/Division-Algorithms ,
   and https://fossies.org/dox/gmp-6.3.0/mpn_2generic_2tdiv__qr_8c_source.html .
   Candidate native kernel for 2-32 64-bit limbs: proven
   3-by-2 quotient estimate/preinverse (including its
   initialization cost), full-width 64x64 -> 128 multiply,
   tight fused mulsub with exact carries and add-back,
   normalization and scratch reused for an entire operation.
   Inspect compiler-generated assembly first. The current
   quotient-digit variable has type WideUInt: explicitly
   verify whether its multiplication by divisor words
   produces unnecessary 128-bit arithmetic rather than
   true 64x64->128 word multiplication. A type fix
   alone is not the architectural goal. Use a C++20
   implementation inside FloatAP, no target-CPU branches.

2. **MPFR-inspired quotient approximation with certified rounding.**
   MPFR div.c has fast approximate-quotient paths with
   proved small error and uses extra exact work when
   the approximation is too close to a rounding boundary.
   Source:
   https://fossies.org/linux/mpfr/src/div.c .
   Research the representation-specific problem
   Q=floor((N << shift)/D) for p-bit FloatAP output:
   compute a candidate Q from high limbs or a reciprocal
   at p+guard-bit precision, with an explicit *integer*
   bound on the error. Then obtain a mathematically
   certified decision using exact residual comparisons
   R = scaled_N - Q*D or a provably sound sticky/guard
   test. Correct the candidate until 0 <= R < D.
   Nearest-even requires comparing 2R to D, including
   exact ties and parity; directed modes require sign
   and R!=0. If a short certification is inconclusive,
   fall back to exact full division; never use a
   probabilistic rounding assertion or change the
   existing numerical oracle. Evaluate the cost of
   rare fallback cases separately. This is a new
   output-specific algorithm, not another local
   Knuth-loop tweak.

3. **Burnikel-Ziegler / divide-and-conquer for larger operands.**
   Dividing 2n by n through recursive blocks reduces
   complexity using Karatsuba or higher multiplication,
   but the GMP manual explicitly notes that the
   crossover is *above* the limb count at which
   Karatsuba helps. Consequently this must not be
   presumed faster for the 16-limb divisor typical
   of 1024-bit FloatAP workloads. References:
   https://gmplib.org/manual/Divide-and-Conquer-Division ;
   Burnikel and Ziegler, Fast Recursive Division (1998),
   https://www.mpi-inf.mpg.de/departments/algorithms-complexity/publications/reports .
   Prototype and measure crossover with operand length
   and quotient length, without architecture-specific
   routing.

4. **Newton reciprocal / blockwise Barrett for high precisions.**
   GMP's blockwise Barrett computes a reciprocal of a
   selected divisor segment and uses multiplication
   to generate quotient blocks. It is reserved for
   its largest divisions, not automatically for
   128-1024 bits. Reference:
   https://gmplib.org/manual/Block_002dWise-Barrett-Division .
   Research an integer fixed-point reciprocal with
   controlled error and final exact quotient/remainder
   correction. Avoid floating heuristics that break
   exact rounding. Prototype at 1024, 2048, 4096,
   8192+ bits to measure when the inverse cost is
   amortized. Use Karatsuba multiplication only when
   the measured operand-length crossover warrants it.

5. **Data movement and FloatAP adapter costs.**
   Isolate time in mpz_export for two operands,
   mpz_import for the result, normalization, reciprocal
   setup, the quotient-digit kernel, remainder
   handling and rounding. A GMP-compatible transient
   read-only representation may remove copies during
   the migration, but cannot be mistaken for the
   eventual GMP-free owning significand.
   Never change the existing FloatAP public type or
   benchmark workload just to improve measured scores.

### Execution protocol

A. Benchmark these as *independent algorithmic candidates*
in an isolated harness implementing the same mathematical
division and rounding contract. Record dense, sparse,
carry-intensive, exact-quotient, half-way, and adversarial
inputs at 53/128/256/1024 and beyond. Include operation
setup/import/export cost, not only kernel throughput.
Check exact integer quotient/remainder against GMP, and
public FloatAP/FloatMP results against all existing
independent mathematical oracles. Sanitizers and platform
CI are mandatory.

B. First prototype the full **GMP-style 3-by-2 preinverted,
fused-mulsub basecase** as the realistic short/medium
competitor. In parallel, evaluate the **certified
approximate-quotient** design from MPFR. Do not spend
further commits on single-digit percentage loop tweaks.
Retain the best exact algorithm per *operand/quotient
limb-length* threshold once measured on Apple ARM64,
with platform-independent behavior. BZ and Barrett are
research candidates for sizes where their setup amortizes.

C. A candidate is not considered successful merely because
a standalone x86-64 microbenchmark improved. Use the
unchanged, rebuilt FloatMP/FloatAP 20000x31 suite on
idle Apple ARM64; use FloatMP and unaffected FloatAP
operations as noise sentinels. A significant widespread
control regression invalidates the performance run and
requires user feedback before further coding. The hard
gate remains *baseline parity or better in absolute
median nanoseconds for construct_rational and all three
division rounding modes at all 53/128/256/1024 bits*.
No further GMP-free substitutions until that gate passes.

**Scope of this note:** research and algorithm selection only.
No production source modifications, claimed ARM64 speedups
or widened benchmark coverage are introduced here.

## G4: fused multiply/subtract carry recurrence (18 October 2026)

**Algorithmic improvement.** After researching GNU MP's
`mpn_submul_1` primitive, replace the two separate
product-carry and subtract-borrow chains in
`FloatAP::divrem_wide` with one exact limb recurrence:
`product = digit * divisor_limb + carry`,
`result_limb = old_limb - low(product)`, and
`carry = high(product) + (old_limb < low(product))`.
The returned carry is directly subtracted from the
top dividend limb, followed by the unchanged Knuth
quotient correction/add-back. This is not a
mere sparse-data shortcut: every quotient-digit
multiply/subtract on every wide dense divisor benefits.
The mathematical invariant is that the incoming carry
equals the prior high product word plus the prior
subtraction borrow. No information is lost by merging
them. The carry fits a 64-bit word: its maximum high
product word B-1 is attainable only with product low
word zero, which precludes a subtraction borrow.
Sparse and dense versions remain separate compile-time
specializations inside the existing private FloatAP
method; in sparse mode a zero divisor limb is skipped
only if the *combined* carry equals zero.
No architecture-specific code, new production
files, changed public API, or GMP calls are added.

**Independent pre-integration experiments.**
A standalone division kernel implemented the identical
combined carry recurrence and passed 100000
deterministic randomized exact quotient/remainder
comparisons with GMP `mpz_tdiv_qr`, including
dense and sparse divisors at varying limb lengths.
The same kernel passed 100000 randomized divisions
under Clang AddressSanitizer + UndefinedBehaviorSanitizer.
Kernel-only dense 1024-bit x86-64 GCC -O3
experiment: approximately 1225.71 ns before
and 498.47 ns after (2.46x faster).
Clang -O3: approximately 774.40 ns before
and 552.57 ns after (1.40x faster).
These isolated results are **not** the public
FloatAP benchmark nor evidence of Apple ARM64
parity. The existing exact tests and checks
remain, plus additional FloatAP/FloatMP
differential ratios across 65-1536 bits,
both signs, inverse operands, and all rounding
modes, emphasizing limb/carry boundaries.
See GNU MP 6.3.0 Low-level Functions manual:
https://gmplib.org/manual/Low_002dlevel-Functions ,
especially `mpn_submul_1`, whose documented
output is the product high word plus
borrow-out from the subtraction.

**Acceptance gate.** Check Linux/macOS/Windows CI;
run unchanged 20000x31 benchmark on idle Mac
ARM64. Compare absolute medians against the
beginning-of-branch frozen 20000x31 baseline
for construct_rational and div_nearest/down/up
at 53/128/256/1024 bits. Examine FloatMP
and untouched FloatAP controls first, flag
wide environmental regressions and await
user feedback before further code changes.
Do not consider GMP-free migration at parity
until *all* these cases meet baseline.

### G4 GCC build fix (18 October 2026)

The initial G4 commit `28f1068` passed Clang Linux Release,
AppleClang macOS Release/Coverage, Clang ThreadSanitizer and
Windows Release but failed both GCC 13 Release and GCC 13
Coverage at the build step with warnings-as-errors enabled.
The local fused-submul lambda had introduced a
`carry` variable inside `divrem_wide` where an outer
`carry` variable already existed. GCC with `-Wshadow`
diagnoses this pattern; rename the inner accumulator to
`submul_carry` consistently. This is an identifier-only
correction and leaves all arithmetic and tests unchanged.
CI must confirm the diagnosis; do not claim GCC passing
before its check succeeds.

## G5: end-to-end division-stage feasibility study (18 October 2026)

**Status: diagnostic experiment, not a new algorithmic change in
production.** Following the 20000x31 ARM64 result at `efc8293`,
where 1024-bit absolute FloatAP medians remained 949.810 ns
(`construct_rational`), 881.790 ns (`div_nearest`),
803.321 ns (`div_down`), 842.252 ns (`div_up`) against
the original 508.450, 292.404, 255.475, and 260.046 ns,
stop loop-level micro-optimization. Compiler- and
end-to-end-aware algorithm selection is now required.

A standalone C++20/GMP simulation on x86-64 reproduces
the general *shape* of the benchmark's 64 deterministic
rational input pairs and paired significand divisions.
It implements fused native 64-bit Knuth limbs, tests
integer quotient and remainder exactly against GMP
for every case, and compares nine medians of 10000
operations. At 1024 bits, GCC 14.2 kernel-only division
takes 400.9 ns for rational inputs and 418.7 ns for
significand inputs; full simulated native operations
take 541.6 and 517.5 ns, vs. GMP simulations
393.7 and 293.4 ns. Clang 17 kernel-only results
are 355.9 and 424.6 ns; full native results
448.4 and 509.7 ns, vs. GMP 395.6 and
302.1 ns. Each operand pair's GMP export cost is
about 48-54 ns, considerably less than the cost of
native division itself in these simulations.
These experiments **do not measure built FloatAP**,
have an extra materialized numerator shift compared
with `divrem_wide`, and are on x86-64, not Apple
ARM64. Do not compare their absolute ns to
the frozen 20000x31 Mac baseline. Stage measurements
overlap and are not additive.

**Literature-based algorithmic candidate.** MPFR
`src/mulders.c` implements short approximate
division with proven output error bounds: its
documented high-quotient bound is
`-2(n-1) < N/D - Q <= 4` under the algorithm's
normalization/units. MPFR `src/div.c` attempts
this approach only above `MPFR_DIV_THRESHOLD`,
then falls back to full integer division when
certification is inconclusive. Sources:
https://fossies.org/linux/mpfr/src/mulders.c
and https://fossies.org/linux/mpfr/src/div.c .
Also see Harvey and Zimmermann,
*Short Division of Long Integers*, ARITH-20 2011.

**Next experiment (separate from production):**
build a portable C++20 short-division prototype
with rigorously derived integer bounds on its own
approximate quotient. Form an exact enclosing
interval for the quotient at the target precision
plus guard bits and use the specified rounding mode
on *both* endpoints; accept the result only when
they agree. For nearest-even, preserve parity and
exact-halfway semantics; for directed modes,
preserve inexactness and sign. Fall back to
existing exact `divrem_wide` otherwise.
Record certificate success fraction, exact fallback
rate, end-to-end time including preparation and
correction, and exact oracle conformance on dense,
sparse, and adversarial cases. Short division is a
candidate, not an assumed win at 128-1024 bits.
Separately study tuned basecase quotient-digit
estimation, and measure Burnikel-Ziegler/Newton
crossovers only when justified by operand size.
No new GMP paths, platform-specific implementation
or independently owned FloatAP code.

**Gate:** keep `efc8293` as production
reference until an independently validated
algorithmic improvement appears. No speculative
production edits or changed benchmark inputs.
Only a trustworthy idle Apple ARM64 run at
20000 iterations / 31 samples can establish
baseline recovery. Flag widespread controls
regressions and stop for user feedback.

## G6: exact-residual truncated-limb division feasibility (18 October 2026)

**Outcome: reject for the original 128–1024-bit performance gate;
no production changes.** The standalone
`floatap_exact_residual_short_division_probe.cpp` (research
artifact retained outside the repository) uses the
identity `A=q0*D+E` after truncating `t` low
base-2^64 limbs from both operands. Let
`A=Ah*B^t+Al`, `D=Dh*B^t+Dl`,
`(q0,rh)=divrem(Ah,Dh)`. Then reconstruct
`E=rh*B^t+Al-q0*Dl` with native limbs,
and recover an exact quotient and remainder:
for E>=0 add `floor(E/D)` to q0; for E<0
subtract `ceil((-E)/D)` from q0, with
the exact appropriately adjusted remainder.
All `A=QD+R`, `0<=R<D` invariants hold,
so existing exact rounding can be reused.
GMP is an isolated test oracle only.

**Verification.** 10000 deterministic random GMP
quotient/remainder comparisons on GCC -O3,
10000 on Clang -O3, and 20000 under Clang
ASan/UBSan (two cut policies per case): PASS.
All 64 benchmark-shaped rational and 64
significand-divisor inputs at each of 128,
256, 1024, 2048, 4096, 8192 bits agreed
exactly with the existing native divider
before timing. No claim is made of completed
public FloatAP tests for this prototype.

**Measurements** (standalone x86-64, medians of
7 samples x 5000, not Apple ARM64 FloatAP;
speedup = original fused limb divider / new
adaptive exact-residual divider):
at 1024 bits, rational construction
GCC **0.68x** / Clang **0.73x**;
significand division GCC **0.45x** /
Clang **0.42x**. All are significant
regressions and *not* attributable to
host CPU noise because comparisons were
repeated under both compilers and the
operation-specific pattern is consistent.
Sparse rational construction at 4096 bits
GCC **2.58x** / Clang **2.23x** and
8192 bits GCC **5.94x** / Clang **3.73x**.
This suggests a possible high-precision
sparse-data crossover, but does NOT
justify integration into FloatAP now.

**Literature reconciliation.** MPFR *Algorithms
and Proofs*, section 3.5.2,
https://www.mpfr.org/algorithms.pdf ,
derives `q1-2 < A/D < q1+1` only under
its exact normalization, which keeps
the approximate quotient and high divisor
at the same working precision. Cutting
the divisor below the p-bit target
invalidates that small-error assumption.
MPFR section 3.5.4 and
https://fossies.org/linux/mpfr/src/div.c
select Mulders short division only at
sufficient working size and retain
exact fallback. A faithful
certified-rounding implementation must
carry its own proof and boundaries.

**Next technical priority.** To address
128–1024-bit dense significand division,
prototype the full GMP-style 3-by-2
preinverse quotient-digit estimator and
fused `submul` with benchmarked scratch
and conversion overhead. For precision
far above 1024, separately prototype
Harvey–Zimmermann/Mulders short division.
Avoid implementing the rejected low-limb
truncation in production. Keep the
original frozen 20000x31 ARM64 baseline
and CPU-noise sentinel gates unchanged.

## Architectural reset: independently implemented algorithms, not engines (18 October 2026)

**User constraint.** No integration or source-code copying from external
multiprecision engines (GMP/MPFR/Arb/FLINT or equivalents) as replacements
for FloatAP. Mathematical algorithms and published techniques may be
studied and independently implemented, subject to appropriate checks for
any patent or other applicable intellectual-property constraints.
A software license attached to an existing implementation must not
be mistaken for a restriction on a mathematical idea, nor does citing
a paper authorize verbatim copying of its pseudocode or source.
Preserve FloatMP as reference, not an implementation dependency for FloatAP.
All FloatAP-specific implementation remains inside its existing class
and source files; no new engines or independently owned helpers.

**Reject incremental search.** The 1024-bit FloatAP division baseline
gap still exceeds 3x after several standalone kernel trials. Small
(roughly 3-10 percent) or compiler-specific changes to quotient-digit
estimation, including the independently tried 3-by-2 preinverse
candidate, are not the acceptance target. Do not keep running tiny
GCC/Clang optimization comparisons when they cannot plausibly
recover the required full-operation baseline. No G7 arithmetic
production integration.

**Promising GENERAL-PURPOSE mathematics from published literature.**
Fredrik Johansson, *Efficient implementation of elementary functions
in the medium-precision range*, ARITH 22 (2015),
https://arxiv.org/abs/1410.7176 ,
studies 53 through about 4096-bit exp, sin, cos, log and atan using
low-overhead *fixed-point* operations, functional-equation / table
argument reduction and Smith-style **rectangular splitting**, with
multiple Taylor denominators combined into single-word operations.
The paper's approximately 3x (cos) to 30x (atan) maximum
speedups are versus MPFR **3.1.2 in 2015**, and the implementation
uses GMP mpn; neither these old benchmark numbers nor its
license-covered source code may be transplanted into FloatAP.
This demonstrates a portable *algorithm family* worth independently
investigating, not a predicted gain for our existing atan which
already has input-centered special cases and certified bounds.

Rectangular splitting reduces the number of full-precision
multiplications needed for an N-term suitable series from
O(N) to O(sqrt(N)), exchanging them for cheaper
single-limb operations. For example, at a reduced
atan argument |x|<=1/16 and about 1024 output bits,
a direct alternating Taylor series may require
roughly 128 terms; rectangular splitting reduces
the order of magnitude of full multiplications,
although initialization, argument reduction,
tables, and rigorous error bounds must be counted.
The same mathematical core can serve atan, exp,
log via atanh, sin and cos. Larger precision may
justify **binary splitting and bit-burst** with
functional-equation argument decomposition,
while modern computer arithmetic provides
Karatsuba/Toom/reciprocal methods for big operands.
Source: Brent & Zimmermann, *Modern Computer
Arithmetic*, Cambridge 2010,
https://arxiv.org/abs/1004.4710 .
For bit-burst atan see
https://fredrikj.net/blog/2014/12/faster-high-precision-arctangents/ .

**Where FloatAP is today.** Its general
`FloatAPTrigonometricDetail::small_atan` already
computes an outwardly certified fixed-point
Taylor sum, but obtains limbs and operations
through GMP `mpz_t` objects and updates
powers/terms successively. `Fraction` uses
`mpq_t`, and the interval/refinement wrappers
are already mathematically rigorous. The project
should not discard that correctness logic or
assume published Arb-to-MPFR speedups apply
to FloatAP's already optimized centered atan
case. The major architectural candidate is a
**single native, size-adaptive fixed-point limb
substrate** (multiplication, squaring, truncation,
scaled add/sub, single-limb scaling, exact
outward error accounting) internal to FloatAP,
plus a reusable certified rectangular-splitting
series evaluator. This targets multiple
functions together instead of individually
patching division, atan and exp.

**Acceptance before production integration.**
Derive fresh, independent formulas and error
bounds; prototype one representative general
(not specially centered) atan/exp workload
at 53/128/256/1024 bits, with full
argument reduction + evaluation + round
certification + allocation included.
Compare against existing FloatAP public
results, FloatMP differential controls
and independent rational interval oracles,
including all rounding modes and hard-to-round
boundary cases. Only integrate if improvement
is **substantial and repeatable on full
operations** while preserving correctness;
otherwise leave production code unchanged.
Table contents, if introduced, must be
independently generated and certified, not
copied. The existing mandatory
20,000-iteration / 31-sample idle-ARM64
baseline gate on construct_rational,
div_nearest/down/up at 53/128/256/1024
bits remains in force before other GMP
substitutions. If unaffected controls
show significant regressions, stop for
user feedback before coding or committing.

## G8: independent, certified rectangular-splitting atan prototype (18 October 2026)

**Result: mathematically verified, but REJECTED for production.**
An independent stand-alone C++20 implementation was constructed
without copying or integrating a third-party multiprecision engine.
The candidate implements base-2^32 native unsigned integers,
fixed-point intervals with outward floor/ceiling arithmetic, and
rectangular splitting for general `atan(x)` where
`0 < |x| <= 1/2`. GMP and MPFR are used only inside an
explicitly separate reference/test harness; they are not called by
the candidate algorithm. No FloatAP production code was changed.

**Certificate.** With `t=x^2`, `y=x^4`,
`P_N = sum(y^k/(4k+1))`, `M_N = sum(y^k/(4k+3))` for
`0<=k<N`,
`atan(x) = x*(P_N-t*M_N)+R` with
`0<=R<=x^(4N+1)/(4N+1)` for nonnegative x up to 1/2.
All fixed-point interval multiplications, scalar divisions
and subtractions round outward; `4N*r>=W+8`
for an integer `r>=1` satisfying `x<=2^-r` bounds
the omitted tail by one unit of the `W`-bit fixed-point
scale. At a requested precision and rounding mode, the
result is accepted only if both interval endpoints round
to the same significand/exponent; otherwise the working
precision doubles. Oddness handles negative arguments.
Rectangular grouping by roughly `sqrt(N)` reduces
full multiprecision interval multiplications from 544
to 54 at 1024 bits for input 1/3; scalar divisions
remain 540 per evaluation.

**Correctness:** 112 deterministic inputs (100 randomized
fractions), 53/128/256/1024 bits, positive and
negative signs, all four rounding modes:
**3584 certified and exactly matching MPFR reference
results**, with the native direct-series implementation
independently agreeing. Clang AddressSanitizer and
UndefinedBehaviorSanitizer passed the same 3584 cases.
Scope is restricted to finite nonzero arguments of
magnitude at most 1/2; this is not a FloatAP public
test-suite replacement.

**Paired cost, including input scaling, fixed-point
interval propagation, allocations, and certification.**
A stand-alone GCC -O3 x86-64 run used 31 interleaved
and order-rotated median samples of 32 calls per
candidate. At p=1024 and x=1/3:
native direct Taylor 1222.6 microseconds;
native rectangular splitting 339.5 microseconds;
a GMP fixed-point Taylor comparator mirroring the
existing `FloatAP::small_atan` recurrence 64.1
microseconds. At 53,128,256 bits for the same
argument, native rectangular vs existing-style
comparator was respectively 10.9/1.57,
18.6/2.87 and 33.9/5.86 microseconds.
For p=1024, x=1/2 it was 360.1/110.5
microseconds, while existing FloatAP additionally
has a centered cached atan(1/2) fast path.
The standalone GMP comparator is **not**
a built FloatAP end-to-end benchmark; no raw
timings here may be compared with the frozen
Apple ARM64 20000x31 performance baseline.

**Conclusion.** The 2-4x relative improvement
against our *own naive native direct-series kernel*
is misleading: this independently implemented
rectangular approach is roughly 3-7x slower
than the existing-style fixed-point comparator
at the precisions of interest. The bottleneck
is the unspecialized heap-allocated native
multiprecision limb substrate, not merely
the number of Taylor terms. Johansson's
large historical performance gains combine
argument-reduction tables, efficient low-level
integer multiplication and rectangular splitting;
rectangular splitting alone does not reproduce
those wins. Reject this production integration
and avoid further 5-10 percent micro-tuning.
A different, substantially more efficient common
native arithmetic architecture would need to
prove broad **actual FloatAP operation-level**
wins before any implementation attempt proceeds.

The source and detailed repeatable report were
delivered as conversation artifacts
`floatap_rectangular_atan_prototype_20261018.cpp`
and `floatap_rectangular_atan_report_20261018.md`;
no new production files are introduced. The
previous 20000x31 ARM64 gate for
construct_rational and the three division
rounding modes remains binding.

## Architecture-specific limb kernels: mandatory engineering contract
### Decision and scope (18 October 2026)

The user explicitly authorized separately tuned **ARM64/AArch64**
and **x86-64** implementations for the native multiprecision
substrate, with especially clear documentation. This supersedes
earlier implementation notes that prohibited CPU-specific
code; it **does not** authorize an external arithmetic engine
or weaken any correctness or performance requirement.

**Separation of layers.** Public FloatAP representation,
special values, signed exponent handling, exception behavior,
IEEE-directed rounding and exact-certification rules remain
architecture-independent. Only private low-level limb
operations may differ by ISA. Candidate operations include
full-width 64x64->128 multiplication, multiply-add/subtract
with incoming/outgoing carries, add-with-carry,
subtract-with-borrow, shifts, normalization, short/long
multiplication, squaring, and later division and
fixed-point interval operations. These are prospective
candidates, not implemented features.

**Portable reference.** Each optimized primitive MUST have
a portable C++20 implementation with specified operand
limb ordering, widths, integer overflow/wrap semantics,
carry/borrow outputs and normalized-result invariants.
Retain the portable implementation as an always-buildable,
directly testable baseline. An ARM64/x86 optimized kernel
must implement the **identical mathematical contract**,
including aliasing restrictions, zero lengths, exact
overlap behavior, scratch ownership, and error cases.
Never create two floating-point rounding specifications.

**ISA-specific exploration, not premature assumptions:**
- ARM64/AArch64: evaluate `MUL`/`UMULH` for the low/high
  halves of a 64x64 product and `ADCS`/`SBCS`
  where carrying arithmetic benefits. Verify exact
  instruction semantics and emitted assembly against
  authoritative ISA documentation before enabling.
- x86-64: evaluate baseline integer `MUL`/`ADC`/
  `SBB` and optional BMI2 `MULX` and ADX
  `ADCX`/`ADOX`. BMI2 and ADX are separate
  CPU capabilities and MUST NOT be assumed for
  generic x86-64 or based solely on CPU vendor/name.
- SIMD/NEON/AVX, inline assembly and compiler intrinsics
  are possible only where a real end-to-end advantage
  is demonstrated and the contract remains rigorous.
  Do not assume vectorization speeds carry-dependent
  arithmetic. Review compiler support and possible
  calling-convention constraints separately for
  GCC, Clang/AppleClang and MSVC.
- Other CPUs, unsupported compilers and unsupported
  feature combinations MUST select the portable path.

### Architecture dispatch and portability rules

**Document every supported feature and selection mechanism.**
Prefer compile-time dispatch where the target flags and ABI
guarantee the instructions are available. Optional x86-64
features require reliable feature detection or a
properly isolated feature-target build plus compatible
runtime selection; mere architecture detection is not
feature detection. If runtime selection is used, cache
the immutable decision rather than adding repeated
CPUID checks in inner loops. Do not introduce virtual
backend interfaces or per-scalar dynamic dispatch.
Provide a test-controllable way to force the portable
path and individual supported optimized paths.
No unsupported ISA instruction may execute on a
fallback CPU. Explain minimum supported macOS/Windows/
Linux configurations, compiler flags and link/build
constraints, and how unsupported targets behave.

**One shared arithmetic ownership layer.** Place reusable
unsigned magnitudes, owning storage, limb views, scratch,
portable arithmetic and the optional ARM64/x86-64
primitives in the **internal shared numeric core**,
for eventual consumption by `Integer`, `Rational`
and `FloatAP`. Such genuinely generic code is
**not FloatAP-specific**, and may use new internal
source/header files following Ariadne conventions.
Do not create multiple private multiprecision engines.
Keep FloatAP-specific representation, IEEE rounding,
dyadic interval certification and transcendental
algorithms inside `include/numeric/floatap.hpp`
and/or `src/floatap.cpp`, within the existing
class/private implementation. Avoid the generic
`Rational` canonicalization cost on FloatAP hot paths.
`FloatMP` remains untouched as the reference/control.
The current `mpz_t` owning significand remains
until a separately validated migration decision.

### REQUIRED per-kernel documentation record

Each architecture-specific kernel MUST have a
clear section **in this journal**, a matching
nearby source comment explaining assumptions,
and a concise table comparing it with the
portable reference. Record **all** of the following:

1. **Kernel identifier and use sites:** semantic
   contract, calls, exact inputs/outputs, little-
   endian *limb order* independent of byte order,
   limb width and edge cases.
2. **Architectures and features:** ARM64 or x86-64,
   baseline/optional instruction sets, required
   compiler versions if any, intrinsics or assembly,
   flags and runtime detection. Explain why every
   selected instruction is legal on the target.
3. **Proof obligations:** carry/borrow recurrence,
   maximum accumulator bounds, overflow behavior,
   exact quotient/remainder or product invariants,
   conditions for addback, rounding certification,
   permitted source/destination aliasing, scratch
   initialization and lifetime.
4. **Storage behavior:** temporary buffer layout,
   stack/heap ownership, alignment, maximum sizes,
   in-place cases, failure handling, thread safety,
   and exception guarantees.
5. **Dispatch:** exact selection predicate and
   fallback path, what happens on unsupported
   instruction sets and how to force each path
   in testing; no hidden performance switches
   based on benchmark operand identities.
6. **Correctness evidence:** independent arithmetic
   oracle, exhaustive small cases where feasible,
   randomized/adversarial carry/borrow and
   limb-boundary tests, every FloatAP rounding
   mode and pathological arguments, ASan/UBSan,
   relevant Windows/Linux/macOS CI configurations.
   Verify bit-identical results across paths.
7. **Performance evidence:** host CPU model,
   OS and compiler/version/flags, selected ISA
   path, warmed-up and otherwise controlled run,
   sample counts, absolute ns (not only ratios),
   setup/import/export/certification overhead,
   sizes and data distributions, portable-vs-
   optimized AND original frozen ARM64 baseline.
   Explain crossover thresholds and regressions
   rather than hiding weak cases in a geomean.
8. **Maintenance:** limitations, fallback
   guarantees, tests that fail when variants
   disagree, links to authoritative ISA/algorithm
   sources, journal section and commit that
   introduced the implementation. Do not copy
   external library source or instruction
   sequences blindly; independently derive
   and explain the implementation.

### Implementation sequence and acceptance gates

Build the common limb representation, allocation
and scratch-buffer discipline **first** and retain
a portable, provable arithmetic reference.
Implement/test one end-to-end meaningful
ARM64 and x86-64 kernel family, starting with
high-product + fused multiply/subtract and
small/medium multiplication. Only then assess
adaptive multiplication (e.g. Karatsuba) and
division/rounded fixed-point operations using
complete FloatAP workloads. Choose different
crossover thresholds per architecture **only if
documented and backed by repeatable measurements**.
Do not multiply code paths without substantial,
full-operation benefit.

The strict existing benchmark gate is unchanged:
on an idle Apple ARM64 host, 20000 iterations,
31 samples, verify stable FloatMP and untouched
FloatAP controls before comparing **absolute**
FloatAP median ns for `construct_rational`,
`div_nearest`, `div_down`, `div_up` at
53, 128, 256, 1024 bits to the original baseline.
Do not accept small improvements that leave
the large baseline gap untouched. Significant
broad control regressions require STOP and
user input before any further code changes.
No further GMP operation removal until that
gate passes. Once the architecture-aware
substrate is validated, transcendental
algorithms can reuse it without introducing
a second multiprecision engine.

**Documentation rule is a merge/commit gate:**
no ARM64/x86-64 optimized production code may
be published without updating this document
in the same change and supplying its complete
per-kernel record. The present commit changes
**documentation only**: no architecture-specific
implementation is claimed yet.

### Shared-core clarification (18 October 2026)

The user's new instruction explicitly requires that ARM64/x86-64
kernels and the common multiprecision storage be reusable for
replacing both `mpz_t` and `mpq_t`. The active architecture is

```text
             native unsigned arithmetic core
        (portable + ARM64/x86-64 limb implementations)
                   /       |       \
            Integer     Rational    FloatAP
            signed      exact       significand /
            integer     fraction    certified dyadics
```

`Rational` uses signed integers, exact numerator/denominator,
GCD and canonicalization; `FloatAP` uses the low-level
unsigned operations directly with its own exponent and
rounding rules. SIMD/assembly tuning does not alter
signed-integer division conventions, exact fraction
normalization or FloatAP rounding. This clarification
**supersedes** the earlier restriction confining all
primitive development to the FloatAP source pair.
That restriction continues to apply to
**FloatAP-specific** code only. The user has not
requested immediate production migration of
`Integer`/`Rational`; the initial performance
and correctness target remains FloatAP. Shared-core
changes update this journal and
`integer_rational_design.md` in the same commit.

## G9: restore beginning-of-`remove-gmp` ratio division (18 October 2026)

**Decision.** Before implementing the jointly owned and
architecture-optimized numeric limb substrate, restore exactly
the original wide-ratio/division execution path from
`9823a3b127712ee39e44b634dbaba49e50ea0054`,
the `remove-mpfr` baseline from which `remove-gmp`
began. The source of truth is that historical
`FloatAP::from_ratio` function, not a rollback of
the entire branch. `FloatAP::div` still calls
`from_ratio` with positive significands and sign metadata.

**Active division policy.** `denominator == 1` continues
to use the optimized in-place integer rounding.
For `precision <= 64` and sufficiently small numerator
and denominator (up to 64 bits), retain the
original stack-only `SmallUnsigned` divide-and-round
fast path. For all other nontrivial ratios, use the
historical signed-magnitude-aware **GMP
`mpz_tdiv_qr` quotient/remainder** path with exactly
the original magnitude classification, binary scaling,
nearest-even tie handling, directed rounding by sign,
and carry/exponent overflow handling. In particular,
`construct_rational` and all three division rounding
modes at 128/256/1024 bits now use the historical
wide GMP algorithm. 53-bit inputs may use either the
unchanged small fast path or GMP fallback according
to their actual magnitudes, as they did originally.

**Code scope and cleanup.** Replace only the
`FloatAP::from_ratio` implementation with the
historical function, adapting
`FloatAPSmallDivision::Unsigned` to the current
nested private `SmallUnsigned` spelling (the
initial experimental header has already been
removed on this branch). Remove the now-unused
FloatAP-private large `Unsigned`, `RoundedRatio`,
`WideScratch`, native long division and
GMP import/export adapters and declarations,
rather than leaving dead alternate kernels.
Do not change `FloatAP::div`, special-value
semantics, `FloatMP`, benchmark code,
`Integer`/`Rational`, tests, or other FloatAP
arithmetic/transcendental paths. Existing public
rounding tests exercising those operations remain
active. The previous experiments and their
measurements are retained as historical journal
entries, **not** descriptions of current code.

**Next phase remains unchanged.** A generic shared
portable unsigned core with separately documented
and independently implemented ARM64/x86-64
specializations will eventually serve
`FloatAP`, `Integer`, and `Rational`.
The restored GMP division is temporary; this
commit makes **no claim of completed GMP removal**.

**Verification gate.** Check all GCC/Clang/
AppleClang/Windows jobs and numerical oracles,
then rebuild and measure the unchanged ARM64
20,000 x 31 suite. A previous reference report
from the beginning of the branch contained
1024-bit FloatAP median ns values
`construct_rational=508.450`,
`div_nearest=292.404`, `div_down=255.475`,
`div_up=260.046`. These are historical
values, not measurements of this rollback.
If previously untouched FloatMP or FloatAP
controls regress materially, STOP and request
user review before further implementation.

## BigUnsigned B0 cross-reference (18 October 2026)

The common portable `BigUnsigned` class and first
`compare`/`add`/`sub`/shift methods have been added
in `include/numeric/big_unsigned.hpp` and
`src/big_unsigned.cpp` with independent CTest
and GMP test oracle. This is the shared
`Integer`/`Rational`/`FloatAP` substrate,
not an additional FloatAP-specific engine.
**No FloatAP code paths changed.** Its current
GMP wide-ratio division and the original
20,000x31 ARM64 performance acceptance gate
are unchanged; ARM64 and x86-64 optimized
kernels have **not** been implemented yet.
The authoritative storage and arithmetic
contracts for phase B0 are recorded in
`integer_rational_design.md`.

### Approved conditional architecture files (18 October 2026)

The user approves, **when separate implementations
become warranted**, the exact filenames
`src/big_unsigned_arm64.cpp` and
`src/big_unsigned_x86-64.cpp`. The portable
reference stays in `src/big_unsigned.cpp`,
with the common `BigUnsigned` API in
`include/numeric/big_unsigned.hpp`.
Neither CPU-specific file is currently created.
Their eventual introduction requires the
full per-kernel documentation, independent
testing and meaningful end-to-end performance
evidence specified in the contract above.
The same generic kernels remain available
to `FloatAP`, `Integer` and `Rational`;
no parallel FloatAP-only engine is authorized.
This naming decision makes no functional change.

### B0 portability correction (18 October 2026)

The portable `BigUnsigned` inline-limb member
is renamed from `_inline` to `_inline_limbs`
to address MSVC 14.51 parsing failures.
Whitespace in the new files is brought closer
to Ariadne style. This is not a new arithmetic
algorithm or performance claim; see
`integer_rational_design.md` for diagnostics.

### Shared BigUnsigned B1 multiply contract (18 October 2026)

Portable `BigUnsigned::mul`, `mul_low`, and `square` are available to the eventual FloatAP migration without linking external GMP in this shared arithmetic module. `mul_low` returns the exact low `k` 64-bit limbs, not an approximation or a rounded FloatAP result. `square` currently reuses `mul`. No FloatAP production code or frozen benchmark inputs have changed. The carry bound, independent GMP test oracle and B1 testing are documented in `integer_rational_design.md`.

### BigUnsigned B2 low-level reusable buffers (18 October 2026)

The shared core now defines non-owning read-only/writable 64-bit limb views, per-instance reusable scratch, and `add_into`, `sub_into`, `mul_low_into` kernels for disjoint buffer inputs/outputs. `assign` supports self-referential limb subspans with canonical normalization. Scratch acquisitions invalidate earlier views and zero the acquired range. FloatAP code and its existing GMP wide division are unchanged; the detailed contracts and GMP-based tests live in `integer_rational_design.md`. No ISA specializations are implemented in B2.


### BigUnsigned B3: optional x86-64 wide products (18 October 2026)

The shared unsigned limb layer adds full-product `mul_into` and explicitly callable portable `mul_into_portable`/`mul_low_into_portable` entry points. A separately compiled `src/big_unsigned_x86-64.cpp` may use baseline x86-64 wide multiplication (conditional GNU/Clang 128-bit unsigned type or MSVC `_umul128`), while `src/big_unsigned.cpp` retains the correct portable four-32x32 reference. Both compute the same exact schoolbook result, support modular low truncation, prohibit output/source overlap and allocate no workspace. No BMI2/ADX assumption or runtime ISA probing is introduced. The x86 code is compiled only for supported CPU/compiler pairs; other targets automatically use portable C++20. Explicit portable entry points make forced-path differential tests possible, without modifying FloatAP. No ARM64 kernel or FloatAP end-to-end improvement is claimed. Design proof, testing record, historic preallocated-kernel timings and validation gate are in `integer_rational_design.md`.


### BigUnsigned B4 optional Comba multiplication experiment (18 October 2026)

Experimental `mul_low_into_comba`/`mul_into_comba` entry points expose a column-wise unsigned multiprecision multiplication kernel on x86-64; other architectures retain the portable schoolbook fallback. The ISA-specific source remains `src/big_unsigned_x86-64.cpp`, using baseline wide integer multiplication only, not BMI2/ADX. All three-word Comba accumulator/carry bounds, no-aliasing, exact low/full output contracts, no-allocation behavior, differential GMP checks and benchmark CLI are recorded in `integer_rational_design.md`. The existing B3 selected production path, public `BigUnsigned::mul` and all FloatAP/GMP-wide operations are unchanged. Optional `--include-comba` measures the alternative only after the four frozen B3 benchmark modes; it is not evidence of any end-to-end FloatAP speedup. Full CI and comparable same-host performance evidence remain prerequisites before introducing any new dispatch threshold.


### BigUnsigned B5 forced fused-product diagnostic (18 October 2026)

The shared BigUnsigned arithmetic core now offers `mul_into_fused` and `mul_low_into_fused` as **experimental explicit-only** kernels on x86-64, for measuring a schoolbook inner step expressed as one exact unsigned 128-bit product-plus-destination-plus-carry. The bound `(2^64-1)^2+2(2^64-1)<2^128` guarantees mathematical exactness. GCC/Clang (including AppleClang) can compile the optional conditional `__int128` expression; MSVC x64 uses its existing B3 native kernel, and unsupported architectures use the portable C++20 reference. No new ISA extensions or changes to production dispatch, FloatAP's GMP wide division, operand pools, rounding contracts, or performance baselines. The optional `--include-fused` benchmark mode runs after the established four B3 modes and compares against GMP. Isolated GCC/Clang oracle and performance diagnostics show compiler-dependent results; detailed data and required acceptance gates are in `integer_rational_design.md`.


### BigUnsigned B6 compiler-specific selected multiply (18 October 2026)

Following the B5 user-reported integrated Ubuntu results, the shared BigUnsigned selected production `mul_low_into` uses the exact fused 64x64->128 product-plus-accumulator kernel on **GCC x86-64 only**, with `NUMERIC_BIG_UNSIGNED_GCC_FUSED_MULTIPLY=ON` (default); setting the option `OFF` restores B3 selected schoolbook. Clang/AppleClang/MSVC x86-64 retain B3 selected, and ARM64/other platforms retain portable C++20. No runtime dispatch or nonbaseline x86 ISA assumption. Forced B3 entry points `mul_into_schoolbook`/`mul_low_into_schoolbook` and `--include-schoolbook` in the benchmark allow precise same-build comparisons; explicit Comba/fused tests remain. The immutable input/no-output-alias/full-or-low/zeroed-output/no-internal-allocation contract, exact carry bound and GMP differential oracle continue to apply. User-reported 1024-bit selected B3/fused preallocated medians were 323.861/205.060 ns, compared to 136.262 ns for direct GMP; **this is not an end-to-end FloatAP speedup**. No FloatAP production code, rounding or historical benchmark workload changed. Further gates and full numerical table live in `integer_rational_design.md`.


### BigUnsigned B8 experimental portable Karatsuba (18 October 2026)

A separate portable C++20 `src/big_unsigned_karatsuba.cpp` provides an **explicit-only** Karatsuba full multiply with caller-owned output and caller-owned reusable scratch; it does not replace GCC x86-64 fused production schoolbook or change any FloatAP code. For `min(n,m)>8` limbs and both upper halves nonempty, split at `floor(max(n,m)/2)`, calculate z0/z2 directly into disjoint result slices and use scratch for two exact half-sums and the cross term `(a0+a1)(b0+b1)-z0-z2`. The scratch recurrence and bounds, no-aliasing/no-internal-allocation contracts, GMP differential tests and benchmark controls are in `integer_rational_design.md`. The optional `--include-karatsuba` reports `prealloc_karatsuba`; isolated GCC/Clang probes found the implementation slower than fused schoolbook through 4096 bits and modestly faster near 8192 bits. No measured performance improvement for FloatAP, Integer or Rational may be inferred. No production selection threshold is enabled.


### B8.1: Karatsuba cutoff diagnostic (18 October 2026)

Integrated Ubuntu B8 results were slower than GCC fused schoolbook at 1024 through 8192 bits, with Karatsuba/selected ratios 1.680x, 1.576x, 1.413x and 1.079x. B8.1 adds opt-in recursion cutoff selection (minimum 8, unchanged default 8) and matching workspace sizing, plus --karatsuba-cutoff N with --include-karatsuba. The production selector and FloatAP remain unchanged. See integer_rational_design.md for evidence and test gates.


### BigUnsigned B8.1 Ubuntu cutoff measurements (18 October 2026)

The complete user-run GCC x86-64 balanced-operand sweep of Karatsuba recursion cutoffs 8/16/32/64 found **32 limbs** to be the strongest of these candidates: 2629.425 ns versus selected fused 2971.200 ns at 4096 bits, and 8461.684 ns versus selected fused 12257.133 ns at 8192 bits, measured within the same respective runs. At 2048 bits cutoff 32 gave 815.782 versus selected fused 805.034 ns, so a practical crossover is not yet demonstrated below 4096 bits. The original fused production path stays selected, and no FloatAP performance claim or dispatch change follows. Full four-cutoff table, measurements, CI status and limitations are documented in `integer_rational_design.md`.


### B9 ARM64 fused kernel experiment (18 October 2026)

ARM64 receives an explicit native fused `BigUnsigned::mul_into_fused` kernel, built in `src/big_unsigned_arm64.cpp` using baseline AArch64 multiply-high and fused 128-bit C++ arithmetic. The native ARM64 kernel is initially benchmark-only: `NUMERIC_BIG_UNSIGNED_ARM64_FUSED_MULTIPLY` defaults OFF until real macOS ARM64 benchmark evidence; ON explicitly selects it for production evaluation. Existing GMP differential multiplication tests exercise the native forced fused path, `--include-fused` is accepted on ARM64, x86-only candidates remain x86-only. GCC x86-64 production selection and FloatAP are unchanged. Details and validation conditions are in `integer_rational_design.md`.


### B9 Apple Silicon ARM64 integrated measurement (18 October 2026)

The user measured the native ARM64 fused kernel with production dispatch still **portable** (CMake option OFF) and 20,000x31 original benchmarks. At 256/1024/2048 bits `prealloc_selected` took 26.977/349.917/1379.546 ns; explicit `prealloc_fused` took 16.592/180.819/802.565 ns (1.63x/1.94x/1.72x faster). GMP reused-output took 13.463/95.004/312.606 ns. The 128-bit GMP reused-output measurement appeared unusually high and merits rerun. This is a real user macOS integrated BigUnsigned measurement, **not** proof of FloatAP performance. The next experiment is a correctness-tested ARM64 CMake-option-ON build to measure actual selected owning behavior. See `integer_rational_design.md` for full rows, caveats and gates. No ARM64 default or FloatAP code change yet.


### B9 ARM64 fused dispatch ON: measured owning gains (18 October 2026)

With `NUMERIC_BIG_UNSIGNED_ARM64_FUSED_MULTIPLY=ON`, the user verified that `selected kernel=fused` is reported on Apple Silicon, and that the owning public `new_selected` operation improves relative to the earlier OFF build: 256 bits 44.527 to 32.946 ns, 1024 bits 368.960 to 202.990 ns, 2048 bits 1396.992 to 836.858 ns. Forced/preallocated fused measured 16.592/180.544/802.896 ns at those sizes, and selected preallocated 16.600/180.769/805.240 ns, corroborating the dispatch. The 2048-bit owning GMP `new_mpz` remains faster (336.356 ns). This confirms a shared BigUnsigned improvement on Apple Silicon, **not** a FloatAP end-to-end speedup. macOS ARM64 CTest verification in the option-ON build is still pending; the ARM64 default remains OFF. The full five-size measurements and comparative limitations are in `integer_rational_design.md`.


### B9.1: native ARM64 fused kernel is the default without a toggle (18 October 2026)

Following user-confirmed green macOS ARM64 project tests and integrated owning/preallocated speedups, the ARM64 fused BigUnsigned multiply becomes the unconditional native default on supported GNU/Clang/AppleClang builds. The original `NUMERIC_BIG_UNSIGNED_ARM64_FUSED_MULTIPLY` CMake option and `ARIADNE_BIG_UNSIGNED_ARM64_FUSED` selection macro are removed entirely. The historical B9 OFF/ON sections above are superseded, but the explicitly callable portable reference remains. x86-64 dispatch and FloatAP source code are unchanged; the shared BigUnsigned improvement does not establish a FloatAP end-to-end performance win.


### B10 ARM64 carry-chain candidate (18 October 2026)

BigUnsigned adds an experimental explicit carry-chain ARM64 kernel, exposed to the integrated benchmark with `--include-arm64-carry`, leaving selected fused ARM64 and FloatAP code unchanged. GMP oracle coverage extended; no production promotion or FloatAP speedup is claimed. Details in `integer_rational_design.md`.


### B10.1: shared BigUnsigned source layout (18 October 2026)

Karatsuba (including recursion, workspaces and cutoff APIs) is moved into `src/big_unsigned.cpp`; the separate `src/big_unsigned_karatsuba.cpp` is removed from disk/build. Fused GCC/Clang `unsigned __int128` schoolbook code formerly duplicated in both x86-64 and ARM64 files is now implemented once in `big_unsigned.cpp`. Only genuinely architecture-specific kernels stay in `big_unsigned_x86-64.cpp` (B3/Comba, MSVC native fallback) and `big_unsigned_arm64.cpp` (explicit carry-chain assembly candidate). This refactor does not intentionally alter production kernel selection, Integer/Rational/FloatAP integration, tests or benchmark semantics; CI and native macOS rerun are still required and performance equivalence must be measured rather than assumed.


### B10.1 Apple Silicon ARM64 explicit carry results (18 October 2026)

In the user-run integrated benchmark following source-layout consolidation, the ARM64 carry-chain candidate (`--include-arm64-carry`) measured 154.744 ns at 1024 bits and 571.048 ns at 2048 bits, versus selected fused 182.898 ns and 806.990 ns within the same run (1.18x and 1.41x throughput improvement). Reused-result GMP measured 95.062/320.725 ns. The selected owning BigUnsigned path is still fused; these numbers do **not** measure the effect of switching the public owning path or FloatAP end-to-end performance. Native macOS differential CTest on B10.1 is pending; the ARM64 dispatch remains unchanged. The exact data and acceptance gate are in `integer_rational_design.md`.


### Second B10.1 Apple Silicon carry benchmark (18 October 2026)

In the user-run 256/1024/2048/4096-bit integrated BigUnsigned benchmark, the explicit ARM64 carry kernel took 17.831/153.998/570.777/2263.667 ns versus fused selected 22.617/179.512/796.533/3673.671 ns in the same run. The 1024/2048-bit speedups were reproduced across two user runs; the 256-bit result varied, and 4096 bits newly showed a 1.62x throughput gain versus the currently selected fused kernel. GMP reused-destination at 4096 bits took 986.669 ns. These results do not imply FloatAP end-to-end speedup and do not yet justify dispatch promotion without a reported native macOS ARM64 GMP-oracle CTest. No FloatAP source or production dispatch was modified. See `integer_rational_design.md` for the full measurements.


### B10.2: default ARM64 BigUnsigned multiplication selects explicit carry chain (18 October 2026)

Following two integrated Apple Silicon measurements and user-confirmed native GMP-differential `ctest` green in the prior build, the ARM64 selected BigUnsigned multiply now calls the existing architecture-specific `MUL/UMULH/ADDS/ADC` carry kernel. The shared `__int128` fused implementation and portable reference remain explicitly callable; no CMake flag or extra source file is introduced. x86-64 behavior and FloatAP code stay unchanged. A newly built ARM64 selected-path `ctest` and owning benchmark comparison remain pending; do not claim an end-to-end FloatAP speedup. See `integer_rational_design.md`.


### B10.2 measured selected ARM64 carry in public owning BigUnsigned multiply (18 October 2026)

The user benchmarked a newly compiled B10.2 selected ARM64 carry build. Compared with the previous native ARM64 fused dispatch, owning `new_selected` at 1024 bits improved 202.990 to 172.594 ns (~1.18x) and at 2048 bits improved 836.858 to 605.933 ns (~1.38x); the 53/128/256-bit ownership times did not show clear changes. The new owning BigUnsigned/GMP ratios were 1.580x (1024) and 1.800x (2048) with matched result lifecycle. Forced and selected ARM64 carry matched at 256+ bits, supporting the dispatch; unusually different timings at 53/128 bits merit checking benchmark effects. User native macOS GMP differential CTest was green **before** the selected-kernel switch, with carry forced, but has not been reported anew after promotion. This is BigUnsigned-only benchmark evidence and does not demonstrate FloatAP end-to-end performance. Detailed measurements are in `integer_rational_design.md`.


### B10.2 selected ARM64 carry validation closed (18 October 2026)

The user confirmed the updated **post-dispatch-change** macOS ARM64 `ctest --output-on-failure` fully green. GitHub Actions Unix/Windows/Coverage on the functional ARM64 carry default commit `d01c505` also completed successfully. These results close the previously pending validation gate for production BigUnsigned ARM64 selected carry. The earlier native owning BigUnsigned benchmark and remaining GMP gap are recorded in `integer_rational_design.md`; no FloatAP end-to-end speedup is claimed.


### B11 x86-64 explicit carry-chain experiment (18 October 2026)

BigUnsigned adds a **non-production** x86-64 GNU/Clang `MULQ/ADDQ/ADCQ` schoolbook specialization in the existing architecture-specific source file, with a GMP-differential test and opt-in `--include-x86-64-carry` benchmark. MSVC keeps the native B3 fallback and no new file or CMake flag is introduced. GCC x86-64 selected fused, ARM64 selected carry and FloatAP code are unchanged. An isolated GCC/Clang x86-64 correctness check passed 1.5 million random/edge tuples; project-level CI and user Ubuntu integrated measurements are still pending. See `integer_rational_design.md` for algorithm, constraints and acceptance criteria.


### B11.1 portable BigUnsigned fallback in dedicated backend source (18 October 2026)

Source layout now uses `big_unsigned_fallback.cpp` for the portable C++20 reference kernel; `big_unsigned_x86-64.cpp` and `big_unsigned_arm64.cpp` remain the sole architecture-specialized sources. `big_unsigned.cpp` retains common BigUnsigned methods, Karatsuba, public dispatch, and the unduplicated GCC/Clang `__int128` fused algorithm shared across native targets. The fallback translation unit is **always compiled** for GMP-differential control and benchmarking, while CMake selects at most one native architecture file. This is strictly organizational: selected x86-64 and ARM64 multiplication dispatch and FloatAP operations are unchanged, no algorithm-specific extra .cpp files or CMake flags. CI and host benchmarks must confirm regression-free behavior.


### B11 x86-64 integrated carry benchmark outcome (19 October 2026)

The user benchmarked the B11 opt-in x86-64 `MULQ/ADDQ/ADCQ` multiply against the selected GCC fused kernel on the same Ubuntu machine at 53/128/256/1024/2048/4096 bits (20,000 iterations, 31 samples). Explicit assembly was slower at five widths (the remaining 128-bit case effectively tied), including 939.067 vs 761.449 ns at 2048 bits (+23.3%) and 3156.871 vs 2975.414 ns at 4096 bits (+6.1%). Therefore keep GCC fused selected and explicit x86 carry experimental, without modifying FloatAP. Both the B11 experiment commit and B11.1 three-backend source split passed GitHub Unix/Windows/Coverage CI. Full measurements and caution about owning GMP comparisons are in `integer_rational_design.md`.


### B12 experimental Karatsuba write optimization (19 October 2026)

The existing common-source opt-in Karatsuba kernel was adjusted to avoid redundant output clearing and high-word scanning on shifted cross-term accumulation. Exact arithmetic and public API contracts are unchanged, and new GMP-oracle tests exercise prefilled outputs, odd/asymmetric recursive splits, workspace sentinel guards and cutoffs 8/16/32/64/128. GCC x86-64 fused and native ARM64 carry remain the selected production kernels; FloatAP code and design are not modified. New B12 test/CI and integrated performance measurements remain pending. Full rationale and suggested crossover experiment are in `integer_rational_design.md`.


### B12 Ubuntu integrated cutoff-32 Karatsuba results (19 October 2026)

The user's post-B12 balanced-operand x86-64 benchmark found opt-in preallocated Karatsuba with cutoff 32 faster than selected GCC fused at 2560/3072/3584/4096/8192 bits, but approximately tied at 2048. At 4096 bits Karatsuba measured 2582.760 ns vs fused 3023.645 ns (14.6% less time), and at 8192 bits 8264.658 vs 12387.829 ns (33.3% less time). CI (Unix, Windows, Coverage) for functional B12 commit `014eedd` passed. The opt-in Karatsuba benchmark excludes scratch allocation/sizing and does not model owning public `mul`; no default dispatch change or end-to-end FloatAP speedup can yet be justified. See `integer_rational_design.md` for complete data and fairness qualification.


### B13 shared unsigned primitive benchmark (19 October 2026)

Introduced a benchmark-only `benchmark_big_unsigned_arithmetic` executable (no new arithmetic/backend source): add/sub with both preallocated and owning destinations, compare as nonallocating read-only, and two left/right shifts as owning operations, all against matching GMP lifetimes and independently checked exact results. Production `BigUnsigned`, `FloatAP`, `Integer` and ISA dispatch are unchanged. This prepares measurement before implementing the common unsigned division/remainder kernel. See `integer_rational_design.md` for the benchmark contract and CLI; CI and native measurements pending.


### B14 single-limb division and fallback-source rule (19 October 2026)

The user clarified that generic arithmetic with **no native x86-64/ARM64 specialized counterpart** belongs in `src/big_unsigned.cpp`. Only portable alternatives to ISA-specialized algorithms belong in `big_unsigned_fallback.cpp`. Thus the pure C++20 unsigned single-limb quotient/remainder implementation and tests are in the common source; `big_unsigned_fallback.cpp` is unchanged. GMP differential tests cover 6000 cases, normalization and exceptions. A general multiprecision-divisor division is still needed before native Integer/Rational migration. Production multiplication dispatch and FloatAP remain unchanged.


### B14.1 Windows GMP 64-bit remainder oracle correction (19 October 2026)

The B14 Windows Release build succeeded, but the single-limb division
test compared a 64-bit native remainder against `mpz_get_ui`, whose
`unsigned long` result is only 32 bits on Windows x64. The test now
uses full-width `mpz_import` plus `mpz_cmp` and explicitly covers a
remainder above 32 bits; all existing randomized GMP-differential
checks remain unchanged. Production `BigUnsigned`, `FloatAP`,
architecture dispatch and backend source placement remain unchanged.
Windows and other CI results for this correction are pending. See
`integer_rational_design.md` for the detailed B14.1 record.


### B15 multi-limb shared BigUnsigned division (19 October 2026)

Implemented exact owning and caller-buffer unsigned `divrem` in shared `src/big_unsigned.cpp`, using normalized radix-2^64 quotient-digit estimation with full-product correction, B14 single-limb delegation, and explicit scratch sizing. No change to selected architecture multiplication, portable ISA fallback or FloatAP. Six thousand GMP-differential project test cases and error/buffer regression checks are added. A standalone Linux prototype passed 90,000 GMP comparisons; Windows/Unix/Coverage CI for the integrated changes are pending. See `integer_rational_design.md` B15 for detailed contracts and validation limits.


### B15.1 compile fix: reuse shared checked subtraction (19 October 2026)

B15 introduced a division-local `subtract_in_place` with `noexcept`
that collided with the existing throwing Karatsuba helper in the same
anonymous namespace. AppleClang ARM64 diagnosed the mismatched exception
specifications under `-Werror`. Remove the duplicate and forward-declare
the existing checked helper for division. Its final-borrow guard remains
active for both consumers. No arithmetic, FloatAP, ISA dispatch, public
API, tests or build flags change. B15.1 cross-platform CI is pending;
see `integer_rational_design.md` for details.


### B16 shared BigUnsigned GCD and decimal IO (19 October 2026)

B15.1 cross-platform CI was green. Shared `src/big_unsigned.cpp` gains Euclidean unsigned `gcd`, strict ASCII decimal parsing and canonical unsigned decimal formatting, using existing division and portable 32-bit half-product operations without GMP dependencies. No change to FloatAP, ISA dispatch or fallback source placement. An isolated warning-clean Clang harness passed 5000 limb-format GMP comparisons, 5000 decimal input GMP comparisons and 1100 GCD GMP comparisons. Project CTest and CI for B16 are pending. Consult `integer_rational_design.md` for the contracts and migration sequencing.


### B17 isolated signed-magnitude Integer substrate (19 October 2026)

Following B16 CI green, the native `SignedInteger` class is added
as a separate shared-core consumer of `BigUnsigned`. It implements
canonical sign/magnitude, decimal conversions, signed arithmetic,
truncation-toward-zero divrem with dividend-sign remainder, powers,
shifts and checked 64-bit casts, and has its own GMP-differential
tests. Legacy `Integer::_mpz` and its GMP-dependent consumers remain
unchanged while integration is staged. See `integer_rational_design.md`
for the migration boundaries; B17 CI is pending.


### B17.1 clarification: transient SignedInteger test stage (19 October 2026)

Keep B17 `SignedInteger` for cross-platform GMP-differential verification; do not revert it before transferring validated operations to the existing `Integer`. After the direct `Integer` storage and dependency migration, delete the temporary `SignedInteger` header, source, CTest target and dedicated test (migrating coverage into `test_integer.cpp`). The final architecture must expose only the existing `Integer` as the signed multiprecision type; `BigUnsigned` remains its GMP-free magnitude engine. See `integer_rational_design.md` B17.1 for mandatory completion gates.


### B18 native Integer replacement without modifying FloatAP arithmetic (19 October 2026)

The existing `Integer` class has been refactored in place to canonical `BigUnsigned` magnitude plus sign, removing its permanent GMP storage. `SignedInteger` stays only as temporary GMP-differential verification and must be deleted once tests/benchmarks validate the migrated Integer. The few dependencies on `Integer::_mpz`/`get_mpz()` in GMP-backed Rational/Dyadic/FloatMP/FloatAP use a narrow internal import/export adapter; these other numeric types retain their own representations and core algorithms. The benchmark's signed quotient/remainder GMP control now uses truncation toward zero to match the coherent Integer division identity. Integrated Unix/Windows/Coverage and benchmarks for B18 are not yet verified. See `integer_rational_design.md` B18 for implementation and testing details.


### B19 shared decimal throughput and Integer benchmark correction (19 October 2026)

B18 native ARM64 and x86-64 measurements showed `to_decimal` dominating total time. In the shared `BigUnsigned` decimal conversion only, divide each radix-`2^64` word by `10^9` through two 32-bit halves using portable 64-bit arithmetic, avoiding 64-step bitwise division. Preserve `BigUnsigned::divrem_limb` and general `divrem` unchanged. The Integer-vs-GMP benchmark now compares one native `cmp` against one GMP `mpz_cmp` and adds separate nontrivial signed quotient/remainder rows without changing the original workload. All new values are verified exactly against GMP outside timing. A standalone 800,000-case mathematical oracle passed, but integrated CI and performance outcomes are pending. 4096-bit measurement remains available, though the next developer runs can temporarily specify `--bits 53,128,256,1024,2048`. See integer_rational_design.md B19.


### B19.1 decimal type alias fix (19 October 2026)

Within `namespace Ariadne`, the free function `to_decimal` must refer to the nested limb alias as `BigUnsigned::Limb`. Corrected two declarations/expressions; algorithms and benchmarks unchanged. CI pending.


### B19.2 BigUnsigned comparison access and benchmark count correction (19 October 2026)

Following B19.1 green cross-platform CI and ARM64/x86-64 user
benchmarks, the shared `BigUnsigned::compare` now chooses its
inline/heap limb views once and compares direct digits instead of
repeated `limb()` calls, retaining its exact value ordering and
`noexcept`. A 64-limb equal-high-prefix regression was added.
The benchmark's extra 32-pair division-input check was moved outside
the original 32-element loop so the reported number of distinct
inputs is accurate; all prior comparisons and timed operations
remain intact. No FloatAP algorithm or architecture-specific kernel
changes. Cross-platform CI and measured speedups pending; see
`integer_rational_design.md` B19.2.


### B20 shared 128/64 quotient-digit acceleration (19 October 2026)

From the B19.2 ARM64/x86-64 six-size GMP-differential benchmarks,
nontrivial division is now the dominant common slowdown. B20
modifies the shared `BigUnsigned` 128/64 helper only: use plain
64-bit division for `high==0`; direct `unsigned __int128`
division for GCC/Clang; the MSVC x64 `_udiv128` intrinsic
with `<immintrin.h>`; and retain the pure C++20 bitwise
fallback elsewhere. The helper is confined to `namespace
Ariadne` and serves single-limb division and long-division
quotient estimation. Keep long division, Integer signed semantics,
FloatAP algorithms and the benchmark corpus unchanged. Add
focused GMP oracle boundaries for maximum quotients and all
64 divisor normalization shifts without removing prior cases.
An independent GCC/Clang isolated reference check passed
1,200,252 selected-vs-fallback cases on each compiler.
Integrated Unix/Windows/Coverage CI, native performance and
MSVC verification remain pending; consult `integer_rational_design.md`.


### B21 removal of transient signed-magnitude migration scaffold (19 October 2026)

Following successful B20 cross-platform validation, the B17
temporary independent signed integer type is retired. Its header,
implementation and dedicated CTest executable are deleted, with
their respective CMake source/test registrations removed. Retain
only the production `Integer` public signed type owning
`BigUnsigned` magnitude and sign; continue to use the existing
temporary boundary adapter for still-GMP-backed Rational/Dyadic/
FloatMP/FloatAP. The historical B17 2500-case GMP-differential
corpus, boundary and signed shift semantics are moved to the
production `test_integer`, preserving its B18 1700 independent
cases. Signed shifts are checked using the shared native magnitude
shifts plus sign, and GMP; this does not add an `Integer` shift API
or change FloatAP arithmetic. Project cross-platform CI remains
pending for B21; see `integer_rational_design.md`.


### B22 recursive decimal formatting for large Integer magnitudes (19 October 2026)

Based on B21's cross-platform green CI and native ARM64/x86-64
benchmarks, optimize only `BigUnsigned::to_decimal` for large
magnitudes. Retain the existing portable two-64/32 radix-10^9
formatter up to 16 limbs; for larger values recursively divide
by balanced cached exact powers `10^(9*2^i)`, padding each
low-order decimal chunk to its fixed width. All helpers live
inside `namespace Ariadne`. The powers are local to one
conversion (no shared cache), and the recursion reuses shared
native `divrem`. Preserve `Integer` signed-magnitude,
the retired signed-class state, FloatAP and rational internals,
GMP boundary adapters, parsing and benchmark workloads.
Add 50 direct GMP-oracle power-of-ten and zero-run boundary
cases to `test_big_unsigned` without removing any existing
randomized validation. An isolated alternative Boost-based
algorithm passed 6058 GMP differential cases, but integrated
Unix/Windows/Coverage and native ARM64/x86-64 speedup evidence
remain pending. See `integer_rational_design.md` B22.


### B23 single-output division after B22 (19 October 2026)

B23 makes `Integer::quot` and `Integer::rem` use new
single-output native `BigUnsigned` primitives backed by one
templated internal division core, while preserving the public
two-output `divrem_into` contract and the identical long-division
quotient-digit correction. Unrequested quotient/remainder object
allocation, final unnormalization or quotient stores are skipped
as appropriate; single-limb remainder avoids quotient storage.
No FloatAP, Rational, Dyadic, GMP adapter, arithmetic algorithm,
benchmark definition or B22 decimal-formatting change. Existing
GMP-randomized division tests are extended with single-output
comparisons and boundary checks. Cross-platform CI and measured
performance remain pending. Details: `integer_rational_design.md`.


### B23.1 report label clarification (19 October 2026)

Rename original trivial Integer division benchmark results to
`quotient_trivial` and `remainder_trivial`; use the canonical
`quotient` and `remainder` for the additional meaningful
division workload previously suffixed `_nontrivial`. This alters
only reporting labels in `benchmark_exact_arithmetic_workload.hpp`,
not sample inputs, benchmark order, GMP verification, production
arithmetic or FloatAP. When comparing B23.1 with older reports,
map both old names to their new counterparts; see
`integer_rational_design.md`. B23.1 CI remains pending.


### B24 decimal parsing with optional 128-bit accumulation (19 October 2026)

B24 introduces a wide-radix path only in native
`BigUnsigned::from_decimal` for decimal strings exceeding
309 characters on GCC/Clang targets with `unsigned __int128`.
Parse groups of up to 18 decimal digits and propagate full
64-bit limb carries with one exact 128-bit intermediate.
The prior portable 9-digit base converter, input validation,
and <=309-digit path are unchanged, with portable/MSVC
fallback for all sizes. It does not change decimal output,
signed Integer semantics, FloatAP/Rational/Dyadic algorithms,
GMP boundary bridges, or the integer benchmark workload.
Add 45 independent GMP oracle parsing boundaries and four
long malformed-input cases while retaining all existing tests.
An isolated GCC/Clang/GMP probe passed and suggested possible
large-input performance gains; integrated CI and native
two-host benchmarks are pending. See `integer_rational_design.md`.


### B25 optional radix-10^18 leaf decimal output (19 October 2026)

The B25 change applies only to `BigUnsigned::to_decimal`:
on GCC/Clang with exact `unsigned __int128`, leaf values
of 9 or more 64-bit limbs use repeated 128-bit quotient and
remainder by the constant `10^18`, emitting 18-digit padded
blocks. The unmodified portable base-`10^9` implementation
remains for values of at most 8 limbs and on non-`__int128`
toolchains such as MSVC. B22 recursive split powers and cutoff,
B24 decimal parsing, GMP boundary adapters, Integer sign
semantics, Rational/Dyadic/FloatAP, and the benchmark corpus
are untouched. Add 50 direct GMP-oracle decimal formatting
checks across the 8/9-limb selector and large-number boundaries,
preserving all earlier test coverage. An isolated GMP-validated
x86 leaf-conversion probe suggested improvements but integrated
cross-platform CI and native ARM64/x86-64 benchmark results
are pending. See `integer_rational_design.md` B25.


### B25.1 compiler-native ARM64 decimal-format selector, no CMake flags (19 October 2026)

B25 was beneficial for large `to_decimal` on ARM64, but the
radix-10^18 128-bit leaf division regressed substantially
on x86-64. B25.1 gates the B25 implementation using the
compiler predefined `__aarch64__` target macro plus
`__SIZEOF_INT128__` and GCC/Clang, retaining the B19
radix-10^9 path on x86-64, MSVC, other architectures
and small ARM64 values. This requires **no additions or
changes to CMake, flags or build definitions**. No
FloatAP/Rational/Dyadic, GMP adapter, parsing, recursive
decimal-split or benchmark behavior changes. Add 18
independent GMP-oracle radix-power boundary cases and
retain all preceding tests. Cross-platform CI and host
performance benchmarks for B25.1 remain pending;
see `integer_rational_design.md`.


### B26 baseline x86-64 DIVQ for shared unsigned division (19 October 2026)

To address the measured x86-64 general quotient/remainder
slowdown relative to GMP, B26 specializes the internal
128/64 `divide_wide` helper in `src/big_unsigned.cpp`
for GNU/Clang x86-64 using a single baseline `divq`.
The quotient/remainder input/output registers and
hardware nonoverflow precondition are explicit. Keep
the 64-bit direct path, MSVC `_udiv128`, non-x86
`unsigned __int128` and portable fallback unchanged.
Preserve all arithmetic semantics, normalization,
allocation, correction logic, `FloatAP`, `Rational`,
`Dyadic`, decimal formatting/parsing, GMP bridges and
benchmark workloads. Do not alter CMake or build flags.
Add 2048 exact constrained 128/64 GMP-oracle cases
without removing existing division tests. A standalone
GCC/Clang/GMP checker passed 500072 cases against the
wide C++ implementation and checked an independent
GMP sample; GitHub CI and native performance measurements
for B26 remain pending. See the B26 entry in
`integer_rational_design.md`.


### B27 fused long division and smaller scratch (19 October 2026)

The native multiprecision divider now calculates each
quotient digit by a single pass that multiplies the
normalized divisor by the estimated digit and subtracts
the resulting low limb from the working dividend,
propagating the high product plus borrow. A negative
result triggers decrement/add-back corrections (at most
two for normalized long division). This eliminates
the previously allocated per-digit product buffer and
reduces `divrem_scratch_limbs(n,m)` from
`n+2*m+2` to `n+m+1` for n>=m>=2.
The kernel uses GCC/Clang unsigned 128-bit products,
MSVC `_umul128` or a portable 32-bit half-product
fallback as appropriate, with no new CMake flags
or ISA requirements. Preserve the B26 x86-64 DIVQ
estimate, B23 selective outputs, all Integer semantics,
FloatAP/Rational/Dyadic, GMP adapters, and benchmark
operands. Add 2050 exact GMP-oracle regression fixtures
including two independent double-correction cases;
do not remove older cases. Three independent local
GMP stress configurations (GCC wide, GCC portable
and Clang wide), each with 256400 divisions, passed.
Integrated Unix/Windows/Coverage CI and native
two-host benchmarks remain pending. See B27 in
`integer_rational_design.md`.

### B28 equal-leading-limb exact division shortcut (19 October 2026)

The shared `BigUnsigned` divider now has an exact fast path for
equal-length (at least two limbs), normalized operands with equal
leading limbs. Their quotient is mathematically 0 or 1; quotient-only
compares limbs, and remainder/joint modes use one subtraction with
borrow to obtain the exact result. B27's multiply-subtract divider,
other length cases, signed Integer behavior, FloatAP algorithms,
Rational/Dyadic, GMP adapters, selected outputs and B27 scratch
contract are unchanged. No architecture-specific instruction,
feature dispatch, CMake change or new API is introduced.
A standalone GMP differential stress probe passed 25,000 constrained
randomized cases for all output modes under the GCC/Clang paths,
including a Clang ASan/UBSan run. The production test suite gains
2048 additional equal-leading-limb GMP-oracle cases. The near-equal
subset is measurably faster on the isolated x86-64 GCC and Clang
microprobes, but production ARM64/x86-64 benchmark benefits remain
unverified, and full-mix timings can fluctuate. See B28 in
`integer_rational_design.md`; GitHub CI/Codecov checks are deferred
at user request.

### B29 selected sparse-divisor multiply/subtract kernel (19 October 2026)

The reusable `BigUnsigned` normalized division now selects an
optional C++20 compile-time zero-limb specialization only for at
least eight divisor limbs, eight quotient digits, and at least half
zero limbs after normalization. Multiplication by a zero divisor
limb is replaced with an exactly equivalent carry/borrow update.
The dense B27 multiply-subtract kernel and B28 quotient-0/1
shortcut remain unchanged; all three output modes use the same
mathematical division. No extra platform or build requirements,
CMake modifications, ABI/API changes, or FloatAP/Rational/Dyadic
algorithm changes. A 32000-case standalone differential GMP stress
test passed both GCC and Clang ASan/UBSan variants; production
tests add 1024 targeted sparse/dense-mixed oracles without removing
previous fixtures. Isolated x86-64 GCC/Clang profiling indicates
a benefit for sparse half-width divisors; integrated ARM64/x86-64
results remain pending. See B29 in `integer_rational_design.md`.

### B30 x86-64 DIVQ-based decimal output for shared BigUnsigned (19 October 2026)

B30 accelerates only the shared `BigUnsigned::to_decimal`
radix-10^18 leaf formatter on GNU/Clang x86-64, selecting it
for >=9 limbs using B26's preexisting 128/64 hardware `DIVQ`
helper, not B25's regression-prone compiler-generated
`unsigned __int128` division. The 128/64 precondition
follows from the decimal remainder being strictly below
10^18. B25.1 ARM64 wide formatting, B19 portable/MSVC
radix-10^9 formatting and the <=8-limb path remain
unchanged, as do recursive splitting, FloatAP algorithms,
Integer/Rational/Dyadic contracts, B29 division, GMP
adapters, benchmarks, compiler flags and CMake.
Standalone x86-64 GCC/Clang experiments (including
a recursive proxy with GMP internal splitting) indicate a
speedup, but integrated performance is **not** established.
A 30,000-case independent GMP differential formatter
stress test passed GCC and Clang ASan/UBSan; add 128
focused decimal GMP fixtures while preserving previous
tests. See B30 in `integer_rational_design.md`.
GitHub CI and Codecov remain deferred by user request.

### B30.1 restore B29 x86-64 decimal conversion after measured regression (19 October 2026)

B30's official x86-64 native Integer benchmark disproved the
isolated GNU/Clang DIVQ-radix-10^18 win: 1024/2048/4096-bit
`to_decimal` regressed from B29's 2247.362/5736.002/14893.242 ns
to 3715.747/7599.686/18242.587 ns. Numerical verification
passed on both ARM64 and x86-64, and the unchanged ARM64
wide-radix path exhibited only small run-to-run differences.
Because 1024-bit decimal conversion is nonrecursive, the
x86-64 leaf change itself is sufficient to explain the
observed regression. B30.1 restores the entire B29
`src/big_unsigned.cpp` byte-for-byte, thereby restoring
radix-10^9 output on x86-64 while preserving B25.1 ARM64
radix-10^18 output, portable/MSVC compatibility, B27-B29
arithmetic, parsing, APIs and benchmarks. B30's 128 GMP
decimal regression cases are kept and their architecture
comment updated. No FloatAP algorithm or representation
changes. A new integrated B30.1 run is still needed;
the B29 performance data are historical, not new
measurements. See B30.1 in `integer_rational_design.md`.
GitHub CI and Codecov are excluded by user request.

### B31 decimal leaf string rendering without temporary blocks (19 October 2026)

The shared `BigUnsigned` radix-10^9 leaf formatter
now pre-reserves its decimal blocks and assembles
the output into a single exactly sized string,
writing each lower nine-digit block backward as
four fixed-width base-100 digit pairs plus one
high digit from a constexpr 200-character table.
The most significant block retains variable width.
Arithmetic division by 1e9 is **unchanged**,
as are ARM64 B25.1 radix-10^18 leaves, B22
recursive splitting, B27-B29 division, GMP adapters,
all FloatAP algorithms and public numeric semantics.
A standalone 30000-case direct GMP differential
formatter test passed both GCC and Clang
ASan/UBSan; additional GMP fixtures target two-digit
zero padding and radix carry boundaries. GCC/Clang
isolated x86-64 leaf benchmarks favored the new
string assembly at <=16 limbs, but integrated
Integer and ARM64 effects remain unverified.
See B31 in `integer_rational_design.md`.
GitHub CI and Codecov are not consulted by user request.

