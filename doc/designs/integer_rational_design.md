# Integer, Rational, and GMP-independent arithmetic core

**Status:** initial design, not yet implemented.  
**Work branch:** `remove-gmp`.  
**Baseline:** `remove-mpfr` commit `9823a3b127712ee39e44b634dbaba49e50ea0054`, with Unix, Windows and Coverage passing.  
**Companion journal:** [`floatap_design.md`](floatap_design.md).

This is the **living design and implementation journal** for removing GMP
from the numeric library. Update this document for every relevant
arithmetic-core, `Integer`, or `Rational` change; update the FloatAP
journal **in the same commit** when FloatAP behavior or performance
is affected. Avoid recording prospective changes as implemented.

## Goals and non-goals

1. Remove GMP from **production** `Integer`, `Natural`, `Rational`,
   `FloatAP`, and all other numeric paths using GMP. No hidden GMP
   dependency through conversions, helper classes, or cached constants.
2. Share one portable, tested **unsigned limb arithmetic core**.
   Preserve separate algorithms and storage strategies where needed:
   `FloatAP` must not pay for generic rational canonicalization.
3. Preserve the mathematical **correctness** contracts needed by the
   library (exact values, comparisons, rounding where applicable, and
   Rational NaN/positive/negative infinity), but **do not preserve**
   legacy GMP-wrapper implementation details, public GMP fields,
   constructor/accessor ABI, or unnecessary API quirks. Explicitly
   decide and test any change in the observable arithmetic contract
   before adapting the rest of Ariadne.
4. Target C++20, GCC, Clang/AppleClang and MSVC with warnings-as-errors.
   Do not require `__int128`, GCC-only intrinsics, or platform GMP
   limb layout. Accelerated code can be optional with portable fallback.
5. Keep **FloatAP** as the primary regression/performance gate:
   at each of 53, 128, 256, 1024 bits require worst operation
   AP/MPFR slowdown **<2x**, geomean **<1x**, sum of operation times
   **<1x**, as the established targets. Also report absolute times.
   These are objectives, **not** established guarantees.
6. Benchmark future native `Integer` against GMP `mpz_t`, and native
   `Rational` against GMP `mpq_t` (not `mpz_t`, which cannot
   represent a general rational). Correctness checks are mandatory
   outside timing regions.

GMP/MPFR may remain installed **in benchmark and differential-test
targets** during migration. The eventual production link must be
independent of both. This project does not require removing
comparison backends or diagnostic tools before the native types work.

## Existing code and compatibility inventory

At branch creation:
- `include/numeric/integer.hpp`: `Integer::_mpz` is a **public**
  `mpz_t`; `Integer(const mpz_t)` and `get_mpz()` expose GMP.
  `src/integer.cpp` delegates its operations to `mpz_*`.
- `include/numeric/rational.hpp`: `Rational::_mpq` is a **public**
  `mpq_t`; `Rational(const mpq_t)` and `get_mpq()` expose GMP.
  `src/rational.cpp` delegates arithmetic/canonicalization to GMP.
- **Important special-value semantics:** `Rational::inf(+/-)` and
  `Rational::nan()` use a denominator of zero in existing code.
  `ExtensionOperations<Rational>` treats numerator 0/denominator 0
  as NaN and nonzero numerator/denominator 0 as infinity. A new
  fraction type must preserve externally visible handling via an
  explicit kind/tag, **not** assume Rational contains finite
  fractions only. Determine whether signed zeros are observable
  before changing representation.
- `include/numeric/floatap.hpp`: `FloatAP::_significand` is
  `mpz_t`, `significand()` returns `mpz_srcptr`, and conversion
  helpers accept `mpz_srcptr` or `mpq_srcptr`.
  The FloatAP transcendentals in `src/floatap.cpp` use
  `Fraction(mpq_t)`, GMP fixed-point integer intervals, and
  cached outward enclosures.
- `src/floatap_small_division.hpp` already provides a bounded
  32-bit-limb unsigned exact division prototype (capacity 18);
  it is **not** a general unlimited-precision type. The GMP
  import/export adapter is outside this file.
- Also audit `Natural`, `Dyadic`, `Decimal`, floating wrappers,
  `ExactNumber`, Python bindings, downstream modules, serialization,
  libraries and CMake target links for direct/indirect GMP ABI use.
  Do not assume replacing three fields removes GMP from the build.

These exposed `_mpz`/`_mpq` fields and GMP constructors/accessors
are **legacy implementation interfaces, not compatibility requirements**.
The user has explicitly authorized replacing the wrappers rather than
preserving their representation or GMP-specific ABI. Inventory and
migrate callers across Ariadne; a separate, temporary GMP adapter is
permitted for differential testing and transitional integration only.
Preserve meaningful mathematical behavior, especially extended
Rational values, and document changes to division conventions,
errors and other observable semantics explicitly rather than
carrying historical artifacts forward by default.

## Source-file style and documentation contract

All new and substantially modified C++ files must **follow the existing
Ariadne conventions**, using adjacent `numeric` sources as templates,
rather than introducing a separate formatting style. The concrete
reference files are `include/numeric/integer.hpp`,
`include/numeric/rational.hpp`, `src/integer.cpp`,
`src/rational.cpp` and `tests/test_integer.cpp`.

- **File header and license:** reproduce the Ariadne banner with
  `numeric/<filename>` (or the corresponding test filename),
  the full existing GPLv3-or-later notice and correct, non-invented
  attribution. Do not use the short experimental header from
  `src/floatap_small_division.hpp` as the template for a new
  reusable core file.
- **Documentation:** follow the project's `/*! \\file ... * \\brief ... */`
  and `//! \\ingroup NumericModule` / `//! \\brief` / `//! \\sa`
  conventions for public interfaces. Explain arithmetic invariants,
  aliasing, preconditions, quotient/remainder conventions, precision
  and directed-rounding guarantees where necessary, without
  duplicating obvious statements.
- **Layout:** respect Ariadne's include ordering/grouping, header
  guards, `namespace Ariadne`, class access-label indentation,
  declaration layout, blank lines, braces, and whitespace as used
  in the most relevant neighboring files. Use the existing naming
  conventions; do not run a repository-wide reformatter or create
  gratuitous unrelated diffs.
- **Tests:** follow the existing `TestInteger` /`TestRational`
  organization and `ARIADNE_TEST_*` macros. Separate fixed
  independent expected values from secondary GMP differential
  checks; preserve the full existing suites.
- **Review gate:** before each new file is committed, compare its
  banner, includes, documentation, naming, indentation and
  whitespace against the closest existing Ariadne counterpart.
  Treat warnings-as-errors on GCC, Clang/AppleClang and MSVC as
  part of the normal acceptance gate. Do not alter CMake/Windows
  build behavior merely to work around unrelated transient failures.

## Architecture and ownership

```text
            GMP-free unsigned limb arithmetic
             (views, owning storage, scratch)
                          |
        +-----------------+----------------+
        |                 |                |
  signed integer     exact fraction      dyadic fixed-point
        |                 |                |
     Integer            Rational       FloatAP internals
        |                 |                |
        +---------- explicit, exact conversions
```

The **shared layer** implements unsigned magnitude arithmetic and
memory management only; signs, NaN/inf, rounding and floating-point
exponents belong to higher layers. `Integer` is signed magnitude;
`Rational` is signed numerator / positive nonzero denominator for
finite values plus an explicit special-value tag; `FloatAP` stores
a normalized p-bit unsigned magnitude, sign, signed binary exponent
and finite/inf/NaN kind.

**Do not route FloatAP operations through general `Rational`** unless
needed for an explicit exact conversion or non-hot-path fallback.
FloatAP add/sub, multiply, divide, FMA, sqrt and transcendentals
use low-level primitives and outward dyadic intervals directly.
Avoid double rounding: certified endpoints must agree under the
requested mode before returning a FloatAP result.

### Proposed internal APIs (to refine after benchmarks)

- `BigUnsigned`: dynamic arbitrary-size little-endian limb sequence,
  canonical zero, size, capacity, `bit_length`, compare,
  `add`, `sub` (require lhs>=rhs), `mul`, `square`,
  `divrem`, shifts, bit tests, small-divisor `divrem_ui`,
  `gcd`, allocation/resize/swap and exact import/export.
- `LimbView` / `MutableLimbView` / scratch: nonowning spans for
  algorithms that cannot afford copies; documented aliasing,
  failure, and normalization contracts.
- `SignedInteger`: sign and magnitude; define negative/zero
  canonicalization, quotient/remainder conventions separately.
- `FiniteFraction`: numerator sign and positive unsigned
  denominator; canonical numerator/denominator (e.g. 0/1,
  gcd=1), explicit `RationalKind` for ±inf and NaN.
  Benchmark whether cross-cancellation before products is
  faster than unconditional post-operation gcd.
- `BinaryScaledInteger` / outward interval: signed magnitude
  times an exact power of two, two ordered dyadic bounds,
  outward-directed shifts/quotients and exact rounding
  agreement tests. This is **not** just a general Fraction
  with denominator forced to 2^k.
- `GmpAdapter` (transitional and optional): exact
  import/export with `mpz_t` / `mpq_t`, located away from
  the portable core and outside timed native operations.

Naming, object layout, limb width and public headers are
**design candidates**, not frozen decisions.

### Limb width, memory and algorithm policy

Start from the proven portable 32-bit-limb division logic;
benchmark 32-bit and 64-bit storage/primitive strategies
on all target compilers before committing globally.
`uint64_t` provides exact intermediate products of 32-bit
limbs. Support 64-bit wide-limb multiply/divide only
behind a portable implementation or explicit feature checks,
with a strict compiler-independent fallback.

Small-buffer optimization: inline capacity for common 53/128-bit
cases, heap fallback for arbitrary precision, well-defined move
semantics and no fixed-size limit. The current 18-limb prototype
cannot represent arbitrary `Integer` and must be generalized.
Avoid forced heap allocation for 1-4 limbs if benchmarks justify
the inline footprint; report cache/storage costs.

Keep multiple algorithm dispatches by operand length:
single-limb, short multi-limb, medium, large. Measure before
choosing thresholds. Correct asymptotic algorithms matter for
1024+ bits; optimize small sizes without damaging them.
Implement `divrem` once and select rounding at a separate
layer (nearest-even, down, up, toward zero).
Explicitly test exact quotient, remainder, carry, borrow,
limb boundaries, normalization and large shifts.

## Integer semantic contract

The new `Integer` is a fresh signed-magnitude front end backed by
the shared `BigUnsigned`. Keep its useful mathematical capabilities:
construction from standard integers and decimal strings, comparisons,
sign/zero, add/sub/mul, exact square/pow, shifts, increment/decrement,
well-defined quotient/remainder, checked casts and IO. The exact
legacy overload set and GMP-shaped storage/constructors/accessors are
**not requirements**; update internal callers to the new interfaces.
Before deciding on division semantics, inspect current
`src/integer.cpp`:

```cpp
Integer quot(Integer const& z1, Integer const& z2) {
    Integer r;
    mpz_div(r._mpz,z1._mpz,z2._mpz);
    return r;
}
Integer rem(Integer const& z1, Integer const& z2) {
    Integer r;
    mpz_cdiv_r(r._mpz,z1._mpz,z2._mpz);
    return r;
}
```

These use **different GMP division functions**. This is a reason
to define a coherent signed division contract afresh, not a reason
to perpetuate accidental wrapper behavior. Start by documenting
current behavior on all sign combinations; then choose one explicit
quotient/remainder convention and enforce the same identity
`a=q*b+r` for every test. Record intentional migration changes
and update consumers/tests. The low-level unsigned `divrem` remains
exact and independent of the higher-level sign convention.

Check treatment of `Integer` vs `Natural` (including
zero, signs and conversion preconditions), and negative
shifts/powers. Document every intentional breaking
change separately.

## Rational semantic contract

Finite: `n/d` with `d>0`, coprime magnitudes and
canonical zero `0/1`, subject to an explicit decision
on when reduction must be materialized vs lazily
proven. A slow repeated `gcd` would harm `Integer`,
`Rational` and FloatAP performance.

Extended: preserve the currently observable
NaN/positive-infinity/negative-infinity behavior,
including constructors, `is_nan`, `is_inf`,
`is_finite`, comparisons, arithmetic, reciprocal,
zero division, sign and conversions. Record concrete
existing tests for special values before changing
their internal representation. Do not rely on GMP's
undocumented arithmetic on denominator zero:
current `ExtendedOperations` determines the
public behavior; reconstruct it explicitly.

`Rational` is benchmarked against **`mpq_t`**.
Optional comparisons of numerator/denominator primitives
against `mpz_t` are supporting diagnostics, not a
substitute for end-to-end rational benchmarks.

## Validation and performance contracts

**Functional correctness first.** Keep and grow
the existing test suites. Baseline is
`remove-mpfr` at `9823a3b1`: Unix, Windows
and Coverage green at branch creation.

- Integer: fixed exact expected values independent
  of GMP for add/sub/mul/quot/rem/shifts/pow, negative
  operands, zero, sign/magnitude boundaries, 32/64-bit
  limb transitions, 53/128/256/1024-bit and larger
  values; randomized cross-checks with `mpz_t` as
  *secondary* oracle.
- Rational: fixed independent expected fractions
  and special-value outputs; test normalized
  numerator/denominator, arithmetic, compare,
  sign, ±inf/NaN and exception policies; randomized
  `mpq_t` differential testing as secondary oracle.
- FloatAP: preserve exact certified-literal tests
  at four precisions and all rounding modes;
  benchmark continues to compare the 32 paired inputs
  plus 3456 operation results with FloatMP **after**
  timing (not inside `median_ns`). Also test
  non-centered domains and hard rounding cases.
- In all benchmarks, compute the *same* input and
  operation definitions for both backends and verify
  equality **outside timing**. No tolerance or
  double-conversion when comparing exact values.
  Fail with operation, operand sizes/signs, index
  and both results. Compare special-value kinds
  explicitly. The 1000x11 settings and original
  FloatAP operation pool must not be changed to
  manufacture better ratios.

### Planned integer benchmark

Proposed executable: `benchmark_compare_integer_backends`
(name to be finalized with the CMake integration).
Compare native `Integer` vs raw GMP `mpz_t`,
**not** vs another wrapper whose object allocations
are accidentally included for only one side.
Use a common exact operand specification prepared
outside timed regions. Include bit sizes
0/1, 32, 53, 64, 128, 256, 1024, 4096+
and adversarial limb boundaries; multiple signs.
Measure construct/parse, copy/move, compare,
add/sub, mul/square, shift, div/rem, gcd and
integer->text, distinguishing construction
and allocation from the arithmetic kernel itself.
Report per-operation median ns, AP/GMP ratio,
geomean, total-time ratio, worst, sample counts
and absolute times. Add a 1000x11-style default
where sensible, plus a reproducible random seed.
Never change the corpus to favor our backend.

### Planned rational benchmark

Proposed executable: `benchmark_compare_rational_backends`.
Compare native `Rational` vs raw GMP **`mpq_t`**,
with direct exact-rational operand pairs of similar
bit-length, near equal fractions, relatively prime
and large shared factors, dyadic denominators,
cross-cancellation opportunities, zero, and
special values in **separate semantic tests**
when GMP itself does not define the operation.
Measure construct/canonicalize, copy/move,
comparison, +,-,*,/, reciprocal, numerator/denominator
extraction, floor/ceil and exact IO. Report
per-operation and aggregate timing as above.
The `Rational` facade should not gain speed by
skipping a required canonicalization that raw GMP
actually performs; compare equivalent semantics.

### G0b benchmark scaffolding: Integer, Rational, and absolute timings

**Implemented before replacing the GMP wrappers.** Four C++20
executables share one workload header
(`benchmarks/benchmark_exact_arithmetic_workload.hpp`) and one
report/argument-parsing header
(`benchmarks/benchmark_exact_arithmetic_report.hpp`):

| Executable | What is timed | Reference |
| --- | --- | --- |
| `benchmark_compare_integer_backends` | native `Integer` and raw `mpz_t` | GMP C `mpz_*` |
| `benchmark_compare_rational_backends` | native `Rational` and raw `mpq_t` | GMP C `mpq_*` |
| `benchmark_integer` | native `Integer` only | absolute nanoseconds |
| `benchmark_rational` | native `Rational` only | absolute nanoseconds |

They follow the existing FloatAP comparison's defaults:
`--iterations 1000`, `--samples 11`, `--bits
53,128,256,1024` and 32 deterministic input pairs at each
operand size. Other positive bit sizes (e.g. 32, 64,
4096) can be supplied with `--bits`. They support
`--summary-only`, `--details`, standard Ariadne logging,
and comparison-only `--max-slowdown`. The timer uses
`steady_clock`, one unrecorded warm-up sample, an odd-sample
median, and keep-alive barriers, exactly like FloatAP.
Individual rows include operation, size and median ns;
comparisons report geomean, total-time ratio and worst
ratio for each size and all sizes. Absolute benchmarks
report operation times and sums of medians.

`Integer` battery: **12 operations** — decimal
construction, copy, compare, negate, add, sub, mul,
square, quotient, remainder, cube and decimal output.
The raw GMP path uses `mpz_*` directly; until the
Integer replacement, quotient and remainder match
the legacy wrappers' *distinct* GMP quotient modes
for a fair baseline. Document/change the signed
division contract later rather than masking
differences in this benchmark.

`Rational` battery: **15 operations** — fraction
construction (including normalization), copy,
compare, negate, add, sub, mul, div, square,
reciprocal, numerator, denominator, floor, ceil
and conversion to binary64. The GMP reference
uses `mpq_*` and `mpz_*` directly, and includes
`mpq_canonicalize` on construction. Both backends
receive the same exact numerator/denominator
specifications. Alternating operand signs and
shared factors test costs of canonicalization and
signed operations; special `NaN`/±inf cases remain
in separate semantic tests because raw GMP
`mpq_*` does not provide Ariadne's extended
rational arithmetic contract.

Comparison executables rerun the **exact same
shared operations after all timed regions** and
check 32 pairs of inputs and every operation
result at every selected size. The default four
size groups imply 256 paired input comparisons
per executable, and **1536 Integer** or
**1920 Rational** paired operation results.
Exact values are compared as normalized decimal
integer or numerator/denominator strings;
binary64 is compared by its full 64-bit pattern.
Failures report operation, bit size, index, and
the differing values. These are differential
regression checks, **not independent proof** of
correctness; the separately planned fixed
mathematical expected-value test suites remain
essential.

The two benchmark backends currently use GMP
internally on **both sides**: the new comparison
reveals wrapper overhead only, not the future
GMP-free performance. GMP remains a permitted
dependency in these benchmark executables.
The native benchmark functions are explicitly
separate from their reference callbacks; both
include ordinary return-object construction and
destruction as in the FloatAP harness. To assess
pure arithmetic-kernel cost later, add separate
documented in-place microbenchmarks instead of
pretending that the existing facade timings
measure naked GMP kernel calls.

Run four CTest smoke commands with reduced
iterations/samples and the 53-bit group, then
the complete 1000x11 cross-size performance run
on a fixed host. **No performance claims until
the CI and user's measurement are observed.**
The FloatAP benchmark, its input workload,
and its numerical oracle remain unchanged.

### CI and metric gates

Run Unix/GCC+Clang/macOS, Windows/MSVC and Coverage,
and inspect Codecov patch status separately where
available. New helper APIs require unit tests,
algorithm-boundary property tests and, for any
backend dispatch, baseline before/after timings.
Do not lower precision, omit hard input cases,
shorten needed enclosures, weaken rounding or
change tested scenarios simply to pass a
regression or improve metrics. Distinguish
actual source bugs, incorrect expected values
and test-phase dependencies when diagnosing failures.

## Migration phases

**G0 — Inventory/architecture (this commit):**
create `remove-gmp` from the verified `remove-mpfr`
head, retain its history and FloatAP journal, add
this design and the two planned benchmark contracts.
No GMP replacement, benchmark code or numerical
behavior change is claimed.

**G1 — Portable unsigned foundation:** extract
`src/floatap_small_division.hpp` into a reusable,
tested `BigUnsigned` internal core; make capacity
unbounded, add storage and missing operations.
Do not prematurely switch every implementation.

**G2 — Integer replacement:** implement `Integer` and
related `Natural` interfaces on the GMP-free core,
without reproducing GMP-specific fields/accessors.
Migrate internal consumers and run the
`Integer` vs `mpz_t` benchmark. Preserve
mathematical correctness, not wrapper ABI; document
and explicitly test any selected API/semantic changes.

**G3 — Rational replacement:** build the finite and
extended `Rational` type on the new signed/unsigned
integer core, not on `mpq_t`. Preserve NaN, ±infinity,
finite exact arithmetic and the mathematical comparison
contract. Rework callers that relied on GMP fields and
run `Rational` vs `mpq_t` benchmark with equivalent
canonicalization.

**G4 — FloatAP migration:** replace `mpz_t`
significand and GMP `Fraction`/certified interval
hot paths with shared limbs and specialized dyadic
rounding; rerun every precision/mode and existing
performance benchmark. G2-G4 may proceed in
parallel through adapters provided each commit is
fully tested.

**G5 — Project-wide removal:** audit remaining
`mpz`, `mpq`, `mpf`, GMP C++ bindings, Python
bindings, dependencies and CMake/package exports.
Production binaries must build and run with GMP
and MPFR absent; comparison-only benchmarks/tests
may link both explicitly. Do not remove GMP from
build configuration while native code still needs it.

## Decision G1 architecture: shared limb kernels, separate ownership (18 October 2026)

**The primary design constraint is FloatAP performance.** Sharing
the **algorithms and primitive contracts** is encouraged; sharing
the owning `BigUnsigned` *object and allocation policy* is not
required. Do not force a generic integer object, dynamic memory
allocation or rational normalization onto FloatAP hot paths.

Split the implementation into two layers:

1. **Portable low-level limb operations** accepting contiguous
   little-endian unsigned limb views and caller-owned result/
   workspace buffers: comparison, bit length, add/sub/carry/borrow,
   shifts, multiplication, exact divrem and eventually GCD. Define
   normalization, output size, aliasing rules and scratch-size
   requirements precisely. Start with portable `uint32_t` limbs
   and `uint64_t` intermediates; any optional specialization
   must preserve a working GCC/Clang/MSVC fallback. Do not
   expose GMP types or allocator behavior in this layer.
2. **Independent owning representations** using those primitive
   operations: (a) dynamic, arbitrary-size `BigUnsigned`
   with inline storage as justified by measurements, eventually
   backing signed `Integer` and finite/extended `Rational`;
   (b) FloatAP-specific normalized precision-aware significand
   and scratch storage, able to bypass a generic heap-allocating
   object for 53/128-bit fast paths and handle 256/1024+ bits.
   For rigorous transcendental work, allow specialized dyadic
   intervals with outward rounding instead of generic fractions.

It is acceptable for FloatAP to contain *measured and documented*
specialized scalar or short-limb algorithms (including an
optimized 64-bit case) even if those do not invoke the general
`BigUnsigned` facade. Do not duplicate a full independent
arbitrary-size big-integer engine unless cross-backend
measurements clearly demonstrate a need. The existing 18-limb
`FloatAPSmallDivision` supplies useful division algorithm
ideas, but not an owning representation that should simply
be used unchanged by `Integer`.

**First implementation commit (G1a):** introduce an isolated,
GMP-free low-level unsigned-limb kernel with fixed-oracle tests
for zero, carry/borrow, limb boundaries, shifts, comparison
and multiplication, including longer-than-inline inputs and
aliasing semantics. Implement `BigUnsigned` ownership on top
in a second incremental step, without modifying the
current `Integer`/`Rational`/FloatAP implementations.
Build functional tests through CTest; keep **all four
performance benchmarks separate from CTest**, using only
the dedicated `benchmarks` target. Record CI and benchmark
outcomes before any FloatAP backend migration.

**Performance acceptance:** for each proposed change
to FloatAP's storage or arithmetic, compare baseline and
candidate on the same host, with the unchanged 1000x11
workload and the post-timing exact result check. Report
absolute timings and worst/geomean/weighted-total ratios
*separately* at 53/128/256/1024 bits. Retain certified
independent expected-result literals for correctness.
Keep a specialized implementation when the shared
owning type demonstrably harms performance, while
preferring shared low-level algorithms where costs
are equivalent. Do not claim a GMP-free FloatAP simply
because its helper kernel is GMP-free while the
significand or interval representations still use GMP.

## Recommended first implementation milestone: G1a

**Start with an independent unsigned core and its tests, not by
rewriting the public `Integer` or `Rational` wrappers in place.**

1. Freeze a minimal internal contract: canonical zero, unsigned
   magnitude, little-endian limbs, bit-length, compare, add/sub with
   carry/borrow, shifts and multiply/square. Set explicit aliasing,
   overflow/allocation and error contracts; use arithmetic operations
   on caller-provided scratch or output buffers where appropriate.
2. Implement a portable `BigUnsigned` storage type with dynamic
   capacity and a **measured** inline small-value optimization.
   Use standard C++20 32-bit limb arithmetic as the first fully
   portable baseline. Keep limb width and SBO capacity as
   benchmark-guided decisions, not public API commitments.
3. Introduce deterministic fixed expected-value tests independent
   of GMP for 0/1, powers-of-two boundaries, carry/borrow chains,
   aliasing, 32/53/64/128/256/1024+ bits, and stress cases beyond
   inline capacity. Add randomized differential GMP checks in
   test-only code. Follow Ariadne's source-header, GPL, Doxygen,
   indentation and test organization requirements.
4. Add exact unsigned `divrem`, single-limb division and GCD only
   once core add/sub/shift/mul invariants and tests are established.
   Reuse the `FloatAPSmallDivision` long-division *algorithm*
   after reviewing its fixed-capacity edge cases, but do not
   elevate the current 18-limb type into a supposedly arbitrary
   precision implementation.
5. Add the **benchmark infrastructure** for Integer-vs-`mpz_t`
   early, separating end-to-end object costs from kernel costs.
   Exercise the same exact inputs and check correctness after
   timing. A Rational-vs-`mpq_t` benchmark follows the finite
   Rational implementation. FloatAP's unchanged 1000x11
   benchmark is the cross-cutting regression reference.
6. Only then replace `Integer` with a fresh signed-magnitude
   facade and migrate call sites, followed by `Rational`
   (finite fraction plus explicit NaN/±infinity kinds).
   Once the core is mature, migrate FloatAP's unsigned
   significand and fixed-point interval paths, retaining
   exact-rounding certification.

**First commit acceptance:** a self-contained GMP-free unsigned
core plus substantive fixed-oracle tests, Unix/Windows/Coverage
passing, **without** changing current `Integer`, `Rational`,
`FloatAP`, benchmark workloads or build configuration beyond what
is required to add and exercise the new tests. No speculative
performance improvement is claimed before same-host measurements.

## Decision log

### G0 — Branch creation and design baseline (18 October 2026)

Created `remove-gmp` from verified `remove-mpfr`
`9823a3b1`. Chose a **common unsigned arithmetic
core with specialized higher-level implementations**,
not multiple incompatible integer engines and not
a general Rational allocation path for FloatAP.
The first milestone is the living design only.
No new numerical implementation or benchmark
has yet been added. The FloatAP journal remains
the performance reference and is to be updated
in parallel when core changes affect FloatAP.

### G0 clarification — wrapper compatibility is not a goal (18 October 2026)

The user clarified that the current `Integer` and `Rational`
classes are essentially GMP wrappers and **do not need to be
preserved** as implementation classes or as GMP-shaped APIs.
The only explicitly emphasized compatibility requirement is the
handling of `Rational` NaN and ±infinity, in addition to exact
mathematical correctness. Consequently the migration can
replace `mpz_t`/`mpq_t` public fields, GMP constructors and
accessors, and reorganize the classes. Follow a staged caller
migration rather than maintaining a permanent compatibility
facade. Decide quotient/remainder semantics explicitly and
update consumers as needed. Begin with isolated `BigUnsigned`
and independent tests before replacing either type; FloatAP
performance/correct rounding remains the key regression gate.

### G0b clarification — keep benchmarks out of CTest (18 October 2026)

**Current design rule (supersedes previous G0b smoke-test integration):**
CTest and the custom `tests` target contain only functional/unit
tests and their build prerequisites. Benchmarks are separate
executables managed solely by `benchmarks/CMakeLists.txt` and
the `benchmarks` target. Specifically, remove all four
`smoke_benchmark_*` CTest entries, the special
`benchmark_smoke_executables ALL` target, and the benchmark
dependency of the `tests` target. Retain the four named
benchmark executables and their comparison/absolute data paths.

The earlier G0b journal entries below record the historical
smoke-test experiment and its CI workarounds; they are no longer
the intended configuration. A benchmark can still fail by
throwing on mismatched post-timing outputs **when invoked
explicitly**. Such verification is part of the benchmarking
program, not a unit test automatically registered with CTest.
Performance runs use the same-host standard `1000x11`
setup unless the user explicitly specifies otherwise.
No numerical code, input pool, expected result or measured
operation is modified by this separation. Revalidate the
ordinary Unix, Windows and Coverage builds after removal.

### Implementation record template

For every G1+ change record: date, commit,
files/API changed, invariants and exact-rounding
proof/limits, independently expected tests,
differential results, Unix/Windows/Coverage and
Codecov status, absolute timings per precision,
FloatAP worst/geomean/total changes, Integer
and Rational operation timings as available,
known regressions, and next action. Never claim
CI/performance results before observing them.

### G0b implementation log — benchmarks introduced (18 October 2026)

Added four independent executables and shared workload/report
helpers while leaving current Integer, Rational and FloatAP
implementation untouched. The Integer GMP comparator calls
`mpz_*`, and Rational GMP comparator calls `mpq_*`,
not the Ariadne wrapper types. Aggregate and absolute
timings are separated; numerical agreement is verified
after the timing loops. CMake registers all targets and
four lightweight CTest smoke tests. The initial commit
only introduces benchmark infrastructure: green CI,
successful smoke tests and same-host 1000x11 numbers
must be checked and recorded separately.

### G0b follow-up — ensure CTest benchmark executables are built by coverage (18 October 2026)

The first new benchmark commit `f0a60f78` compiled and
passed the Release pipelines on GCC 13, Clang 16,
AppleClang 21 and Windows. Both Coverage jobs, however,
failed in their combined build/test phase, and the
Unix ThreadSanitizer job failed at Test. The connector
exposes failure *steps* but not full job logs, so the
specific failure message is **not yet established**.

The Coverage target builds a project-defined tests
target, which originally depended only on the
unit-test executables. Four new CTest smoke tests invoke
the new benchmark executables, so those executables
must be built whenever the `tests` target is built.
Add four explicit benchmark target dependencies to
`tests`. This does not change any benchmark code
or operation times, and avoids relying on an all-target
build in Coverage. Re-run Unix/Windows/Coverage and
distinguish any remaining sanitizer failures; no
all-green claim is warranted until observed.

### G0b ThreadSanitizer: build smoke targets with default ALL (18 October 2026)

The supplied ThreadSanitizer CTest output from build `build-tsan`
shows **all 19 numerical tests passing**, followed by four
`smoke_benchmark_*` tests marked `Not Run` because their
benchmark executables do not exist. This is definitively a
missing-build-target failure, not a numerical comparison
failure or detected data race. The workflow calls
`cmake --build . --parallel 2` (default `ALL` target)
before `ctest --output-on-failure`, whereas the previous
commit only attached these executables to the *custom*
`tests` target; the two build roots are different.

The solution is a single `benchmark_smoke_executables ALL`
custom CMake target, depending on exactly the four
registered benchmark binaries. The `tests` target also
explicitly depends on it. Hence both TSAN's default build
and Coverage's `tests`-based build produce all executables
that CTest needs. This does **not** remove the smoke tests,
skip TSAN, touch numerical code, change benchmark inputs or
affect the Release performance workload. The preceding
`24e4a290` commit is confirmed Coverage and Windows green,
but Unix remains red from TSAN due to missing binaries.
The current change must pass Unix/GCC/Clang/macOS,
Windows and Coverage before claiming completion.
