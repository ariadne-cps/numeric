# Integer, Rational, and GMP-independent arithmetic core

**Status:** native Integer and Rational implemented; remaining GMP-backed types and performance work are in progress.
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

## Current shared-core decision — reusable native substrate (18 October 2026)

**Active policy; supersedes older FloatAP-first isolation notes
only for genuinely reusable low-level arithmetic.** The user
requires one native multiprecision substrate reusable in
`FloatAP`, in the eventual replacement of `mpz_t` inside
`Integer`, and in the eventual replacement of `mpq_t`
inside `Rational`. Its portable C++20 reference and
optional, clearly documented ARM64/AArch64 and x86-64
kernels are implemented once in the shared internal
numeric layer, not as copies hidden inside FloatAP.
Internal core files are allowed and must follow the
Ariadne banner, style, build and documentation rules
below. This is **not** authorization to integrate
an existing third-party arithmetic engine or its
source. Mathematical ideas may be implemented
independently.

The generic core owns unsigned limb storage/views,
buffer reuse, exact add/sub, shifts, multiply,
square, quotient/remainder and useful fixed-point
building blocks; all ISA-specific dispatch is here.
The three layers above it retain their own semantic
contracts: signed integers and division conventions
for `Integer`, finite exact ratios/GCD and explicit
special values for `Rational`, and normalized
significands/exponents with IEEE-directed certified
rounding for `FloatAP`. FloatAP hot operations
must **not** incur generic fraction canonicalization.

**Migration ordering remains selective:** the
performance gate continues to prioritize FloatAP,
and `Integer`/`Rational` production types are
**not** being migrated by this documentation
decision. Unlike the earlier G1 decision to
isolate FloatAP's first prototype, new
cross-cutting limb primitives should now be
implemented for reuse from the outset.
FloatAP-specific algorithms, rounding and
transcendental code remain confined to
`include/numeric/floatap.hpp` and
`src/floatap.cpp`; only generic arithmetic
may be factored into shared core files.

**Per-ISA documentation is mandatory before implementation:**
every optimized kernel must state its input/output
mathematical contract, tested ISA requirements,
feature detection/dispatch, overflow and carry
proofs, aliasing and scratch layout, portable
fallback, differential tests and sanitizers,
measured absolute performance and compiler/host
details. See the paired
`floatap_design.md` architecture-specific
engineering contract. Documentation for shared
kernels must be updated **here** and in that
companion journal in the same commit. No
unverified improvement claims.

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

### G1 priority change — FloatAP before Integer/Rational (18 October 2026)

The user has explicitly prioritized finishing the GMP-free
implementation of **FloatAP** before doing any `Integer` and
`Rational` replacement or refactoring for code reuse.
The previously recommended shared unsigned foundation
and G1a `BigUnsigned` facade are therefore **not an
immediate milestone**. The `FloatAPNative::Unsigned`
type in `src/floatap_native_unsigned.hpp` is deliberately
private to FloatAP; its correctness and FloatAP performance
are more important than premature reuse. The new
`test_floatap_native_unsigned` is a real functional CTest
test, distinct from the four benchmark executables
which must never be registered with CTest. FloatAP's
public GMP interface, underlying `mpz_t` storage and
certified interval `mpq_t` remain until incremental
migration. This entry describes the initial isolated
kernel only, not completion of GMP removal.
### Migration priority refinement — one FloatAP source pair (18 October 2026)

The user explicitly requires only `floatap.hpp` and `floatap.cpp`
for FloatAP, not a separate header per arithmetic algorithm.
All unsigned kernels and quotient/remainder rounding have been
moved into `floatap.cpp`, with no changes to the public header.
Future GMP removal should directly refactor the existing FloatAP
representation and implementation, retaining FloatMP as a
differential reference and independent fixed-answer rounding
tests. Do not build a temporary compatibility layer or factor
out code for Integer/Rational before measuring an actual need.
Private implementation tests have been replaced with tests of
the FloatAP public mathematical behavior, under CTest.
Benchmark targets remain separate.


### BigUnsigned phase B0 — portable storage and foundation (18 October 2026)

**Implementation introduced**, independently of FloatAP and with no GMP
dependency in `include/numeric/big_unsigned.hpp` or
`src/big_unsigned.cpp`, built into `NUMERIC_SRC`. The
`tests/test_big_unsigned.cpp` CTest target validates it against GMP
as a **test-only reference**. `FloatAP`, `Integer` and `Rational`
retain their current production representations; no operations have
been rewired to use BigUnsigned. FloatMP is unchanged.

**Storage contract.** Nonnegative exact values are stored as a
little-endian sequence of 64-bit unsigned limbs; each represents a
coefficient of `(2^64)^i`. Zero is the empty sequence.
The highest stored limb is always nonzero for nonzero values.
Four limbs (256 bits) fit in the value's inline buffer;
larger values use dynamically sized contiguous storage.
`size()` returns the number of significant limbs,
`limb(i)` returns zero beyond that length,
`limbs()` exposes a read-only, invalidated-on-mutation
span. Constructor normalization strips upper zero limbs.
All normalizing mutations switch back to inline storage
if at most four limbs remain. Copy owns independent storage;
move leaves the source as zero, and `swap` swaps all
representation fields. Array/vector storage is intentionally
independent of GMP limb types and platform byte order.

**Portable operations now implemented:**
`compare` returns -1/0/+1;
`add` and `sub` compute exact nonnegative magnitudes
with carry and borrow limited to 0 or 1;
`sub` throws `std::underflow_error` if the
result would be negative. `shift_left` computes
exact multiplication by `2^k` and checks representable
result-size overflow; `shift_right` is exact
floor division by `2^k`. All operations return
owning values and leave input values unchanged.
No CPU-specific intrinsics, inline assembly or
runtime ISA dispatch are present yet.

**Cost and deliberate limits.** Allocations are avoided for
the common result range of at most four limbs, including
the **normalization of carry-free add/shift results**
that provisionally require an extra limb. Data for
larger values uses standard contiguous vector storage.
Storage operations are foundation code, **not yet
competitive GMP replacement**. Still to be implemented:
writable limb views and external scratch, in-place
arithmetic, multiplication and squaring, exact
division/remainder, small-divisor APIs and GCD;
only then test optional ARM64 and x86-64 variants.
This phase makes no performance claims.

**Verification.** GCC and Clang C++20 with strict
`-Wall -Wextra -Wshadow -Wconversion -Werror -pedantic`
both compile and pass the same deterministic
10,000-case GMP-differential test covering
add/sub/compare, normalization and both shifts
up to 1029 bits on operands with up to 20 limbs;
Clang AddressSanitizer and UndefinedBehaviorSanitizer
pass the same cases. The actual repository CTest
must also pass in CI before calling the change
fully integrated across platforms. These standalone
runs used a minimal local macro shim so they are
not a substitute for repository CMake/CTest.

**Architecture specific code** remains strictly subject
to the mandatory per-kernel journal and ISA requirements
in `floatap_design.md`; there is currently no
specialization to document. Neither benchmarks nor
the frozen ARM64 20,000x31 FloatAP workload was changed.

### Authorized conditional ARM64/x86-64 source split (18 October 2026)

The user explicitly approved the **exact future filenames**
`src/big_unsigned_arm64.cpp` and
`src/big_unsigned_x86-64.cpp` **if and when** distinct
architecture-specific implementations are justified.
These files are **not part of phase B0**, do not exist
yet and must not be created as empty placeholders.

`include/numeric/big_unsigned.hpp` remains the
single architecture-neutral interface.
`src/big_unsigned.cpp` owns the portable C++20
reference implementation, storage, normalization and
any architecture-neutral selection logic.
Genuinely ARM64/AArch64-specific low-level kernels
may eventually live in `src/big_unsigned_arm64.cpp`;
genuinely x86-64-specific low-level kernels may
eventually live in `src/big_unsigned_x86-64.cpp`.
The hyphen in `x86-64` is intentional in this
approved filename. No third-party arithmetic
engine and no duplicated floating-point semantics.

**Split criterion.** Create a dedicated source
file only when a first kernel using native ISA
instructions demonstrates a meaningful advantage
under a reproducible, full-operation workload
and its portable equivalent has passing
bit-exact differential tests. Document separately:
contract and carry bounds, instruction/CPU-feature
requirements, compiler flags, runtime/compile-time
dispatch and portable fallback, aliasing and scratch
lifetime, test forcing of each path, ISA-appropriate
CI and sanitizers, benchmark host/compiler, absolute
timings, and operation-level crossover.
Architecture files must be compiled only on
supported targets; generic builds including
Windows/MSVC must continue to build through
portable fallbacks. No unguarded instruction
assumptions. The per-kernel documentation belongs
in this journal and in `floatap_design.md`
with the same source change.

**Current validation note.** The original B0
commit `be34905` passed Unix and macOS CI jobs
but its Windows Release job failed during **Build**,
before tests. The precise MSVC diagnostic has not
yet been verified; this documentation-only
approval does not resolve the Windows failure.
Do not call phase B0 cross-platform validated
until Windows passes.

### BigUnsigned B0 Windows fix and style (18 October 2026)

The initial `be34905` Windows 2025 Release
job failed in MSVC 14.51.36231 when parsing
`BigUnsigned::_inline` in the class declaration.
The resulting C2226/C2334 and follow-on
C3409/C2143/C2059 diagnostics point to an
identifier collision: `_inline` is interpreted
as a compiler-specific inline keyword/extension.
The member is renamed everywhere to
`_inline_limbs` with no change to the
storage representation, arithmetic algorithms,
tests, or ABI-visible public methods.

The same maintenance commit expands whitespace,
blank-line separation and control-flow formatting
in `src/big_unsigned.cpp` and the dedicated test,
following neighboring Ariadne numeric sources.
Existing differential and adversarial tests are
retained unaltered in their mathematical
coverage. MSVC/Windows CI must pass before
claiming the fix is verified; GCC/Clang
regressions must also be checked. No
architecture-specific kernels or FloatAP
arithmetic changes are introduced.

### BigUnsigned B0 Windows GMP-test bit-count fix (18 October 2026)

MSVC 14.51.36231, compiled with `/W4 /WX`, exposed
C4267 errors in `test_big_unsigned.cpp` at the two
test-only GMP calls `mpz_mul_2exp` and
`mpz_fdiv_q_2exp`. Their shift parameter is
`mp_bitcnt_t`, which can be 32-bit on Windows,
whereas the randomized native shift was held in
64-bit `std::size_t`. The fix keeps
`std::size_t` for the `BigUnsigned` API, verifies
at compile time that `mp_bitcnt_t` can represent
the maximum tested shift of 1029 bits, and
makes the conversion to `mp_bitcnt_t` explicit
**only when invoking the GMP test oracle**.
No warning is disabled or demoted, no
test distribution/case coverage changes, and
no production `BigUnsigned` code is modified.
Revalidate the Windows Build and CTest alongside
the existing Linux and macOS checks.

## Baseline: legacy Integer wrapper versus direct GMP (18 October 2026)

**User-supplied same-run comparison, to freeze BEFORE native Integer
migration:** `benchmark_compare_integer_backends`, 20,000 iterations,
31 measured samples, operand sizes 53/128/256/1024 bits. The
`integer` backend in these numbers is the **existing GMP-backed
`Ariadne::Integer`**, not the new `BigUnsigned`. The
`mpz` backend is the test harness's `RawInteger` RAII wrapper
calling GMP directly. Numerical verification: PASS,
256 paired inputs, 1536 paired operation results.

Each table cell is **direct GMP ns / existing Integer ns**.
The benchmark configuration and input distribution are frozen
for subsequent comparisons; do not change them to improve scores.

| Operation | 53 bit | 128 bit | 256 bit | 1024 bit |
|---|---:|---:|---:|---:|
| construct_decimal | 64.402 / 39.773 | 63.800 / 81.585 | 106.315 / 125.775 | 354.738 / 372.910 |
| copy | 18.521 / 16.554 | 15.885 / 17.267 | 18.492 / 19.990 | 19.590 / 18.798 |
| compare | 1.902 / 3.700 | 1.729 / 4.008 | 2.281 / 4.625 | 4.377 / 7.275 |
| negate | 15.412 / 17.985 | 15.110 / 33.092 | 22.996 / 46.098 | 21.767 / 44.823 |
| add | 21.646 / 34.671 | 18.781 / 48.604 | 26.215 / 51.523 | 28.010 / 55.248 |
| sub | 20.075 / 33.792 | 20.475 / 48.558 | 24.996 / 50.490 | 28.727 / 53.817 |
| mul | 17.925 / 32.502 | 21.508 / 36.058 | 26.075 / 40.219 | 100.867 / 114.446 |
| square | 17.231 / 31.229 | 19.392 / 35.998 | 25.925 / 39.806 | 75.463 / 92.058 |
| quotient | 24.127 / 29.171 | 29.017 / 31.158 | 39.427 / 41.583 | 45.383 / 47.773 |
| remainder | 32.365 / 36.679 | 43.362 / 62.083 | 59.237 / 86.723 | 65.858 / 97.388 |
| pow3 | 20.765 / 47.129 | 30.415 / 55.865 | 43.617 / 68.975 | 287.131 / 316.660 |
| to_decimal | 45.352 / 46.569 | 89.977 / 106.252 | 133.769 / 149.383 | 671.954 / 696.202 |

| Bits | Geomean Integer/GMP | Sum-of-medians Integer/GMP | Worst |
|---|---:|---:|---|
| 53 | 1.346x | 1.234x | 2.270x pow3 |
| 128 | 1.664x | 1.517x | 2.588x add |
| 256 | 1.503x | 1.370x | 2.027x compare |
| 1024 | 1.333x | 1.125x | 2.059x negate |
| ALL 48 | **1.455x** | **1.231x** | **2.588x add at 128 bits** |

**Verified code-path interpretation (not performance causal proof):**

- `Integer` routes arithmetic through out-of-line C++ functions
  (`src/integer.cpp`), while the benchmark's direct GMP
  functions call `mpz_*` in the benchmark translation unit.
  The extra call boundary and inhibited inlining are plausible
  fixed-cost overheads, particularly for cheap operations.
- `Integer::Integer()` performs `mpz_init(_mpz)` followed by
  `mpz_set_si(_mpz,0)`: the second step is redundant because
  `mpz_init` already initializes zero. The benchmark's
  `RawInteger()` calls only `mpz_init`.
- **Compare is not operation-for-operation equivalent:**
  `IntegerOperations::compare` uses
  `a<b ? -1 : b<a ? 1 : 0`, potentially invoking
  `mpz_cmp` twice via `lt`; `GmpIntegerOperations`
  invokes `mpz_cmp` exactly once and classifies
  its result. The 1.7-2.3x comparison gap is
  compatible with, but not isolated proof of, this difference.
- `pow3` goes through `pow(a,Nat(3u))` and
  another C++ wrapper layer, while the raw path
  calls `mpz_pow_ui(...,3ul)` directly.
  Decimal output goes through `Integer::literal()`,
  an owned `new[]` intermediate and `String`;
  the direct path uses an owned vector and
  `std::string`. Decimal input is likewise
  wrapped by a `String` construction.
- Each timed iteration constructs and destroys the
  result object in both variants. Since measured
  times include this lifecycle and call overhead,
  they do **not** isolate the speed of the
  underlying multiprecision algorithms.
  The benchmark measures the *whole public operation*,
  which remains appropriate for an end-to-end target.
  Measurements are gathered GMP first and Integer
  second at each size; validate stability/control
  before attributing small residual gaps to any
  individual factor.

**Design implication.** Baseline is 23.1% slower in
the measured total across all operations and sizes
(or 45.5% slower in the equally-weighted geomean),
not uniformly 25% per operation. `BigUnsigned`
must be benchmarked separately against **raw GMP**,
with identical mathematical semantics and operand
distributions, while future native `Integer`
also needs to beat or match the public-wrapper
end-to-end baseline. Optimize algorithmic core
and API/lifetime overhead as separate phenomena.
Do not make performance claims for `BigUnsigned`
from this legacy Integer dataset. No arithmetic
or benchmark implementation changes were made
as part of this documentation entry.

### BigUnsigned B1: portable full and truncated multiplication (18 October 2026)

Implemented in the common numeric substrate, not FloatAP: `mul(a,b)` computes the exact nonnegative product; `square(a)` is currently `mul(a,a)`; `mul_low(a,b,k)` computes the product modulo `2^(64*k)` without calculating higher limbs. Inputs are immutable, output canonical; overlarge sizes raise `std::length_error`. No GMP or `__int128` in production.

**Arithmetic proof:** break 64-bit factors into 32-bit halves, form four 32x32 products, then combine their high/low halves with a bounded middle sum to obtain exact 128-bit `high:low`. In row multiplication with radix `B=2^64`, the next-limb carry arises from the high product plus two overflow bits from accumulating the existing destination limb, low product and prior carry. The full combined sum is bounded by `(B-1)^2 + 2(B-1) = B^2-1`; therefore the next carry fits 64 bits. `mul_low` ignores contributions with `i+j>=k`, so its modular result is exact. The destination is zero-initialized and disjoint from source spans.

**Tests:** 20,000 stand-alone random GMP-differential comparisons, GCC 14.2, Clang 17, and Clang ASan/UBSan passed. The integrated CTest adds 6,000 more deterministic randomized pairs, including all-one and zero limbs, varying sizes up to 34 limbs, low truncation to 0/1/4/5 or random limb lengths, full product, square and commutativity. The standalone sanitizer harness was corrected to call `mpz_clears`, after which leak checks passed. The project CI must still validate MSVC and complete CTest. No performance improvement is claimed; a specialized fast square, in-place views, scratch reuse and ARM64/x86-64 kernels remain future work. The existing Integer and FloatAP 20,000x31 benchmarks are unchanged.

### BigUnsigned B2: limb views, reusable scratch, direct-output kernels (18 October 2026)

Portable shared substrate adds `LimbView` and `MutableLimbView` (non-owning `std::span` aliases), `BigUnsigned::Scratch` with sixteen inline 64-bit limbs and vector-backed reuse for larger requests, plus `BigUnsigned::assign(LimbView)`. `assign` normalizes its result and supports self-referential subspans. `Scratch::acquire(count)` zeroes the returned range, invalidates previous views, and reuses vector capacity when sufficient; scratch belongs to the caller, without globals or concurrency sharing.

The new `add_into` and `sub_into` accept two read-only limb views and one disjoint writable destination whose length must cover the larger input. They write the low output limbs and zero higher available destination limbs, returning the final 0/1 carry or borrow. `sub_into` computes a modular subtraction even when the mathematical difference is negative; the public `sub` keeps its checked nonnegative contract. `mul_low_into` computes the exact product modulo `2^(64*output.size())`, including clearing the destination on zero input. All direct-output kernels forbid input/output overlap and allocate no internal memory; caller-provided scratch can be reused. Existing owning `add`, `sub` and `mul_low` now delegate to these same kernels. The multiplication carry proof from B1 applies unchanged.

Local deterministic differential verification: 5,000 additional GMP-oracle randomized pairs across limb boundaries, zero/all-ones operands, variable output lengths, explicit carry/borrow and result normalization, scratch capacity/reuse and subspan self-assignment. GCC and Clang strict C++20 builds as well as Clang ASan/UBSan passed locally in a shim-assisted standalone build. Full repository CI including Windows/MSVC remains the acceptance gate; no architecture-specific kernel, FloatAP, legacy Integer/Rational or frozen benchmark behavior changes. No speedup is claimed.


### BigUnsigned B3: selectable x86-64 multiplication (18 October 2026)

**Scope and identity.** B3 introduces `mul_into` (full exact product into a disjoint caller-owned buffer of at least `a.size()+b.size()` limbs), `mul_into_portable`, and `mul_low_into_portable`, alongside the existing `mul_low_into` selected kernel. Empty factors zero all output limbs; truncated output computes the exact residue modulo `2^(64*output.size())`; all kernels zero the full output and allocate no scratch. Inputs are immutable, little-endian 64-bit limbs. `mul` uses the direct full-product output kernel. No FloatAP, Integer, Rational, or existing benchmark workload is modified.

**ISA and dispatch.** `src/big_unsigned_x86-64.cpp` is compiled only on x86-64/AMD64 with GCC, Clang, AppleClang, or MSVC. Its 64x64->128 unsigned product uses GCC/Clang unsigned `__int128` in this conditional source or MSVC x64 `_umul128`; neither is required by the C++20 reference in `src/big_unsigned.cpp`. x86-64 baseline integer multiplication is sufficient: no BMI2, ADX, CPUID, runtime dispatch, assembly, or optional instruction flags. The selected path calls the architecture kernel only for two or more limbs in *both* inputs, otherwise falls back to the portable kernel. Other targets always use the reference; `has_native_multiply()` reports whether the optional implementation was built. Explicit portable entry points are available to tests and microbenchmarks. The CMake CPU and compiler conditions are a build-time feature gate.

**Proof and memory contracts.** Represent radix `B=2^64`. For one inner product, `a[i]*b[j]+out[i+j]+carry <= (B-1)^2+2(B-1)=B^2-1`; thus the carry is a 64-bit limb. Adding low-product, old destination and prior carry generates at most two carry-out flags, which are added to the high product. Full and truncated products use identical recurrence and exact modular truncation. Caller owns output, which must not overlap either source (including partial overlaps); no internal allocation, global scratch, or mutable dispatch state. Contract violations for too-small nonempty full output raise `std::length_error` before touching the result; full/low outputs for empty factors are zero. No alignment beyond naturally valid 64-bit spans is required. Each call is independent and thread-safe when buffers are disjoint.

**Verification.** The original B2 GMP-differential suite is preserved. Added 6,000 deterministic GMP-oracle product cases across zero, all-ones, unequal lengths, inline/heap boundaries and oversized destinations. The selected and explicitly portable kernels are compared word-for-word for full and low products; output overwrite, scratch reuse and insufficient full-buffer exception are exercised. A matching isolated x86-64 algorithm probe passed 45,736 edge/random cases with GCC 14.2, Clang 17, and Clang ASan/UBSan. These local checks are **not** a substitute for the project's full CMake/CTest and Windows/MSVC/Apple ARM64 validation.

**Performance and acceptance.** Previously reported x86-64 preallocated-output 20,000x31 microbenchmark medians (GCC 14.2, -O3) were portable/optimized: 53 bit 10.29/10.30 ns, 128 bit 17.97/16.14 ns, 256 bit 57.78/32.27 ns, 1024 bit 877.04/386.57 ns, 2048 bit 3576.64/1558.36 ns; these are **user-provided prior-session measurements**, not new measurements of this integration. In particular, 1024-bit speedup was 2.27x in a kernel-only comparison. The optimized path is not claimed faster end-to-end than GMP, and ARM64 has no optimized implementation or validated timing yet. Before using B3 as a FloatAP replacement, measure same-host absolute arithmetic times, portable/native/GMP and whole public FloatAP operations at 53/128/256/1024 bits; preserve all existing frozen 20,000x31 comparisons and numerical oracles. No adaptive performance crossover is claimed from these results. The source paths are `src/big_unsigned.cpp` and `src/big_unsigned_x86-64.cpp`; these reference kernels must continue to agree bit-for-bit. ISA reference: x86-64 baseline unsigned `MUL` (Intel SDM); MSVC `_umul128` compiler intrinsic contract. No third-party arithmetic source was copied.


### BigUnsigned B3: first integrated Ubuntu x86-64 benchmark and CLI alignment (18 October 2026)

**User-reported result, B3 commit `f8e56dc`.** On Ubuntu x86-64, `benchmark_big_unsigned_multiply` ran with 20,000 iterations and 31 median samples; the selected x86-64 path was active. Timings are nanoseconds per operation, with the originally frozen deterministic operand construction:

| Bits | Portable into | Selected into | Direct GMP `mpz_mul` | Owning `BigUnsigned::mul` | Selected / GMP |
|---|---:|---:|---:|---:|---:|
| 53 | 14.8784 | 15.189 | 7.8989 | 30.5604 | 1.92293x |
| 128 | 24.0364 | 20.1557 | 9.72135 | 35.493 | 2.07334x |
| 256 | 61.653 | 28.863 | 17.3753 | 58.2633 | 1.66116x |
| 1024 | 593.563 | 319.858 | 136.202 | 349.82 | 2.3484x |
| 2048 | 2317.65 | 1239.11 | 437.547 | 1254.69 | 2.83195x |

The selected kernel is faster than the portable reference at 128 bits and above (about 1.86x at 1024 bits) but remains **slower than raw GMP at every measured size** (2.348x at 1024 bits). `BigUnsigned::mul` additionally includes owning-result costs, whereas the direct GMP destination can reuse its allocation. These numbers come from the user's completed benchmark report, not a new assistant-side run, and do not establish ARM64 performance.

**Reporting-only follow-up.** Align the standalone multiplication benchmark with the other Numeric benchmarks' `CommandLineInterface` and logging format: named `--iterations/-i`, `--samples/-n`, `--bits/-b`, `--summary-only`, `--details`, `--max-slowdown` (selected/GMP), standard logging/help switches, and the legacy positional iteration/sample form. Defaults remain 20,000 / 31 and operand sizes 53/128/256/1024/2048 bits. The resulting table reports backend, operation, bits, and median ns; aggregate comparisons report equally weighted geometric mean, sum-of-medians ratio and worst slowdown separately for portable, selected and owning operations. The operation pool, deterministic inputs, warm-up, median definition, measurement order, and before/after GMP correctness checks are unchanged. Only argument handling and presentation differ. Tests and benchmarks remain separate.

### BigUnsigned B3: Ubuntu x86-64 benchmark after CLI/report integration (18 October 2026)

**User-reported same-host rerun** of the updated benchmark at `f87dcd9`, with the unchanged five deterministic operand sizes, 20,000 iterations and 31 median samples. The selected native x86-64 path was active. Units are median ns per public/kernel call; `mpz_mul` reuses an initialized GMP destination and `owned` constructs the result. These are not performance measurements obtained by the assistant.

| Bits | Portable into (ns) | Selected into (ns) | Raw GMP (ns) | Owned `BigUnsigned::mul` (ns) |
|---|---:|---:|---:|---:|
| 53 | 15.190 | 15.189 | 7.899 | 30.526 |
| 128 | 24.302 | 20.049 | 9.722 | 35.729 |
| 256 | 62.008 | 28.861 | 17.377 | 58.594 |
| 1024 | 592.303 | 319.524 | 137.275 | 349.872 |
| 2048 | 2350.366 | 1295.613 | 440.460 | 1254.717 |

Reported aggregates against direct GMP: portable geomean **3.306x**, sum-of-medians **4.968x**, worst **5.336x** at 2048 bits; selected geomean **2.142x**, sum-of-medians **2.741x**, worst **2.942x** at 2048 bits; owned geomean **3.223x**, sum-of-medians **2.823x**, worst **3.865x** at 53 bits. At 1024 bits the selected kernel is about 1.85x faster than portable yet 2.33x slower than raw GMP; at 2048 bits it is about 1.81x faster than portable yet 2.94x slower than raw GMP. Consequently no competitive-GMP claim is justified by B3.

**Interpretation constraints.** The 2048-bit owning operation appears about 3.2% faster than its nominally cheaper selected direct-output kernel despite owning a result; the measurement order, runtime state and compiler may contribute. Do not assume this reflects an allocation benefit. This benchmark uses one deterministic operand pair per size repeatedly: broaden operand patterns in a **separate diagnostic**, without changing or quietly rebasing this existing workload. Compare kernel output/preallocation, equivalent GMP output reuse and end-to-end owning lifecycles separately. Profile actual instruction-level inner loops before attempting an adaptive algorithm or new ISA optimization. Do not alter the mathematical oracle or report microbenchmark gains as full public-operation gains.


### BigUnsigned multiplication benchmark: explicit allocation/dispatch labels (18 October 2026)

Following user review of the B3 report, the benchmark labels its three BigUnsigned cases as `prealloc_portable`, `prealloc_selected`, and `new_selected`. The first calls `BigUnsigned::mul_into_portable` with caller-provided output, the second calls `BigUnsigned::mul_into` with caller-provided output, and the third calls the public `mul(a,b)`, which also dispatches to the selected architecture kernel but owns/normalizes a newly created result. The raw `mpz` baseline remains unchanged. Output column widths are adjusted for the longer labels. Row labels, corresponding diagnostic labels and aggregate names change, **not** operand generation, timing order, tested algorithms, mathematical oracle, warm-up, iterations, sample counts, allocation or measurement semantics. `--max-slowdown` still gates the preallocated selected-path slowdown against raw GMP.


### BigUnsigned B4 experiment: x86-64 column-wise Comba candidate (18 October 2026)

**Experimental only:** `BigUnsigned::mul_low_into_comba` and `mul_into_comba` provide a forced, allocation-free comparison kernel for the same disjoint read-only limb views and caller-owned output views as B3. The selected production `mul_into`, `mul_low_into`, owning `mul` and `square` continue to use exactly the B3 dispatch. The Comba candidate is implemented only in `src/big_unsigned_x86-64.cpp`, compiled behind the existing x86-64 build gate; on other platforms its explicit entry point delegates to the portable reference. No runtime feature probing, optional CPU extensions or architecture-specific floating-point behavior is introduced.

**Arithmetic and carry proof.** Column `k` sums all partial products `a[i]*b[k-i]` where `0<=i<a.size()` and `0<=k-i<b.size()`, plus the carry from the prior column. The sum is held in three unsigned 64-bit words `(top,hi,lo)` (radix `B=2^64`). Each exact product `(product_hi,product_lo)` adds to `(hi,lo)` using 64-bit wrapping additions, with each carry-out propagated to the next word. For at most `min(a.size(),b.size())` products in any column, the three-word sum is sufficient for any representable input spans: the worst case is `min(n,m)*(B-1)^2` plus carry of at most `min(n,m)*(B-1)`, giving a high-word count bounded by the number of products and representable span lengths. After emitting `lo` to output[k], shift `hi` and `top` down to serve as the following column's initial carry, and zero `top`. The algorithm computes the exact full result or its low `output.size()` limbs; output is completely zeroed beforehand and nonoverlapping input/output is required. No extra scratch, heap allocation, or mutable global state; `mul_into_comba` checks full-buffer length before writing and handles empty operands like B3. x86-64 wide products use the same baseline `MUL` lowering via conditional GNU/Clang unsigned `__int128` or MSVC `_umul128`; no additional ISA feature is needed.

**Tests and portability.** The pre-existing 6,000-case deterministic GMP-oracle B3 test now additionally compares forced Comba and the portable reference word-for-word for both full and modular-truncated products, including zero/all-ones operands, uneven limb sizes, oversized outputs, scratch reuse and insufficient full-buffer exceptions. On non-x86-64 this also exercises the documented reference fallback. A separate Linux x86-64 C++20 probe (25,404 edge/random cases with GMP, containing two-row schoolbook variants and Comba) passed with GCC and Clang, including Clang ASan/UBSan; those are *probe* results, not substitute for integrated CMake/CTest/Windows CI.

**Performance and report controls.** The benchmark retains default 20,000 iterations, 31 samples and the existing 53/128/256/1024/2048-bit operand set. By default the four historical modes `prealloc_portable`, `prealloc_selected`, `mpz` and `new_selected` are unchanged in implementation and timing order. Optional `--include-comba` (x86-64 builds only) first verifies the Comba result against GMP, then measures `prealloc_comba` *after* the historical four modes, and includes its own absolute rows and geomean/total_time/worst comparison against GMP. The default benchmark must remain comparable with the B3 user-reported Ubuntu timings. An isolated probe with GCC/Clang showed inconsistent crossovers between Comba and existing schoolbook depending on compiler and system state; do not claim Comba faster from those exploratory times. Require same-host integrated measurements and green cross-platform CI before considering any production dispatch change. User-facing mode labels follow B3 naming.

**Future acceptance gate.** Keep B3 as production until repeated same-host absolute tests show a significant, stable advantage on realistic inputs, compare GMP directly, preserve correctness tests and FloatAP's frozen performance gates, and document the measured crossover if one exists.


### BigUnsigned B4 measurements and B5 fused-accumulation experiment (18 October 2026)

**B4 user-reported Ubuntu x86-64 measurements (20,000 iterations, 31 samples).** All timings in ns per multiplication, using the frozen deterministic inputs. For 53/128/256/1024/2048 bits respectively: `prealloc_selected` 17.010/19.570/29.800/322.086/1273.343; `prealloc_comba` 15.085/21.436/33.308/330.851/1365.606; direct GMP `mpz` 7.595/9.724/17.154/136.058/443.277; `new_selected` 31.895/35.809/57.149/352.058/1254.529; `prealloc_portable` 15.705/22.587/60.751/591.501/2301.793. Relative to GMP, the `prealloc_selected` geomean is 2.214x, total_time 2.707x, worst 2.873x at 2048 bits; Comba geomean 2.295x, total_time 2.878x, worst 3.081x. Comba is faster only at 53 bits and slower at the other four sizes, so **B3 remains the production selector**. Windows, Unix and Coverage CI all passed on B4 commit `47776c9`.

**B5 hypothesis.** A native x86-64 schoolbook inner step can be expressed as a single unsigned 128-bit mathematical sum `total=(Wide)a[i]*b[j]+out[i+j]+carry`; write its low 64 bits to the destination and retain the high 64 bits as the next carry. The bound `(B-1)^2+2(B-1)=B^2-1` for radix `B=2^64` proves that the intermediate is exactly representable in 128 bits (no overflow). This differs from B3's explicitly separated 64-bit low/high product and two overflow flags only in how the compiler sees the accumulation; it does not change the arithmetic algorithm or operation count. Conditional x86-64 GCC/Clang use `unsigned __int128`, while MSVC x64 delegates this optional entry point to its existing, verified native B3 multiply, and non-x86-64 platforms delegate to the portable C++20 reference. No BMI2/ADX, runtime CPU dispatch, global scratch or new allocation is introduced. Both input/output aliasing and full/truncated destination constraints remain those of B3.

**B5 validation.** A standalone C++20 Linux x86-64 diagnostic built with GCC 14.2 and Clang 17 and verified 12,000 deterministic random/adversarial cases (including all-ones, zero, unequal input and output sizes, modular low truncation) against direct GMP. GCC -O3 20,000x31 medians (old/new, ns): 128 bits 13.94/10.92; 256 bits 29.40/19.68; 1024 bits 388.33/249.47; 2048 bits 1559.72/958.70. Clang -O3 medians (old/new, ns): 128 bits 11.51/12.45; 256 bits 25.87/27.60; 1024 bits 283.10/292.31; 2048 bits 1154.16/1114.53. **These are standalone probe measurements, not an integrated project or GMP speedup claim**. Compiler-dependent results motivate a gated benchmark comparison before any production path change.

**B5 implementation and benchmark protocol.** Expose optional `mul_low_into_fused` and `mul_into_fused`, with no changes to `mul_low_into`/public owning `mul`. Extend integrated 6,000-case CTest to compare the B5 forced path against the GMP oracle and the B3 portable reference, including full, low, extra output words and size errors. Add `--include-fused` to the existing exact multiplication benchmark, which produces a separate `prealloc_fused` row and GMP aggregates only when explicitly requested. Keep 20,000x31 defaults and the order/content of the four B3 baseline measurements intact; experimental Comba is likewise optional. A local 12,000-case oracle does not establish full project CI correctness, nor does performance on one compiler establish a dispatch crossover. Gate any future selector change on successful Unix/Windows/Coverage CI and repeated same-host GCC/Clang measurements. Both shared design journals must remain synchronized.


### BigUnsigned B6: GCC x86-64 selected fused multiplication (18 October 2026)

**Decision and measurements.** A user-reported integrated Ubuntu x86-64 20,000-iteration/31-sample B5 benchmark measured `prealloc_selected` vs forced `prealloc_fused` medians (ns): 53 bits 15.442/13.076, 128 bits 18.834/13.169, 256 bits 29.297/20.148, 1024 bits 323.861/205.060, 2048 bits 1263.749/807.183. Raw GMP times were 8.056/9.425/18.946/136.262/437.660 ns. The corresponding selected vs GMP aggregate geomean and sum-of-medians slowdowns were 2.098x/2.705x, while forced fused gave 1.463x/1.734x. Fused was faster at all five recorded sizes, by about 1.58x at 1024 bits and 1.57x at 2048 bits; **GMP remains faster at every tested width**. These are user-reported integrated figures on Ubuntu GCC x86-64, not an independent cross-compiler acceptance claim. Clang probe performance in the B5 journal was inconsistent; no Clang/MSVC/ARM64 speedup may be inferred.

**B6 build and dispatch.** When `CMAKE_SYSTEM_PROCESSOR` identifies x86-64/AMD64 and `CMAKE_CXX_COMPILER_ID` is `GNU`, the CMake option `NUMERIC_BIG_UNSIGNED_GCC_FUSED_MULTIPLY` defaults to `ON` and defines the private compile-time `ARIADNE_BIG_UNSIGNED_GCC_FUSED` switch. `BigUnsigned::mul_low_into` then selects the proven B5 unsigned `__int128` fused inner loop, including one-limb cases. The public owning `mul` and `mul_low` inherit that selection through the existing `mul_into` and `mul_low_into` entry points. GCC x86-64 with the option `OFF`, Clang/AppleClang x86-64 and MSVC x64 retain the B3 selected kernel and its one-limb portable fallback. Other architectures retain the C++20 portable reference. There is no runtime CPUID, nonbaseline instruction, ABI change, shared mutable state, input/output aliasing or scratch allocation. B5's carry bound `(B-1)^2+2(B-1)=B^2-1` applies without modification. The `selected_multiply_uses_fused()` query reports the build-time selection.

**Controlled comparison and regression.** `mul_low_into_schoolbook` and `mul_into_schoolbook` explicitly force the former B3 selected schoolbook (or portable fallback elsewhere), in addition to pre-existing portable, Comba and fused explicit kernels. The existing 6,000-case GMP-differential CTest now verifies forced B3 full and truncated kernels against the same oracle and selected path, including empty/all-ones/variable operands, larger output buffers, and insufficient-full-buffer exception. `--include-schoolbook` adds an opt-in `prealloc_schoolbook` row after the four standard report entries and any other experimental comparisons; normal benchmark defaults, operands, measure order and 20,000x31 counts stay unchanged. `prealloc_selected` and `new_selected` **now legitimately use fused on GCC x86-64 default builds**; the older B3 timing is available via explicit schoolbook and an option-OFF control build. Note this changes the meaning of `selected` relative to previous commits but not what its name promises: selected is always architecture/compiler-dependent. A same-host A/B should compare the exact library operations (including owning paths) and raw GMP. Acceptance of B6 requires complete Unix/Windows/Coverage CI, ideally a matched GCC option-OFF CTest and repeat GCC/Clang/Apple ARM64 benchmark checks before any further promotion. B6 does not migrate Integer/Rational or FloatAP.


### BigUnsigned B6: selected fused Ubuntu x86-64 integrated results (18 October 2026)

**User-provided measurements** for B6 commit `63ce810` (GCC x86-64 with `selected kernel=fused`, 20,000 iterations and 31 median samples). Units: ns per multiplication, fixed operand pool unchanged since B3.

| Bits | prealloc_portable | prealloc_selected | mpz (reused destination) | new_selected | prealloc_fused | prealloc_schoolbook |
|---:|---:|---:|---:|---:|---:|---:|
| 53 | 15.370 | 13.385 | 7.595 | 26.008 | 13.376 | 16.523 |
| 128 | 24.312 | 14.277 | 9.425 | 30.102 | 13.682 | 20.049 |
| 256 | 45.406 | 20.274 | 17.385 | 47.616 | 19.606 | 30.137 |
| 1024 | 591.746 | 204.928 | 136.173 | 239.030 | 203.752 | 320.640 |
| 2048 | 2290.353 | 757.071 | 442.111 | 785.913 | 755.872 | 1270.594 |

Reported prealloc_selected vs reused mpz: geomean **1.517x**, total_time **1.648x**, worst **1.762x at 53 bits**. Reported new_selected vs reused mpz: geomean **2.478x**, total_time **1.842x**, worst **3.424x at 53 bits**. Forced schoolbook vs reused mpz: geomean **2.223x**, total_time **2.706x**, worst **2.874x at 2048 bits**. A B5 user-provided earlier session measured new_selected 350.572 ns at 1024 bits and 1278.520 ns at 2048 bits, against B6 239.030 and 785.913 ns respectively; these are different runs, not a strictly paired old/new owning benchmark. B6 selected and forced fused agree closely, validating routing/performance in the integrated Ubuntu benchmark. No inference for Clang/MSVC/ARM64; GMP remains faster at each size.

**Validation.** GitHub Actions Unix, Windows and Coverage succeeded for `63ce810`. Independent B6 user results, not assistant-side execution.

### BigUnsigned B7: opt-in matched-result-lifetime GMP comparison (18 October 2026)

B6's default `mpz` reuses an already initialized and grown GMP destination across iterations while `new_selected` creates and destroys a fresh owning BigUnsigned. Both statistics are valid but do **not** measure identical result lifetimes. The optional `--include-new-mpz` diagnostic adds a `new_mpz` case, constructing a fresh `mpz_t` via `mpz_init`, multiplying the *same* preconstructed GMP inputs via `mpz_mul`, and destroying the new result via `mpz_clear` on **each timed iteration**. Its nonzero sanity guard parallels the existing owning BigUnsigned case. Before and after timing, the new-result GMP calculation is checked for exact identity against the existing GMP reference, outside the timed loops. This provides a second comparison of `new_selected/new_mpz` with matched high-level result lifetimes, but different storage representations and implementation internals; it is not a claim of perfectly equivalent allocator behavior or GMP performance. The report adds a separate `new_selected / new_mpz` geomean/total_time/worst aggregate when requested, while preserving the historical raw-mpz ratios and `--max-slowdown` gate.

The existing four default modes, their order, operand generation, warm-up, sample counts, reference check, selected-kernel dispatch, and existing optional experiments are unchanged. The `new_mpz` measurement is placed last, so default output is bit-for-bit semantically unchanged except for unavoidable runtime variance. No production code, Rational or FloatAP paths are changed. Validate the optional benchmark on Ubuntu and retain CI as a compile/regression gate.


### BigUnsigned B7: Ubuntu x86-64 result-lifetime benchmark (18 October 2026)

**Measurement source.** User-reported integrated benchmark run on Ubuntu x86-64 with selected GCC fused kernel; 20,000 iterations, 31 samples, unchanged deterministic inputs and original mode timing order; optional `--include-new-mpz --include-schoolbook --include-fused`. Values below are medians in ns per multiplication. This is not an independent assistant-side timing or a measurement of ARM64/Clang/MSVC.

| Bits | prealloc_selected | reused mpz | new_selected | new_mpz | new_selected / new_mpz |
|---:|---:|---:|---:|---:|---:|
| 53 | 13.099 | 7.837 | 27.080 | 28.250 | 0.959x |
| 128 | 14.885 | 10.626 | 31.335 | 28.722 | 1.091x |
| 256 | 21.273 | 17.157 | 47.803 | 29.632 | 1.613x |
| 1024 | 205.278 | 138.559 | 233.324 | 149.986 | 1.556x |
| 2048 | 759.731 | 440.800 | 786.426 | 451.984 | 1.740x |

The additional reported modes (53/128/256/1024/2048 bits) were portable 15.370/23.819/45.232/591.363/2304.662 ns, forced fused 13.081/13.682/19.615/203.619/756.209 ns, and forced B3 schoolbook 17.314/20.246/30.337/320.331/1221.233 ns. Reported aggregates: `prealloc_selected / mpz` geomean **1.493x**, sum-of-medians **1.649x**, worst **1.724x** at 2048 bits; `new_selected / mpz` geomean **2.433x**, sum-of-medians **1.831x**; `new_selected / new_mpz` geomean **1.355x**, sum-of-medians **1.635x**, worst **1.740x** at 2048 bits. All three GitHub Actions workflows (Unix/Windows/Coverage) passed on B7 commit `6722b74`.

**Diagnostic decomposition, not a causal profiling proof.** Subtracting each independently sampled median gives the nominal cost of creating rather than reusing the result: `new_selected-prealloc_selected` is 13.981/16.450/26.530/28.046/26.695 ns across the five sizes; `new_mpz-mpz` is 20.413/18.096/12.475/11.427/11.184 ns. The observed arithmetic/preallocated gap `prealloc_selected-mpz` is 5.262/4.259/4.116/66.719/318.931 ns. Thus at 256 bits lifecycle overhead is comparatively important, but at 1024 and particularly 2048 bits most of the absolute gap is already present with reusable output buffers. These differences of medians are descriptive, not an isolated attribution of allocator vs canonicalization, and paired timing order and compiler optimizations remain potential confounders.

**B8 direction (not yet implemented).** Preserve the fully validated GCC fused B6 selector. Investigate result-lifetime/normalization costs separately from arithmetic on shorter sizes; at 1024–2048 bits explore genuinely faster multiplication algorithms or loop structure against exact GMP oracles. A dedicated operand-distribution diagnostic should span balanced and unbalanced lengths and avoid quietly replacing the five-case frozen benchmark. Do not widen `BigUnsigned` inline storage globally or change public representation solely to optimize the 256-bit benchmark: that affects every `Integer`/`Rational`/future `FloatAP` object. Avoid claiming a GMP competitive implementation before same-host `prealloc_selected/mpz` and `new_selected/new_mpz` gates pass.


### BigUnsigned B8: opt-in allocation-free Karatsuba experiment (18 October 2026)

**Decision and scope.** B7 showed that a remaining 1024/2048-bit arithmetic gap versus raw GMP persists even with reusable output. B8 therefore tests an independently implemented Karatsuba full-product kernel, while **preserving production B6 selected fused GCC x86-64, existing B3/B4/B5 forced kernels and all four default benchmark modes**. The experimental `BigUnsigned::mul_into_karatsuba(a,b,out,scratch)` computes an exact full product into an output view of at least `a.size()+b.size()` limbs. Inputs, output and scratch must be mutually disjoint; all result words including excess output are overwritten; neither this kernel nor recursive calls allocates. `karatsuba_scratch_limbs(a.size(),b.size())` computes a checked upper bound on required reusable 64-bit scratch words. Insufficient output or scratch raises `std::length_error` before modifying the destination (nonempty inputs). An empty factor zeroes the destination. No GMP dependency enters production.

**Algorithm and arithmetic invariant.** Split at `k=floor(max(n,m)/2)`; the fast path applies only when both operands have nonempty upper halves and `min(n,m)>8` limbs. Let `z0=a0*b0`, `z2=a1*b1` and `t=(a0+a1)(b0+b1)-z0-z2=a0*b1+a1*b0`. Recursively write z0 and z2 directly into nonoverlapping output slices. Use two scratch-resident half sums, each with one extra high limb, and one scratch-resident cross-product buffer. Subtractions are exact nonnegative operations and cross-term is added to the shifted destination with carry propagation. The recursion scratch requirement is `2*(sa+sb)+max(work(z0),work(z2),work(sumA*sumB))`, where `sa=max(|a0|,|a1|)+1`, `sb=max(|b0|,|b1|)+1`; recursive calls reuse scratch sequentially. Every sum and carry is computed from unsigned 64-bit limbs with overflow checks; the full output is exact at radix 2^64. The recursion's base case uses the current selected (not necessarily portable) B6 schoolbook multiply. This is a full-product implementation only: no general truncated Karatsuba path or production crossover is claimed. The file `src/big_unsigned_karatsuba.cpp` is portable C++20 and requires no `__int128`, third-party source or CPU-specific instructions.

**Verification.** A local standalone prototype, independently compared against GMP and the selected fused reference, passed **47,116** deterministic edge/random test cases with GCC 14.2 and Clang 17, including zero limbs, all-ones, sparse/zero operands, balanced/asymmetric dimensions up to 128 limbs, and output guards. Clang ASan/UBSan passed **4,616** focused differential cases. Project CTest is extended with 6,000 preexisting randomized GMP comparison cases and 700 wider cases (up to 128 limbs) for the B8 full kernel, plus scratch insufficiency and output-guard tests. Those integrated modifications still require Unix/Windows/Coverage CI acceptance. No mutation of the production selector or test oracle is permitted to rescue an unsuccessful experiment.

**Isolated performance (not project results).** A standalone 20,000-iteration/31-sample GCC x86-64 microbenchmark recorded fused/Karatsuba (ns): 1024 bits 316.92/410.54, 2048 bits 1185.21/1480.69, 4096 bits 4703.88/5086.85, 8192 bits 19268.6/16420.2. Clang on the same probe recorded 1024 bits 290.36/398.26, 2048 bits 1078.98/1535.44, 4096 bits 4306.24/5110.39, 8192 bits 17939/16803.1. These are **exploratory isolated measurements**, not the unchanged integrated Numeric benchmark, and show that this initial implementation is slower through at least 4096 bits; only at 8192 bits did a modest speedup appear. The crossover will depend on compiler, input distribution and runtime environment; there is insufficient evidence to make Karatsuba production-default, particularly for the project's 53/128/256/1024/2048-bit performance gate.

**Benchmark integration.** The optional CLI flag `--include-karatsuba` appends a `prealloc_karatsuba` row and a GMP slowdown aggregate **after** the frozen default measurements and existing optional experiments. Its scratch is sized/acquired *before timing*, then reused. The input operand pair, median protocol, 20,000x31 defaults, raw GMP and owning lifetimes, and `--max-slowdown` semantics remain unchanged. For assessing large-input crossover, pass `--bits 1024,2048,4096,8192` explicitly and check `prealloc_karatsuba` versus `prealloc_selected` on the same Ubuntu host. Do not conflate a high-bit diagnostic with replacement of frozen baselines or a whole FloatAP performance result. B8's purpose is to establish evidence, not enable subquadratic dispatch prematurely.


### BigUnsigned B8 build fix: testing macro shadow (18 October 2026)

Ubuntu GCC found an integrated CTest **test-source compilation error**, not evidence of an arithmetic failure, in the newly added B8 wide-limb loop: `ARIADNE_TEST_ASSERT(...)` internally declares `auto result=(expression)`, so a local scratch view named `result` caused name shadowing and use-before-deduction when passed inside the assertion. Rename that local view consistently to `product_view` in the allocation, Karatsuba call, GMP check, and high-word zero assertion. Preserve unchanged the operand generation, output size, GMP oracle, and all B8 assertions; do not weaken test coverage. The build issue was reported with GCC `-Wshadow`, and the full Unix/Windows/Coverage CI for the corrected commit remains the acceptance gate.


### B8.1: Karatsuba cutoff experiment (18 October 2026)

The user-reported Ubuntu 20,000-iteration, 31-sample timings show Karatsuba cutoff 8 slower than fused selected at 1024 (346.542 vs 206.249 ns), 2048 (1256.520 vs 797.097 ns), 4096 (4203.549 vs 2975.416 ns), and 8192 bits (13583.344 vs 12589.667 ns). This contradicts the isolated-prototype crossover suggestion; no production promotion is justified. B8.1 adds parameterized scratch sizing and Karatsuba multiplication at cutoff >=8, preserving the original default cutoff 8 and unchanged GCC fused dispatch. --include-karatsuba --karatsuba-cutoff N reports the threshold and compares to GMP and fused at fixed sizes. The widened random GMP oracle covers cutoff 8,16,32,64. Same threshold must be used for sizing and multiplication. No runtime workspace allocations, production changes, or FloatAP changes. The CI still needs to be checked for the new commit.


### BigUnsigned B8.1: complete Ubuntu cutoff sweep, 8/16/32/64 (18 October 2026)

**Source:** four user-run integrated Ubuntu GCC x86-64 benchmarks on commit `56453a4`, 20,000 iterations, 31 samples, balanced 1024/2048/4096/8192-bit operands, units ns. Each column is from a *different* benchmark invocation; comparison with selected must use the corresponding run rather than a cross-run baseline.

| Bits | Karatsuba cutoff=8 | cutoff=16 | cutoff=32 | cutoff=64 |
|---:|---:|---:|---:|---:|
| 1024 | 357.162 | 218.163 | 215.290 | 215.269 |
| 2048 | 1293.077 | 881.020 | 815.782 | 829.656 |
| 4096 | 4346.714 | 2967.023 | 2629.425 | 2985.271 |
| 8192 | 13940.751 | 9652.585 | 8461.684 | 9178.001 |

Corresponding preallocated selected fused times by cutoff-run (8/16/32/64, ns): 1024 209.197/205.349/205.269/205.131; 2048 761.314/807.764/805.034/820.501; 4096 2966.787/2969.927/2971.200/2969.405; 8192 12464.999/12257.438/12257.133/12255.704. At cutoff=32, compared *within the same run*, Karatsuba is slower by 4.88% (1024) and 1.34% (2048), and faster by 11.50% (4096) and 30.97% (8192). With cutoff=64 the 4096-bit case falls back to nonrecursive selected multiplication (observed 2985.271 vs selected 2969.405 ns), while 8192 bits still benefit from a top-level split (9178.001 vs selected 12255.704 ns). Cutoff 32 is the best measured candidate among these four; 1024-cutoff 32 and 64 results are indistinguishable at this precision.

**Decision:** retain fused B6 production dispatch, retain Karatsuba experimental-only. The crossover is currently bracketed between 2048 and 4096 bits *for this balanced-input GCC x86-64 measurement configuration*, not established as an architecture-general threshold. Repeat same-host runs with alternating ordering and intermediate operand sizes near 2560/3072/3584 bits, and consider asymmetric operand widths before any production change. The range in selected 2048-bit timings (761.314 to 820.501 ns across runs) is a reminder that a single crossover estimate can be noise-sensitive. No empirical evidence for macOS ARM64, Clang or Windows optimization. GitHub Actions Unix, Windows and Coverage all **passed** for B8.1 commit `56453a4`. Test correctness and production paths remain unchanged.


### B9 ARM64 native fused multiplication candidate (18 October 2026)

An architecture-specific AArch64 multiplication kernel is added in `src/big_unsigned_arm64.cpp`, compiled for aarch64/arm64 with GCC, Clang or AppleClang. Its explicit `mul_into_fused` path uses unsigned 128-bit multiplication and accumulation, lowering to baseline ARM64 `MUL`, `UMULH` and carry operations; compilation to those instructions was checked with Clang 17 targeting AArch64. This confirms code generation, **not performance**. The arithmetic bound `(2^64-1)^2+2(2^64-1)=2^128-1` proves each wide intermediate fits in 128 bits. Truncated and full multiplication, disjoint inputs and outputs, fully overwritten result, and no internal allocations preserve the B5 semantics. The existing GMP differential tests already exercise forced fused full and low paths on ARM64. `has_native_multiply` now means x86-64 OR ARM64 native support, while `has_x86_64_multiply` specifically identifies availability of Comba/B3 schoolbook. The benchmark permits `--include-fused` on ARM64 but correctly restricts x86-specific modes.

**Dispatch discipline:** ARM64 production selection remains the portable kernel by default. The CMake option `NUMERIC_BIG_UNSIGNED_ARM64_FUSED_MULTIPLY=ON` enables fused ARM64 dispatch explicitly for controlled experiments; the default is OFF, as there are no integrated user macOS ARM64 speed measurements. GCC x86-64 fused production selection and Windows/non-ARM platforms remain unchanged. The default four benchmark modes remain exactly the same; on macOS ARM64 run `benchmarks/benchmark_big_unsigned_multiply --include-fused --include-new-mpz -t dark` and compare `prealloc_portable`, `prealloc_selected`, `prealloc_fused`, `mpz`, and owning results. Accept production promotion only after correct macOS CI/local CTest and repeated same-host measurements demonstrate stable benefit. No FloatAP code change or ARM64 performance claim.


### B9 native Apple Silicon ARM64 benchmark, user-reported (18 October 2026)

**Source and configuration.** User ran the integrated `benchmark_big_unsigned_multiply` on macOS ARM64 from commit `3b34604` with `--include-fused --include-new-mpz`, selected kernel reported **portable**, meaning the new native ARM64 fused kernel was only force-invoked and `NUMERIC_BIG_UNSIGNED_ARM64_FUSED_MULTIPLY` was **OFF**. The benchmark used 20,000 iterations, 31 samples, the five frozen operand sizes, and unchanged deterministic inputs. All times in ns per multiplication:

| Bits | prealloc_portable | prealloc_selected | mpz (reused output) | new_selected | prealloc_fused ARM64 | new_mpz |
|---:|---:|---:|---:|---:|---:|---:|
| 53 | 13.735 | 12.242 | 5.802 | 10.842 | 10.131 | 22.408 |
| 128 | 13.973 | 13.096 | 25.671 | 12.069 | 10.569 | 22.846 |
| 256 | 26.977 | 26.977 | 13.463 | 44.527 | 16.592 | 26.798 |
| 1024 | 348.967 | 349.917 | 95.004 | 368.960 | 180.819 | 109.246 |
| 2048 | 1380.987 | 1379.546 | 312.606 | 1396.992 | 802.565 | 337.758 |

Same-run portable-selected to forced-fused speedups: **1.208x**, **1.239x**, **1.626x**, **1.935x**, **1.719x** at the respective sizes, equivalently -17.2%, -19.3%, -38.5%, -48.3%, -41.8% execution time. Thus the explicit ARM64 `__int128` kernel has an important integrated speedup over the baseline portable path on the user's Apple Silicon machine, especially at 1024/2048 bits. Reported per-mode GMP aggregates: selected/MPZ geomean 2.037x, total_time 3.937x; forced fused/MPZ geomean 1.341x, total_time 2.255x; new_selected/new_mpz geomean 1.428x, total_time 3.532x. **Caveats:** GMP's reused-output 128-bit 25.671 ns exceeds its newly owned 22.846 ns, and at 53/128 bits `new_selected` apparently outpaces `prealloc_selected`. Treat these cases as possibly affected by compiler optimization, timing noise or benchmark design; do not assume isolated cause or universal advantage. The fused ARM64 1024/2048-bit reference ratios were about 1.90x/2.57x against reused-output GMP; GMP still outperforms at those larger widths. The commit `3b34604` passed Unix/Windows/Coverage CI, but those workflow results alone do not establish the new kernel was exercised on macOS ARM64.

**Next gate, not yet a default-dispatch change.** Repeat the same user macOS ARM64 benchmark using a separate build configured with `-DNUMERIC_BIG_UNSIGNED_ARM64_FUSED_MULTIPLY=ON` and run project correctness tests on that build. Confirm the benchmark banner `selected kernel=fused`, compare `prealloc_selected` to forced `prealloc_fused` and the newly selected owning `new_selected` to `new_mpz`, retaining `--include-fused --include-new-mpz` with original five sizes. Prefer multiple repeated paired runs before enabling ARM64 fused by default. The present OFF build and all non-ARM64 production code remain unchanged.


### B9 Apple Silicon ARM64 fused selector ON: user-reported owning results (18 October 2026)

**Source:** User reran the integrated benchmark on macOS Apple Silicon, 20,000 iterations, 31 median samples, five original operand sizes, after setting `NUMERIC_BIG_UNSIGNED_ARM64_FUSED_MULTIPLY=ON` and rebuilding. Banner explicitly reported `native available=true | selected kernel=fused`. Units ns; matched rows within this run:

| Bits | prealloc_portable | prealloc_selected (fused) | mpz (reused) | new_selected | prealloc_fused (forced) | new_mpz |
|---:|---:|---:|---:|---:|---:|---:|
| 53 | 9.927 | 9.552 | 5.235 | 11.065 | 9.546 | 24.262 |
| 128 | 15.379 | 12.554 | 10.354 | 11.508 | 11.496 | 24.556 |
| 256 | 27.723 | 16.600 | 13.217 | 32.946 | 16.592 | 26.817 |
| 1024 | 350.102 | 180.769 | 95.121 | 202.990 | 180.544 | 109.433 |
| 2048 | 1393.152 | 805.240 | 312.863 | 836.858 | 802.896 | 336.356 |

**Comparison against previous OFF run on the same user's Mac** (separate benchmark invocations, not paired trials): `new_selected` OFF=10.842/12.069/44.527/368.960/1396.992 ns versus ON=11.065/11.508/32.946/202.990/836.858 ns at 53/128/256/1024/2048 bits. ON/OFF speedup is approximately 0.98x/1.05x/1.35x/1.82x/1.67x. Small-input 53- and 128-bit owning differences are not reliable signals; the 256-bit and especially the 1024/2048-bit gains are large. Within the ON run, `prealloc_selected` and `prealloc_fused` match closely at 53,256,1024,2048 bits; the ~9.2% difference at 128 bits suggests microbenchmark variability or call-path effects, not a distinct algorithm. The earlier OFF-run GMP 128-bit anomaly (25.671 ns raw vs 22.846 ns new) did not recur: the ON-run raw GMP was 10.354 ns. These discrepancies underscore the need for repeated measurements and ensuring compiler optimizations do not render fixed-operand benchmarks unrepresentative.

**Remaining GMP gap and acceptance.** User-reported ON-run geomean/total_time: selected preallocated vs reused GMP 1.685x/2.346x, new_selected vs reused GMP 2.017x/2.508x, and new_selected vs new_mpz 1.039x/2.101x, worst 2.488x (2048 bits). At 2048 bits, even an owning-to-owning comparison is ~2.49x slower than GMP. These are `BigUnsigned` measurements, **not** FloatAP end-to-end benchmarks. GitHub Actions Unix/Windows/Coverage are green on B9 `3b34604` and documentation commit `e827a57`; macOS ARM64 project CTest in the selected-fused ON build has not been reported and therefore must not be claimed. Keep ARM64's default option OFF pending a clean macOS ARM64 test run and ideally repeated same-host measurements. Candidate for a future default change is now supported by *real* user integrated timings, not by a standalone prediction.


### B9.1: promote native ARM64 fused multiplication without a CMake toggle (18 October 2026)

**Acceptance:** The user reported the macOS ARM64 project tests green with fused dispatch active. Paired selected/forced-fused integrated measurements on Apple Silicon were approximately 16.600/16.592 ns at 256 bits, 180.769/180.544 ns at 1024 bits and 805.240/802.896 ns at 2048 bits; the owning public API also improved materially compared with the earlier portable selection. The earlier B9 notes describe the historical OFF/ON experiments; those settings are now **retired**.

**Final production design:** For supported native ARM64 GNU/Clang/AppleClang builds, `ARIADNE_BIG_UNSIGNED_ARM64` alone enables compilation *and* selection of the fused ARM64 kernel. Remove both the public CMake option `NUMERIC_BIG_UNSIGNED_ARM64_FUSED_MULTIPLY` and the redundant private macro `ARIADNE_BIG_UNSIGNED_ARM64_FUSED`; no transitional flag, runtime toggle, or compatibility branch is retained. The separately exposed `mul_into_portable` and `mul_low_into_portable` remain sufficient for reference comparisons. A direct ARM64 test assertion now checks that the native fused path is selected. x86-64 GCC retains its independent selector/option, while other x86 compilers keep B3 and unsupported platforms remain portable. No changes to FloatAP code or claims about FloatAP performance. New revision requires Unix/Windows/Coverage CI verification; the reported macOS tests covered the prior option-enabled but arithmetically identical implementation.


### B10 ARM64 explicit carry experiment (18 October 2026)

New B10 experimental AArch64 row-schoolbook multiplication `mul_low_into_arm64_carry` / `mul_into_arm64_carry`: baseline instruction sequence MUL, UMULH, ADDS, ADC with two carry propagations per partial product, preserving the 128-bit bound and allocation-free full/truncated semantics. Clang 17 assembled an isolated instruction probe for both AArch64 Linux ELF and Apple arm64 Mach-O; **no native Apple Silicon timing has been obtained**. The production fused ARM64 selector is unchanged. The forced kernel uses portable fallback off ARM64, the integrated `--include-arm64-carry` benchmark mode requires ARM64 and appends an optional comparison. Existing randomized GMP-differential full/truncated, overflow and buffer-size tests are extended to this mode. Any promotion requires correctness checks and same-host user benchmark. From the macOS build directory: `benchmarks/benchmark_big_unsigned_multiply --include-arm64-carry --include-new-mpz --bits 53,128,256,1024,2048,4096 -t dark`. Compare `prealloc_selected` against `prealloc_arm64_carry`. No performance improvement is asserted.


### B10.1: source layout consolidation (18 October 2026)

Per repository architecture rule, only architecture-specific implementations have dedicated translation units. The experimental allocation-free Karatsuba implementation, its cutoff selection, sizing, recursive helpers and public method definitions have been moved **without algorithmic changes** into `src/big_unsigned.cpp`; `src/big_unsigned_karatsuba.cpp` is deleted and removed from `src/CMakeLists.txt`. Its GMP-differential tests and integrated cutoff benchmark are unchanged. GCC/Clang's identical `unsigned __int128` fused accumulation implementation is now maintained **once**, as `multiply_limbs_fused` in the shared `src/big_unsigned.cpp`. The duplicate bodies and entry points in `src/big_unsigned_x86-64.cpp` and `src/big_unsigned_arm64.cpp` are removed. ARM64 source retains only its explicit native `MUL/UMULH/ADDS/ADC` experiment; x86 source retains B3 schoolbook/Comba and the MSVC `_umul128` implementation. The shared fused helper keeps the x86-64 MSVC fallback to B3, the forced fused wrapper's architecture availability, and GCC x86-64/ARM64 default dispatch unchanged. On toolchains without `__int128`, the helper's fallback remains portable; unsupported platforms select the existing portable reference. This is a structural refactor, **not** a claimed speedup: recheck actual compiler output and same-host performance after consolidation, as moving a kernel between translation units can affect optimizer behavior. The user requires no nonarchitecture extra .cpp files for multiplication algorithms; do not recreate a standalone Karatsuba source. CI and macOS ARM64 test/benchmark results for this new layout are not yet verified.


### B10.1 Apple Silicon ARM64 carry-chain user benchmark (18 October 2026)

The user ran the integrated benchmark from the B10.1 consolidated source layout with `--include-arm64-carry --include-new-mpz -t dark`, 20,000 iterations, 31 samples, default operand sizes, selected kernel=fused. All figures are median ns:

| Bits | prealloc_portable | selected fused | GMP reused | new_selected | GMP new | explicit ARM64 carry |
|---:|---:|---:|---:|---:|---:|---:|
| 53 | 18.656 | 13.277 | 5.746 | 11.344 | 26.040 | 9.929 |
| 128 | 14.598 | 12.008 | 9.208 | 11.346 | 24.629 | 11.173 |
| 256 | 26.994 | 16.585 | 13.146 | 32.823 | 26.579 | 16.617 |
| 1024 | 350.810 | 182.898 | 95.062 | 202.892 | 109.094 | 154.744 |
| 2048 | 1383.787 | 806.990 | 320.725 | 844.354 | 335.694 | 571.048 |

Same-invocation fused-selected / explicit carry speedup = 1.337x/1.075x/0.998x/1.182x/1.413x for 53/128/256/1024/2048 bits. Larger balanced 1024- and 2048-bit cases show 15.39% and 29.24% lower time, while 256-bit times are effectively identical. Versus reused GMP, the carry kernel is 1.728x/1.213x/1.264x/1.628x/1.780x slower across these cases, whereas fused-selected is 2.310x/1.304x/1.262x/1.924x/2.516x. The new native assembly offers actual user measured benefit at 1024/2048 bits. The selected fused kernel itself remains near previous 180.769/805.240 ns measurements at 1024/2048 bits, giving a coarse (not conclusive) cross-run indication of no consolidation regression.

**Decision:** keep the current ARM64 production fused selector unchanged until the user runs macOS ARM64 `ctest --output-on-failure` on the newly consolidated and carry-enabled source tree. The 6000-case GMP differential tests exercise the explicit ARM64 carry path only when run on ARM64. Unix CI alone cannot certify native ARM64 inline assembly correctness. Rerun representative timings on the same machine, preferably reversed bit-size ordering and including 4096 bits, before selecting the assembly by default. The benchmark uses fixed operands and fixed per-backend measurement ordering, so it is not a universal performance conclusion and notably tiny-operand owning results can be affected by optimization and run-to-run noise. Crucially the current `new_selected` still runs the fused path: a potential future promotion must measure the *owning* public `mul` path against `new_mpz`. Do not add a CMake experimentation flag, and do not create extra nonarchitecture source files. The two architecture files remain the only exceptional compilation units; the fused and Karatsuba algorithms are kept in the shared file.


### B10.1 second ARM64 carry-chain benchmark, including 4096 bits (18 October 2026)

The user ran `benchmarks/benchmark_big_unsigned_multiply --include-arm64-carry --bits 256,1024,2048,4096 -t dark` from `build-release` on Apple Silicon. Selected remained `fused` (not the carry candidate); 20,000 iterations, 31 samples. Same-run median ns:

| Bits | prealloc_portable | selected fused | GMP mpz reused | new_selected fused | ARM64 carry |
|---:|---:|---:|---:|---:|---:|
| 256 | 43.019 | 22.617 | 16.048 | 36.623 | 17.831 |
| 1024 | 348.498 | 179.512 | 94.696 | 203.654 | 153.998 |
| 2048 | 1371.981 | 796.533 | 311.896 | 829.667 | 570.777 |
| 4096 | 5354.337 | 3673.671 | 986.669 | 3738.015 | 2263.667 |

Carry vs selected fused same-run speedups: 1.268x (256), 1.166x (1024), 1.395x (2048), 1.623x (4096). Carry execution-time savings: 21.2%, 14.2%, 28.3%, 38.4%. Compared with reused GMP, carry/MPZ slowdown is 1.111x/1.626x/1.830x/2.294x; selected fused/MPZ slowdown 1.409x/1.896x/2.554x/3.723x. Across this four-size subset the user benchmark reported selected/MPZ geomean 2.245x, carry/MPZ geomean 1.660x. These aggregate ratios are **not** directly comparable to earlier five-size aggregates.

**Repeatability and qualification:** B10.1's first Apple Silicon benchmark showed carry 16.617/154.744/571.048 ns at 256/1024/2048 bits versus selected fused 16.585/182.898/806.990 ns. The second run reproduces the strong carry benefits at 1024 and 2048 bits, adds a large 4096-bit benefit, but shows variable relative behavior at 256 bits (first run tied, second run faster). One interpretation is timing variability at short durations; do not infer a reliable 256-bit threshold or claim universal speedup. Because carry remains opt-in, user `new_selected` still measures fused and does not yet represent owning carry performance. The user has **not yet reported** native macOS ARM64 `ctest` results for the consolidated B10.1 code, so do not change production dispatch at this stage; retain current fused ARM64 selected kernel. GitHub Actions for the B10.1 source-layout refactor `5d63d7e` passed Unix/Windows/Coverage; those jobs are not proof of ARM64 inline assembly execution. Any carry promotion requires native GMP oracle tests green on Mac and subsequent selected owning A/B measurements. Architecture-only extra translation-unit rule remains in force.


### B10.2: promote ARM64 carry-chain to production (18 October 2026)

**Acceptance evidence:** The user reported `ctest --output-on-failure` **green** on the latest consolidated Apple Silicon B10.1 build; those native tests invoke the full randomized GMP differential oracle for forced ARM64 carry, including full/truncated results, all-one/zero/asymmetric operands and oversized buffers. Two user-run integrated macOS benchmark invocations measured real carry-chain wins over selected fused at 1024 and 2048 bits; the second also measured 4096 bits: selected fused 3673.671 ns, carry 2263.667 ns (1.62x faster); at 2048 bits fused 796.533 ns, carry 570.777 ns (1.40x faster). At 256 bits results varied across runs and no size-specific threshold is asserted. This is measured integrated evidence, not an isolated prototype prediction.

**Change:** On all supported ARM64 GCC/Clang/AppleClang native builds, `BigUnsigned::mul_low_into` and therefore `mul_into` and owning `mul` select the existing explicit ARMv8-A `MUL/UMULH/ADDS/ADC` carry-chain implementation unconditionally. No CMake options or extra translation units are created. The shared portable `mul_into_portable` and common `mul_into_fused` remain available for comparisons. `selected_multiply_uses_fused()` is now **false on ARM64**, correctly indicating that selected is no longer the shared fused kernel; benchmark output reads `selected kernel=ARM64 carry` based on native architecture availability. The ARM64 regression test asserts selection metadata and the full existing GMP comparisons exercise the selected carry implementation. The GCC x86-64 fused selection, Clang/MSVC x86-64 B3 fallback, and other architecture portable fallback are unchanged.

**Remaining validation:** GitHub Actions for this *new selector commit* must complete; reported macOS `ctest` passed on the previous implementation where carry was explicitly forced, so the selected-path switch itself still needs a fresh local `ctest` and owning `new_selected` benchmark against `new_mpz` after the new commit is rebuilt. Performance is not assumed to carry over mechanically from preallocated forced calls; compare to same-host historical baselines. Full end-to-end FloatAP performance remains unmeasured.


### B10.2 Apple Silicon production carry dispatch benchmark (18 October 2026)

**Source.** The user compiled the B10.2 ARM64 selected-carry build and ran `benchmark_big_unsigned_multiply --include-arm64-carry --include-fused --include-new-mpz -t dark` (20,000 iterations, 31 samples, original five sizes). The banner confirmed `native available=true | selected kernel=ARM64 carry`. User-reported median times in nanoseconds:

| Bits | portable | selected ARM64 carry | GMP reused | new_selected (owning) | common fused forced | new_mpz owning | ARM64 carry forced |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 53 | 18.619 | 13.735 | 7.156 | 11.325 | 10.842 | 22.762 | 9.423 |
| 128 | 13.971 | 12.144 | 9.467 | 11.037 | 11.029 | 23.485 | 10.846 |
| 256 | 26.979 | 16.617 | 13.146 | 33.354 | 16.594 | 26.758 | 16.621 |
| 1024 | 350.731 | 154.312 | 94.854 | 172.594 | 183.508 | 109.221 | 154.731 |
| 2048 | 1376.983 | 574.690 | 313.325 | 605.933 | 802.894 | 336.690 | 571.248 |

**Owning speedups verified against the previous fused-selected user run (separate invocations):** at 1024 bits, 202.990 -> 172.594 ns (1.176x speedup, ~15.0% less time); at 2048 bits, 836.858 -> 605.933 ns (1.381x, ~27.6% less time). The tiny widths 53/128/256 bits have no robust owning performance trend. With equal result lifecycle, `new_selected / new_mpz` is 0.498x, 0.470x, 1.247x, 1.580x and 1.800x at the five sizes, respectively; reported geometric mean 0.963x, but sum-of-times ratio 1.608x (large costs matter). Reused-output selected / reused GMP geomean 1.562x; forced ARM64 carry / reused GMP geomean 1.415x. At 256, 1024 and 2048 bits, selected and forced carry are nearly identical, verifying the new dispatch. At 53 and 128 bits, selected appears slower than forced carry by about 46% and 12%; because their arithmetic implementation is the same this is a call-path, optimizer, benchmark-order or timing issue to investigate, **not** evidence of two different algorithms. Forced fused remains available as a performance control.

**Validation state:** The user earlier reported macOS ARM64 `ctest` green with the carry kernel explicitly forced, prior to the B10.2 production selector change. No separate local macOS CTest result has yet been reported for this post-promotion build; benchmark's correctness oracle is not a replacement for full CTest. Unix/Windows/Coverage CI statuses should be inspected for this new commit. Both design journals are updated for the measured owning promotion. No extra translation units or performance CMake toggles. Does not demonstrate an end-to-end FloatAP speedup.


### B10.2 validation closed: selected ARM64 carry, post-promotion (18 October 2026)

**User confirmation:** `ctest --output-on-failure` on the updated macOS ARM64 build was reported **fully green after** the B10.2 production dispatch change. This supersedes the earlier pending native-CTest qualification in the B10.2 benchmark record. Native tests now exercise `mul_low_into`/public `mul` with carry selected as well as the forced carry path, including GMP-differential checks; the previously documented owning and forced-vs-selected Apple Silicon measurements remain the performance evidence. **GitHub Actions on functional dispatch commit `d01c505`: Unix, Windows, Coverage all completed successfully.** The subsequent journal-only `4bf1886` CI jobs were still running when checked and should not be confused with the completed functional-gate CI. The native ARM64 carry kernel is now the accepted default. Shared fused and portable kernels remain explicit controls. No new CMake flags or nonarchitecture source files are added; no FloatAP end-to-end speed claim.


### B11 experimental baseline x86-64 carry-chain assembly candidate (18 October 2026)

In the existing architecture-only `src/big_unsigned_x86-64.cpp`, B11 adds `mul_low_into_carry_x86_64`, an **opt-in** row-schoolbook with explicit GNU/Clang x86-64 `MULQ`, `ADDQ`, `ADCQ`: each product computes the exact 128-bit total `a[i]*b[j]+output[i+j]+carry`, using RDX:RAX for the product, then two ADD/ADC pairs to propagate the two carries into the high word. GCC/Clang extended asm ties the low input/output to `+&a` (RAX), the high output to `=&d` (RDX), keeps other operands in nonoverlapping registers, and clobbers condition codes. No BMI2, ADX, CPU feature probes, scratch allocations or new translation units. MSVC retains the existing B3 native schoolbook fallback; non-x86 architectures retain portable fallback. This candidate does not modify any production dispatch, including ARM64 carry and GCC x86-64 fused.

**Correctness approach:** Locally on Linux x86-64, the exact isolated carry instruction sequence was compiled with GCC 14.2 and Clang and matched unsigned `__int128` for 1,500,000 random/edge tuples, including all-ones. This isolated check is **not** a substitute for project tests. The 6,000-case GMP-differential oracle is extended for both full and low truncated x86-64 carry results, output overwrite, asymmetric sizes, and an undersized-full-result exception, with portable fallback on unsupported architectures. The integrated benchmark adds `--include-x86-64-carry`, which prints optional `prealloc_x86_64_carry` after existing modes; it is accepted only on GNU/Clang x86-64 builds, preserving frozen defaults and measurement ordering.

**Performance gate:** No integrated x86-64 measurements of this kernel are yet available and no claim of speedup is made. After a green full project test suite on Ubuntu x86-64, compare `prealloc_selected` (GCC fused), `prealloc_x86_64_carry`, and `mpz` within the same benchmark invocation. Suggested one-line command from build directory: `benchmarks/benchmark_big_unsigned_multiply --include-x86-64-carry --include-new-mpz --bits 53,128,256,1024,2048,4096 -t dark`. The x86-64 selected production path is unchanged pending repeated same-host evidence, including public owning results if ever promoted. Maintain architecture-only additional .cpp file policy.


### B11.1 separate portable fallback translation unit (18 October 2026)

Per user-requested three-backend source organization, the **portable C++20 limb multiplication** (`multiply_wide` 32-bit decomposition and `multiply_limbs` row-schoolbook) was moved **without changing the arithmetic** from `src/big_unsigned.cpp` to `src/big_unsigned_fallback.cpp`. This new file contains the public `BigUnsigned::mul_low_into_portable` definition and is built **unconditionally** on every platform, retaining the independently available portable GMP-differential oracle and `prealloc_portable` A/B benchmark even on native architectures. `src/big_unsigned_x86-64.cpp` still holds all x86-64-native B3/Comba/experimental `MULQ/ADCQ` kernels and `src/big_unsigned_arm64.cpp` the selected ARM64 `MUL/UMULH/ADDS/ADC` carry kernel. These are the three backend-specific .cpp files. `src/big_unsigned.cpp` stays shared for representation, additions/subtraction, full product size checks, public dispatch, owning helpers, Karatsuba, and the one shared GCC/Clang `__int128` fused implementation (which must **not** be duplicated between x86 and ARM as noted in B10.1). Its no-`__int128` path now calls `mul_low_into_portable` in the fallback TU. CMake always adds `big_unsigned_fallback.cpp`, and uses one `if()/elseif()` branch to include **at most one** native backend and corresponding macro for supported architecture/compiler. Unsupported targets use only the fallback. Neither production selectors nor benchmark/test APIs nor data allocation/lifetime change. Existing x86 GCC fused switch continues unchanged; no new flags.

This is a **structural refactor**, not a measured performance improvement. Compilation/linking and CI on supported toolchains, independent full and truncated GMP tests, and same-host benchmark confirmation are required; moving the portable implementation to another TU can change optimization and invocation overhead. Do not create any algorithm-specific extra TU, including Karatsuba: only `big_unsigned.cpp`, `big_unsigned_fallback.cpp`, `big_unsigned_x86-64.cpp`, and `big_unsigned_arm64.cpp` are expected for BigUnsigned implementations.


### B11 Ubuntu x86-64 integrated benchmark: MULQ/ADCQ loses to GCC fused (19 October 2026)

The user ran the existing integrated BigUnsigned multiplication benchmark with `--include-x86-64-carry --include-new-mpz --bits 53,128,256,1024,2048,4096 -t dark` on their native x86-64 build, reporting `selected kernel=fused`, 20,000 iterations, 31 samples. The precise same-invocation medians (ns) are:

| Bits | prealloc portable | prealloc GCC fused selected | GMP prealloc | owning selected | owning GMP | opt-in x86-64 MULQ/ADCQ |
|---:|---:|---:|---:|---:|---:|---:|
| 53 | 16.156 | 12.893 | 7.595 | 27.487 | 28.359 | 14.277 |
| 128 | 25.667 | 17.124 | 9.723 | 32.173 | 29.768 | 17.011 |
| 256 | 46.884 | 22.368 | 17.154 | 49.970 | 30.293 | 23.842 |
| 1024 | 593.858 | 205.685 | 135.967 | 234.659 | 148.040 | 221.831 |
| 2048 | 2302.470 | 761.449 | 439.639 | 786.487 | 453.858 | 939.067 |
| 4096 | 9022.733 | 2975.414 | 1389.409 | 3105.734 | 1389.616 | 3156.871 |

The explicit x86 carry kernel time relative to GCC fused selected is +10.7%, -0.7%, +6.6%, +7.8%, +23.3%, +6.1% across these widths (positive means **slower**). 128 bits is effectively tied and no general speedup exists. Reported geometric-mean slowdown against GMP reused-output is 1.672x for selected fused, 1.819x for x86 carry, and 3.628x portable; selected owning/new GMP geometric mean is 1.482x and total-time ratio 2.037x. At 4096 bits, GMP reused is 1389.409 ns against GCC fused selected 2975.414 ns (~2.14x). This suggests that higher-level multiplication algorithms/crossover (already-experimental Karatsuba) may be a more promising direction at large sizes than forcing a baseline MULQ/ADCQ instruction sequence, but requires new integrated data.

**Decision:** B11 explicit x86-64 `MULQ/ADDQ/ADCQ` stays **opt-in only**; no promotion to GCC x86-64 production dispatch. ARM64 production carry remains selected and unchanged. Possible causes include GNU/Clang inline-asm register restrictions and compiler optimization/scheduling differences; no causal claim has been proved from instruction-level profiling. CI for B11 kernel commit `7c3f9a7` and subsequent B11.1 fallback-refactor `d86c877` is fully green (Unix, Windows, Coverage each). The user supplied performance numbers, but did not explicitly report a separate local full `ctest` outcome for this new x86-64 build; do not infer full local test status from the benchmark alone. Keep the recorded B11 opt-in integrated GMP check and test oracle independent of the performance decision. No new source files, ISA flags, CMake switches, or production dispatch changes are justified by these measurements.


### B12 experimental Karatsuba write-volume optimization (19 October 2026)

Existing B8.1 Ubuntu GCC x86-64 integrated cutoff-sweep data show cutoff=32 was the most promising measured candidate (Karatsuba ~2629 ns at 4096 bits and ~8462 ns at 8192 bits), with an estimated balanced-operand crossover between 2048 and 4096 bits; the observation is historical, not an independently reproduced result for this change. This B12 experiment modifies **only** the experimental Karatsuba implementation in the shared `src/big_unsigned.cpp`: (1) `add_shifted` adds the actual cross-term length and propagates any remaining carry only until it clears, rather than scanning every high destination word; it still validates all discarded high cross-term words, preserving overflow checks; (2) recursive internal nodes no longer zero the output before two child products, since `should_split` implies both low halves have exactly `half` limbs and `z0` and `z2` entirely cover disjoint portions of the exact output; (3) the public full-product Karatsuba clears only extra caller destination limbs, leaving all exact product limbs to the recursive writer. Base cases still use the existing architecture-selected low multiplication. The algorithm, recursion splits, minimum cutoff 8, workspace sizing, signatures, and production dispatch do **not** change.

The existing 700 GMP-differential wide cases with cutoffs 8/16/32/64 are retained. Additional deterministic all-ones/randomized limb patterns test odd and asymmetric sizes including 17x17, 31x33, 63x67, 79x80, 127x129, 255x257, and the non-splitting 1x180 and 8x8 cases; output is prefilled nonzero and oversized, and scratch is bounded by sentinel limbs on either side. This is correctness coverage, not artificial coverage metric inflation. Recheck source build and Unix/Windows/Coverage CI and the full test suite; project CTest has not yet been run on a user's host for this experiment. Neither correctness nor performance improvement is asserted before validations. No new `.cpp` files, CMake flags, or other architecture changes. For performance assessment, use integrated opt-in `--include-karatsuba --karatsuba-cutoff 32 --bits 2048,2560,3072,3584,4096,8192 -t dark` alongside selected fused and GMP, with repeated run ordering to determine crossover. Do not select Karatsuba automatically from unverified new timings.


### B12 Ubuntu GCC x86-64 integrated Karatsuba cutoff=32 benchmark (19 October 2026)

The user ran the post-B12 `benchmark_big_unsigned_multiply --include-karatsuba --karatsuba-cutoff 32 --bits 2048,2560,3072,3584,4096,8192 -t dark` in their native Ubuntu x86-64 release build. Selected production kernel remains fused GCC. The integrated protocol was 20,000 iterations x 31 samples and GMP oracle validation. Median ns for the **same invocation**:

| Operand bits | portable prealloc | fused selected prealloc | GMP reused output | fused owning (`new_selected`) | Karatsuba cutoff 32 prealloc |
|---:|---:|---:|---:|---:|---:|
| 2048 | 2305.408 | 758.911 | 441.952 | 788.671 | 764.395 |
| 2560 | 3718.500 | 1387.758 | 670.152 | 1417.699 | 1141.891 |
| 3072 | 5095.817 | 1686.657 | 955.902 | 1780.820 | 1519.572 |
| 3584 | 6915.939 | 2279.586 | 1091.484 | 2318.154 | 2003.538 |
| 4096 | 9027.228 | 3023.645 | 1404.778 | 2999.035 | 2582.760 |
| 8192 | 36197.528 | 12387.829 | 4456.663 | 12489.829 | 8264.658 |

Relative to same-run fused preallocated, Karatsuba time is +0.72%, -17.72%, -9.91%, -12.11%, -14.58%, and -33.29% respectively (negative is lower/better). Observed crossover for these **balanced** operands on this one host lies **between 2048 and 2560 bits**, but is not established for asymmetric operands, other targets, or end-to-end owning `mul`. Aggregate slowdown geomean vs reused-output GMP is 2.069x for selected fused versus 1.756x for Karatsuba; cost-weighted ratio 2.386x vs 1.804x, and at 8192 bits Karatsuba remains 1.85x slower than GMP. B12's B8.1 historical cutoff=32 measurements for 2048/4096/8192 bits were 815.782/2629.425/8461.684 ns (different invocations). Post-B12 values 764.395/2582.760/8264.658 are somewhat lower, but the differences are small enough compared with cross-run noise that **B12's own optimization effect is not proved** without controlled repeated A/B measurements.

**Validation:** B12 functional commit `014eedd` completed GitHub Unix, Windows and Coverage CI successfully; its new and existing Karatsuba GMP-differential tests therefore pass on CI. User has not separately reported native `ctest` results for this commit. The current opt-in `prealloc_karatsuba` benchmark preallocates reusable output **and scratch** outside timing and precomputes scratch sizing, while public owning `mul` allocates its own result and does not currently select Karatsuba. **Do not compare preallocated Karatsuba directly with owning `new_selected` or promote it yet.** Before production promotion, measure an exact owning Karatsuba path including realistic scratch sizing/acquisition and result ownership against `new_selected`/`new_mpz`, test asymmetric operands, and repeat around the crossover. Preserve ARM64 carry and x86 fused production dispatch and existing source structure; no additional algorithm-specific .cpp files.


### B13 BigUnsigned multi-operation integrated benchmark (19 October 2026)

Per the agreed sequencing, stop optimizing multiplication and measure already-implemented unsigned operations before designing division. Added `benchmarks/benchmark_big_unsigned_arithmetic.cpp` to the **benchmark-only** executable list in `benchmarks/CMakeLists.txt`, leaving production sources, runtime dispatch, `BigUnsigned` APIs, CMake architecture selection, and all old multiplication benchmark cases unchanged. This is an executable benchmark, **not** an extra algorithm/backend source file. Default sizes 53/128/256/1024/2048 bits, 20,000 iterations, 31 samples; `--bits`, `--iterations`, `--samples`, `--summary-only` and Ariadne logging controls mirror the established integrated benchmark CLI. Deterministic equal-bit-size operands with nonzero high bits; each pair is ordered so unsigned subtraction cannot underflow. All values are checked exactly against GMP **outside timing** (including after timing). Comparison uses sign-normalized `mpz_cmp`.

Per operation, measure `add` and `sub` with preallocated destination buffers vs GMP's reused `mpz_t`, and separately owning `BigUnsigned` return values vs fresh `mpz_init`/operation/`mpz_clear`. Compare two existing owning shifts left/right by 1 and 65 bits with equivalent owning GMP shifts, excluding shift primitive prealloc claims because the current BigUnsigned API has no preallocated shift entry points. Compare `compare` vs GMP's `mpz_cmp` as a non-allocating read-only primitive. This produces 7 operation rows per operand size (add/sub/compare/shl_1/shr_1/shl_65/shr_65); aggregate preallocated and owning ratios separately, not mixing unmatched lifecycle semantics. Timings use a warm-up, steady_clock, a median per sample and compiler keepalive barriers, but are **not** end-to-end Integer performance because `Integer` is still GMP-based. As with the multiplication benchmark, fixed deterministic inputs and benchmark order should be interpreted cautiously; later add asymmetrical/zero/equal/carry-saturated cases for performance modeling and richer workloads as needed. Build/CI and native host measurements remain pending.

**Next implementation step only after baseline measurement:** define and test generic unsigned `divrem` in shared `BigUnsigned`; start with divisor one limb and progress to multiprecision, preserving exact quotient/remainder identity and explicit zero-divisor errors. Then `gcd` and decimal parsing/output, then new signed-magnitude `Integer` integration. Avoid guessing ratios or changing production behavior in this commit.


### B14 initial single-limb unsigned division in shared core (19 October 2026)

**Source-location rule clarified by user:** `big_unsigned.cpp` implements generic algorithms for which there is **no native ISA-specialized version**. `big_unsigned_fallback.cpp` contains the generic fallback **only for an algorithm that also has an x86-64 and/or ARM64 specialized implementation** (currently portable schoolbook multiplication). Native `big_unsigned_x86-64.cpp` and `big_unsigned_arm64.cpp` contain their own specialized kernels. The experimental single-limb division was therefore placed exclusively in the shared `src/big_unsigned.cpp` instead of `big_unsigned_fallback.cpp`, and the fallback file was **left unchanged**. This rule supersedes any prior ambiguous statement that all portable primitives belong in fallback; it also applies to prospective generic `divrem`, `gcd`, and conversion functions until genuine architecture specializations exist.

B14 adds exact arbitrary-length unsigned division by a **single nonzero 64-bit limb**: `BigUnsigned::divrem_limb_into(dividend_view,divisor,quotient_view)` returns a remainder and fully overwrites the quotient view (zeroes extra words), with no heap allocation and a disjoint-views contract. The owning free function `divrem_limb(BigUnsigned const&,Limb)` returns a normalized quotient/remainder pair. Invalid divisor 0 raises `std::domain_error`; too-short quotient raises `std::length_error` before mutating output. The implementation uses a pure C++20 radix-2 bitwise division of the 128-bit remainder/digit pair (no `__int128` or compiler-ISA dependencies), maintaining remainder < divisor each bit, including when a logical 65th bit is set. This is a portable **correctness-first baseline**, not a claimed high-speed alternative to GMP. General multilimbed-divisor division is **not** yet implemented.

Tests: independent exact cases for radix boundary division and zero, 6000 deterministic randomized GMP-differential cases (zero/saturated/powers-of-radix/64-limb dividends, divisor 1, 2, max, high-bit, random), both owning and preallocated outputs, normalization, extra-word zeroing, quotient*divisor+remainder identity, zero-divisor errors and too-small output without modification. The B13 add/sub/compare/shift measured baselines remain as supplied by user; there is no justified microoptimization of them in this commit. After CI and native project tests pass, continue with multi-limb `divrem`, not premature Integer migration.


### B14.1 Windows single-limb GMP-oracle width fix (19 October 2026)

The B14 Windows 2025 Release job compiled successfully but `test_big_unsigned`
failed repeatedly in `test_division_limb` at the oracle assertion
`rem==mpz_get_ui(r.value)`. On Windows x64 (LLP64), GMP's
`mpz_get_ui` returns a 32-bit `unsigned long`, whereas the native
remainder and the single-limb divisor are 64-bit. This was a **test
oracle truncation**, not evidence of an arithmetic kernel defect.

The test now imports the full native 64-bit remainder using `mpz_import`
and checks exact GMP equality with `mpz_cmp`. A deterministic case
also explicitly checks a nonzero remainder above 32 bits. All 6000
GMP-differential randomized cases, quotient comparisons, output guards,
normalization, identity checks, and exception cases remain intact.
No production implementation, kernel selection, API, build flags, or
benchmarks change. The isolated full-width GMP oracle operation was
compiled and checked on Linux; this is **not** a Windows or project-level
CTest result. Re-run Unix, Windows and Coverage CI and report their
actual outcomes before closing the Windows validation gate.


### B15 common multi-limb unsigned division (19 October 2026)

After the user-confirmed B14.1 Windows green CI, added exact unsigned `divrem` to the shared `src/big_unsigned.cpp`, leaving native multiplication and portable ISA fallback unchanged. Owning `divrem(a,b)` returns canonical nonnegative quotient and remainder satisfying `a=b*q+r`, `r<b`. The no-allocation `divrem_into` accepts pairwise-disjoint canonical limb inputs, output buffers and scratch; it fully overwrites outputs, zeroes extra limbs, validates sizes before mutation, and throws on a zero divisor. For normalized lengths n>=m>=2, `divrem_scratch_limbs` returns n+2*m+2 with checked size arithmetic; other cases need zero scratch.

Algorithm: normalize the multi-limb divisor to set its highest bit, estimate each base-2^64 quotient digit using the existing portable (high:low)/limb division, cap at one limb, compute the full product with the existing selected multiplication kernel, correct the estimate downward while its product exceeds the current dividend window, then subtract and unnormalize the final remainder. For single-limb divisors reuse B14; for smaller dividends return zero quotient and input remainder. No new backend source or ISA assumption is introduced. This is a correctness-first reference, not a performance-parity claim against GMP.

A standalone Linux C++20 algorithm prototype passed 90,000 randomized/adversarial GMP comparisons. Production regression now includes 6000 deterministic GMP quotient/remainder differential cases (owning and caller-buffer paths), fixed zero/small/equal/radix boundary checks, exact identity, nonzero scratch guards, zeroed output tails, zero-divisor and short-buffer failures without mutation, plus scratch-size overflow. **Project CTest and Unix/Windows/Coverage CI are pending** for this commit; the standalone test is not a substitute. Future work: gcd and decimal import/export, then signed-magnitude Integer, only after production validation.


### B15.1 shared subtract helper collision fix (19 October 2026)

B15 failed to compile on AppleClang ARM64 with `-Werror`: the new
division-local `subtract_in_place(Write,Read) noexcept` conflicted with
the existing later Karatsuba `subtract_in_place(Write,Read)`, which can
throw `std::logic_error` if the subtraction's nonnegative precondition
is violated. Both helpers were in the same unnamed namespace and had
the same argument types; the exception-specification mismatch caused
the diagnostic. This is a compile-time regression, **not** evidence
that the division algorithm returned an incorrect value.

Remove the redundant division-local implementation, declare the shared
helper before the division routine, and reuse the existing checked
implementation for both quotient-product subtraction and Karatsuba.
This preserves division's exact correction loop and its nonnegative
precondition while keeping the existing borrow guard. No changes to
API, buffers, ISA selection, CMake flags, GMP oracle, or test coverage.
The previously described standalone prototype checks do **not** verify
the integrated B15.1 fix. Re-run Unix/Windows/Coverage CI and host
CTest before claiming integrated correctness.


### B16 BigUnsigned GCD and strict decimal conversion (19 October 2026)

B15.1 Unix, Windows, and Coverage workflows passed at `feb15cc`. B16 adds generic GMP-free `gcd(BigUnsigned,BigUnsigned)`, `from_decimal(std::string_view)` and `to_decimal(BigUnsigned const&)` in the existing shared `src/big_unsigned.cpp`. GCD uses Euclidean remainder with canonical zero rules `gcd(0,0)=0`, `gcd(x,0)=x`. Parsing accepts only a nonempty ASCII unsigned decimal sequence, including leading zeroes; signs, whitespace, prefixes and nondigits raise `std::invalid_argument`. Output always uses canonical digits, with zero as `"0"`.

The decimal parser consumes up to nine digits per block and multiplies existing radix-2^64 limbs by at most 10^9 using checked 32-bit half-product arithmetic, without a required 128-bit compiler extension. Decimal output repeatedly divides a temporary mutable limb copy by 10^9 using the existing portable wide-limb divider, accumulating and zero-padding nine-digit blocks. Neither operation uses GMP in production or alters architecture dispatch. These are correctness-first baselines, not competitive performance claims.

An isolated Clang C++20 harness with warnings-as-errors passed 5000 native limb-to-decimal GMP round trips, 5000 independent decimal-string GMP comparisons (including 1900-digit cases), 1100 randomized GMP GCD comparisons and rejected malformed ASCII inputs. Added project regression cases against GMP and fixed boundary checks. **Integrated Unix/Windows/Coverage CI remains pending** on this new commit; the isolated harness is not a project-level result. This concludes the planned prerequisite primitive set for the initial signed-magnitude Integer design, but does not constitute Integer migration.


### B17 native signed-magnitude foundation, before Integer facade migration (19 October 2026)

B16 CI `3ea38f4` is green on Unix, Windows and Coverage. Directly
replacing `Integer::_mpz` is not yet isolated: live consumers in
`rational.cpp`, `dyadic.cpp`, and `floatmp.cpp` use exposed GMP
members/accessors, while the existing Integer header and constructors
expose GMP types. A one-commit replacement could break unrelated
numeric semantics and leave test failures hard to localize.

B17 therefore introduces a separate `SignedInteger` native
signed-magnitude class in `include/numeric/signed_integer.hpp` and
`src/signed_integer.cpp`. It contains a canonical `BigUnsigned`
magnitude and a sign flag, normalizing negative zero. Supported
operations: signed/unsigned 64-bit construction (including INT64_MIN),
strict decimal parse/output, sign, comparison, negate, add/sub/mul,
square, exponentiation by squaring, unsigned-magnitude signed divrem,
shifts and checked 64-bit extraction. Division is **truncation toward
zero**, with remainder having the dividend's sign: `a=q*b+r` and
`abs(r)<abs(b)`, including all sign combinations. Right shift
truncates toward zero; `pow(0,0)=1`. These are explicit new
contracts, not accidental preservation of the old mismatched
`mpz_div` quotient and `mpz_cdiv_r` remainder.
No implicit GMP storage, arithmetic, or conversions are used.

A separate GMP-differential test target checks independent sign
boundaries, exact signed identities and random wide operands, including
the 64-bit minimum, canonical zero, and invalid inputs. Existing
`Integer`, `Natural`, `Rational`, `FloatAP` and legacy APIs
are left intact in B17; migrating consumers and replacing the public
facade are **remaining work**, not completed by introducing this core.
Full B17 cross-platform project CI remains pending.


### B17.1 decision: SignedInteger is temporary verification scaffolding (19 October 2026)

The user clarified that B17's `SignedInteger` **must remain temporarily** to independently validate sign-magnitude behavior against GMP before the in-place public `Integer` migration. This corrects the previous inclination to immediately revert B17; **no revert was applied**. There must be **only one signed arbitrary-precision type in the final architecture: `Integer`**.

Acceptance sequence: (1) validate the standalone B17 `SignedInteger` implementation and its GMP-differential tests across Unix, Windows and Coverage; (2) move its verified signed-magnitude invariants and operations into existing `include/numeric/integer.hpp` and `src/integer.cpp`, using `BigUnsigned` directly, without a permanent wrapper; (3) migrate dependent `Rational`, `Dyadic`, `FloatMP` and other consumers of `_mpz`/`get_mpz` using explicit conversions or narrow integration adapters, retaining existing meaningful semantics and the agreed signed division convention; (4) transfer and expand verification into existing `tests/test_integer.cpp`; (5) **delete** `include/numeric/signed_integer.hpp`, `src/signed_integer.cpp`, `tests/test_signed_integer.cpp` and their CMake target/source entries. Check that no code or build reference to `SignedInteger` remains. This deletion is a required completion gate for the migration, not optional cleanup.

B17 remains an intermediate development and test stage; it does not establish a second production API. This documentation-only clarification does not alter B17 code or its CI results.


### B18 direct native Integer storage and arithmetic (19 October 2026)

B17 `SignedInteger` Unix, Windows and Coverage CI succeeded at `31c81e1`; the temporary signed verification class and its existing tests remain in the tree for side-by-side validation. **B18 migrates the actual `Integer` class in-place**: `include/numeric/integer.hpp` now owns a canonical `BigUnsigned` magnitude and a sign flag, without a persistent GMP `_mpz` or `get_mpz()`; `src/integer.cpp` uses the shared BigUnsigned arithmetic for signed integer construction, comparison, add/sub/mul, square, unsigned-magnitude quotient/remainder, powers, decimal IO, increment/decrement, absolute value, bit-length and checked casts. Move operations reset the moved-from sign to avoid negative zero. The legacy GMP constructor is removed; the previous direct-access tests use the external adapter.

Signed quotient/remainder are deliberately coherent: division truncates toward zero, the remainder has the dividend's sign, and `a=b*q+r` with `abs(r)<abs(b)`; this replaces the previous mixed `mpz_div` quotient and `mpz_cdiv_r` remainder behavior. The existing Integer-vs-GMP benchmark retains its corpus, timing controls and exact operations, but its direct GMP quotient/remainder oracle now uses `mpz_tdiv_q` and `mpz_tdiv_r` to compare mathematically identical semantics. Checked extraction is explicit on overflow; it no longer silently truncates via `mpz_get_si`.

An **internal, temporary boundary adapter** `include/numeric/detail/integer_gmp_bridge.hpp` converts signed `Integer` values to/from GMP with exact `mpz_export`/`mpz_import` on 64-bit limbs and canonical sign. The few existing interface sites in `rational.cpp`, `dyadic.cpp`, `floatmp.cpp` and `floatap.cpp` that previously accessed `Integer::_mpz`/`get_mpz` now use the adapter; these classes **continue to store GMP** and retain their arithmetic, except that rational round/ceil/floor use native Integer operations after extraction rather than abusing Integer as mutable GMP scratch. No broad rewrite of these classes or FloatAP algorithm is intended. The bridge is not Integer storage and is excluded from native Integer arithmetic hot paths.

Extend `test_integer.cpp` with independent signed boundary assertions, all sign combinations for truncating quotient/remainder, malformed decimal rejection, 64-bit min/max, checked conversions, canonical zero, powers and 1700 GMP-differential randomized cases over signed wide decimal operands. Keep all existing `test_integer` coverage and B17 standalone `test_signed_integer` tests. **Cross-platform CI / project CTest and native benchmark performance evidence are pending; no performance claim.** After successful B18 migration regression and benchmarking, remove the temporary `SignedInteger` source/header/CTest target and transfer any unique validation into `test_integer`. The internal GMP adapter may remain temporarily while Rational/Dyadic/FloatMP/FloatAP themselves still use GMP; it must not reinstate GMP-dependent fields in Integer.


### B19 decimal format speed and representative division benchmark (19 October 2026)

Native B18 benchmark measurements on macOS ARM64 and x86-64 (20,000 iterations, 31 samples, 53/128/256/1024/2048/4096 bits) verified 384 paired input values and 2304 paired operation results. For 4096 bits, `to_decimal` was 115.476x GMP on ARM64 and 209.944x on x86-64, dominating total-time ratios. The benchmark's legacy quotient inputs always had `abs(a)<abs(b)`, hence only zero quotient; retain them for historical comparability, but do not misrepresent them as general division performance. The comparison wrapper also performed two native less-than tests versus one GMP compare. These were measurement issues, not demonstrated arithmetic defects.

B19 makes the **targeted common-only** decimal-format change in `src/big_unsigned.cpp`: while dividing the copied limb vector by `10^9`, split each 64-bit limb into high and low 32-bit halves and perform two ordinary 64-bit quotient/remainder steps, maintaining the invariant `0<=remainder<10^9<2^32`. Each `(remainder<<32)|half` is bounded below `2^64`; this is portable C++20 on Windows LLP64 and AppleClang/GCC/Clang without `__int128` or new ISA/backend sources. The existing arbitrary-divisor `divrem_limb` and `divrem` implementations remain untouched. The exact fixed and GMP-differential `test_gcd_decimal` suite remains intact. An isolated Clang warnings-as-errors differential check compared 800,000 high/low inputs and partial remainders against a 128-bit reference and passed; this does not establish an integrated build or performance improvement.

The Integer benchmark compares using **one** `cmp` call versus one `mpz_cmp`. Existing original quotient/remainder rows and the original 32-value pools remain unchanged, while two additional `quotient_nontrivial`/`remainder_nontrivial` rows use a separate deterministic 32-value signed division pool with single-limb, half-width, near-equal and small divisors. Both backend pools are constructed outside timing. The benchmark's exhaustive cross-backend oracle now checks the additional input pairs and two operations at each width, with a fixed expected count of 14 instead of 12; historical 12-op per-size aggregation must not be directly compared to new 14-op aggregate. No changes to Rational benchmark paths, CI settings, timed sample counts, or data exclusions. **4096 bits remains supported and untouched**; for the next provisional benchmark invocation select `--bits 53,128,256,1024,2048` explicitly for speed, then restore 4096 after measuring `to_decimal` improvements.

Integrated Unix/Windows/Coverage CI and native ARM64/x86-64 benchmark results for B19 are pending. Keep `SignedInteger` temporary and migrate its unique tests before final deletion, independently of benchmark throughput.


### B19.1 scoped type fix (19 October 2026)

AppleClang reported two unqualified `Limb` uses in free function `Ariadne::to_decimal`. The function was already in `namespace Ariadne`; `Limb` is a nested alias of `BigUnsigned`. Qualify both as `BigUnsigned::Limb`. No other source, test, benchmark, algorithm or namespace is changed. Integrated CI remains pending.


### B19.2 equal-prefix compare improvement and accurate benchmark verification count (19 October 2026)

The user supplied B19.1 native ARM64 and x86-64 benchmark runs (20,000
iterations, 31 samples). The shared `to_decimal` optimization reduced
2048-bit conversion from 111251 to 5113 ns on ARM64 (21.76x) and from
240108 to 7589 ns on x86-64 (31.64x). The x86-64 run also included
4096 bits: conversion dropped from 952785 to 27951 ns (34.09x).
These are same-user before/after runs, not simultaneous A/B trials.

The most striking remaining platform difference is `compare`:
at 2048 bits native/GMP is approximately 1.02x on ARM64 but 16.05x
on x86-64. `BigUnsigned::compare` repeatedly invoked `limb()`
for equal high prefixes, including repeated bounds and inline-vs-heap
selection. B19.2 changes only this shared implementation to acquire
both `LimbView` spans once and compare their limb values directly,
preserving length-first comparison, signed ordering, canonical rules
and `noexcept`. No ISA-specific source or kernel selection changes.
Existing GMP-differential randomized compare coverage remains intact;
added a deterministic 64-limb equal-prefix case and high-limb mismatch
case. An isolated C++20 compare-view implementation passed 250000
deterministic reference comparisons when compiled with both Clang and
GCC, with warnings-as-errors. That does not replace project CTest or
establish an observed speedup.

The new benchmark division-input oracle in B19 checked all 32
additional input pairs inside the pre-existing 32-element input loop,
incorrectly reporting 10560 distinct inputs for five widths or 12672
for six. B19.2 moves that validation loop outside the old input loop,
preserving exactly the same comparisons, with no weakening and no
change in timed measurements. The accurate counts are 640 distinct
input values for five widths or 768 for six, with unchanged 2240 or
2688 per-operation result checks respectively. This is not a
performance change.

The B19.1 Unix/Windows/Coverage CI was green at `a7b7bd1`, and its
native benchmark outputs passed numerical verification. The B19.2
full project and cross-platform CI plus new ARM64/x86-64 performance
results are **pending**. Next investigate nontrivial division using
representative quotient sizes, without changing or reducing its GMP
differential coverage. Keep 4096 available but use
`--bits 53,128,256,1024,2048` for interim repeats.


### B20 selected native 128/64 unsigned division and GMP boundary tests (19 October 2026)

The user's B19.2 native ARM64 and x86-64 six-size Integer-vs-GMP
benchmarks (53, 128, 256, 1024, 2048, 4096 bits; 20,000 iterations,
31 samples) both passed exact numerical verification: 768 distinct
input values and 2688 result values. At 4096 bits, the user measured
`Integer/GMP` for nontrivial quotient as about 15.56x on ARM64 and
13.16x on x86-64; the x86-64 `compare` regression was largely removed
by B19.2. These ratios are empirical baselines, not forecasts.

B20 changes **only** the shared 128/64 digit division selector inside
`src/big_unsigned.cpp`. The invariant for
`divide_wide(high,low,divisor)` is `divisor!=0 && high<divisor`,
so both quotient and remainder fit in 64 bits and
`0<=remainder<divisor`. If `high==0`, plain C++ 64/64 quotient and
remainder suffice. With GNU/Clang `__SIZEOF_INT128__` the nontrivial
path computes `(uint128(high)<<64)|low` and divides directly,
suppressing only the accepted compiler-extension `-Wpedantic`
diagnostic; with genuine MSVC x64 it uses the Visual Studio 2019+
`_udiv128` intrinsic from `<immintrin.h>`. Other compilers and
architectures retain the unchanged 64-iteration pure-C++20
`divide_wide_portable` fallback, annotated `[[maybe_unused]]`
to avoid spurious warnings when native operations are compiled.
There is no new backend translation unit or ISA assembly, no changes
to long-division correction / remainder algorithms, no new public
API and no functions outside `namespace Ariadne`. Both
`divrem_limb_into` and `divrem_into` call this one selected helper.
The existing `to_decimal` fast radix-`10^9` path remains unchanged.
Quotient and remainder at `Integer` level continue to share
`divrem` and the existing signed truncation convention.

Outside the project, the selected GNU/Clang implementation was
compiled with GCC and Clang at `-std=c++20 -O3 -Wall -Wextra
-Wpedantic -Werror` and compared against the original bitwise
fallback over **1,200,252** exact 128/64 cases per compiler
(1,200,000 deterministic random cases plus 252 boundary cases),
all passing. This is an isolated check, not evidence of project
CI, of MSVC intrinsic correctness or of speedups.

Extend the existing `test_big_unsigned` GMP-differential suite
without relaxing any previous test: cover eight divisor
extremes, all-high/low quotient boundaries (including
`high=divisor-1,low=UINT64_MAX`), 64-bit exact remainders
on LLP64, plus **all 64 long-division normalization shifts**
with ten saturated/near-boundary dividend patterns each.
The previous 6000+6000 random oracle cases and checked
preallocated-buffer invariants remain untouched. No benchmark
corpus, timed loops, sizes or special mode switches are changed.
For the next six-size ARM64 and x86-64 comparisons, use the
existing `benchmark_compare_integer_backends` executable with
`--bits 53,128,256,1024,2048,4096 --iterations 20000
--samples 31`. **Unix, Windows, Coverage CI and native benchmark
measurements are pending for B20**, so no new performance ratio
is claimed.


### B21 complete retirement of transient signed class (19 October 2026)

After successful B20 native ARM64 and x86-64 benchmark verification and
green Unix/Windows/Coverage CI for `c767b45`, retire the independent
B17 migration scaffold. The production `Integer` class continues to
own canonical `BigUnsigned` magnitude and sign. Delete exactly
`include/numeric/signed_integer.hpp`,
`src/signed_integer.cpp`, and `tests/test_signed_integer.cpp`;
remove its source from `src/CMakeLists.txt` and its test executable
from `tests/CMakeLists.txt`. Do not replace the scaffold with another
signed numeric wrapper. No changes to production arithmetic,
benchmarks, ISA dispatch, FloatAP implementation or GMP boundary
adapter are necessary.

Before removal, transfer the complete independent 2500-case B17
GMP-differential random corpus (original seed 0x51C0FFEE, operand
construction, sign patterns, zeros, add/sub/mul, truncating quotient
and remainder, remainder sign and algebraic identity) into a new
`TestInteger::test_signed_regression` test section. Preserve B18's
existing distinct 1700-case random differential corpus. Check the
fixed sign/negative-zero, signed extrema, parsing failures,
checked conversions, zero-power, negative square, signed 20-case
division matrix, division-by-zero, and 130-shift pattern from B17.
Since public `Integer` does not expose a direct signed shift API,
test the same right-shift truncation-to-zero semantics by combining
the native `BigUnsigned::shift_right` with `Integer::from_magnitude`,
and independently through `quot(x,2^bits)`, comparing both with
GMP `mpz_tdiv_q_2exp`. Similarly check native magnitude left shift
with restored sign against GMP `mpz_mul_2exp`. Keep exact checked
GMP oracle conversions and canonical nonnegative zero assertions.
Do not add a new permanent signed class or public API.

The deletion is an architecture completion gate, not merely a
build-time convenience. All production code continues to reside
under `namespace Ariadne`. Full Unix/Windows/Coverage and CTest
results for this migration are pending; do not assert complete
validation until green. The benchmark remains unchanged (all six
sizes, fourteen operations) and can be rerun without migration
specific modifications.


### B22 balanced large-integer decimal formatting (19 October 2026)

After completion of the transient signed-class retirement in B21,
the user provided native ARM64 and x86-64 six-size Integer/GMP
benchmarks: both passed 768-input / 2688-result verification, and
B21 Unix, Windows and Coverage CI were green at `7faa6da`.
The 4096-bit native `to_decimal` cost remained approximately
19.6 microseconds on ARM64 and 27.9 microseconds on x86-64,
about 5.3x and 6.3x GMP respectively. These are B21 baselines.

B22 optimizes only `BigUnsigned::to_decimal` in the shared
`src/big_unsigned.cpp`. Preserve B19's exact radix-`10^9`,
two-64/32-division formatter as `to_decimal_blocks` and select it
unchanged for all values of at most **16 64-bit limbs** (including
zero, signed-magnitude single limb, the inline/heap boundary,
and every 1024-bit input). For larger inputs, build a private
per-call table `P_i=10^(9*2^i)` by exact `BigUnsigned::square`,
caching the exact decimal widths `9*2^i`. Avoid constructing the
expensive final unnecessary square when its minimum bit length
already exceeds half the original input's bit length.
Recursively split a positive value at the largest `P_i` whose
bit length is at most half the current value's bit length,
using the existing exact `BigUnsigned::divrem(value,P_i)`.
Convert the high and low parts recursively, left-padding the low
decimal string with exactly `9*2^i` digits before concatenation.
Because `0 <= low < P_i`, the suffix never exceeds the fixed
width; leading high zeros are absent, and low zeros/long internal
zero runs are exact. Leaves of at most 16 limbs use
`to_decimal_blocks`. Zero is returned as `"0"`. No
floating-point digit estimate, no global or shared mutable cache,
no change to parsing/from_decimal, arithmetic, sign rules,
benchmarks, public API, target-specific dispatch or GMP boundaries.
The helper functions are in an anonymous namespace **inside**
`namespace Ariadne`, never in the global namespace.
This algorithm reuses native divrem and should be evaluated on
both architectures for speed rather than assuming a crossover.

Extend the existing `test_big_unsigned::test_gcd_decimal` with
50 exact GMP-oracle decimal cases: powers `10^(9*2^k)` for
`k=0..9`, each at its boundary, one below, one above,
nine above, and with a long zero-filled run plus a nonzero suffix.
Construct the independent native reference via `mpz_export`
and compare output **exactly** to GMP's `mpz_get_str`;
also verify `from_decimal(to_decimal(x)) == x`.
This supplements, never replaces, the pre-existing 2800
randomized GMP decimal checks, 1800 randomized parse checks
and signed Integer regression suites. A separate Boost
multiprecision implementation of the same balanced split
(not the repository's BigUnsigned class) passed **6058** GMP
differential decimal checks up to approximately 25000 bits.
Its timings are a *proxy only*, not proof of a native performance
improvement. B22 integrated cross-platform CI and the actual
ARM64/x86-64 performance benchmarks are pending.

Measure B22 with the unchanged full 53,128,256,1024,2048,4096-bit
Integer-vs-GMP benchmark (20000 iterations, 31 samples), comparing
B22 to B21 separately for each platform and especially inspecting
1024-bit stability, 2048-/4096-bit `to_decimal`, overall slowdown
and any regression in `from_decimal` and other operations.
Do not optimize or alter the benchmark corpus to manufacture
improvements. Keep 4096 bits included.


### B23 specialized single-output unsigned division (19 October 2026)

Starting with the B22 native decimal formatting baseline, split
only the *result ownership* of unsigned division while retaining the
exact normalized radix-2^64 quotient estimation and correction algorithm.
Inside `BigUnsigned`, a private compile-time-selected
`divrem_selected_into<WantQuotient,WantRemainder>` shares the same
normalization, quotient-digit correction, and working scratch. The
legacy public `divrem_into` and owning `divrem` preserve both output
buffers, checked public buffer sizes, fail-before-mutation semantics,
and normalized quotient/remainder.

New native magnitude `quot(BigUnsigned,BigUnsigned)` and
`rem(BigUnsigned,BigUnsigned)` produce only the requested result.
The quotient buffer has only `n-m+1` limbs where `n>=m`, rather
than `n`; remainder-only allocates no quotient object and skips
all quotient output stores, while quotient-only allocates no
remainder object and skips final remainder unnormalization.
For a one-limb divisor, remainder-only consumes digit divisions
without allocating or storing quotient limbs. A dividend shorter
than its divisor yields canonical zero quotient or the untouched
canonical dividend as remainder. Zero divisor still throws.
All helpers/functions remain inside `namespace Ariadne`.
`Integer::quot` and `Integer::rem` now call the corresponding
single-result unsigned primitive, preserving truncation toward zero,
dividend-signed remainder, and canonical nonnegative zero.

Extend the existing 6000-case random multiprecision GMP oracle by
verifying both independent owning single-result paths against GMP
and the original `divrem`, retaining all old tests and adding
zero/short/single-limb/large-boundary cases and both exception paths.
Benchmark corpus, public signed Integer behavior, GMP adapters,
Rational/Dyadic/FloatAP, architecture multiplication kernels and
B22 decimal formatting are unchanged. This is a targeted isolation
of output work, not a new division algorithm. Full Unix, Windows,
Coverage CI and ARM64/x86-64 performance measurements remain pending;
compare B23 against B22 at all 53..4096-bit widths with 20,000
iterations and 31 samples. Inspect all 4 quotient/remainder workloads.


### B23.1 benchmark operation naming clarification (19 October 2026)

Following the user-requested clarification, rename the four integer
division benchmark *labels only*: historical `quotient` to
`quotient_trivial`, historical `remainder` to
`remainder_trivial`, `quotient_nontrivial` to `quotient`, and
`remainder_nontrivial` to `remainder`. The original 32-pair
corpus satisfies `|a|<|b|`, hence its quotient is zero and remainder
equals the dividend; the additional 32-pair division corpus includes
small divisors, half-width divisors and near-equal magnitudes,
and produces nonzero quotients. Keep every operand, iteration, sample,
operation order, verified result, operation count and timed routine
unchanged: only report names differ. Old and new benchmark output
rows must be matched by this four-way mapping before comparisons.
B23 native results from the user were GMP verified on ARM64/x86-64
(768 pairs and 2688 results on both); its Unix/Windows/Coverage
workflow suite succeeded at commit `16609ab`. The current B23.1
naming-only commit requires fresh CI validation, but no numerical
performance gains are claimed from its label change.


### B24 wide-radix decimal construction (19 October 2026)

After B23.1's report-label-only update, GitHub Unix, Windows,
and Coverage CI all succeeded on `68bc5e9`. B23's last
ARM64 and x86-64 GMP differential benchmarks each passed
768 paired inputs and 2688 paired operation results. As shown by
the previous six-size native benchmark, large decimal construction
remains expensive relative to GMP (roughly 7.6 microseconds
ARM64 and 6.5 microseconds x86-64 at 4096 bits in B23).

B24 changes only native `BigUnsigned::from_decimal` in
`src/big_unsigned.cpp`. Preserve the complete B19 base-`10^9`
split-32-bit implementation and the shared ASCII digit validation
for all strings of **309 characters or fewer**, and for all
compiler targets that do not support `unsigned __int128`.
For longer strings on GCC/Clang targets defining
`__SIZEOF_INT128__`, accumulate decimal chunks of up to **18
digits** (radix up to `10^18`) using one exact 128-bit intermediate
per 64-bit limb: `product=(Wide)digit*radix+carry`,
`digit=(Limb)product`, `carry=(Limb)(product>>64)`.
Because `radix<=10^18<2^64` and `carry<radix`, no product
overflows 128 bits. Both the first possibly short chunk and
subsequent full chunks preserve leading zeros and all limb carries.
The fallback has no ISA or GMP dependencies and remains available
on Windows/MSVC. Helpers and all modified production functions
are inside `namespace Ariadne`. The 309-digit cutoff leaves
the complete <=1024-bit benchmark pool on the old implementation.

A separate C++20/GMP single-process probe independently
compared the existing 9-digit algorithm, new 18-digit algorithm,
and `mpz_set_str` on deterministic strings up to 7000 digits.
It passed all checked values on GCC and Clang. Its timings,
particularly ~1.8x around 617 digits and ~2.7-3.0x around
1234 digits, are isolated **proxies**, not project-integrated
ARM64/x86-64 speedups. Keep the benchmark corpora, 14 operations,
sample count and architecture test fixtures unchanged; measure
actual improvement on each host.

Extend `test_big_unsigned::test_gcd_decimal` with 45 exact
GMP-oracle parsings covering lengths 308, 309, 310, 311, 312,
616, 617, 1234 and 1900, with all-nine carry chains,
long zero runs, leading zeros, mixed digits and alternating
digits. Check directly against GMP values parsed from the
*original input* and against exact GMP decimal output.
Add four malformed >309-digit cases spanning first/middle/last
invalid character positions. Existing 2800 randomized magnitude
oracle cases, 1800 parsing oracle cases, 50 B22 split boundaries
and signed Integer tests are untouched. Unix/Windows/Coverage
CI and native ARM64/x86-64 performance results for B24 remain
pending until the commit is validated. Compare B24 vs B23
using unchanged six sizes, 20,000 iterations, 31 samples;
pay special attention to `construct_decimal` at 2048 and
4096 bits and ensure no regressions in `to_decimal` or division.


### B25 wide-radix large-leaf decimal formatting (19 October 2026)

Building on B22 balanced decimal formatting and B24's successful
wide-radix *parsing*, B25 optimizes the leaf formatter inside
`BigUnsigned::to_decimal`, without changing the B22 split-power
recursion, quotient/remainder, signed Integer behavior, or benchmarks.
For magnitudes of **9 limbs or more** on GCC/Clang targets with
`__SIZEOF_INT128__`, use an unsigned-128-bit 128/64 division by
the constant `10^18` to process one 64-bit input limb per decimal
block pass. A 64-bit running remainder is always below `10^18`
and the combined `(Wide(remainder)<<64) | limb` fits in 128 bits;
each quotient limb is below `2^64`. Collect 18-digit base-`10^18`
remainders and pad all nonleading blocks to 18 digits.
This should reduce scans relative to B19's 9-digit method.
The original two-64/32-division base-`10^9` algorithm is
unchanged for **up to 8 limbs** and on compilers without
native unsigned 128-bit integers, including MSVC.
The B22 recursive cutoff remains 16 limbs; large values
therefore gain the wide leaf inside the existing split recursion.
No public API changes, external dependencies, global caches,
floating-point estimates, or new architecture-specific assembly.

The implementation's fast-path selection occurs in the existing
private `to_decimal_blocks` helper; `to_decimal_blocks_wide`
is an internal helper in the anonymous namespace nested inside
`namespace Ariadne`. Validate with 50 additional GMP-oracle
cases at limb sizes 7,8,9,10,15,16,17,32,64,128 including
binary radix powers, all-ones values, zero gaps, and random
limb patterns. Each case compares the complete decimal string
directly against GMP `mpz_get_str` and roundtrips through
`from_decimal`. Preserve B22's 50 exact decimal-boundary tests,
B24's 45 long-input parse checks, existing GMP random corpora,
and the signed Integer suites.

A *separate* C++20 + GMP local probe ran old 9-digit and new
18-digit leaf loops on 32 exact independent random operands
for 512 through 8192 bits. Both returned identical decimal
strings to `mpz_get_str` on all checked inputs. Its x86
per-value median of nine local repetitions was roughly
2.74us (9-digit) versus 1.95us (18-digit) at 1024 bits,
10.37us versus 5.23us at 2048 bits, and 47.26us
versus 18.32us at 4096 bits; these are **isolated leaf-only
proxies**, not timings of the repository's optimized B22
recursive path and not representative of ARM64. They justify
trying B25 but do not establish actual acceleration.
The output fast path differs from B24's parsing fast path.
Retain the exact six-size Integer/GMP benchmark with 20,000
iterations and 31 samples. Compare B25 with B24 separately
on native ARM64 and x86-64, particularly `to_decimal` at
1024/2048/4096 bits, other operations, and aggregation.
Run Unix/Windows/Coverage CI; all B25 integrated validation
and production speedups are pending until verified.


### B25.1 ARM64-only wide-radix output, no CMake flags (19 October 2026)

B25's full native benchmarks verified 768 independent paired
inputs and 2688 paired operation results on both ARM64 and
x86-64. However, B25's base-10^18 decimal output on x86-64
regressed strongly relative to B24: native `to_decimal` rose
from approximately 2239 to 4212 ns at 1024 bits, 6725 to
9126 ns at 2048 bits, and 19731 to 23495 ns at 4096 bits.
On ARM64, B25 improved the same three sizes from roughly
1590 to 1346 ns, 4556 to 4186 ns and 13138 to 12596 ns.
These measurements justify an architecture-specific
compile-time selection; do not use the isolated B25 leaf
microbenchmark as a production speedup claim.

B25.1 retains the exact B25 base-10^18 unsigned 128-bit leaf
implementation **only when compiler target macro
`__aarch64__` is predefined** and the compiler defines
`__SIZEOF_INT128__` and is GCC/Clang. It continues to use
the B19 base-10^9 leaf for ARM64 magnitudes of <=8 limbs.
For **x86-64 regardless of `__int128` support**, MSVC,
and other targets, the existing portable radix-10^9 path is
always selected. Both compile-time guards in the same
`src/big_unsigned.cpp` implementation are updated.
No CMakeLists.txt file is modified, no user-visible CMake
option, compile definition, new feature flag or build-setting
dependency is introduced: selection uses the compiler's
standard target macro directly. The B22 recursive output
algorithm, B24 parsing, B23 division, benchmark operands and
operation names, test fixtures, and GMP bridges are unchanged.

Keep all B25 50 GMP-differential limb-size cases and add
18 exact power-of-10^(18*k) boundary checks for
k=8,9,10,16,17,18 (exact power, minus one, plus one)
to exercise wide-radix carries and low-block padding on
ARM64 and the portable fallback on x86-64/Windows.
Do not remove any existing regression tests.

B25.1 Unix/Windows/Coverage CI results and native host
benchmarks are pending. Expect the B25 to_decimal results
to be retained on ARM64 and the B24 results to be restored
on x86-64; these are expectations, not assertions of
measured performance. Compare all six established bit
sizes with the same 20,000 iterations, 31 samples and
unchanged operation pool. Do not modify CMake for
architecture selection.


### B26 x86-64 baseline DIVQ for exact 128/64 division (19 October 2026)

After B25.1, the native x86-64 Integer/GMP slowdowns at 4096 bits
were approximately 5.48x for the general `quotient` benchmark
(1982 vs 362 ns) and 5.35x for `remainder` (1947 vs 364 ns).
ARM64 at the same size was approximately 2.93x and 2.88x.
In contrast, `pow3` was approximately 2.64x on x86-64
versus 2.68x on ARM64, so B26 targets division rather than
modifying multiplication or exponentiation.

B26 changes *only* the internal `divide_wide` helper in
`src/big_unsigned.cpp` for baseline GNU/Clang x86-64
(`__x86_64__`, not `_MSC_VER`). Under the already-established
precondition `divisor!=0 && high<divisor`, place the two-limb
dividend in `RDX:RAX` and execute a single `divq` with an
explicit register input for the divisor. The quotient returns
in RAX and remainder in RDX. Use `+a` and `+d` read/write
constraints, plus the `cc` clobber; there is no memory
access inside the assembly. The high=0 shortcut remains
the normal 64-bit division. The precondition ensures that
DIVQ cannot overflow its 64-bit quotient. This baseline
instruction requires no BMI2, ADX or runtime CPU feature
detection. MSVC x64 retains its existing `_udiv128` intrinsic;
non-x86 platforms retain their `unsigned __int128` or
bitwise portable fallback. The kernel's public behavior,
error checks, memory allocation, long-division digit
correction, and GMP bridge remain unchanged.

The change affects the expensive 128/64 quotient estimates
and single-limb division, not trivial division cases where
`|a|<|b|` and the helper is not invoked. It does not
change `Integer::pow`, the multiplication dispatch, decimal
input/output, benchmarks, build flags, or CMake in any way.
No new CMake option, preprocessor build definition, instruction
set requirement or ABI dependency is introduced.

A separate local C++20 GNU inline-assembly checker compiled
with GCC and Clang and passed **500072** exact constrained
128/64 divisions against the arithmetic `unsigned __int128`
implementation, verifying approximately one in 256 against
an independent GMP `mpz_tdiv_qr` oracle. The GCC generated
instruction stream for the isolated helper contains `div`
rather than the wide-division support call. These are
local isolated *correctness and code-generation checks*;
they do not establish performance on the repository's
actual production benchmark.

In the existing `test_division_limb`, retain all B20
6000 random GMP oracle cases, all 128/64 edge cases,
and the division exception/guard tests, and add **2048**
deterministic GMP-oracle 128/64 checks with high=0,
high=divisor-1, random constrained high, saturated low
words, and divisors across the full 64-bit range. The
existing 6000 multiprecision differential cases continue
to exercise all quotient-estimate correction paths.
All tests are compiler-portable; no inline assembly in
test sources.

B26 requires Unix, Windows and Coverage GitHub CI
verification and native x86-64 plus ARM64 Integer/GMP
benchmarks with the same six widths (53,128,256,1024,
2048,4096), 20000 iterations, 31 samples and 14
unchanged operations. Prioritize actual x86-64
`quotient` and `remainder` speedups, report changes
in `quotient_trivial` separately, and confirm ARM64
stability. CI and integrated performance are pending
at commit preparation.


### B27 fused multiprecision quotient-step subtraction (19 October 2026)

B27 replaces the normalized multi-limb division's separate
`mul_into(v,qhat,product)`, full product-vs-window comparison,
product correction and subtraction passes with a **single
multiply-subtract pass** updating the normalized dividend
in place. For each divisor limb calculate the exact 64x64
product plus the incoming high carry, subtract its low limb
from the partial dividend and accumulate the subtraction
borrow into the high carry. Subtract the final carry from
the extra high dividend limb. If this indicates the initial
quotient estimate was too large, decrement it and add the
divisor back, propagating the addition carry; repeat only
when the high limb still denotes a negative result.
The proven 64-bit normalized quotient estimator can
overestimate by up to two. No underestimated digit is
introduced. This is the exact radix-2^64 long division
recurrence, not a new approximation.

The product/borrow operations use `unsigned __int128`
on GCC/Clang, MSVC x64 `_umul128` and explicit low-word
carry handling, and a C++20 half-limb 32-bit product
fallback on other toolchains. B26's x86-64 `divq` quotient
estimator and its existing platform fallbacks remain
unchanged. The architecture-independent fused division
step uses **no inline assembly** and needs no additional
compiler flags, CMake edits, ISA extensions or dispatch.
Quotient-only, remainder-only and joint-output modes
share the same implementation; B23 ownership and
signed Integer truncation semantics are unchanged.

The exact public scratch requirement for n dividend limbs
and m divisor limbs, n>=m>=2, shrinks from
`n+2*m+2` to **`n+m+1`**, retaining only the
normalized dividend (n+1) and divisor (m). There is no
temporary product slice. Required scratch for single-limb
divisors or dividends smaller than the divisor stays zero.
All pre-mutation output and scratch checks remain in place.
Document the new scratch contract in the public header.

A standalone GMP-oracle C++20 division probe checked
256400 constrained randomized/boundary multi-limb
divisions using GCC `__int128`; the identical cases
passed using the pure 32-bit half-product fallback,
and using Clang `__int128`. Cases covered magnitudes
up to 128 limbs, 64 distinct normalization shifts,
saturated limbs, high-zero inputs, near-equal operands,
and quotients requiring one or two add-back corrections.
The oracle recorded **7930 two-correction events**
and no case needing more than two. This validates
the independent local core algorithm, not yet the
whole repository build.

Preserve all previous GMP differential tests and add
two fixed double-correction fixtures discovered by the
independent stress probe and 2048 deterministic random
GMP tests exercising the public joint, quotient-only
and remainder-only division operations, including the
exact new scratch requirement. The existing 6000
randomized multiprecision and 64*10 normalization
cases are left intact. Benchmark corpus, 14 operation
labels, six widths, 20000 iterations and 31 samples
must not change. B27 GitHub Unix/Windows/Coverage CI
and native ARM64/x86-64 integer benchmarks remain
pending until run, and no integrated performance
improvement is claimed yet. Compare B27 with B26
separately on each architecture, especially general
`quotient` and `remainder`, and check
`quotient_trivial` for run-to-run noise.

### B28 exact equal-leading-limb quotient 0/1 fast path (19 October 2026)

The B27 benchmark's 32 meaningful division operands consist of 16
single-limb small divisors, eight approximately half-width divisors and
eight near-equal dividend/divisor pairs. Therefore the x86-64
`quotient`/GMP ratio cannot be attributed solely to the normalized
multi-limb multiply-subtract loop. An independent GCC x86-64 probe
also found that replacing one baseline `divq` with two 32-bit
`divl` instructions for small divisors was slower; that candidate
was **not** integrated.

For normalized n-limb magnitudes a,b with n>=2 and the same most
significant limb h, write a=h*B^(n-1)+a0 and b=h*B^(n-1)+b0, where
B=2^64 and 0<=a0,b0<B^(n-1). Since h>=1, a<2*b.
Consequently floor(a/b) is exactly 0 or 1. B28 adds an early return
inside `BigUnsigned::divrem_selected_into` after all existing
pre-mutation output and scratch checks, and after the single-limb
case. Quotient-only scans from the highest limb and returns 0 for
a<b or 1 otherwise without subtraction. For remainder-only and
joint-output, one unsigned multi-limb subtraction obtains a-b and
its final borrow. A zero borrow produces quotient 1 and remainder
a-b; a nonzero borrow produces quotient 0 and remainder a. All
output spans remain zero-initialized as before. The shared B27
normalization/estimation/fused subtraction and correction path is
unchanged for other operands. No ISA feature, external dependency,
CMake option or public API is added, and the exact scratch
requirement n+m+1 is deliberately unchanged even though this fast
path does not touch it.

An independent faithful B27-core reconstruction plus the B28
shortcut passed 25,000 randomized equal-high-limb cases against GMP
for joint, quotient-only and remainder-only outputs; Clang
AddressSanitizer and UndefinedBehaviorSanitizer also passed.
Add 2048 deterministic equal-leading-limb GMP-oracle cases to
`test_big_unsigned`, with exact equality, both orderings,
borrow/carry boundaries and up to 128 limbs. Verify all output
modes and check that public output tails are cleared and scratch
guards are unmodified. Retain all previous B27 and earlier
regressions; benchmark corpus, operations, iterations and samples
must remain unchanged.

**Performance evidence is isolated, not integrated:** an x86-64
AMD EPYC 9V74 environment using GCC -O3, fixed operand pools and
preallocated output/scratch views measured the near-equal 4096-bit
subset at approximately 113.3 -> 32.3 ns (quotient only) and
133.3 -> 84.1 ns (remainder only). In the same local probe the full
32-case mix was approximately 581.4 -> 562.9 ns (quotient) and
607.0 -> 573.0 ns (remainder). Clang -O3 likewise improved the
near-equal subset but full-mix remainder timings were variable
(599.5 -> 618.9 ns in one run). The isolated test excludes owning
Integer allocations, signed wrappers and the exact repository
benchmark harness; its absolute numbers must not be compared
against native user ARM64/x86-64 benchmark values. No integrated
B28 improvement or ARM64 result is claimed. Repeat native six-size,
14-operation Integer/GMP benchmarks and validate correctness;
GitHub CI and Codecov are explicitly deferred by user request
during the temporary GitHub service difficulties.

### B29 selective zero-limb fused subtraction for sparse divisors (19 October 2026)

The B28 six-size ARM64/x86-64 Integer-vs-GMP benchmark confirms
substantial reductions in trivial quotient/remainder, but general
x86-64 quotient remains about 3.60x GMP at 4096 bits. Disaggregating
the **unchanged** 32-pair general division corpus reveals that eight
half-width divisors of the form `2^(bits/2)+1+i` contain long
zero-limb runs: their quotient/remainder steps dominate the isolated
native cost despite accounting for only one quarter of operands.
This workload is mathematically representative of sparse bignums,
not a new fixture contrived to score better.

B29 keeps B27's fused normalized long-division recurrence and
B28's exact equal-leading-limb quotient-0/1 shortcut unchanged.
The shared product/subtraction helper now has two compile-time
instantiations. For a zero limb of the normalized divisor, the
sparse specialization subtracts only the incoming carry/borrow
from the dividend limb, then propagates the resulting borrow:
`window[i] -= carry; carry = old < carry`. This is algebraically
identical to multiplying that zero limb by the current estimated
quotient digit, including the high-product zero and carry rules.
The dense specialization retains the exact original multiplication
code on GCC/Clang `unsigned __int128`, MSVC `_umul128`, and the
portable 32-bit half-limb fallback.

**Selection occurs once before the quotient-digit loop**, after
B27 normalization: use the sparse specialization only if the
divisor has at least eight limbs, there are at least eight quotient
digits, and at least half its normalized limbs are zero. Instantiate
the two quotient-digit loops using a C++20 templated lambda; do not
add a per-digit mode branch to the dense path. This cutoff avoids
replacing the dense kernel with a branch on each multiplication,
which performed worse on small or dense divisors in probes.
The preexisting output and scratch size checks, correction bound,
public `n+m+1` scratch allocation, signing/rounding boundaries,
benchmark input pool, six sizes, 14 operation labels, 20000
iterations and 31 samples are all preserved. No build/CMake flags,
optional CPU instruction or new public API are introduced.

Exploratory standalone x86-64 AMD EPYC 9V74 C++20 probes replicated
the B28 arithmetic body with exact GMP reference outputs. On the
32-pair corpus at 4096 bits, selecting a sparse multiply-subtract
kernel reduced measured quotient/remainder work by roughly
10-30% across GCC 14 and Clang 17 runs, while isolated half-width
sparse-divisor cases often improved by approximately 20-40%.
Absolute isolated kernel timings are not the official owning
`Integer` benchmark and can be sensitive to compilation/code
layout; **no integrated B29 speedup is claimed**. The
single-limb-preinverse alternative showed benefits only in some
4096-bit GCC runs and was not incorporated.

An independent standalone differential stress test exercised
32000 randomized dense and zero-run divisors up to 128 dividend
limbs, comparing B28 versus the selective B29 path and exact GMP
for quotient-only, remainder-only and joint operations. GCC -O3
and Clang 17 AddressSanitizer/UndefinedBehaviorSanitizer runs
both passed. Preserve all earlier differential suites and add
1024 deterministic long-zero-run and mixed-density GMP oracle
cases in the production `test_big_unsigned` with the same
three output variants. The integrated production test suite and
native ARM64/x86-64 six-size performance gates remain outstanding;
CI and Codecov are not being consulted at user request.

### B30 x86-64 radix-10^18 decimal leaves via B26 DIVQ (19 October 2026)

B29's native x86-64 benchmark reports `Integer::to_decimal` at
approximately 14893 ns versus GMP 4431 ns at 4096 bits (3.36x).
B25's attempt to divide an unsigned 128-bit numerator by 10^18
regressed on x86-64 in the integrated benchmark and was disabled
by B25.1. B30 does **not** restore that implementation: it reuses
the exact B26 hardware 128/64 `divide_wide` helper, which already
executes baseline x86-64 `divq` under GNU/Clang. Each decimal
radix-10^18 quotient step has numerator
`(remainder << 64) | limb`, where `0 <= remainder < 10^18`;
hence the hardware `RDX < divisor` precondition holds and both
the quotient and next remainder are exact 64-bit values.

The leaf formatter `to_decimal_blocks_wide` now selects this
existing helper only for GNU/Clang x86-64 with at least nine
64-bit limbs. Preserve B25.1's original `unsigned __int128`
implementation on ARM64, B19's radix-10^9 implementation
for up to eight limbs and for MSVC/other architectures, the
B22 recursive split powers and cutoff of sixteen limbs,
B24 decimal parsing, B27/B29 division, B28 quotient-0/1,
all output contracts and all Integer/FloatAP/Rational/Dyadic
semantics. No additional assembly, new ISA feature, CMake
change, build definition or public API is introduced.

An isolated x86-64 AMD EPYC 9V74 probe measured
median leaf times for 16-limb operands of approximately
2075 -> 1081 ns under GCC 14.2 and 2260 -> 1152 ns under
Clang 17. A recursive decimal conversion proxy with the
same B22 power-selection rule but **GMP** for internal
split divisions measured GCC 4096-bit values around
10384 -> 6909 ns and Clang values around 11811 -> 6638 ns.
These are algorithmic proxies, **not** the repository's
native `Integer::to_decimal` production timing or proof of
cross-host improvements. They address the prior B25
leaf-only-versus-integrated discrepancy more thoroughly
but cannot settle it; require the real six-size 20,000
iterations/31 samples ARM64 and x86-64 benchmark before
claiming an integrated benefit. The native comparator and
benchmark corpus must remain unchanged.

An independent C++20 formatter probe passed 30,000
direct GMP-oracle cases spanning 1..128 64-bit limbs,
radix powers, saturated inputs, long runs of zeros,
and random limb words, with GCC and Clang AddressSanitizer/
UndefinedBehaviorSanitizer. Preserve all existing decimal
and integer regression suites and add 128 deterministic
GMP-oracle cases at 8,9,10,12,16,17,32,64 limbs, including
the formatter selection boundary and recursive split.
Integrated cross-platform test execution and ARM64/x86-64
production timings remain pending. GitHub CI and Codecov
are not being consulted by user request.

### B30.1 revert x86-64 wide decimal leaf after integrated regression (19 October 2026)

The official, **unchanged** six-size Integer-vs-direct-GMP benchmark
has now been run for B30 on ARM64 and x86-64 with 20000 iterations,
31 samples and all fourteen original operations. Numerical
verification passed on both hosts (768 paired inputs, 2688 paired
results per run). It falsifies the promising isolated B30 proxy
performance result for GNU/Clang x86-64: native
`Integer::to_decimal` times at 1024, 2048 and 4096 bits worsened
respectively from B29's **2247.362, 5736.002, 14893.242 ns** to
B30's **3715.747, 7599.686, 18242.587 ns** (+65.3%, +32.5%,
+22.5%). The all-operation total_time/GMP ratio worsened from
2.477x to 2.781x. All sizes up to 256 bits use the old formatter
and were not expected to change. On ARM64, where the code path
was unchanged, corresponding formatting values were
**1407.481 -> 1346.456**, **3565.394 -> 3501.423** and
**9505.252 -> 9339.346 ns**: the small variations do not
establish a causal speedup.

A key diagnostic is that 1024 bits fits the B22 16-limb
nonrecursive cutoff, so the pronounced regression is already
present in B30's x86-64 decimal **leaf** and cannot be blamed
on recursive division or power construction. The native B29
base-10^9 path divides 32-bit halves by a fixed radix, which a
compiler can strength-reduce; B30's `DIVQ` on 10^18 is not
necessarily cheaper on the user's x86-64 CPU despite halving
the number of radix blocks. CPU differences and generated code
may explain the discrepancy with the independent AMD EPYC
proxy, but the real project benchmark is decisive.

**Decision:** B30.1 restores `src/big_unsigned.cpp` byte-for-byte
from the proven B29 commit `3d425d8`. Thus on x86-64 the
original 10^9 block formatter is selected again. ARM64 retains
B25.1's 10^18 wide-division formatter, and MSVC/portable targets
retain their original B19 path. All B29 sparse division,
B28 exact quotient-0/1, B27 fused arithmetic, parsing,
recursive decimal splitting, public APIs, build flags and
benchmark inputs remain unchanged. B30's 128 deterministic
GMP decimal regressions are **retained** because they validate
both radix-formatting implementations; their comment now
reflects the actual architecture selection. No performance
special casing of benchmark inputs is introduced.

The B30 integrated functional benchmark passed prior to
the revert; B30.1 restores the earlier B29 code, but a new
integrated run of the B30.1 tree and its unit tests remains
outstanding. Do not describe the B29 benchmark timings as
new B30.1 measurements. Future x86-64 decimal experiments
should profile the actual native formatter's constant-radix
machine code and the allocation/copying costs, and verify
performance with the **real** BigUnsigned/Integer end-to-end
benchmark rather than a GMP-based recursion proxy. CI and
Codecov remain deliberately unconsulted at user request.

### B31 pairwise decimal string assembly and allocation reduction (19 October 2026)

B30.1 re-established the reliable B29 x86-64 radix-10^9
leaf after B30's official native 1024/2048/4096-bit
regression. Profiling the reconstructed **exact B30.1
leaf** on an x86-64 AMD EPYC 9V74 with GCC 14 and
Clang 17 showed that GCC strength-reduces the constant
1e9 divisions to integer multiplication and shifts,
rather than executing a runtime hardware DIVQ.
Replacing the divisor is not the current target.

B31 retains that B19 arithmetic, B22 recursive
power/split logic, B25.1 ARM64 1e18 wide-radix path,
B27-B29 division and all public APIs and signed
semantics. It changes **only the rendering and
allocation strategy** of `to_decimal_blocks` with
radix 1e9: pre-reserve an upper bound of three
base-1e9 blocks per 64-bit limb; render the highest
block once via `std::to_chars`; allocate an exactly
sized result string; write the lower blocks backward
from least to most significant using four base-100
digit pairs and one leading digit. A local constexpr
200-character lookup maps 0..99 to fixed two-digit
strings. This preserves exactly nine characters per
lower block (including leading zeros) without the
previous repeated `std::to_string`, `append`, and
temporary-string operations. The result-length
calculation is overflow-checked; the constexpr
table and C++17 `std::to_chars` require no new
architecture support, build flag, or CMake change.
The ARM64 wide leaf is untouched; on ARM64 the
radix-1e9 renderer is still used through eight limbs.

An exploratory alternating-order microbenchmark on
the exact prior 1e9 extraction loop and the new
renderer, with both random and the official
near-power-of-two `2^(bits-1)+c` operand shapes,
measured 16-limb GCC baseline ~2198-2206 ns versus
candidate ~1769-1798 ns (roughly 18-20% faster).
Clang ranged from ~2054-2131 ns to ~1719-1788 ns
in short replicated runs (roughly 13-19% faster).
Eight- and nine-limb leaves often improved ~25-34%;
32-limb stand-alone leaves improved ~6-11%, though
production B22 uses the recursive path beyond
16 limbs. These figures isolate a leaf implementation,
not production `Integer::to_decimal`: **no integrated
B31 speedup has yet been shown**. Their limited
transferability is especially relevant following
B25/B30's misleading standalone wins.

The isolated candidate passed 30000 independent
GMP-oracle values (random 1..128 limb words, saturated
limbs, powers of two, long internal zeros) under
GCC -O2 and Clang AddressSanitizer/UBSan. Extend
`test_big_unsigned` with targeted GMP checks of
powers of ten +/- 0,1,9,10,11,99,100,101,999,
including exact nine-digit zero padding and
recursive split boundaries through 576 digits.
All previous tests and official benchmark inputs,
14 operation definitions, iteration/sample counts,
ARM64 selection, GMP boundary rules and other
native arithmetic are unchanged. Full repository
tests and the ARM64/x86-64 integrated benchmark
remain pending; GitHub CI and Codecov are deliberately
not consulted at user request.

### B31.1 shadow warning fix (19 October 2026)

Rename the local counter constructing B31's constexpr digit-pair table
from `value` to `pair_value` to avoid shadowing the
`BigUnsigned const& value` argument in `to_decimal_blocks`
under `-Wshadow`. Only an identifier changes: the 100 table
entries, output bytes, arithmetic, benchmark workload and APIs
remain identical. No new performance or test results are claimed;
GitHub CI and Codecov remain unconsulted by user request.

### B32 selected GCC x86-64 symmetric squaring and generic pow identity elision (19 October 2026)

B31.1's repeated native ARM64 benchmarks confirmed its decimal
conversion improvement and showed the first 256-bit `pow3` spike
was measurement noise (181.7 ns falling to 124.5 ns). The remaining
major 4096-bit costs are `pow3`, `mul` and `square`.
Inspection of `Integer::pow` found that the first set exponent
bit multiplied the initial accumulator `1` by a power of the
base even though assignment is sufficient. The exponentiation
schedule now copies that first factor and executes all subsequent
squarings/multiplications as before. This works for every
nonnegative Nat exponent, including exponent zero, zero base,
negative inputs and even/odd sign; it does not special-case
the benchmark's exponent three.

`BigUnsigned::square` previously called `mul(a,a)`, which
multiplies every pair of limbs twice. B32 introduces a
**GCC x86-64 fused-build-only**, allocation-free symmetric
square at eight or more limbs. The result's existing 2n
limbs hold the sum of products with i<j, calculated once
for each unordered pair; one left shift doubles all
cross terms; adding each a[i]^2 at positions 2i and
2i+1 completes the square with explicit carry
propagation. The cross-term value cannot overflow
2n limbs because it is less than a^2. Per-row
carry placement at i+n is safe: previous rows
reach at most i+n-1. The full 2n-limb result
cannot overflow, including saturated-limb inputs.
GCC x86-64 chooses this path only when its existing
`ARIADNE_BIG_UNSIGNED_GCC_FUSED` selection is
enabled; turning off the existing CMake option
still exercises the unchanged reference schoolbook
kernel. ARM64 retains its verified assembly
carry-chain multiplier, Clang its selected B3
schoolbook, and MSVC/portable their original paths.
No new public API, scratch, ISA assumption,
CMake option, benchmark workload or native
dependency is introduced.

An independently built x86-64 GCC 14.2 -O3
probe compared the existing fused schoolbook
recurrence against the symmetric implementation
with fixed-size operand pools: at 8/16/32/64
limbs the times were approximately
67.7->50.6 / 232.2->148.6 / 929.2->531.5 /
3855.9->2044.8 ns. These correspond to
25/36/43/47% lower isolated kernel time.
The same reference prototype passed 2500
randomized/edge cases at each of 14 sizes
(1..128 limbs), against the original kernel
and independent GMP, under GCC and Clang;
Clang AddressSanitizer/UndefinedBehaviorSanitizer
also passed. **Integrated B32 performance
has not yet been measured** and the isolated
x86 numbers must not be compared directly
with the six-size Integer benchmark. This
does not establish any ARM64 speedup from
the symmetric algorithm, which is not selected
there.

Extend `test_big_unsigned` with 2112 GMP
differential square and mul(x,x) cases, spanning
the 7/8-limb selection boundary through
128 limbs, saturated words, sparsity and zeros.
Extend `test_integer` with 130 independent
GMP pow exponent/negative-sign/zero oracles.
All existing architecture/kernel selection,
test scenarios and official Integer-vs-GMP
workload remain unchanged. The integrated
production test suite and user-host ARM64/
x86-64 benchmarks remain outstanding;
GitHub CI and Codecov are intentionally
not consulted.

### B33 alias-aware Integer self-multiplication dispatch to B32 square (19 October 2026)

The official Integer-vs-GMP B32 x86-64 run confirms a meaningful
`pow3` improvement from B31.1 (4096-bit native ~9509.5 to
~8195.6 ns), but reports `mul` and `square` almost identical
(4096-bit ~3049.2 and ~3046.4 ns respectively). Source-level
inspection explains the latter: the benchmark's unchanged
`IntegerOperations::square(a)` is defined as `a*a`; the
inline `Integer::operator*` forwards both references directly to
`Integer mul(Integer const&,Integer const&)`, which until B33
always called generic `BigUnsigned::mul`. The B32 symmetric
`BigUnsigned::square` was not reached for that ordinary
self-multiplication use case. This is a genuine production API
opportunity, not a request to alter the benchmark workload.

B33 detects **object identity only** (`&z1==&z2`) in
`Integer::mul`, and uses the already verified B32
`BigUnsigned::square` when the magnitude has >=8 limbs,
only in GCC x86-64 builds with the existing
`ARIADNE_BIG_UNSIGNED_GCC_FUSED` definition.
Distinct-but-equal objects continue through generic
`BigUnsigned::mul`; different values are unaffected.
All smaller GCC inputs, Clang, MSVC, ARM64 Apple M2/M4
and portable configurations retain their previous
binary-multiplication path. The square of a negative
Integer is positive, so the alias case constructs an
unsigned magnitude with sign=false; zero stays canonical.
This also correctly implements `x*=x`: the original
input is read completely before the assignment of
the resulting Integer. The check is O(1), based on
semantic aliasing of operand references, and does not
recognize benchmark-generated values or inspect
their numeric pattern. B32's existing full-product
size guard and kernel dispatch remain authoritative.

A fresh stand-alone x86-64 GCC 14.2 isolated
fixed-pool test of the exact B32 limb recurrences
measured schoolbook -> symmetric at 8,16,32,64 limbs
around 66.4->49.5 / 230.8->147.3 /
903.6->515.8 / 3625.6->1922.1 ns
(approximately 25/36/43/47% less kernel time).
Clang's isolated prototype also favored symmetry,
but this does **not** change its production selector,
which remains unmodified. The previous independent
GMP-oracle kernel differential and sanitizer tests
continue to support mathematical correctness; B33
adds 144 deterministic native Integer self-alias,
distinct-but-equal and `x*=x` GMP-oracle fixtures at
1,7,8,9,16,17,32,64,128 limbs and varying signs/zeros.
These new production tests and an integrated B33
binary have **not yet been run** by the assistant.
Thus no integrated x86-64 `square` gain may be claimed
before repeating the user's 84-operation benchmark
on the same x86-64 host, ideally with a GCC fused
build and optionally its preexisting option-OFF
control. Preserve the official square definition
`a*a`, all inputs, 20000 iterations/31 samples,
all existing tests and GMP boundaries unchanged.

ARM64 B32 measurements were recorded on Apple
**M2**, whereas preceding B31.1 values came from
Apple **M4**; absolute cross-host differences
must not be interpreted as B32 code regressions.
No ARM64 code path is changed by B33. Neither CI
nor Codecov was consulted at the user's request.

### B34 alias-dispatch investigation: targeted real-production benchmark (19 October 2026)

The B33 full x86-64 benchmark verified 1024/2048/4096-bit native
square improvements of approximately 32/39/40% over B32. Its
generic `mul` results at 53/128/256 bits, however, increased by
approximately 19/13/10% even though B33 changed only a pointer
identity guard in `Integer::mul` and retained the exact pre-B33
multiplication body for distinct operands.

Tested B33's current identity-first guard against size-first
reordering, an explicit unlikely branch hint and both combined,
using alternating-order GCC 14.2 x86-64 stand-alone probes with
64-bit inline capacity four, the exact fused schoolbook arithmetic
and the B32 symmetric square recurrence. Measured native
1/2/4-limb distinct-operation times clustered around
20/23/44-47 ns; comparing a generic call without an alias guard
to a call with the B33 guard generally yielded differences
around 0-1 ns with variable sign. At 8-64 limbs no alternative
had consistent performance wins across replicates. These
isolated probes are *not* the complete production Integer
implementation, so cannot conclusively explain the user's
5-6 ns small-`mul` shifts. They do show that blindly
reordering the guard or attaching a branch hint lacks a
credible performance justification. Binary layout or
inter-run variability remain hypotheses, not established causes.

**Decision:** Preserve B33's current production `Integer::mul`
and symmetric kernel rather than force an unproven B34 dispatch.
Add a separate, opt-in `benchmark_integer_alias_dispatch`
executable, built against the **same actual NUMERIC_SRC** as
the official Integer benchmark. The new diagnostic records
median ns for public `a*b` on distinct operands, distinct
objects with identical values, and same-object `a*a`.
Its references call `BigUnsigned::mul` via
`Integer::from_magnitude`, and use the actual B32
`BigUnsigned::square` as a third self-square control.
Each width uses 32 deterministic operand pairs and signs,
pre-verifies public/native/reference results against GMP,
and alternates seven operation-order slots across each
sample to limit drift. Report the library's actual
`selected_multiply_uses_fused()` selection. This is
a diagnostic tool **only**; do not alter the established
14-operation Integer/GMP benchmark or its corpus, timing,
iteration counts, fairness checks or output.

To run from the same build on x86-64:
`./benchmarks/benchmark_integer_alias_dispatch --bits 53,128,256,512,1024,2048,4096 --iterations 20000 --samples 31`.
Interpret public/reference differences cautiously because
they include call structure as well as the alias guard;
run B32 and B33 builds on the **same CPU**, alternating
order if possible, for a causal performance conclusion.
The C++ source and new CMake target were added but a
full production build and this actual diagnostic have
not yet been executed by the assistant. All Integer,
BigUnsigned, Rational, FloatAP implementation paths,
existing tests and official benchmarks remain untouched.
B32 ARM64 numbers are Apple M2, B31.1 ARM64 numbers were
Apple M4; never use cross-machine absolute times as
evidence of a code regression. CI and Codecov remain
unconsulted by user request.

### B34.1 diagnostic warning fix and ARM64/x86-64 results (19 October 2026)

The user built and executed B34's production-linked diagnostic
on both Apple ARM64 (reported `fused=false`) and GCC x86-64
(reported `fused=true`), with 20,000 iterations, 31 samples.
Both runs ended with `GMP verification PASS`. Clang on ARM64
reported `-Wswitch-default` because the diagnostic's enum
switch did not contain an explicit `default` case. B34.1 moves
the preexisting `std::logic_error("Unknown diagnostic mode")`
throw into the `default` label of that switch, removing only
the warning. It changes none of the seven measured cases,
their order, timing, GMP verification or arithmetic behavior.

The x86-64 GCC results demonstrate B33's genuine benefit
on aliased self-products by comparing its production
`a*a` with the diagnostic's explicit generic multiply:
512 bits: 94.453 vs 102.012 ns (7.4% faster);
1024 bits: 181.105 vs 253.586 ns (28.6%);
2048 bits: 503.167 vs 809.378 ns (37.8%);
4096 bits: 1834.987 vs 3047.787 ns (39.8%).
The direct `BigUnsigned::square` control at 4096 bits
was 1827.186 ns and confirms that the alias path is
reaching the specialized square rather than generic mul.

For distinct operands at 53,128,256 bits, the same
x86-64 run yielded native/generic pairs 51.296/49.827,
42.475/40.296, and 69.793/67.834 ns, respectively,
corresponding to gaps of about 1.5, 2.2 and 2.0 ns.
Those differences include compiler code generation
and call/dispatch structure, so **must not be treated
as the pure alias-check cost**. The 5-6 ns B33/B32
small-multiplication regression has not been causally
established. For 1024,2048,4096 bits, distinct
native/generic were 257.940/257.924,
821.206/819.870 and 3063.295/3063.707 ns,
indicating near parity for large inputs.

ARM64 native/generic measurements were very close,
as expected for a platform that does not enable
B33's GCC fused alias specialization (4096-bit
distinct 2359.202/2355.177 ns, aliased
2362.954/2361.898 ns). These diagnostic numbers
identify an ARM64 host only, not whether it was
M2 or M4: never infer its specific CPU without
explicit user confirmation. The earlier official
ARM64 B32/B33 runs were Apple M2, whereas
B31.1 was measured on Apple M4.

**Decision:** keep B33 unchanged, neither reorder its
alias guard nor add a branch hint. B34.1 contains no
numeric implementation change and does not alter
the official Integer-vs-GMP corpus. The user has
executed the B34 diagnostic, but the B34.1 warning-only
commit has not undergone a new full build by the
assistant. GitHub CI and Codecov were not consulted.

### B35 selective existing Karatsuba for large balanced GCC x86-64 multiplication (19 October 2026)

B34.1 preserves the successful B33 alias-to-symmetric-square
dispatch. The remaining large x86-64 cost is generic distinct
`Integer::mul` (B33 around 3066 ns vs GMP 1425 ns at 4096 bits).
Before introducing any new assembly kernel, B35 measured the
existing B8 `BigUnsigned::mul_into_karatsuba` recurrence against
the GCC-selected B6 fused schoolbook. The existing exact
Karatsuba recursion uses caller-owned scratch, recursively
computed sums, and the same architecture-selected schoolbook
base cases; no new numeric algorithm is added.

An isolated x86-64 GCC 14.2 -O3 A/B probe reproducing the
production Karatsuba workspace, cross-term correction and
fused schoolbook recurrence used rotating, pre-generated
32-pair pools with independent GMP result verification.
At 32 limbs, cutoff 32 offered no useful speedup; at 64
limbs, cutoff 32 improved a dense-mixed preallocated
kernel by approximately 12-15% and a benchmark-shaped
near-power-of-two pool by approximately 16%. In repeated
**owning** tests that allocate output and Karatsuba
workspace each iteration, 64-limb sparse operands improved
from approximately 3730 to 3148 ns (15.6%).
At 96/128 limbs, sparse owning improvements were
approximately 24%/30%. The exact local numbers
depend on the environment and are **not** official
native Integer-vs-GMP B35 timings. An 8000-case
independent differential test compared cutoff-32
Karatsuba with both the fused schoolbook and GMP on
48..137 limbs, including random, saturated, powers
of two, alternating-zero limbs and unequal operand
lengths. GCC -O3 and Clang ASan/UBSan both passed.

B35 therefore changes only **owning** `BigUnsigned::mul`
on GCC x86-64 builds already selecting the existing
`ARIADNE_BIG_UNSIGNED_GCC_FUSED` kernel, and only
when the *smaller* operand has >=64 normalized limbs
and the operand lengths are within 25% of the larger.
The comparison uses `smaller >= larger - larger/4`
to avoid integer overflow. Compute the exact scratch
size via the existing
`karatsuba_scratch_limbs(n,m,32)`, allocate one private
vector, and call the existing exact
`mul_into_karatsuba(...,scratch,32)` into the usual
owned result. No global scratch, changes in output
ownership, new ISA, CMake setting, public API,
or benchmark-specific input tests. The original fused
dispatch remains selected for 53/128/256/1024/2048-bit
benchmarks and strongly unbalanced operands; 4096
bits corresponds to 64 limbs and can use Karatsuba.
All ARM64 Apple M2/M4, Clang, MSVC, and portable
architectures keep B33's original multiplication;
GCC control builds with the fused option OFF also
retain the previous kernel. B32 symmetric
`BigUnsigned::square` and B33 aliased Integer
self-products are unchanged; preallocated
`BigUnsigned::mul_into` stays fused.

Add 448 focused GMP-oracle owning `mul(lhs,rhs)`
and `mul(rhs,lhs)` regression scenarios spanning
63/64/65, 80/96/128 limbs, balanced and skewed
lengths, dense, all-ones and sparse patterns,
with canonical high nonzero limbs. The original
6000 randomized product tests, B32 symmetric
square cases, Karatsuba scratch/guards and
other benchmarks remain unchanged.
Full integrated tests and actual x86-64/ARM64
six-size Integer-vs-GMP benchmarks for B35
remain outstanding; do not claim an integrated
speedup from this kernel proxy. GitHub CI and
Codecov were not consulted by user request.

### B35.1 Clang -Werror regression fix and integrated measurements (19 October 2026)

The project's Clang 16 production test build with
`-Werror -Wrange-loop-construct` rejected the B35 GMP-oracle
regression corpus in `tests/test_big_unsigned.cpp`: its
`for (const auto [n,m]: {std::pair<...>, ...})`
structured binding makes an unnecessary copy of each constant
`std::pair`. B35.1 changes **only** that B35 loop to
`for (const auto& [n,m]: ...)`. The sequence of 14
size pairs, all 32 variants per pair, 448 GMP products
in each multiplication direction, expected results,
arithmetic selection and official benchmarks remain
unchanged. The user's reported Clang diagnostic is the
reason for the change; successful CI compilation remains
to be verified on the new commit.

Two native x86-64 integrated B35 benchmark runs, with
20,000 iterations and 31 samples, measured 4096-bit
`Integer::mul` at **2752.212** and **2750.796 ns**.
B33 on the same architecture had measured
**3065.983 ns**, so B35's 4096-bit multiplication
benefit was approximately **10.2-10.3%**, and the
two B35 measurements differ by 0.05%. B35's
4096-bit `square` (1830.722/1832.047 ns) and
`pow3` (8142.558/8146.656 ns) remained stable.
Both x86-64 runs and the ARM64 B35 benchmark
reported numerical verification PASS with 768
paired inputs and 2688 paired operation results.
Do not infer causality for timing drift in unrelated
operations. Earlier ARM64 B35 host was Apple M2;
never compare it directly with the Apple M4 B31.1
baseline.

At the user's explicit request of 19 October 2026,
**GitHub Actions CI may be inspected again starting
with this B35.1 commit**. Prior exclusions of CI
checks applied to older commits; this permission
does not itself establish a passing build. Existing
numeric test coverage and correctness gates are
preserved, and Codecov status should be described
separately from compilation and functional CI.

### B36 exact high-limb quotient-0/1 extension (19 October 2026)

After B35.1 cleared Unix and Coverage CI (including Clang 16,
AppleClang 21 and ThreadSanitizer), the initial Windows run was
still active and later completed successfully. The user
explicitly reauthorized checking GitHub Actions CI.

The remaining x86-64 division gap (B35 repeat: 4096-bit
`quotient` 1243.979 ns vs GMP 362.277 ns and
`remainder` 1116.712 ns vs GMP 364.261 ns) motivated
a correctness-first examination of the normalized long-
division loop. Skipping multiply/subtract for zero quotient
digits is algebraically valid, but the 4096-bit half-width
sparse-divisor corpus has no zero digits in its quotient.
Do not select this candidate to manufacture benchmark speedups.

B36 instead extends the architecture-independent B28
exact 0-or-1 quotient shortcut for any nonzero canonical
`n==m>=2` magnitudes with different high limbs:
if `a.back()<b.back()` then `a<b` and quotient 0,
remainder a. If `a.back()>b.back()` and
`a.back()-b.back()<b.back()`, then
`b<=a<2*b`, giving quotient exactly 1 and remainder
one full n-limb subtraction. The subtraction's final
borrow is zero in this branch. The subtraction-free
quotient-only case directly writes 1; remainder-only
does not create any quotient. Equal high limbs use
the existing B28 scan/borrow algorithm with unchanged
semantics, including a<b determined from lower limbs.
All other cases use the existing B27/B29 normalized
digit-estimation, sparse multiply-subtract and correction
algorithm. The high-limb comparison avoids overflow
from evaluating `2*b.back()`; no new CPU features,
allocations, scratch, APIs, build flags or CMake changes
are required. Shared logic preserves the exact full
public scratch/output prechecks and tail initialization
on all code paths.

An independent exact-arithmetic Python stress probe
verified the high-limb classification against true
arbitrary-precision integer division on 70,000 random
equal-length inputs (2..128 limbs): 34,845 zero-quotient
cases, 17,690 strict high-dominance one-quotient cases,
and 17,465 fallback cases. This validates the shortcut's
arithmetic classification, not compiled C++ integration
or any runtime performance improvement.

Extend production `test_big_unsigned` with 1,280
deterministic independent GMP-oracle cases, each checked
via combined, quotient-only and remainder-only paths
for below/equal/above high limbs, near doubling
boundaries, all-ones, sparse and arbitrary low words.
All B28/B29, B35 and earlier cases are retained.
The official 14-operation benchmark and all its
operand pools remain unchanged. The half-width divisor
inputs and B28 equal-high near-equal cases are already
handled by their previous paths, so do **not** predict
a substantial benchmark benefit from B36 without
same-machine integrated measurements. Verify no
regressions on Apple ARM64 M2 vs M2 and GCC/Clang
x86-64. GitHub Actions may be checked; Codecov is
not required to establish numeric correctness.

### B37 native ARM64 symmetric square for Integer self-products and powers (19 October 2026)

Apple M2 B36 integrated timings reveal a disproportionate ARM64
square cost relative to GMP: at 1024/2048/4096 bits,
native `square` was 220.983/656.942/2298.379 ns versus
GMP 88.683/268.323/772.415 ns (2.49x/2.45x/2.98x).
The corresponding `pow3` operation was
590.819/1892.640/7303.648 ns versus
352.962/921.167/2818.185 ns (1.67x/2.05x/2.59x).
The official benchmark's `square` implementation is `a*a`,
not the named `sqr` API. The power algorithm itself is
already a binary method with the B32 first-identity-multiply
elision; `pow3` requires one `square(base)` and one
general product. Prioritize the underlying square rather
than a benchmark-specific `pow3` shortcut.

Until B37, ARM64 `BigUnsigned::square(a)` called generic
`mul(a,a)`, performing all n*n products; B32's specialized
symmetric square dispatch was restricted to GCC x86-64.
B37 adds `BigUnsignedDetail::square_into_symmetric_arm64`
in the existing ARM64 TU. It computes only off-diagonal
`a[i]*a[j]` for i<j using exactly the selected baseline
`MUL/UMULH/ADDS/ADC` carry-chain assembly, doubles the
cross-term limb array, and adds all diagonal squares with
exact 128-bit unsigned accumulation and carry propagation.
The sum is the mathematical square and fits exactly in
2*n limbs. There are n*(n+1)/2 wide products versus n*n
for general multiplication. Inputs/output must be
disjoint; the own-result wrapper provides the exact
output size. No runtime CPUID, non-baseline ARM ISA,
global storage, new API, CMake setting, ABI change
or additional heap scratch is introduced.

Dispatch ARM64 AppleClang/GCC/Clang with selected native
`ARIADNE_BIG_UNSIGNED_ARM64` only for magnitudes
>=8 limbs; preserve the previously selected carry-chain
schoolbook for smaller operands. On ARM64, extend the
existing B33 self-alias `Integer::mul(z,z)` to call
`BigUnsigned::square` for >=8-limb same-object operands,
with a size-first branch; numerically equal but separate
objects retain generic multiplication, including the
usual sign semantics. All paths on x86-64 GCC/Clang,
Windows/MSVC and portable targets are unchanged.
The new kernel also accelerates `sqr` and internal
`Integer::pow` squaring automatically. Do not alter
the official 14-operation Integer-vs-GMP benchmark,
its six widths, GMP oracle, iterations or samples.

An independently compiled C++20 reconstruction of
the exact cross/doubling/diagonal recurrence passed
**8160** GMP-oracle cases with Clang 17
AddressSanitizer/UndefinedBehaviorSanitizer, including
saturated, zero, sparse, alternating-zero and carry
boundaries at 1..128 limbs. The six-instruction
AArch64 multiplication-accumulation assembly block
also cross-compiled with Clang 17's
`--target=aarch64-none-elf`, yielding the intended
MUL, UMULH, ADDS and ADC instructions.
These are local algorithmic and assembly-syntax
checks, **not** a native Apple M2 performance
measurement or complete project build.

Add 1,248 GMP-oracle unsigned square fixtures at
1,4,7,8,9,15,16,17,32,33,64,65,128 limbs,
comparing both `square(x)` and generic `mul(x,x)`,
plus 216 GMP-oracle signed large `Integer::pow`
cases at 7..64 limbs for exponents 2/3/5.
Preserve B32/B33 and all previous tests.
The full project test suite, GitHub Actions for B37
and same-M2 31-sample integrated performance
comparison against B36 remain pending until run.
Specifically gate on `square` at 1024,2048,4096
bits, `pow3` at all widths, `mul` as a control,
and confirm no regression at 53/128/256 bits.
No integrated speedup is claimed before those
measurements. Keep M2 and M4 baseline results
separate. The user has reauthorized GitHub CI
inspections; report third-party Codecov statuses
separately.

### B38 skip exact-zero rows in native ARM64 multiplication and symmetric square (19 October 2026)

The B37 Apple M2 integrated benchmark improved 4096-bit
`square` from B36's 2298.379 to 1476.506 ns (35.8%) and
`pow3` from 7303.648 to 6380.048 ns (12.6%). At
1024/2048 bits, `pow3` improved to 506.440/1667.752 ns.
All B37 Unix, Windows and Coverage GitHub Actions later
completed successfully (including AppleClang, Clang 16,
GCC and ThreadSanitizer). The generic ARM64 product in
`pow(x,3)` is an unbalanced product of a roughly n-limb
operand and its roughly 2n-limb square. For values
containing long zero-limb runs, both the native generic
schoolbook kernel and B37's symmetric square were
needlessly executing wide multiplications by exact zero.
This is a general arithmetic property, not a `pow3`
or benchmark-operand-specific recognition.

B38 preserves B37's six-instruction
MUL/UMULH/ADDS/ADC accumulation for each **nonzero**
outer-limb row but skips the entire row when `a[i]==0`
in `mul_low_into_carry_arm64`. Since the output is
zero-filled up front, a skipped row adds exactly zero;
its final carry slot was never touched by earlier
rows. The same exact-zero row elimination applies to
`square_into_symmetric_arm64` for the off-diagonal
upper triangle. It is safe because slot `output[i+n]`
is not touched before row i. Doubling and diagonal
square accumulation still process **all** limbs,
including those skipped in the off-diagonal step.
There is no division, dynamic input classification,
new allocation or new ISA. This also preserves
full, truncated, oversized-output and asymmetric
`BigUnsigned::mul_low_into` contracts, including
the reversed order of the operands; no special
cases of `Nat(3)` or of the official operand pool
were introduced. The generic multiplication API,
`Integer::pow`, x86-64, AppleClang selection,
CMake build flags and official benchmarks remain
unchanged. A dense input incurs only an outer-row
zero comparison, so check actual Apple M2 small and
dense cases before accepting a performance claim.

A local x86-64 portable C++20 reconstruction of
the ARM64 generic row-schoolbook recurrence compared
the unmodified and zero-row-skipping algorithms
against independent GMP over **16,000** full/low/
oversized output cases with operands 1..140 limbs,
random, saturated, alternating sparse and zero.
All passed. On the same proxy, 64x128-limb sparse
products reduced the *inner schoolbook work* from
64 rows to 2 and measured about 6471 -> 249 ns,
while dense 64x128 results varied by about 1%.
These are **x86 proxy** results, not actual
ARM64 or production integrated `pow3` timings.
A separate 16,000-case GMP check of the exact
B37 symmetric square recurrence with zero-row
skipping passed under Clang ASan and UBSan.

The production GMP regression corpus now additionally
covers 336 mixed dense/sparse and balanced/unbalanced
unsigned products (both orders), each with six low
truncation sizes, and 120 signed GMP-oracle `pow3`
cases on 8..128 limbs. All previously tested
saturated, no-zero, high/low-boundary, x*=x and
power cases remain intact. The main Apple M2
benchmark remains the unchanged 84-pair
Integer/GMP 31-sample workload. The complete
production build/CI and B38 same-host integrated
benchmarks are **pending**; do not attribute
the x86 proxy's large speedups to Apple Silicon.
Compare B38 only against B37 on the same Apple M2,
and use B37 x86-64 as an unchanged-backend control.

### B39 adaptive ARM64 dense/sparse schoolbook and selective owning Karatsuba (19 October 2026)

The same-M2 dense benchmark with cutoff-32 Karatsuba
compared B37 and B38 directly. On 4096-bit dense
operands, preallocated ARM64 selected schoolbook
regressed from B37's **2284.963 ns** to B38's
**2548.485 ns** (+11.5%), whereas the preallocated
Karatsuba control changed only 1964.965 -> 1972.463 ns
(+0.4%) and GMP 993.337 -> 992.431 ns. At
1024/2048 bits, B38's selected times were
154.600/576.502 ns versus B37 154.896/573.425 ns.
Karatsuba cutoff 32 was slightly slower than
selected schoolbook at 1024/2048 bits, but at
4096 bits it improved B38 preallocated time
2548.485 -> 1972.463 ns (22.6%) and B37
2284.963 -> 1964.965 ns (14.0%). The Karatsuba
control reuses caller-owned scratch, so those
timings do **not** establish an equally large
speedup when the owning `mul` allocates scratch.

B39 addresses this with two exact compile-time
instantiations of the **existing** native ARM64
MUL/UMULH/ADDS/ADC row-schoolbook recurrence.
`mul_low_into_carry_arm64_impl<false>` has no
zero-row guard, matching B37's dense inner loop.
`...<true>` retains B38's exact-zero-row skip.
A cheap linear prescan selects the sparse version
only when at least `max(2,floor(active_rows/8))`
of the active outer rows are zero; the scan stops
once that count is reached. For active row counts
<8 or multipliers <8 limbs, retain the B38 sparse
implementation without any prescan. In full and
truncated cases active rows are
`min(a.size(),output.size())`. No reordered
operand assumptions, allocations or altered
alias/scratch/output contracts are introduced.

The B37 symmetric square gets the same template
separation: the cross-product triangle is dense
and branch-free unless the zero-row prescan
selects B38's sparse version; doubling and the
full diagonal update remain identical on both
paths. The last input limb is excluded from the
cross-row density count because it is never an
outer cross row. The branch only chooses an
equivalent exact algorithm, not a different
mathematical result.

On ARM64 only, **owning** `BigUnsigned::mul`
now dispatches to the existing generic B8
Karatsuba recurrence (cutoff 32) if both input
magnitudes have >=64 limbs, are balanced within
25% of the longer magnitude, and each contains
no more than floor(length/8) zero limbs. The
test is early-exit for sparse operands, so the
two-nonzero-limb 4096-bit official benchmark
still uses sparse B38, and the approximately
n x 2n general product inside `pow3` remains
on the original ARM64 schoolbook path.
Only this owning path allocates the exact private
scratch; the preallocated `mul_into` path and
the public scratch/output contracts do not
change. Non-ARM64 implementations including
the existing GCC x86-64 B35 Karatsuba dispatch
remain unchanged. No specific exponent,
benchmark input shape, new ISA, CMake option,
compiler flags, GMP adapter or API changes.

Add 380 GMP-oracle test shapes over 19
dimension pairs (7..128 limbs), including
63/64, 64/85 and 64/86 boundaries,
both multiplication orders, exact zero-row
density threshold crossings, balanced,
unbalanced and saturated-high operands,
plus square and low-truncation verification.
Preserve all B38 cases, B37 power tests,
independent GMP-oracle comparison and the
official 84-result benchmark unmodified.
Actual native M2 timings of the adaptive
dense path and of owning Karatsuba, along
with integrated pow3 results and CI,
remain **unmeasured on B39**. Compare
B39 dense 4096-bit `prealloc_selected`
against B37 2284.963 ns, and new_selected
against B37 2318.360 ns; do not conflate
preallocated Karatsuba and owning timing.
Confirm the sparse 4096-bit `mul`,
`square` and `pow3` do not regress
relative to B38 158.777/228.442/565.898 ns.

### B40 eliminate medium-size density prescans and Karatsuba leaf overhead (19 October 2026)

The same-M2 B39 follow-up of the 31-sample official
benchmark preserved B38's main sparse results:
4096-bit `mul` 150.442 ns, `square` 237.315 ns,
`pow3` 563.371 ns (versus B38
158.777/228.442/565.898 ns); the total-time
Integer/GMP ratio stayed 1.316x. All three B39
GitHub Actions workflows completed successfully.

A separate same-M2 dense multiplication benchmark
identified *new* B39 dispatch costs at moderate
sizes. Compared with B38, dense preallocated
selected products rose 154.600 -> **166.696 ns**
at 1024 bits (+7.8%) and 576.502 -> **600.090 ns**
at 2048 bits (+4.1%). For the 4096-bit dense
case B39 recovered much of B38's regression:
2548.485 -> 2343.858 ns, while large owning
products used Karatsuba and improved
2575.535 -> 2265.110 ns. Repeated B39 density
classification in cutoff-32 Karatsuba's
32/33-limb recursive leaf products is also
consistent with the measured preallocated
Karatsuba changes from B38 to B39:
158.329 -> 170.012 ns (1024),
580.327 -> 603.987 ns (2048),
1972.463 -> 2065.662 ns (4096).
These observations motivate but do not
individually prove the source of all timing drift.

B40 conservatively **bypasses the B39 linear
density prescan** in the ARM64 carry-chain
multiplier whenever the active outer row count
`min(a.size(),output.size())<64` **or** the
multiplier contains fewer than 64 limbs.
Those cases dispatch directly to the already
validated B38 exact-zero-row implementation;
the dense B37-like branch-free specialization
is still available at >=64x64 limbs after
B39's exact zero-density scan. This restores
the B38 dispatch behavior in 1024/2048-bit
dense full products, cutoff-32 Karatsuba's
small recursive base cases, and ordinary
low-truncated multiplication. At high
sizes no density threshold, algorithm,
instruction set or output contract changes.

For ARM64 symmetric squaring, there are only
`n-1` off-diagonal outer rows; B40 bypasses
classification at `n<=64`, invoking the
B38 symmetric kernel directly at 4096 bits
and below, including cases with nonzero
dense diagonals. At larger sizes B39's
adaptive choice stays in place. There are
still two compile-time native
MUL/UMULH/ADDS/ADC variants, no new
allocations, no scratch changes, no new
features or runtime flags. Dense balanced
>=64-limb *owning* multiplication retains
B39's cutoff-32 Karatsuba selection and
density gates. Small, sparse or imbalanced
owning products continue down the selected
ARM64 schoolbook. The Karatsuba recursion
itself is unchanged; only the base-case
multiplier ceases to repeat an O(n) scan
for n=32/33. Caution: this should not be
described as a measured B40 speedup until
native M2 data are collected.

Add **288** deterministic GMP-oracle test
shapes at eighteen balanced and unbalanced
size pairs spanning 7,8,16,31,32,33,
63,64,65 and 128 limbs. Check owning
multiplication in both orders, low products
in both orders, the public symmetric square,
and for large shapes also explicit cutoff-32
Karatsuba with caller-owned scratch.
Cover dense, all-ones, moderately sparse,
threshold-zero, and highly sparse limb
patterns. Keep all earlier tests,
architecture-independent GMP reference,
both design journals and the official
84-operation benchmark unchanged.
Full production C++ CI and same-M2
dense/sparse benchmarks for B40 are
pending; do not compare ARM M2 absolute
timings to any other CPU.

### B41 avoid temporary negation and trim common high limbs in subtraction (19 October 2026)

In the same-M2 B40 31-sample official Integer/GMP benchmark,
`sub` at 1024/2048/4096 bits measured
71.156/97.917/159.254 ns versus GMP
33.467/41.702/57.587 ns, i.e. 2.13x/2.35x/2.77x.
The signed benchmark operand corpus has
`a=+/- (2^(bits-1)+13+2*i)` and
`b=+/- (2^(bits-1)+31+5*i)`, with
32 deterministic sign combinations. Half of
the subtractions have opposite signs and
therefore require a true n-limb **addition**
of unsigned magnitudes. The other half
subtract nearly identical n-limb magnitudes
and produce a small signed result. Treating
every timed `sub` as a long unsigned borrow
chain would be technically incorrect.

Before B41, `Integer::sub(z1,z2)` invoked
`add(z1,neg(z2))`. `neg` produces a new
owning `Integer` and **copies the entire
`BigUnsigned` magnitude**, even though
the signed subtraction only requires
reading the sign and magnitude of z2.
B41 performs direct signed-magnitude
dispatch without materializing this negated
copy. If signs differ, add magnitudes and
keep z1's sign; if signs match, compare
magnitudes, return canonical zero when
equal, otherwise subtract larger-minus-
smaller and select z1's sign for larger z1,
its inverse for larger z2. All four sign
combinations, zero, x-x and x-=x preserve
the exact original result and canonical sign.

Unsigned `BigUnsigned::sub(a,b)` formerly
pre-scanned the full numbers with `compare`,
allocated `a.size()` output limbs, subtracted
all limbs using `sub_into`, then normalized
the entire result. B41 first rejects a shorter
minuend without allocation. For equal-width
canonical inputs, it scans from the high limb
to find the highest **different** position k.
If no position differs, return zero directly;
if a[k]<b[k], throw `underflow_error`
without modifying an operand. Otherwise the
high equal prefix cancels, and the entire
nonnegative difference is represented by
the low k+1 limbs. Allocate precisely that
result width and call the existing exact
`sub_into` on only those limbs, then
normalize. For a strictly longer minuend
there is no high equal-width prefix, and
the full-width subtraction remains exact.
This avoids duplicate whole-number
comparison within `BigUnsigned::sub`
and reduces work/allocations for all
equal-prefix values, not only the benchmark
pool. Full low-range borrow propagation,
underflow behavior, canonical-zero and
output-normalization invariants remain
preserved. The lower-level public
`sub_into` signature, bounds checks,
output-tail clearing and borrow-out behavior
are not modified. No new ISA,
platform-specific selection, CMake flags,
benchmark operands, library dependencies,
scratch changes or API changes are introduced.

A standalone GMP differential probe of the
common-prefix subtraction recurrence
passed **18,000** random and adversarial
cases on x86-64 under GCC ASan/UBSan,
including 7,399 correct unsigned underflow
exceptions. Production tests add 1,536
unsigned GMP comparisons including reverse
order and underflow, plus 768 signed GMP
subtractions across all four sign pairs,
with `-=` and same-object `x-=x`
checks. All existing architecture-wide
and GMP-backed regression scenarios
remain intact. Note: the production
tests and full CI are **pending** on
the B41 commit. The same-M2 integrated
benchmark must measure B41 versus B40,
especially `sub` at 1024/2048/4096
bits, while using `add` and `negate`
as controls. Avoid claiming a speedup
before same-host timing data; inspect
Windows/macOS/Linux CI separately from
external Codecov.

### B41.1 Clang 16 compatibility of same-object subtraction regression (19 October 2026)

The native Apple M2 B41 full 31-sample official
benchmark confirmed large speedups without changing
the official workload. B40 -> B41 `sub` at
1024/2048/4096 bits improved
71.156 -> 27.744 ns (61.0%),
97.917 -> 44.462 ns (54.6%) and
159.254 -> 73.335 ns (54.0%).
The 84-case geomean Integer/GMP improved
0.841x -> 0.785x and total_time
1.312x -> 1.293x; GMP differential
verification passed (768 paired inputs,
2688 paired operation results).
The new largest relative bottleneck is
`to_decimal` at 4096 bits (9434.652 ns
versus GMP 4069.252 ns, 2.319x).

B41's Unix CI failed during **compilation**
on Ubuntu Clang 16 Release and
ThreadSanitizer at
`tests/test_integer.cpp:463`, not
during an arithmetic test:
`-Werror,-Wself-assign-overloaded`
diagnosed the intentional
`aliased-=aliased` regression
statement. Both macOS AppleClang,
Ubuntu GCC, GCC/AppleClang Coverage,
and Windows workflows passed B41.
This compiler-specific diagnostic
does not by itself establish an
arithmetic or sanitizer defect.

B41.1 changes **only** the intentional
same-object test's spelling to:
`Integer const& same_object=aliased;`
followed by `aliased-=same_object;`.
Both names still reference the exact
same `Integer` object; the test
retains its aliasing, canonical-zero
and nonnegative-sign assertions.
No source arithmetic, API,
compiler warning settings,
test scenarios, workload, oracle,
scratch contracts, or other
platform paths change. A local
Clang 17 `-Wself-assign-overloaded
-Werror` compilation of the
same reference-alias idiom passed;
full Clang 16 release/TSan and
all CI remain pending until the
B41.1 workflows finish.

### B42 single-allocation ARM64 radix-10^18 decimal formatting (19 October 2026)

B41.1 passes all Unix, Windows and Coverage GitHub Actions, including Clang 16 Release and ThreadSanitizer. B41 Apple M2 measurements of `to_decimal` at 1024/2048/4096 bits were 1514.429/3675.169/9434.652 ns; GMP at 4096 bits was 4069.252 ns (2.319x). B42 retains the exact existing recursive splitting and wide 128/64 quotient/remainder algorithm but optimizes formatting of the ARM64 radix-10^18 leaves: the decimal block vector reserves at most two blocks per input limb, then one exactly sized result string is filled using `to_chars` for the leading block and constexpr digit-pair lookup for all nine pairs in each remaining 18-digit block. No per-block `std::to_string` or temporary allocation is used. Other platforms' radix-10^9 path, all arithmetic kernels and the official benchmark are unchanged. A standalone C++20 reconstruction of the wide division and formatter passed 1800 GMP `mpz_get_str` comparisons with Clang ASan/UBSan over 0..256 limbs, including sparse, saturated and zero values. Production regression adds 208 independent GMP decimal cases at 8..256 limbs, including 8/9 and 16/17 dispatch boundaries. B42 production CI and Apple M2 timings are pending.

### B43 append-only recursive decimal output and one root string (19 October 2026)

B42 on Apple M2 reduced native `to_decimal` versus B41.1 from
1514.429 -> **1178.808 ns** at 1024 bits (-22.2%),
3675.169 -> **3085.894 ns** at 2048 bits (-16.0%) and
9434.652 -> **8564.325 ns** at 4096 bits (-9.2%).
The corresponding GMP B42 4096-bit conversion was 4077.340 ns,
leaving a 2.100x slowdown. The B42 official benchmark verified 768
input pairs and 2688 operation results. All Unix, Windows and
Coverage GitHub Actions passed (external Codecov statuses remained
separate). The residual is not proven to lie in any single stage:
decimal-power construction, multiprecision `divrem`, and leaf
conversion were not isolated on the user's Apple M2.

B43 targets a separate, demonstrable source of allocation/copy
overhead in the **shared** B22 recursive string assembly. Previously
`to_decimal_recursive(value,powers)` returned a new string at
every recursive node, constructing separate high and low strings,
padding the low string to the exact divisor width and concatenating
both. B43 keeps the same power table, selected splitting divisor,
exact `divrem` quotient and remainder, 16-limb recursion cutoff
and architecture-selected leaf conversion, but threads **one**
`std::string& output` down the recursion. At each high branch it
passes only the inherited minimum decimal field width not already
owned by the lower branch; each low branch receives exactly its
divisor's decimal width. Only leaf nodes format or emit leading
zeroes, so root output has no leading zeros and all lower fields are
exactly zero-padded. This handles recursively split *remainders*
that themselves have many digits, including zero, small residuals
and long internal decimal-zero runs. The root pre-reserves
`bit_length/3+1` characters as a conservative, overflow-safe
upper bound on decimal length (`log10(2)<1/3`), and shares its
already computed bit length with power-table creation.
There are no caches, GMP calls, altered benchmarks, new public
APIs, new arithmetic kernels, ISA requirements or division changes.

Unlike B42, this string-assembly optimization is **shared**
across ARM64 and x86-64, using each architecture's existing
radix-10^18 / radix-10^9 leaf formatter; it does not reenable the
historically slower 128/64 decimal division on x86-64.
Do not assume a cross-platform speedup until independently
measured on each host. The change removes intermediate recursive
string allocations, not the leaf string creation or the
multiprecision division cost; the empirical 7000 ns (4096-bit M2)
goal remains unverified.

Before integration, a standalone C++20 reconstruction of both
the original recursive string assembly and B43's exact minimum-
width append recurrence was run with independent GMP quotient/
remainder and GMP decimal conversion for the leaves. All
**880** randomized/adversarial 0..256-limb scenarios passed under
Clang AddressSanitizer and UndefinedBehaviorSanitizer,
including powers of two, zeros, saturated operands, random
dense inputs and distant nonzero limbs. Production GMP decimal
tests add **80** deliberately nested decimal-split and
zero-padding cases across 144,288,576,1152,2304 decimal
exponents, with independent `mpz_get_str` oracles.
All earlier B22/B25/B31/B42 decimal scenarios and the
official 84-operation benchmark are retained.
B43 CI and same-M2 performance measurements are pending.



### B44 standalone production-decimal stage diagnostic (19 October 2026)

B43's Apple M2 native 4096-bit to_decimal changed from 8564.325 to
8512.402 ns (-0.61%), while GMP changed from 4077.340 to 4116.558 ns.
This is not conclusive evidence of a speedup. All B43 Unix, Windows
and Coverage GitHub Actions passed; external Codecov patch/project
checks failed separately. Keep the official 84-operation Integer/GMP
workload, operands and reporting unchanged.

B44 introduces only benchmarks/benchmark_integer_decimal_profile.cpp
and its separate CMake benchmark target, without modifying production
arithmetic, ISA dispatch, public APIs, or any existing benchmark file.
It reuses integer_specifications(bits), the exact 32-case official
signed Integer input pool, and profiles the corresponding magnitudes.
Both production to_decimal and an independent B43-style split/append
reconstruction are compared with GMP's decimal oracle before timing.
The default sizes are 1024/2048/4096 bits, with 20000 iterations and
31 samples. Example:
./build/benchmarks/benchmark_integer_decimal_profile --bits 1024,2048,4096 --iterations 20000 --samples 31

The diagnostic measures complete native to_decimal, and independently
replays B22 power-table construction via production square, recursion
divisions via production divrem, leaf conversions via production
to_decimal (its <=16-limb leaf dispatch), and append/padding with
precomputed exact leaf strings. The 1024-bit size bypasses recursive
splits. The recursions are reconstructed outside all timed regions.
These are same-input non-additive stage PROXIES, not exact in-situ
call-stack timings; do not sum their fractions or claim causal
bottlenecks without repeatable same-host evidence. Apple M2 B44
measurements remain pending.


### B44.1 logging-compatible output and user diagnostic measurements (19 October 2026)

The B44.1 standalone benchmark fixes an output-only bug: although
`logging_command_line()` already parses `-t/--theme`, `-v/--verbosity`
and `-s/--scheduler`, B44 emitted its benchmark report through
`std::cout`, outside the Ariadne logger. The B44.1 benchmark now
uses `LOGGING_PRINTLN` for all normal report lines, including component
measurements and the final verification banner. This restores normal
logger theme/verbosity behavior, following all other benchmark
executables. It does not alter the sample loops, clock, input pool,
GMP oracle, arithmetic, B43 recursion or any other benchmark.
The standalone diagnostic's `--summary-only` and `--details` options
are shared parser options and do not change its stage measurements.
Fatal errors still go to `std::cerr`.

The user-supplied B44 diagnostic (20000 iterations, 31 samples, 32
official-input magnitudes per size) passed GMP verification, with:
- 1024 bits: full production decimal 1073.046 ns, no recursion;
- 2048 bits: full 2942.338 ns; power-table replay 261.742 ns;
  recursive divrem replay 960.735 ns; production leaf conversion
  1649.181 ns; append/padding replay 40.525 ns;
- 4096 bits: full 8236.562 ns; power-table replay 361.692 ns;
  recursive divrem replay 4440.488 ns; production leaf conversion
  3387.277 ns; append/padding replay 50.379 ns.
These are independent, non-additive proxies: no individual replay
percentage is proven to be exact production call-stack attribution.
They suggest leaf conversion dominates the measured replays at 2048
bits and recursive division at 4096 bits, while append/padding is
small at both. Measurements are user-provided, not independently
reproduced; native performance improvement is not claimed.
B44's Unix, Windows and Coverage GitHub Actions passed; the separate
Codecov patch/project checks were failing.


### B45 remove exact whole-zero divisor limb prefixes (19 October 2026)

The user-supplied B44 4096-bit Apple M2 diagnostic measured production
to_decimal 8236.562 ns; independently replayed recursive divrem
4440.488 ns, leaf conversion 3387.277 ns, decimal-power construction
361.692 ns, and append/padding 50.379 ns. Stage timings are separate
non-additive proxies, not exact production CPU attribution. At 2048
bits leaf conversion was the largest measured replay (1649.181 ns).

B45 optimizes general shared BigUnsigned long division for a divisor
with k>=1 complete zero low limbs. Let B=2^64, a=a_hi*B^k+a_lo
with 0<=a_lo<B^k, and b=b_hi*B^k. Then the quotient equals
floor(a_hi/b_hi), and the remainder equals
(a_hi mod b_hi)*B^k+a_lo. After the original argument-size and
scratch validations and output initialization, recurse into the
same selected kernel on the canonical high input limb subspans,
reusing the existing scratch area. Copy the original low dividend
limbs into the remainder only when requested. The original
joint-output, quotient-only and remainder-only results, public
scratch contract, zero-filled surplus output limbs and no-alias
requirements are all preserved. This exact optimization skips
normalization, quotient-digit products and correction work for
provably zero divisor limbs. Decimal powers 10^(9*2^j) contain
many such low limbs because of their factor 2^(9*2^j).

Add 504 deterministic GMP-differential shifted-divisor cases
across sizes 2,3,5,8,16,32,64 limbs and zero-prefix selections,
checking all three public result modes and guarded output buffers.
Add 40 GMP division cases using decimal powers 10^72,10^144,
10^288,10^576 and magnitudes 2^(bits-1)+offset for bits
2048 and 4096, matching the official benchmark operand family.
No original tests or official benchmark workloads are changed.
Both design journals describe the general arithmetic proof.
Full CI and same-host Apple M2 benchmarking are pending.
Do not claim an actual B45 speedup until those measurements.


### B46 ARM64 radix-10^18 exact reciprocal leaf division (19 October 2026)

B45's independent Apple M2 stage benchmark showed a reduction of the
4096-bit decimal recursion's divrem replay from 4440.488 to 3678.875
ns and full to_decimal from 8236.562 to 7520.004 ns. The unchanged
official Integer/GMP benchmark measured native to_decimal at 4096 bits
falling from B43 8512.402 ns to B45 7818.717 ns; native/GMP geomean
changed 0.782x to 0.775x and total_time 1.197x to 1.151x. Correctness
verification passed 768 paired inputs and 2688 operation results.
All B45 Unix, Windows and Coverage Actions succeeded; external
Codecov patch/project statuses remained failing. The remaining B45
stage replays at 4096 bits were divrem 3678.875 ns, leaf conversion
3464.458 ns, power construction 380.217 ns and append 51.048 ns;
these separate timings are non-additive. At 2048 bits the leaf replay
was 1674.200 ns and the full production conversion 2813.100 ns.

B46 changes ONLY the existing native ARM64 wide-leaf formatter
to_decimal_blocks_wide in src/big_unsigned.cpp. It keeps the radix
10^18, >=9-limb ARM64 dispatch, work-block extraction, B42
pre-reserved blocks and exactly-sized two-digit output, B43 recursion,
B45 division, x86-64/Windows base-10^9 fallback, and the unchanged
84-case official benchmark. ARM64 Clang currently lowers the wide
128/64 constant division of each leaf work limb to __udivti3.
B46 replaces that division with an exact reciprocal calculation using
64x64->128 products and at most one quotient adjustment.

Let B=2^64 and d=10^18 (<B/2), and take R=floor(B^2/d), which equals
floor((B^2-1)/d) because d does not divide B^2. For each leaf
division step, N=h*B+l with 0<=h<d and 0<=l<B. Define
q0=floor(N*R/B^2). Since
0<=N/d-N*R/B^2 < d/B < 1,
q0 is either floor(N/d) or one less, never greater. Compute q0
using the high 128 bits of the product of two 128-bit values:
R=r_hi*B+r_lo, N=h*B+l,
q0=h*r_hi + floor((h*r_lo+l*r_hi+floor(l*r_lo/B))/B).
The cross term fits in unsigned 128 bits with this d. The exact
provisional remainder N-q0*d is in [0,2*d), strictly below B;
thus the 64-bit modular subtraction l-q0*d is already exact.
If provisional remainder>=d, increment q0 and subtract d once.
No floating-point estimate, external arithmetic dependency, new
public API, global cache, or ISA-specific assembly is introduced.

Before integration, a C++20 proof-of-implementation probe checked
2,643,208 independent hi/lo pairs against the compiler's exact
unsigned 128-bit quotient and remainder under Clang ASan/UBSan;
26,303 cases required the single correction. An independent
standalone full base-10^18 leaf formatter using the same recurrence
matched GMP mpz_get_str for 9,000 randomized/adversarial operands
at 9,10,11,12,15,16,17,32,64 limbs under Clang ASan/UBSan.
An ARM64 Clang 17 assembly probe compiled the new recurrence into
inline integer instructions without any __udivti3 call, whereas
the previous implementation emitted such a call. The assembler
probe does NOT establish Apple M2 runtime speedups.

Production tests add 352 GMP-oracle decimal boundary cases at
8,9,10,11,12,13,14,15,16,17,32 limbs, including all-ones, zero
runs and 10^18-related repeating words. All older test corpora and
benchmarks are preserved. Run standard Unix/Windows/Coverage CI and
both 20,000-iteration/31-sample B44 diagnostic and the official
Integer/GMP comparison on the same Apple M2 before accepting B46
as a performance improvement. A regression would require revisiting
or reverting the reciprocal path; no runtime speedup is claimed yet.


### B47 lower wide decimal parsing threshold and align input chunks (19 October 2026)

B46's Apple M2 official benchmark passed all 768-input/2688-operation
comparisons; construct_decimal at 1024 bits was 819.235 ns versus
GMP's 405.552 ns, giving a 2.020x worst-case slowdown, even though
to_decimal improved to 609.552 ns at 1024, 2273.581 ns at 2048,
and 6528.560 ns at 4096. The direct cause is visible in
BigUnsigned::from_decimal: the B24 unsigned-128 radix-10^18 parser
was enabled only when input.size()>309, while the benchmark's
2^1023+offset values have exactly 308 decimal digits. They therefore
used the B19 slower radix-10^9 two-half-products-per-limb path.
This behavior is architectural, not a benchmark anomaly.

B47 makes three isolated changes to the existing C++20 parsing path
on GCC/Clang compilers supporting unsigned __int128:
(1) select the 18-digit wide radix when the input has >=32 digits;
(2) pre-reserve ceil(input.size()/18) 64-bit limbs, a conservative
bound because 10^18<2^64, with overflow-safe size arithmetic;
(3) consume the first possibly short decimal group first, and then
use the compile-time constant multiplier 10^18 on every remaining
group. Carry accumulation stays exact in unsigned 128-bit, without
floating point, GMP, new dependencies, persistent caches, changes
to the public API or changes to the benchmark corpus. Inputs shorter
than 32 digits still use the existing radix-10^9 portable path.
MSVC and compilers without unsigned __int128 continue to use that
portable path at all lengths, preserving cross-platform semantics.
The ASCII validation is unchanged and precedes either implementation.

The isolated Clang x86-64 C++20 proxy (not an M2 or full Integer
benchmark) measured the 308-digit, 1024-bit input parser at about
573.9 ns for the old 9-digit path and 263.9 ns for the aligned,
reserved 18-digit path. At 128/256 bits the corresponding proxy
was 48.8/100.8 ns (old) versus 33.8/55.9 ns (new). These are
independent algorithm measurements and do not establish native
production speedups. A separate standalone C++20 implementation
of the exact aligned algorithm passed 8400 GMP-differential
comparisons under Clang ASan/UBSan, spanning 1..4096 decimal digits,
threshold boundaries, random and adversarial leading zeros, and
long zero runs. Production regression adds 181 independent GMP
parse comparisons: 26 decimal lengths x 6 scenarios and 5 bit
sizes x 5 official-shape values, including lengths 31/32/33 and
307/308/309/310. Preserve all historical decimal, parse, and
signed Integer suites; no benchmark source is changed.

As for all prior performance changes, run all Unix/Windows/Coverage
CI, the official 84-operation Integer/GMP comparison at
53,128,256,1024,2048,4096 bits, and the B44/B46 diagnostic with
20000 iterations and 31 samples before claiming a B47 speedup.
The primary target is construct_decimal at 1024 bits, not
to_decimal; verify no regression at short sizes and on ARM64,
x86-64, and portable toolchains. B47 integrated CI and native
Apple M2 measurements are pending.


### B48 ARM64 small-divisor fast path in shared divide_wide (19 October 2026)

B47's unchanged official 84-operation Apple M2 Integer/GMP
benchmark measured native construct_decimal at 1024 bits as
415.533 ns versus GMP 405.544 ns (1.025x), solving the longstanding
worst case. The overall total_time changed to 0.981x and geomean
to 0.755x; correctness passed 768 paired inputs and 2688 operation
results. The new worst case is quotient at 2048 bits, 271.283 ns
versus GMP 139.525 ns (1.944x), while remainder at 2048 bits was
293.979 ns versus GMP 153.952 ns (1.909x).
The original division-only benchmark pool has 32 fixed pairs:
16 with small, single-limb divisors, eight with half-width
divisors and eight near-equal magnitudes. The benchmark and its
operand pool must remain unchanged.

B48 optimizes the shared low-level divide_wide(high,low,divisor)
when the existing precondition high<divisor holds and
divisor<=2^32-1, on native ARM64 GCC/Clang only. Previously the
unsigned-128 numerator quotient/remainder may call the compiler's
__udivti3 helper; x86-64 already uses B26's hardware DIVQ.
Write N=high*2^64+low. Divide
U=high*2^32+(low>>32) by divisor to obtain the high 32 bits of
the quotient and remainder rU. Divide
L=rU*2^32+(low&0xffffffff) to obtain the low 32 quotient bits
and final remainder. The invariants high<divisor<2^32 and
rU<divisor imply U<2^64, L<2^64, and both quotient halves fit
32 bits, so both stages use exact native 64-bit division, never
a 128-bit arithmetic runtime helper. The high==0 fast path,
ARM64 larger-divisor existing implementation, x86-64 assembly
DIVQ, Windows MSVC _udiv128 and portable fallback are unchanged.
The improvement feeds all three joint, quotient-only and
remainder-only BigUnsigned division modes, including the ordinary
one-limb divisors used by Integer's official benchmark. Public
contracts, scratch, GMP independence and APIs are unchanged.

Standalone C++20 Clang ASan/UBSan checked 1,200,156 exact small
quotient/remainder pairs against independent unsigned __int128
division, including saturated quotient, high=0 and divisor 1.
A Clang 17 AArch64 assembly probe confirmed the two-stage narrow
implementation emits native UDIV/MSUB and no __udivti3 call;
the previous wide implementation emitted __udivti3. This is
assembly verification, not an Apple M2 timing measurement.
Integrated GMP regression adds 560 divrem_limb/into cases over
10 dividend sizes and 14 divisors (including the 2^32
boundary and all-ones), four operand patterns, and another
336 GMP-differential quotient/remainder/joint cases over eight
sizes, 14 divisors and three patterns. All original coverage
and the official 84-operation workload are preserved. Benchmark
on the same M2 before claiming a B48 speedup; check in particular
quotient/remainder at 1024, 2048 and 4096 bits, as well as
to_decimal to detect any unintended ARM64 regression.


### B49 ARM64 pre-inverted narrow quotient loop (20 October 2026)

Both B47 and B48 benchmarks were measured on the SAME Apple M4
ARM64, so they support a valid within-host comparison. B48 versus
B47 reduced native remainder timings at 1024/2048/4096 bits from
121.002/257.183/546.194 ns to 103.625/221.458/467.669 ns
(-14.36%/-13.89%/-14.38%). Quotient timings changed more modestly
from 110.192/236.183/489.544 to 108.296/226.177/484.854 ns
(-1.72%/-4.24%/-0.96%). On M4, the B48 official worst relevant
large-size case remains quotient 2048 bits at 226.177 ns versus
GMP 119.260 ns (1.896x). B48 84-operation geomean and total_time
were 0.758x and 1.021x, compared with B47 M4 0.766x and 1.019x;
the original 2.305x B47 worst compare 4096 timing appears noisy
and must not be ascribed to B48. Both benchmarks retained GMP
verification on 768 input pairs and 2688 operation results.
All B48 Unix/Windows/Coverage Actions passed; external Codecov
patch/project checks failed separately.

B49 targets quotient generation for single-limb divisors 2..2^32-1
on native ARM64 GCC/Clang. These account for 16 of the 32 official
division-only workload cases (small divisors 3..35), so the goal
is to accelerate the quotient-only path while RETAINING B48's
remainder-only small-divisor fast path. The existing shared
BigUnsigned::divrem_limb_into function used divide_wide per limb;
B48's divide_wide requires two ARM64 UDIV operations per 64-bit
limb for a 32-bit divisor. B49 selects a new exact pre-inverted
path when the dividend contains >=8 limbs, the divisor is between
2 and 2^32-1 inclusive, and native ARM64 unsigned 128-bit
multiplication is available. Outside that case the old loop,
MSVC, x86-64 and portable paths are unchanged. The optimization
is general for arbitrary 64-bit limb inputs, NOT specialized to
the benchmark's operand values.

Let d be the divisor, B=2^32, and R=floor(2^64/d). The reciprocal
is computed once using a 64-bit quotient/remainder operation:
R=floor((2^64-1)/d) + [((2^64-1)%d)==d-1].
For each 64-bit dividend word, split its high/low 32-bit halves;
each 64-bit partial numerator X is strictly below d*B<=2^64.
Set q0=floor(X*R/2^64) from the high half of a 64x64 product.
Since 0<=X/d-X*R/2^64<X/2^64<1, q0 is either the exact
floor(X/d) or one smaller. The provisional residual X-q0*d
is in [0,2*d), so one compare/subtraction (and increment of
q0) yields an exact 32-bit quotient half and remainder. Each
full quotient word combines both halves. The approach replaces
repeated hardware divisions with UMULH/MSUB-style operations
and requires no new public API, dependency or global table.
The B48 remainder-only path and shared long-division digit
estimator remain unchanged; the large-divisor paths are untouched.
The returned full remainder, output zero-filling, public failure
semantics and no-alias requirement remain unchanged.

Before integration an independent C++20 implementation completed
1,200,000 exact 128/32 quotient/remainder comparisons and 5000
multi-limb GMP mpz_tdiv_qr comparisons under Clang ASan/UBSan.
An AArch64 Clang 17 assembly probe for the core quotient step
showed native UMULH and MSUB, without division helper calls.
Those checks establish arithmetic correctness for the isolated
implementation, NOT performance on Apple M4. Integrated tests add
784 deterministic GMP-differential (input size, divisor, pattern)
cases across 7/8/9/16/32/64/96/128 limb lengths and divisors
spanning the 2^32 selector, checking owning quotient/remainder
and direct output with two zeroed surplus limbs. All existing
benchmarks, earlier GMP tests, B45 divrem and B46/B47 decimal
paths are retained unchanged.

Measure the existing official 84-operation benchmark on M4,
20,000 iterations and 31 samples, comparing B48/M4 with B49/M4.
Primary outcomes: quotient 1024/2048/4096 ns and GMP-relative
ratio. Guardrails: remainder improvements from B48 should remain,
total_time should not regress materially, and all platform CI
jobs should pass. A performance improvement is NOT claimed
until same-host measurements are supplied.


### B50 wider 19-digit input radix for construct_decimal (20 October 2026)

The B49 official 84-case benchmark is measured on the SAME Apple M4
ARM64 as B47 and B48. For B49, native construct_decimal at 2048 bits
was 1009.473 ns versus GMP 744.775 ns (1.355x); at 4096 bits,
2883.946 ns versus GMP 1726.542 ns (1.670x). The overall 84-case
native/GMP total_time was 1.016x and geomean 0.747x, with successful
verification of 768 input pairs and 2688 paired operation results.
Although B49 reduced native quotient at 2048 bits to 198.350 ns,
construct_decimal remains a large contributor to total time.

B50 is a narrowly scoped change to the wide 128-bit-supported
BigUnsigned::from_decimal path that B47 selected for >=32-digit
inputs. Each radix-10^18 block currently contains at most 18
decimal digits in one 64-bit limb. B50 uses the larger exact
radix 10^19=10000000000000000000, still strictly below 2^64,
to consume 19 digits per pass and reduce full multiprecision
Horner scans. The first possibly shorter group is still consumed
first, all subsequent groups use the constant radix, and the
pre-reservation changes from ceil(decimal_digits/18) to
ceil(decimal_digits/19). Each wide product has the form
digit*radix+carry, where digit<2^64, radix<=10^19<2^64 and
carry<radix. Hence the maximum product is
(2^64-1)*radix+(radix-1)=2^64*radix-1<2^128,
so the existing unsigned __int128 calculation stays exact.
No divide-and-conquer algorithm, caches, new public API,
dependency, GMP production call or benchmark alterations.
The historical base-10^9 portable path on Windows/MSVC and
other compilers without unsigned 128-bit remains unchanged;
the unsigned-wide path is still selected at >=32 digits.
ASCII validation and signed Integer construction are unchanged.

A preliminary isolated Clang 17 x86-64 C++20 comparison (NOT
an Apple M4 result) of otherwise identical Horner parsers
showed radix-10^19 close to radix-10^18 at 1024 bits, roughly
3-5% faster at 2048 bits, and roughly 1-3% faster at 4096
bits. A prototype balanced parser and a 38-digit two-limb
radix were slower in that environment and were NOT integrated.
An independent standalone test verified 12000 adversarial/random
decimal inputs against GMP and the prior radix-10^18 parser
under Clang AddressSanitizer and UndefinedBehaviorSanitizer,
including thresholds 31/32/33, exact radix boundaries,
leading zeros, saturated digit strings, 308 digits and
2048/4096-bit-sized inputs. Production regression adds 234
GMP-oracle parsing cases (39 lengths x 6 patterns) spanning
19-digit group boundaries, 1024-bit and 2048/4096-bit
construction shapes; previously existing tests remain unchanged.
The numbers above are input-parser proxy measurements only,
not independent Apple M4 evidence.

Acceptance requires the unchanged official 84-operation
benchmark on Apple M4 at 20000 iterations and 31 samples:
compare B49/M4 with B50/M4 construct_decimal at 128/256/
1024/2048/4096 bits, overall total_time/geomean, and
other operations for regressions. All standard CI matrices
must pass. Keep B49 as the rollback baseline if the 19-digit
wide path does not show a repeatable same-host benefit.


### B51 first x86-64 performance campaign: revisit B46 reciprocal leaves (20 October 2026)

The user provided the baseline *x86-64* 84-operation Integer/GMP
benchmark (20,000 iterations, 31 samples, numerical verification
PASS on 768 input pairs and 2688 operation results). At 4096 bits
native versus GMP were: construct_decimal 2637.401 versus 2347.551
ns, mul 3131.380 versus 1409.130 ns, square 1823.646 versus
960.669 ns, quotient 1215.565 versus 361.724 ns, remainder
1175.650 versus 365.284 ns, pow3 8192.439 versus 3650.179
ns, and to_decimal 11715.378 versus 4431.388 ns. All-operation
geomean was 1.459x and total_time 2.052x. This x86-64 hardware
is not identified by model/compiler flags in the submitted log;
M4 ARM64 wall-clock timings cannot establish same-CPU benefits.
The M4 B50 total_time was 1.015x, not a legitimate cross-machine
speedup baseline for individual changes.

The code review found 1) x86-64 production multiply dispatch
depends on compiler: GCC enables the fused kernel and associated
large balanced Karatsuba and symmetric-square specializations,
while Clang generally chooses the B3 row-schoolbook path;
2) portable and x86-64 division kernels already include x86
DIVQ and several B28/B36 quotient fast paths; and 3) x86-64
decimal leaves remained on radix-10^9 in B50. Historical B25/B30
radix-10^18 experiments REGRESSED in the official x86-64
benchmark. B30 used a real 128/64 DIVQ loop, which is not the
same algorithm as B46's exact pre-inverted radix-10^18, so
these old measurements do not settle the current hypothesis.

B51 changes only the existing compile-time guards around
to_decimal_blocks_wide and its 9+-limb selection in
src/big_unsigned.cpp: native x86-64 GNU/Clang targets providing
unsigned __int128 now reuse the *unchanged B46* quotient/remainder
reciprocal and B42 18-digit pair formatter, just as ARM64 does.
The B46 quotient uses 64x64->128 products and at most one
remainder correction; it requires neither a runtime 128-bit
divide nor an x86-specific instruction. The full B19 base-10^9
fallback remains selected for <=8-limb magnitudes, MSVC and
other architectures; ARM64 behavior is unchanged. The shared
recursive split, GMP independence, 84-case official workloads,
public APIs, x86-64 multiply/divide kernels and their dispatch
remain unchanged. No arbitrary handwritten assembly is
introduced before the available competing mulq/adcq, Comba,
fused and Karatsuba kernels are benchmarked on the actual
target host. This is a **controlled first x86-64 experiment**.

Before production integration, a local standalone x86-64 Linux
probe on Intel Xeon Platinum 8573C (not the user's target CPU)
reconstructed the B31 10^9 and B42/B46 10^18 formatting
algorithms and checked each result with independent GMP
mpz_get_str for 256 dense/sparse/saturated operands at each
of 8,9,12,16,32 limbs (1280 GMP comparisons per variant).
Clang 17 -O3 baseline-x86 medians (ns):
- 8 limbs: narrow 482.030, reciprocal wide 304.750;
- 9 limbs: narrow 548.132, reciprocal wide 336.512;
- 12 limbs: narrow 956.700, reciprocal wide 524.170;
- 16 limbs: narrow 1697.650, reciprocal wide 847.686;
- 32 limbs: narrow 7771.950, reciprocal wide 3069.380.
GCC 14 -O3 baseline-x86 medians (ns):
- 8 limbs: narrow 551.198, reciprocal wide 346.259;
- 9 limbs: narrow 645.396, reciprocal wide 366.258;
- 12 limbs: narrow 1093.210, reciprocal wide 553.028;
- 16 limbs: narrow 1998.440, reciprocal wide 797.383;
- 32 limbs: narrow 8642.830, reciprocal wide 2632.700.
These are **isolated local CPU timings**, not B51 official
results and not a same-host comparison with the uploaded
x86-64 baseline. Clang 17 AddressSanitizer/UndefinedBehavior
Sanitizer and LeakSanitizer also passed the exact paired
formatter checks against GMP. The production test suite adds
288 deterministic GMP-oracle decimal values at 7/8/9, 16/17
and larger size boundaries, preserving all earlier tests.
Run CI across GCC, Clang, AppleClang, Windows/MSVC and
ThreadSanitizer, and compare the official x86-64 benchmark
on the user's machine; revert the x86 selector if the
integrated experiment repeats B30's earlier regression.

**Campaign next step**: on the same x86-64 CPU/compiler,
run benchmark_big_unsigned_multiply with --include-fused,
--include-x86-64-carry, --include-comba, --include-schoolbook
and --include-karatsuba, plus --karatsuba-cutoff 32 to profile
native selected, alternative preallocated kernels and GMP.
Use diagnostic data to decide whether MULQ/ADCQ, BMI2/ADX
gated kernels, square specialization or Karatsuba thresholds
are justified for that target, rather than assuming M4 ARM64
instruction-level tuning transfers directly.


### B52 x86-64 GCC 32-limb assembly dispatch informed by native benchmarks (20 October 2026)

The user's official B51 x86-64 benchmark (20,000 iterations,
31 samples, unchanged 84 cases) demonstrated that enabling the
B46 exact reciprocal base-10^18 leaf under x86-64 significantly
accelerated Integer::to_decimal versus B50 on the SAME host:
1024 bits 1793.826->720.012 ns (-59.86%); 2048 bits
4526.997->3280.192 ns (-27.54%); 4096 bits
11715.378->9396.203 ns (-19.80%). Overall geomean improved
1.459x->1.370x, total_time improved 2.052x->1.842x,
and worst slowdown improved 3.400x->3.257x. GMP verification
passed 768 paired input fixtures and 2688 result fixtures.
This CONFIRMS that the B46 reciprocal implementation, unlike
the historical pre-B46 DIVQ implementation rejected in B30.1,
can beat the old x86 radix-10^9 output path in production.

The user also ran the existing benchmark_big_unsigned_multiply on
that SAME x86-64 host, requesting all experimental variants,
20,000 iterations and 31 samples with Karatsuba cutoff 32.
The selected production kernel was reported as fused (x86-64
GCC). For preallocated dense general products:
bits         fused       MULQ/ADCQ      relative assembly delta
128          16.409      16.709         +1.83%
256          22.367      23.956         +7.10%
1024         205.249     220.956        +7.66%
2048         954.267     810.044        -15.11%
4096         2964.676    3153.710       +6.38%
At 2048 bits the historical portable and Comba alternatives
were much worse; Karatsuba cutoff 32 was 961.873 ns, close
to fused 954.267. At 4096 bits Karatsuba cutoff 32 was
2994.552 ns, again near fused; changing every general product
to MULQ/ADCQ would REGRESS the observed 4096-bit result.
The x86-64 hardware CPU model and target compiler flags are not
supplied in the log: do not extrapolate these absolute values
to a different CPU, and do not claim BMI2/ADX availability.

B52 selects the already GMP-tested no-extension MULQ/ADDQ/ADCQ
kernel *only* in BigUnsigned::mul_low_into on builds with
ARIADNE_BIG_UNSIGNED_GCC_FUSED (native GNU x86-64) and
a.size()==b.size()==32 limbs and output.size()>=64 limbs.
Other sizes, every truncated output below 64 limbs and
Clang/AppleClang/MSVC/ARM64 all retain their previous selectors.
The existing assembly has been benchmarked on the user's host
and exercises the same full-product recurrence; no new assembly,
CPU feature probing, public API or benchmark code is introduced.
No value-dependent benchmark specialization exists: the 32-limb
dispatch is valid for arbitrary limb values and applies equally
to standalone owning products and Karatsuba base cases. The
owning GCC Karatsuba 64-limb path uses 32x32 z0/z2 leaves,
so it may also benefit; its 33x33 cross term remains fused.
Only a new on-target integrated run can confirm either gain.

The independent GMP test additions exercise 9 pairs of lengths
31/32/33 x 31/32/33 limbs and 32 randomized/adversarial
patterns each (288 distinct source pairs), comparing the
owning output and exact caller-owned result against GMP,
the selected kernel against the original fused kernel
word for word, oversized destination zero-filling, and
truncated low-product arithmetic against GMP. The earlier
source-level and integrated GMP coverage, official Integer
benchmarks and multiplication diagnostic workloads remain
unchanged. Run all CI configurations; compare B51/B52
official benchmarks on the same x86 host, especially
2048/4096-bit mul/pow3 and total_time/geomean, and rerun
the existing isolated kernel benchmark: at 2048 bits
prealloc_selected should now track prealloc_x86_64_carry
rather than prealloc_fused, while selected at 1024/4096
should stay close to fused. If this prediction fails,
inspect build flags / path before adjusting dispatch.

Subsequent campaign stages should profile the existing DIVQ
quotient estimator, reciprocal based small-divisor paths,
and optional BMI2/ADX kernels only after checking actual
target ISA and paired benchmark results; do not silently
expand the narrow 32-limb selector without evidence.


### B53 port ARM64 B38/B40 sparse-row multiplication to x86-64 (20 October 2026)

On the user's SAME x86-64 host, the B52 native multiplication
diagnostic confirmed the 32x32-limb MULQ/ADDQ/ADCQ dispatch:
2048-bit preallocated selected 953.096 ns (B51) to 853.847
ns (B52), 10.41% faster; at 4096 bits cutoff-32 Karatsuba
2994.552->2675.843 ns, 10.64% faster. In the official
84-operation benchmark, 2048-bit Integer mul fell
1030.988->861.190 ns (-16.47%), and 4096-bit Integer
mul 3117.157->2845.357 ns (-8.72%). At 4096 bits
pow3 remained approximately unchanged
8363.683->8362.002 ns, which exposes the cost of
executing full schoolbook rows for sparse powers.
B52 official total_time/GMP 1.843x, geomean 1.435x
and numerical verification passed (768 paired inputs,
2688 operation results). Across B51/B52 GMP and
microbenchmark medians move as well; do not confuse
cross-run ratio noise with code regressions.

B53 transfers the *general mathematical* optimizations
of ARM64 B38/B40 to native x86-64 multiplication, without
special-casing the benchmark's pow3 operands. The x86
GCC selected fused unsigned __int128 row-schoolbook
now has two templated implementations: sparse skips any
exact-zero outer limb row after initially zero-filling
the output; dense has no zero check inside the outer
loop. The existing x86 GNU/Clang MULQ/ADDQ/ADCQ kernel
gets the same two specializations, preserving instruction
syntax and carry propagation. If a row multiplier ai==0,
every partial product, carry and output update from that
row is identically zero, including the final carry slot
which earlier rows cannot touch. Thus skipping a row is
exact for all full/truncated and noncanonical raw inputs.

As in ARM64 B40, when active outer rows =
min(a.size(),output.size()) is <64 or b.size()<64,
select the sparse variant without scanning; larger
cases first scan up to a threshold
max(2,floor(active_rows/8)) exact zero limbs and
select sparse only upon reaching it, otherwise dense.
The prescan exits early when sparse. This transfers
ARM64's protection against costly dense-data branches
at large sizes. The GCC 32x32 B52 assembly selector
remains unchanged; it now benefits from zero-row
elision too. The owning GCC Karatsuba size/balance
selector is unchanged (a separate candidate for B54),
and existing x86 symmetric square remains unchanged.
Clang's default B3 schoolbook and MSVC are unchanged;
only the explicit x86 carry benchmark on Clang uses
the new template selection. ARM64 and portable
kernels and all public exception, view ownership,
no-alias, exactness and scratch contracts are unchanged.
Neither BMI2 nor ADX is required.

Before integration, a standalone GMP-based exact
128-bit/assembly test of these identical recurrences
verified 8000 full/truncated/oversized comparisons
under GCC 14 and Clang 17 on a local AMD EPYC proxy,
and again 8000 comparisons under Clang
AddressSanitizer/UndefinedBehaviorSanitizer. The
proxy is NOT the user's x86 target and cannot
establish a performance claim. Production test
corpus adds 16 operand-length pairs x24 adversarial/
random patterns = 384 source pairs. Each compares
six destination sizes (0,1,half,full-1,full,full+2)
against independent GMP modulo 2^(64*length), with
byte-exact comparison of the selected, forced fused
and forced x86 carry results. Coverage spans the
B52 31/32/33-limb selector, B40 63/64/65 adaptive
prescan, 128-limb unbalanced operands, saturated,
zero, alternating, sparse and dense patterns.

Run the unchanged official x86-64 Integer/GMP benchmark
(20000 iterations, 31 samples) and the existing
benchmark_big_unsigned_multiply alternatives; primary
target is reducing pow3 at 1024/2048/4096 bits and
preserving B52's 2048/4096 mul gains. Inspect dense
1024/2048/4096 preallocated multiply measurements
for regressions and avoid treating Apple M4 ARM64
absolute timings as same-host comparisons. B53 has
no on-target performance measurements yet. All
standard CI checks should pass before acceptance.


### B53.1 fix doubled x86-64 inline-assembly escape sequences (20 October 2026)

The user's native x86-64 GCC -O3 build failed in
src/big_unsigned_x86-64.cpp with assembler
"invalid character '\\' in operand 1" near line 200.
Root cause: B53 generated the five consecutive
MULQ/ADDQ/ADCQ assembly string fragments with C++ source
spelling `\\\\n\\\\t` (double backslashes), which
produces literal backslash characters instead of
the required newline/tab separators. B53.1 changes
exactly these five source occurrences to `\\n\\t`.
No instruction, register constraint, algorithm,
dispatch, public API, test corpus, benchmark, or
platform logic is modified. The fixed sequence was
compiled locally on x86-64 with GCC 14.2 and Clang 17
(-O3 -std=gnu++20 -Wall -Wextra -Wpedantic -Werror)
and each compiler passed 250000 exact products
against unsigned-128 arithmetic. This isolated
assembly smoke check does not replace full project
CI or the user's native B53 benchmark. B53 official
same-host measurements remain pending.


### B54 x86-64 GCC adapt square and Karatsuba to operand density; retune 32-limb dispatch (20 October 2026)

The user's B53.1 native x86-64 benchmark (20,000 iterations,
31 samples; independent GMP verification PASS for 768 paired
inputs and 2688 results) confirms the dominant impact of
B38/B40-style exact-zero row skipping. Versus B52 on the
same x86-64 CPU, 1024/2048/4096-bit Integer::mul changed
252.938/861.190/2845.357 to 79.513/114.049/701.262 ns;
pow3 changed 738.724/2193.092/8362.002 to
310.053/686.622/2176.326 ns. Full benchmark total_time/GMP
fell 1.843x to 1.361x; geomean 1.435x to 1.284x. Dense
preallocated selected multiplication was essentially stable
at 1024/2048/4096 bits (204.242/853.847/2969.094 ns B52
versus 205.667/861.181/2934.027 ns B53). The B53.1
build escape fix restored correct MULQ/ADCQ asm newlines;
all seven GitHub Actions eventually passed (external
Codecov patch/project failures reported independently).

B54 addresses three remaining architecture-general
opportunities, without operand-value-specific
benchmark special cases:

1. **Retune the B52 32x32-limb dispatch**: on the SAME
   B53 x86 machine and compilation, explicit dense
   2048-bit `prealloc_fused` measured 801.142 ns,
   whereas forced `prealloc_x86_64_carry` measured
   858.091 ns and selected 861.181 ns. These new
   B53 fused measurements supersede B52's earlier
   preference for assembly (B52 fused 1008.084 ns
   versus assembly 854.808 ns). B54 removes the
   special-case 32x32 selected assembly routing,
   using B53's exact sparse-aware fused recurrence
   at all GCC input dimensions. The explicit ASM
   entrypoint remains available for comparison and
   optional testing. Sparse 32x32 performance is
   NOT YET established on the user's machine: only
   B54's official 84-case benchmark can validate
   that changing the 32x32 selection does not undo
   sparse gains.

2. **Transfer ARM64 B38/B40 adaptive symmetric square**:
   the existing B32 GCC x86-64 square forms the upper
   off-diagonal triangle and doubles it, then adds
   the diagonal products. B54 introduces sparse/dense
   compile-time instantiations. If an outer limb is
   zero, omit its entire off-diagonal row. Earlier
   triangular rows never write the row's final
   output slot `output[i+n]`, and the output is
   initialized to zero. All diagonal products and
   carry passes are preserved. For <=64 limbs,
   select the sparse variant without a prescan;
   for >64 limbs classify the first n-1 active
   outer rows with B53's exact-zero threshold
   max(2,floor((n-1)/8)) and select dense without
   per-row checks when the threshold is not reached.
   The public square() size >=8 dispatch, original
   arithmetic, allocation and output contracts remain.
   This should especially reduce cost for sparse
   1024/2048/4096-bit squares, previously roughly
   166.670/492.103/1821.579 ns on x86 B53.

3. **Transfer ARM64 B39 owning Karatsuba density gate**:
   the B35 GCC x86 owning mul() used cutoff-32
   Karatsuba for all >=64-limb nearly balanced
   products, including sparse powers where B53's
   zero-row schoolbook is now cheap. B54 requires
   both operands to have <=floor(size/8) exact-zero
   limbs before choosing this same balanced
   Karatsuba path. Otherwise fall through to the
   current sparse-aware row-schoolbook selected
   kernel. Inputs of insufficient size or balance
   bypass density scanning as before; no new scratch
   allocation occurs for rejected Karatsuba.
   Dense products retain cutoff-32 Karatsuba.
   The same size/balance and zero-density criteria
   were independently used on ARM64 in B39.

This is a *controlled new B54 variant*, not a measured
native on-target speedup claim. A standalone x86-64
square kernel reproducer compiled with GCC and with
Clang ASan/UBSan verified 13312 exact GMP squares
and equality to B32's original algorithm per build
across 8..128 limbs and saturated, zero, sparse,
dense and threshold-boundary inputs. Production test
addition supplies 17 length pairs x32 patterns =
544 general multiplication cases (both operand orders
and forced fused reference), plus 13 square sizes
x64 patterns =832 GMP square/self-product pairs:
1376 additional randomized/adversarial fixtures
across B52/B53 32x32 and B35/B39 63/64/65
crossover boundaries. All earlier test cases are
preserved; the official benchmark and its 84 original
workloads remain byte-for-byte unchanged.

Acceptance: run full CI including GCC, Clang,
AppleClang, Windows/MSVC, Coverage and sanitizer.
Repeat the unchanged 31-sample, 20,000-iteration
official x86-64 benchmark and the explicit
preallocated dense multiplication diagnostic with
fused/carry/Karatsuba alternatives. Compare B53/B54
on the same host for 2048-bit sparse mul,
1024/2048/4096-bit square, 4096-bit mul/pow3,
dense 2048/4096-bit multiplication and the
full geomean/total_time. If B54 harms sparse 2048
multiplication or dense 4096 Karatsuba, revert or
retune only the responsible selection without
changing any benchmark input or test expectations.


### B55 x86-64 GMP-proved MULQ/ADCQ multiply-subtract division step (20 October 2026)

The user's B54 84-operation x86-64 benchmark (20,000 iterations,
31 samples, 768 paired input fixtures / 2688 operation results
PASS) reduced the native/GMP total_time from B53's 1.361x to
1.153x, geomean 1.284x to 1.181x, by transferring ARM64 sparse
square and Karatsuba density selection. The remaining worst case
was quotient at 2048 bits: 558.359 ns native vs 163.808 ns GMP,
3.409x; at 4096 bits quotient 1209.668 vs 361.524 ns (3.346x),
remainder 1120.816 vs 364.121 ns (3.078x).
These workloads include both narrow single-limb divisors and
middle-size multiprecision denominators: do not assume their
aggregated median isolates any one subpath.

B55 examined the first hypothesis, directly copying B49's
preinverted narrow divisor strategy from ARM64 to x86-64.
A standalone x86-64 AMD EPYC 9V74 GCC 14.2 and Clang 17
proxy compared exact scalar DIVQ 128/64 division with
B49's two reciprocal multiply-high steps per 64-bit quotient
digit, and full remainder-only equivalents, all using the
unchanged divisor-independent input recurrence. All 101376
deterministic/random/adversarial GMP quotients and remainders
passed per compiler, for 1..128-limb dividends and 18 divisors
spanning 2, powers of two, small primes and the 2^32 boundary.
On this proxy B49 reciprocal was *consistently slower*:
GCC 32 limbs divisor=35 gave 119.167 ns for DIVQ and
181.581 ns for reciprocal quotient; the rem-only path
127.614 versus 206.651 ns. Clang 32 limbs divisor=35
gave 119.395 versus 221.048 ns for quotient. The exact
CPU differs from the user's unknown x86 model, but the
results disprove automatic portability and justify retaining
B26's existing native DIVQ for narrow divisors on x86.
The ARM64 B49 code and all related selection stay UNCHANGED.

The second hypothesis targets the genuine multi-limb divisor
path, which normalizes the dividend and denominator, estimates
each quotient digit with divide_wide, and subtracts
digit*divisor from an m+1 limb window. B27 already fused
the multiply and subtract in C++, but on x86-64 its
128-bit wide arithmetic still incurs separate borrow
comparisons. B55 adds a baseline native
MULQ -> ADDQ/ADCQ -> SUBQ/ADCQ *single carry chain*
for the subtraction recurrence, selected only when
m>=8 limbs on GNU/Clang x86-64, after all output and
workspace contracts have been validated. Each step
computes high:low = divisor[i]*digit+carry and then
window[i]-=low while incorporating the borrow into
the next carry. The exact algebra is unchanged.
SparseZero==true retains B29's handling of normalized
zero divisor limbs by subtracting only incoming carry.
Final high-limb subtraction and underflow flag are
identical. No BMI2, ADX, new runtime CPU gating,
API, heap allocation, GMP production usage or
benchmark modifications; ARM64, MSVC, non-x86
targets and all denominators with 1..7 limbs
retain the original C++ code.

A local AMD EPYC proxy of the exact native recurrence
was compiled with -O3 -std=gnu++20 -Wpedantic -Werror
using GCC 14.2 and Clang 17, plus Clang ASan/UBSan.
For each of 1,2,4,8,16,32,64,128 limbs it verified
4096 random/adversarial cases against independent
GMP integer subtraction modulo 2^(64*(m+1)) and
the negative-result flag, and required byte-exact
equality against the original B27 fused recurrence:
32768 successful GMP checks for each build. In the
GCC local timing proxy, the 16-limb kernel fell
48.10 -> 43.17 ns (~10.2%) and the 64-limb kernel
78.99 -> 71.40 ns (~9.6%); Clang 16-limb
19.31 -> 16.44 ns (~14.9%). Small 1/2 limb cases
were similar/slower, motivating B55's >=8-limb
selector. These are isolated local proxy timings,
NOT evidence of the integrated user-host speedup.

The production regression adds 6 denominator sizes
(7,8,9,16,32,64) x 3 additional dividend lengths
(+2,+9,+32) x 24 adversarial/random patterns =
432 GMP-oracle numerator/denominator pairs, validating
joint divrem, quotient-only, remainder-only, preallocated
outputs, oversized zero-filled destination words and
scratch guard sentinels. This augments rather than
weakens existing B27/B29/B36/B45 and B48/B49 tests.
No adjustments to the official 84-case workload
or benchmark expectation. CI must pass on GCC,
Clang, macOS AppleClang, Windows MSVC and all
sanitizer/coverage configurations. Run the original
official user-host x86 benchmark at 20k iterations
and 31 samples, measuring quotient/remainder at
1024/2048/4096 bits, aggregate total_time/geomean
and verifying all operation results. Preserve B54
as the rollback baseline if division gets slower
or no repeatable native improvement appears.


### B56 category-resolved x86-64 division profiling (20 October 2026)

A second on-target x86-64 B55 official 84-case run (20,000
iterations, 31 samples, GMP verification PASS for 768 paired
input fixtures and 2688 results) reproduced the first:
- 1024-bit quotient 284.508 -> 283.950 ns (-0.20%);
- 2048-bit quotient 530.541 -> 531.040 ns (+0.09%);
- 4096-bit quotient 1103.046 -> 1106.378 ns (+0.30%);
- 2048-bit remainder 537.519 -> 540.614 ns (+0.58%);
- 4096-bit remainder 1114.359 -> 1114.473 ns (+0.01%);
- 4096-bit to_decimal 8672.970 -> 8703.651 ns (+0.35%).
Total_time/GMP 1.117x -> 1.115x and geomean
1.150x -> 1.148x. The B55 benefits were repeatable
at the granularity of these measurements; B55 remains
the accepted production arithmetic baseline.

Before further instruction-level changes, distinguish
which *algorithmic division category* dominates. The
unchanged official workload mixes small single-limb
divisors, half-width multiprecision denominators and
near-equal operands. Its overall quotient/remainder
medians cannot isolate the expensive subpath, nor can
their medians be decomposed into additive stage costs.
B56 intentionally makes NO PRODUCTION ARITHMETIC CHANGE.

Add a separate optional executable,
benchmarks/benchmark_big_unsigned_division_profile.cpp,
wired into benchmarks/CMakeLists.txt, with:
- exact deterministic 32-pair GMP-verified pools for
  each (size,category), default 1024/2048/4096 bits;
- categories limb_small, limb_wide, short_2, short_8,
  quarter_dense, half_dense, half_sparse, long_dense,
  near_less_q0, near_greater_q1, shifted_half and
  power_of_two;
- separately timed owning BigUnsigned quot/rem/divrem
  versus matching GMP mpz_tdiv_q/r/qr with fresh
  owning outputs created and released inside timed loops;
- separately timed caller-owned BigUnsigned::divrem_into
  versus GMP mpz_tdiv_qr with preallocated reusable
  output buffers; result paths are labelled nonadditive;
- GMP oracles for each public result and the preallocated
  result, before and after the timed loops. The official
  benchmark and all production tests are unchanged;
  user-specified --bits, --iterations and --samples
  use the established CLI parser (whole 64-bit multiples
  at >=512 bits; default 20k iterations/31 samples);
- a by-operation report of worst native/GMP case so
  dense long division can be distinguished from a
  one-limb DIVQ or B36 near-equal fast path.

This is **diagnostic instrumentation only**, with no
numerical or speedup claims before the executable
passes CI and is run on the user's same x86-64 host.
A recommended initial diagnostic run is:
  cmake --build build --target benchmark_big_unsigned_division_profile -j8
  ./build/benchmarks/benchmark_big_unsigned_division_profile -t dark \
    --bits 1024,2048,4096 --iterations 4000 --samples 15
Repeat slowest categories at 20000 iterations/31 samples
before choosing B57. Inspect whether multiprecision
quotient/remainder paths still dominate, or if overhead
for small/near-equal denominators is more consequential.
Avoid mistaking preallocated and owning cost differences
for arithmetic speedups. Do not change oracle coverage
or adjust benchmark inputs to make a faster path win.


### B56.1 fix move-only GMP fixture pool construction (20 October 2026)

Before native testing, isolated GCC 14.2 C++20 compilation
found that the newly introduced B56 diagnostic
`make_pool()` returned `std::array<Fixture,32>` by
value even though Fixture contains the explicitly
noncopyable/nonmovable Gmp RAII wrapper. Named-value
return optimization is optional and the compiler
must find a legal move constructor: GCC rejected it
with "use of deleted function std::array<Fixture,32>::array".
B56.1 adds noexcept move constructor and assignment
to the diagnostic-only Gmp wrapper, using mpz_init
and mpz_swap rather than shallow copying opaque
GMP storage. Both live objects remain valid for
independent destruction; self move-assignment is
guarded. A standalone exact reproducer compiled
with GCC 14.2 and Clang 17 under -O3 -std=gnu++20
-Wall -Wextra -Wpedantic -Werror, linked GMP,
returned a 32-element array by value, moved it
again and checked the expected GMP integer; both
executables passed. There is no modification to
any production library code, test source, official
benchmark workload or diagnostic measurement logic.
Run the complete benchmark and CI to establish
integration and native timing; that has NOT yet
been confirmed.


### B56.2 ensure exhaustive shape switch has required default (20 October 2026)

The B56.1 diagnostic CI recorded GCC 13 Release and
AppleClang 21 Release build failures, while Clang 16
Release/TSan passed. CI enables WARNINGS_AS_ERRORS=ON.
Review of the new diagnostic found a 12-case switch
that intentionally enumerated every Shape but omitted
the explicit default label required by the project's
-Wswitch-default policy. B56.2 moves the existing
logic_error into a switch default case without changing
which inputs or cases can execute. An isolated C++20
reproducer confirmed -Wswitch-default -Werror rejects
the old fully-enumerated switch on both GCC 14 and
Clang 17. This is a warning-policy compliance fix,
not a change to the arithmetic, oracles, measured
operations, input fixtures or official benchmark.
The exact CI error annotations were not accessible
through the GitHub connector; confirmation requires
new CI results. No user-host B56 timings have
been recorded.


### B57 x86-64 long single-limb preinverse and universal power-of-two division (20 October 2026)

The completed B56 user-host diagnostic, with 32 GMP-verified
pairs per category at 1024/2048/4096 bits, reported severe
slowness in single-limb division at 4096 bits: limb_small
quot 1281.487 vs GMP 200.639 ns (6.387x), limb_wide
quot 1349.135 vs 196.434 ns (6.868x), and limb_wide
rem 1322.548 vs 196.262 ns (6.739x). The preallocated
comparison likewise showed high arithmetic costs.
GMP verification passed. The unchanged official
Integer benchmark's nontrivial quotient/remainder
corpus has 16/32 one-limb divisors (50% of inputs),
8/32 ~half-width divisors and 8/32 near-equal
operands, as verified in
benchmarks/benchmark_exact_arithmetic_workload.hpp.
B56 CLI and production arithmetic remain intact.

B57 evaluated a **different 64/64 normalized pi1
reciprocal** from ARM64 B49's two 32-bit stages.
For B=2^64 and normalized d in [B/2,B), compute
inv=floor((B^2-1)/d)-B. This fits 64 bits. Its
exact construction is a SINGLE B26 DIVQ computing
(hi=(~d),lo=UINT64_MAX)/d; the high limb is
strictly <d, so DIVQ cannot overflow. Normalize
the entire quotient stream virtually, without
allocations: shift the divisor by countl_zero,
start the remainder at the top discarded bits,
and produce each normalized input limb from
adjacent input limbs. The pi1 recurrence uses
one 64x64->128 product plus modulo-2^64
subtraction/correction to give the exact digit
and next normalized remainder, then shifts the
final remainder down. The result is mathematically
identical to division using DIVQ. GMP's published
preinversion formulation provides useful
independent methodological context; the new
recurrence is implemented in GMP-free production
using native C++ unsigned __int128 products.

A standalone Linux AMD EPYC x86-64 reproduction
checked 32768 quotient/remainder input pairs
independently against mpz_tdiv_qr, including
1..128-limb inputs, small and full-width
divisors, all-ones, zero runs and powers of B,
with -O3 -Wpedantic -Werror GCC 14 and Clang 17.
The revised single-DIVQ inverse setup was also
verified with both compilers over the same cases.
The median proxy times varied by operand size
and compiler: in representative GCC runs at
64 limbs divisor=35, DIVQ quotient ~272 ns,
normalized reciprocal quotient ~229 ns;
at 32 limbs the reciprocal was not better,
and for remainder-only it was often worse.
Clang's 32/64-limb normalized reciprocal was
generally not faster than DIVQ. These tests run
on a proxy, not the user's x86 CPU, and are NOT
integrated B57 performance claims. Therefore
the reciprocal is selected **only** for >=64
raw quotient limbs on GNU x86-64 GCC, not Clang,
not MSVC, not ARM64 and not remainder-only
(which remains on B26 DIVQ); existing AMD
and other x86 microarchitectures may differ.
Small quotient loops retain B26 and ARM64 B49
paths as appropriate. No BMI2/ADX required.

B57 also recognizes exact *power-of-two divisors*
and performs the quotient/remainder by shift and
mask. BigUnsigned::divrem_limb_into uses the
fast exact 2^s path, including 2^0=1; the
multi-limb divrem_selected_into handles
b=B^k*2^s and all three compile-time
result modes. Its fast path is after all
existing error, output-size and scratch-size
checks, so it maintains their fail-before-
mutation contracts; the results are fully
zero-filled as before. For arbitrary dividends
including a<b or noncanonical raw limb inputs
on the single-limb API, the exact mathematical
identity holds. The selected multi-limb
algorithm does not allocate intermediate buffers
or change public exception behavior; ARM64,
Windows/MSVC and other architectures can
benefit from this ISA-independent power-of-two
optimization. The B45 trailing zero-limb and
B36 near-equal optimizations continue to handle
all non-power-of-two inputs unchanged.

A standalone 6534-case GMP differential test
for multi-limb powers of two passed both GCC 14
and Clang 17 under ASan/UBSan, covering
1..128-limb dividends, 1..128-limb divisors,
s=0/1/17/31/32/63, dense/sparse/saturated
values and n<m. The production suite adds
1440 raw-limb single-divisor pairs across
0/1/7/8/16/31/32/63/64/65/96/128
limbs, 15 different divisors and eight
randomized/adversarial patterns, plus 2160
power-of-two multi-limb cases across nine
dividend sizes, eight divisor sizes, five
bit positions and six patterns. Each
independently verifies GMP quotient and
remainder and checks extra result words,
every output mode and scratch sentinels.
The official 84-case benchmark and the B56
diagnostic benchmark are UNCHANGED.

The main same-host acceptance test is to repeat
the original 84-case benchmark at 20k iterations
and 31 samples and the existing B56 profile at
4k iterations and 15 samples, checking separate
limb_small, limb_wide, power_of_two and
short_2 categories at 1024/2048/4096 bits
for wins and regressions. Since GCC 64-limb
reciprocal behavior is experimental on the
user's unknown x86 CPU, retain B56/B55 as
rollback baseline if results regress.


### B57.1 restore warnings-as-errors clean bit-count conversions (20 October 2026)

The user's native x86-64 GCC build of B57 succeeded but warned at
three new conversion sites in src/big_unsigned.cpp, lines 502, 525
and 819: std::countr_zero/std::countl_zero return int whereas the
destination variable is unsigned. All seven B57 GitHub Actions
checks failed at their build steps with WARNINGS_AS_ERRORS=ON,
including GCC, Clang, AppleClang and Windows. B57.1 makes these
three conversions explicit via static_cast<unsigned>; the
mathematical count range 0..64 guarantees they are nonnegative.
GCC and Clang compiled isolated C++20 -Wsign-conversion -Werror
reproducers. No algorithm, thresholds, test corpus or APIs change.

The user's first on-target B57 official benchmark (GMP verification
PASS) demonstrated a performance tradeoff. At 4096 bits native
quotient fell from B55 1106.378 ns to 707.983 ns (-36.0%), while
4096-bit construct_decimal rose from 2675.133 to 3175.472 ns;
total_time/GMP rose 1.115x to 1.149x and geomean 1.148x to
1.163x. Re-run under controlled conditions before attributing
unrelated decimal construction timing changes to B57.

B56/B57 division profile at 4096 bits: limb_small quot
1281.487->314.578 ns (-75.5%), limb_wide quot
1349.135->312.612 ns (-76.8%), power_of_two rem
749.174->87.441 ns (-88.3%). Generic single-limb remainder
limb_wide rem 1322.548->1306.280 ns, nearly unchanged,
and short_2 quot 1665.143->1706.745 ns (+2.5%).
The 84-case official benchmark and existing tests remain unchanged.
Full CI retest required.

### B57.2 fix AppleClang B57 single-limb test shadow warning (20 October 2026)

AppleClang arm64 Release CI failed to compile tests/test_big_unsigned.cpp
with -Wshadow -Werror: the B57 test loop introduced a local
`remainder` while test_division_limb() already held the structured-
binding variable `remainder` from the fixed radix-boundary case.
B57.2 renames only the nested B57 local to `limb_remainder` and
updates the corresponding GMP 64-bit import and owning-result
comparison. This preserves the original GMP oracle, complete test
fixtures, assertions, numerical semantics, public APIs, production
kernels and benchmark workloads. It is solely a warnings-as-errors
source compatibility fix; full AppleClang and cross-platform CI
validation is still pending.

### B58 diagnostic: decimal cutoff-16 versus cutoff-32 and recursive split costs (20 October 2026)

On-target B57.2 official x86-64 benchmark (20000 iterations, 31 samples,
GMP verification PASS) measured total_time/GMP 1.146x and geomean
1.176x. At 2048/4096 bits production to_decimal was
3159.524/8717.960 ns, versus GMP 1861.035/4453.706 ns.
The independent B44 diagnostic measured, at 2048 bits, full
3115.608 ns, division replay 1212.801 ns and leaf replay
1411.688 ns; at 4096 bits, full 8595.959 ns, division replay
4789.633 ns and leaf replay 3056.365 ns. Replay fractions are
not additive and are not exact in-situ attribution.

B58 is deliberately measurement-only, not a production threshold
change. The existing production formatter, 16-limb cutoff,
arithmetic, official 84-operation benchmark and GMP workload
remain intact. An internal benchmark-only symbol in
src/big_unsigned.cpp duplicates the B22/B43 split recursion with
a 32-limb leaf cutoff while invoking the SAME production
to_decimal_blocks, square and divrem implementations. The B44
diagnostic now times production cutoff-16 against cutoff-32 on
the same 32 official-input magnitudes at 1024/2048/4096 bits.
Both outputs must match independent GMP before timing and the
experimental result is rechecked afterward. The diagnostic
also captures the exact production division operands at
each recursive depth, verifies quotient and remainder
against an independent mpz_tdiv_qr oracle, and times owning
native divrem against owning GMP mpz_tdiv_qr for identical
operand pairs grouped by depth. The split results, like
the original stage replays, are independent non-additive
proxies, not literal in-situ call-stack measurements.

Run on target:
  cmake --build build --target benchmark_integer_decimal_profile -j8
  ./build/benchmarks/benchmark_integer_decimal_profile --bits 1024,2048,4096 --iterations 20000 --samples 31
Promote cutoff or division changes only after repeatable on-host
GMP-verified gains in the unchanged 84-operation official benchmark
without regressions at other sizes or on ARM64. B58 source
compilation, integrated execution, CI and user-host timings
remain to be verified.

### B59 promote measured GCC x86-64 decimal cutoff-32 (20 October 2026)

B58 user-host x86-64 GCC benchmark, with 20000 iterations,
31 samples and GMP verification PASS on all 32 official inputs,
showed production cutoff-16 versus experimental cutoff-32:
- 1024 bits: 679.306 versus 679.785 ns (ratio 1.001x);
- 2048 bits: 3115.386 versus 2022.839 ns (ratio 0.649x,
  35.1% faster);
- 4096 bits: 8575.616 versus 6692.417 ns (ratio 0.780x,
  22.0% faster).
The original B58 recursive-division depth diagnostic also
confirmed native owning divrem / GMP mpz_tdiv_qr ratios:
- 2048 bits: depths 0/1 = 2.424x/2.819x;
- 4096 bits: depths 0/1/2 = 2.516x/2.818x/3.011x.
These independently timed depth replays are non-additive and
do not assign exact production call-stack costs.

B59 promotes the measured 32-limb leaf cutoff ONLY for
__x86_64__ GCC (not Clang or MSVC) with unsigned __int128.
ARM64, x86 Clang, Windows and other targets retain cutoff-16.
The existing duplicated B58 cutoff-16/32 recursion and top-level
formatter are consolidated into a single template parameterized
by a compile-time cutoff: production selects the architecture
constant while both diagnostic versions remain callable. No
division, decimal leaf, multiplication, GMP boundary or
public arithmetic API changes. The 84-operation official
benchmark and original 32-input benchmark workload are unchanged.
The standalone profile compares production against explicit
diagnostic cutoff-16 AND cutoff-32, preserving repeatable
same-host pre-promotion comparisons rather than letting the
experimental comparator become identical to production.
Both cutoff variants and each recursive split retain independent
GMP validation. Add 216 independent GMP binary-export/oracle
formatter checks at 15/16/17, 31/32/33 and 63/64/65 limbs,
covering exact powers, saturated, sparse and random inputs.

The B58 diagnosis suggests optimizing multiprecision division
for decimal-power divisors after this promotion; B59 itself
does not alter division. Projected improvements to the overall
time/GMP figure are NOT established by a stage proxy. Repeat
the official 84-operation benchmark at 20000 iterations
and 31 samples on the same host; validate all release,
sanitizer, coverage and Windows CI before treating B59 as
cross-platform verified. This documentation records the
user's B58 results, not an independently executed B59
integration or on-target timing.

### B60 diagnose GCC x86-64 single-limb remainder-only pi1 (20 October 2026)

B59's 84-case same-host x86-64 GMP-verified official benchmark
(20k iterations, 31 samples, 768 operand pairs/2688 operation
results) reported total_time/GMP=0.977x and geomean=1.162x,
versus the B57.2 1.146x and 1.176x respectively. Thus
the primary weighted-total parity objective is now met,
though a single successful full run does not guarantee
cross-run stability. At 2048 bits total_time/GMP=0.945x;
at 4096 bits 0.922x. Production 2048/4096 to_decimal
was 2068.069/6834.012 ns, versus 3159.524/8717.960
ns in B57.2. Unrelated construct_decimal changes should
not be attributed mechanically to the B59 selector.
All B59 GitHub Actions workflows (Unix, Windows, Coverage)
have completed successfully as checked on 20 October 2026.

The B59 production-versus-cutoff diagnostic independently
reported at 1024 bits 684.641/684.808/685.323 ns
(production/16/32); at 2048 bits 2036.237/3111.808/
2031.133 ns; at 4096 bits 6685.362/8586.050/
6686.932 ns. All full conversions and recursive split
results matched the independent GMP oracle.
The old cutoff-16 production split replay remains
distinct from B59's actual cutoff-32 call graph:
2048 depth-0/1 native/GMP=2.435x/2.862x and
4096 depth-0/1/2=2.521x/2.844x/3.007x.
These non-additive diagnostics must NOT be treated
as exact new production stack attribution.

B60 does NOT change the successful B59 production
division/decimal selector or the immutable 84-operation
benchmark. The B56 division profiler receives one
additional owning-result `rem_pi1` row on GCC x86-64
for its existing exact 32-pair `limb_small` and
`limb_wide` categories. The diagnostic-only function
in src/big_unsigned.cpp reuses the B57 exact normalized
preinverse digit recurrence, discards quotient digits,
and returns an owning BigUnsigned remainder. The
baseline `rem` row and GMP `mpz_tdiv_r` owning oracle
remain intact, all on identical original inputs;
the new candidate's results are independently checked
against GMP before and after timings, including
large and small divisors. Other hosts retain the
existing production rem code; no public header/API
or benchmark input pool has been changed.

Before integration, an isolated exact C++20 GMP-linked
reproducer tested 42000 deterministic cases per
compiler under GCC and Clang, with input lengths
8/16/32/63/64/65/128 limbs, zero/saturated/random
numerators and varied small/large odd divisors.
Both compilers passed -O3 -Wall -Wextra
-Wsign-conversion -Werror and all GMP oracle checks
after removing a test-harness-only unused helper.
This proves the isolated recurrence, not an
integrated build or speedup. No on-target B60
measurements or full CI results have yet been
verified.

For same-host selection use:
  cmake --build build --target benchmark_big_unsigned_division_profile -j8
  ./build/benchmarks/benchmark_big_unsigned_division_profile --bits 1024,2048,4096 --iterations 4000 --samples 15
If `rem_pi1` beats the production `rem` row
repeat at 20000 iterations/31 samples and then
consider a GCC-only measured production selector,
with independent 84-case aggregate and geomean
gates. Reject the candidate if it loses, regardless
of earlier isolated proxy assumptions. Keep B59
as rollback baseline for weighted-total performance.

### Post-B59 performance acceptance gates: per-size parity and worst-case cap (20 October 2026)

User-defined target for the x86-64 Integer vs GMP performance campaign
is now strictly stronger than the prior aggregate-total target:

1. Each of the six official sizes (53, 128, 256, 1024, 2048 and
   4096 bits) must individually satisfy total_time/GMP < 1.000x.
2. Every official (operation, size) pair must satisfy
   Integer/GMP < 2.000x; equivalently ALL worst < 2.000x.
3. Both gates must hold in the SAME unmodified official 84-operation
   benchmark, with 20000 iterations, 31 samples and independent
   GMP correctness PASS. Track the global total and geomean as
   diagnostics, but neither substitutes for either acceptance gate.
   Preserve the production semantics, benchmark input corpus and
   result-lifetime matching. Recheck near-threshold results with
   repeated same-host runs; do not claim durable success from
   noise-level wins.

The latest user-reported x86-64 B59 official baseline is:
  bits    total_time/GMP   per-size worst/GMP
    53        1.295x            2.296x (pow3)
   128        1.493x            2.518x (pow3)
   256        1.423x            2.084x (quotient)
  1024        1.074x            3.238x (quotient)
  2048        0.945x            3.319x (quotient)
  4096        0.922x            3.073x (remainder)
   ALL        0.977x            3.319x (quotient, 2048)
Global geomean remains 1.162x but is not a pass/fail condition.

The user also supplied the independent ARM64 reference:
  bits    total_time/GMP   per-size worst/GMP
    53        0.587x            1.016x (to_decimal)
   128        0.989x            1.801x (pow3)
   256        0.999x            1.682x (pow3)
  1024        0.856x            1.813x (quotient)
  2048        0.973x            1.944x (quotient)
  4096        1.018x            1.860x (remainder)
   ALL        0.981x            1.944x (quotient, 2048)
ARM64 meets the worst-case gate, but its 4096-bit total is
still above 1.000x and 128/256-bit totals have little margin.
Avoid regressing ARM64 while tuning x86; if the two acceptance
gates are applied cross-platform, ARM64's 4096-bit aggregate
also needs improvement.

Next work priority is determined by the FAILED gates, not just
the largest absolute timing or a global ratio:
- Improve 53/128/256-bit total-time ratios without sacrificing
  large-bit wins; diagnose short-operand owning-operation overhead
  and 53/128-bit pow3 before promoting microkernels.
- Reduce 1024/2048-bit quotient worst cases to <2x while
  preserving their established GMP-oracle coverage.
- Reduce 4096-bit remainder worst below 2x; B60's rem_pi1
  is only a diagnostic candidate, not a production speedup.
- Maintain the already passing 2048/4096-bit x86 total-time
  gates and the ALL total-time/GMP 0.977x baseline.

The B60 diagnostic commit 351ff6761711 leaves the B59
production arithmetic unchanged. B60 CI and same-host
rem_pi1 timing are not yet confirmed. A change may improve
one failure yet regress another: require an updated complete
84-case comparison before promotion, with a recorded
rollback baseline. No algorithm or benchmark change is
made by this acceptance-policy journal entry.

### B60.1 fix namespace of diagnostic declaration under GCC (20 October 2026)

The user's Ubuntu x86-64 GCC B60 build failed before running the division
profile. The new forward declaration for
`rem_limb_preinverse_diagnostic` in
benchmarks/benchmark_big_unsigned_division_profile.cpp was mistakenly
written inside the file's anonymous namespace, creating a *nested*,
distinct `Ariadne` namespace. `BigUnsigned` belongs to the global
`::Ariadne` namespace, so the compiler rejected the declaration and
subsequent calls. B60.1 moves only that forward declaration into
global `::Ariadne`, before opening the anonymous namespace. The
benchmark cases, diagnostic algorithm, GMP verification, official
benchmark, production selectors, tests and public headers are all
unchanged. Native integration/CI must still confirm the fix; no
performance result is claimed for B60.1.

### B61 promote measured single-limb remainder-only pi1 on GCC x86-64 (20 October 2026)

On the user's same x86-64 GCC host, the B60.1 12-shape B56
division profiler (4000 iterations, 15 samples, 32 deterministic
independently GMP-verified inputs per category) measured exact
owning `rem` versus `rem_pi1` times:
  bits   divisor     rem ns     pi1 ns    rem/pi1
  1024   limb_small  439.277    106.466    4.13x
  1024   limb_wide   335.457     77.273    4.34x
  2048   limb_small  631.118    131.504    4.80x
  2048   limb_wide   655.604    129.613    5.06x
  4096   limb_small 1240.534    238.935    5.19x
  4096   limb_wide  1310.074    237.269    5.52x
The diagnostic worst `rem_pi1`/GMP ratio was 1.206x
(limb_wide, 4096 bits) versus `rem`/GMP 6.680x.
These are diagnostic B56 categories, NOT the official 84-case
worst ratios. The user's report passed its independent GMP oracle.
At 4096 bits the short_2 multi-limb rem still took 1804.791 ns
vs GMP 280.796 ns (6.427x); this is NOT addressed by B61.

B61 accepts only the experimentally supported GCC x86-64
single-limb remainder path for >=8 canonical dividend limbs
and a non-power-of-two divisor. The implementation extracts
the validated B60 exact normalized pi1 recurrence into one
shared internal helper, called by production rem-only selected
division and by its diagnostic wrapper. No quotient buffer or
scratch allocation is introduced. The B57 power-of-two fast
path executes before the candidate. Quotient-only and joint
division preserve B57's existing selection. Other compilers,
short inputs, ARM64 and Windows retain their old paths.

The B56 diagnostic retains both the original `rem` row
(now selecting production B61) and `rem_pi1` candidate,
and additionally times GMP-oracle-verified `rem_divq`,
explicitly preserving the older owning DIVQ recurrence
for direct A/B measurement after promotion. The 32-pair
pool, 12 shapes and official 84-case benchmark are untouched.
Add GMP-oracle owning rem and joint divrem cases covering
zero, 7/8/9, 63/64/65 and 128 limbs, varied odd/wide/power
divisors and dense/sparse/saturated operands. Numerical
comparisons remain independent; do not change test
expectations to make the optimization pass.

The performance acceptance gates remain all six per-size
total_time/GMP < 1.000x and all 84 operation/size worst
ratios < 2.000x, with no significant ARM64 regressions.
The B60.1 diagnostic result alone does not establish that
B61 meets either integrated gate. CI and full user-host
benchmark remain to be checked before calling B61 an
integrated success. Roll back to B59/B60.1 if the new
production selector regresses independent correctness
or aggregate performance.

### B62 diagnostic: quotient-only single-limb pi1 threshold at 16/32 limbs (20 October 2026)

The user's B61 integrated official 84-operation x86-64 benchmark
(20000 iterations, 31 samples; 768 GMP-verified paired inputs
and 2688 correct operation results) measured per-size totals
53:1.452x, 128:1.538x, 256:1.486x, 1024:1.092x,
2048:0.962x, 4096:0.941x and ALL:0.999x. The global
worst is now quotient at 2048 bits, 3.386x, followed
by quotient at 1024 bits, 3.282x. The 4096-bit
official remainder improved from B59 1120.469 ns to
B61 716.668 ns, moving its official ratio from
3.073x to 1.968x. At 1024 and 2048 bits the remainder
fell from 286.099/539.469 ns to 176.753/320.884 ns.
The aggregate total worsened from B59 0.977x to 0.999x
despite the remainder gain: construction and other
operations moved as well, and the short sizes cannot
invoke the new >=8-limb remainder kernel. Do not
attribute unrelated shifts to B61; repeat full runs
under comparable host conditions to estimate noise.

B57 already selects exact quotient-producing normalized
pi1 on GCC x86-64 when the input has >=64 limbs.
At 1024/2048 bits (16/32 limbs), it still uses
per-digit native DIVQ. B62 adds **diagnostic-only**
caller-owned quotient comparators for GCC x86-64
on the B56 limb_small/limb_wide 32-pair exact GMP
pools. The pi1 comparator replays B57's normalized
inverse estimation and digit correction for any
input length, including 16 and 32 limbs. The DIVQ
comparator uses the old exact recurrence independently
of current production size selection. Both fill the
same reusable quotient output span and return a 64-bit
remainder. Each quotient AND remainder is checked
against the independent GMP oracle before and after
the timed phases; the GMP timed counterpart uses
reusable `mpz_tdiv_q` output for fair caller-owned
comparison. This kernel-only preallocated timing is
not an owning `quot()` speedup claim and cannot be
added to the original owning or divrem rows.

B62 changes NO production quotient/remainder dispatch,
arithmetic, official benchmark corpus, test cases,
GMP oracle or FloatAP behavior. It adds two B56 rows:
`quot_pi1_prealloc` and `quot_divq_prealloc` for
both single-limb divisor categories, and records
their matched GMP preallocated quotient results.
After CI, measure on the SAME host using:
  ./benchmarks/benchmark_big_unsigned_division_profile -t dark --bits 1024,2048,4096 --iterations 4000 --samples 15
Only if 16/32-limb pi1 wins robustly in matched kernel
measurements should we experiment with a lower GCC-only
production quotient threshold, followed by GMP
verification and the untouched 84-case acceptance
gates (all six per-size totals <1x; worst <2x).
The new diagnostic is not yet compilation/CI-verified
or independently measured on the user's host.

### B63 x86-64 GCC quotient-only normalized pi1 threshold 16 limbs (20 October 2026)

The user's B62 x86-64 GCC B56 diagnostic (4000 iterations,
15 samples, original 32-pair GMP-verified operand pools) measured
independent caller-owned `quot_pi1_prealloc` against the exact
`quot_divq_prealloc` recurrence. At 1024 bits:
  small 90.369 vs 300.657 ns (pi1 3.33x faster),
  wide  83.630 vs 318.015 ns (pi1 3.80x faster).
At 2048 bits:
  small 147.365 vs 614.124 ns (4.17x),
  wide  144.513 vs 638.188 ns (4.42x).
At 4096 bits:
  small 271.682 vs 1213.953 ns (4.47x),
  wide  266.891 vs 1293.786 ns (4.85x).
The diagnostic worst pi1/GMP was 1.541x (1024-bit wide);
old DIVQ worst was 7.052x (4096-bit wide).
All quotient AND remainder results matched the independent GMP oracle.
These are preallocated kernel figures, not owning `quot()` or
the official 84-case aggregate; do not equate their ratios.

Crucially, the unchanged official 32-case integer division workload
in benchmark_exact_arithmetic_workload.hpp uses single-limb
divisors for i%4==0 and i%4==3, i.e. 16 of 32 inputs.
For i%4==1 and i%4==2 it uses half-width and near-equal
divisors. Consequently B63 can improve the real official
`quotient` benchmark at 1024/2048 bits without modifying
the workload, but cannot address the multi-limb half of
those operation cases. It is not yet proven that the worst
official ratios fall under 2x.

B63 promotes the measured `>=16`-limb pi1 threshold in
`BigUnsigned::divrem_limb_into` ONLY for GCC x86-64
with __int128, replacing its earlier `>=64` threshold.
It does not change the exact preinverse algorithm, the
power-of-two fast path, the B61 rem-only >=8-limb
selection, multi-limb Knuth division, ABI or public
headers. Clang, ARM64, MSVC and short quotient
inputs retain their existing selection.
The B62 pi1 versus DIVQ diagnostic controls and their
32-input GMP oracle stay unchanged as an A/B fallback.
Extend the B57 test count list to 15/16/17 and
31/32/33 limb sizes in addition to existing 63/64/65
cases; quotient and 64-bit residual remain GMP-verified,
including dense, sparse, saturated and zero inputs.
No test expectations or official benchmark inputs are changed.

Acceptance remains ALL six official per-size total_time/GMP
strictly <1x and all 84 per-operation slowdown ratios
strictly <2x, checked in one GMP-verified official run
with 20000 iterations and 31 samples. Repeated on-target
runs should confirm borderline wins. B61 is a baseline,
not proof of B63 outcomes. All B62 Unix/Windows/Coverage
GitHub Actions workflows completed successfully; B63
integrated CI and performance validation are pending.

### B64 short-operand pow/square ownership diagnostic (20 October 2026)

B63 on-target 20,000x31 GCC x86-64 official benchmark verified
768 paired GMP inputs and all 2688 paired operation outputs.
Per-size total_time/GMP ratios were:
53=1.434x, 128=1.468x, 256=1.423x,
1024=0.977x, 2048=0.863x, 4096=0.895x; ALL=0.933x.
Per-size worst operation ratios were:
53 pow3=2.366x, 128 pow3=2.476x,
256 square=2.087x, 1024 quotient=1.969x,
2048 compare=1.879x, 4096 compare=2.067x.
The overall worst is 2.476x (pow3, 128 bits).
B63 reduced the 1024/2048 quotient to 176.383/306.194 ns,
from B61 292.537/554.274 ns. The 53/128/256 short-size
total-time gates are still far from <1x; the 4096-bit
worst compare margin is small and could be sensitive
to repeat-run noise.

B64 adds a separate GMP-verified benchmark target named
benchmark_integer_small_power_profile. It uses the EXACT
same 32 signed integer_specifications(bits) inputs as the
official benchmark and does not modify the official
84-case test corpus or arithmetic semantics.
It times whole owning Integer results and owning GMP
counterparts: production pow(z,Nat(3u)), a generic
exponentiation-by-squaring alternative which defers
copying the original base until an accumulator factor
is actually needed, and a direct cube chaining
diagnostic as a possible ownership lower-bound control.
The direct cube is NOT to be promoted solely because
the official workload includes exponent 3.
It also measures production alias z*z versus an
owning BigUnsigned::square-based signed Integer
construction, retaining the exact output semantics.

All alternatives must pass an independent mpz_pow_ui/
mpz_mul oracle before and after the timed rows, at
exponents 0/1/2/3/4/5/7/8/9 and the original signed
32-input official corpus, with additional zero and
signed-unit identities. The two timing groups
measure whole owned result objects, not non-additive
kernel stages. Different implementations may
allocate internally; this is precisely one of the
effects under investigation and must not be claimed
as a proven bottleneck until measured. As in the
official benchmark, each call creates its result
inside the timed loop and GMP uses fresh owning
mpz_t results. The explicit generic lazy algorithm
is a diagnostic candidate, not a B64 production
change. No change is made to any production path,
public API, existing tests or FloatAP; this new
target is appended to benchmarks/CMakeLists.txt.

Same-host command, from build and with user's preferred
theme and single-line formatting:
  ./benchmarks/benchmark_integer_small_power_profile -t dark --bits 53,128,256,1024,2048,4096 --iterations 20000 --samples 31
The result will determine whether deferred-copy generic
exponentiation can reduce the 53/128-bit pow3 worst
without the unsupported conclusion that one micro
optimization will bring all three small-size totals
under 1x. B64 integrated compilation/CI and actual
GCC user-host timings are not yet confirmed.

### B65 assess ARM64 B37-B40 transfers at short widths and MSB-first power (20 October 2026)

Following the user's question about lessons from the earlier ARM64
square/power campaign, inspect both B37-B40 and x86 B32/B53/B54:
- ARM64 B37's triangular/symmetric square was already present on
  GCC x86-64 since B32, with the SAME >=8-limb production gate.
- B38's skip-zero outer rows and B39/B40's adaptive dense/sparse
  protection were already transferred to x86 by B53 and B54.
- The existing balanced/dense Karatsuba and zero-row dispatch
  serves large inputs; it does not address 1/2/4-limb work.
There is no missing ARM64 assembly kernel that can be copied
verbatim to x86, and the x86 production square under eight limbs
still uses generic multiply.

The user's B64 whole-operation GMP-verified (20000 iterations,
31 samples, same 32 official signed input pairs per size) measured:
  bits   pow3 production    lazy     direct chain      square prod/direct
    53       98.435       92.252        56.310          31.173/30.067
   128       95.375       88.665        78.088          32.736/31.369
   256      119.706      115.919       105.947          52.842/52.384
  1024      232.114      209.079       178.515          94.957/90.719
  2048      342.477      319.100       282.644         150.606/146.884
  4096      638.445      616.891       555.272         274.935/272.190
GMP oracle PASS for production, lazy, direct cube and square,
including powers 0/1/2/3/4/5/7/8/9 and signed zero/unit cases.
The lazy 3-10% gain is much smaller than the 12-43% direct
control's gain. A cube-only production shortcut would be a
fragile and workload-specific response, not a general solution.

B65 adds two diagnostic comparisons to the EXISTING B64
benchmark, without modifying production arithmetic or its
84-operation GMP workload:
(1) Generic *most-significant-bit-first* exponentiation by
squaring, borrowing the original base until the first
square is needed, and then owning intermediate results.
For exponent three, the generic recurrence naturally
computes square(base) then mul(square,base), matching
the direct cube control, but works for any unsigned
exponent (including exponent zero and signed inputs).
The GMP oracle covers exponents 0/1/2/3/4/5/7/8/9/12.
(2) A benchmark-only forced invocation of the existing
GCC x86-64 symmetric-square kernel at 1..7 limbs,
with a bounded 14-limb local output before constructing
the same owning BigUnsigned. Production still dispatches
to generic multiplication under eight limbs; ARM64,
Clang, MSVC and other architectures remain unchanged.
The forced path is additionally checked against independent
GMP for 0..8-limb zero/saturated/sparse/random values
and official signed inputs. The forced comparison is
reported at 53/128/256 bits; 1024+ retains its original
selected square path and has no new short-kernel row.

Both variants return normally owned Integer objects.
MSB-first exponentiation changes the intermediate
multiplication order at general exponents, so it must
not be promoted from a cube-only measurement without
wider correctness and same-host performance coverage
and ARM64 regression checks. Likewise the symmetric
short square should only be promoted if it improves
the *whole owning* square result at all relevant small
sizes, not merely its preallocated arithmetic.
The larger objective is six total_time/GMP <1x and
all 84 worst ratios <2x; optimizing power alone cannot
be assumed to close the ~1.4x 53/128/256 total gap.
B65 native timings and integrated CI are pending.
Command from build (one line; user-preferred -t dark):
  ./benchmarks/benchmark_integer_small_power_profile -t dark --bits 53,128,256,1024,2048,4096 --iterations 20000 --samples 31

### B66 promote GMP-verified MSB-first signed power on GCC x86-64 (20 October 2026)

The user ran B65 on the original 32 signed GMP-verified official
inputs (20,000 iterations, 31 samples), obtaining these
production -> generic MSB-first pow3 results:
  53 bits:   100.293 -> 65.688 ns (0.655x, 1.558x GMP)
 128 bits:    97.224 -> 87.200 ns (0.897x, 2.294x GMP)
 256 bits:   124.522 ->114.540 ns (0.920x, 1.940x GMP)
1024 bits:   248.337 ->202.110 ns (0.814x, 0.526x GMP)
2048 bits:   369.999 ->322.656 ns (0.872x, 0.261x GMP)
4096 bits:   651.561 ->588.354 ns (0.903x, 0.159x GMP)
GMP verification PASS for signed powers 0/1/2/3/4/5/7/8/9/12
and signed zero/unit controls. The diagnostic direct-cube
sequence was faster, but it is not a general justification
to special-case exponent 3. The B65 short symmetric square
control improved 53 bits by only 5.5%, changed 128 bits
by -0.9% and regressed 256 bits by 0.2%; DO NOT lower
the >=8-limb square production gate based on these data.

B66 replaces ONLY GCC x86-64 fused-build Integer::pow's
right-to-left B32 exponentiation with the already
GMP-tested general MSB-first left-to-right binary algorithm.
It starts from the highest set exponent bit, borrows the
immutable input magnitude instead of copying it,
and materializes each required square and multiplication
into an owned temporary before assignment. It handles
exponents 0 and 1 explicitly, preserving 0^0=1,
negative-sign parity and zero normalization. The
exponent is the general Nat argument, not benchmark-only
exponent three. The original B32 implementation remains
selected on ARM64, AppleClang/Clang, MSVC and other builds.
No public headers, FloatAP logic, GMP dependency, square
or multiplication kernel, or the original 84-operation
benchmark and its 32-input corpora are changed.

Extend existing integrated GMP-oracle power tests with
60 signed/unsigned values at canonical limb counts
0/1/2/3/4 and 18 exponents each, including
31/32/33/63/64/65 boundaries and sparse/saturated/random
limbs (1080 oracle comparisons), alongside all previous
power and arithmetic tests. As B65 timings are diagnostics
rather than new integrated production B66 results, require
the original same-host 84-case GMP-verified official benchmark
before claiming B66 achieved a given slowdown reduction.
All six size-specific total_time/GMP must be <1.000x,
and every one of the 84 worst ratios <2.000x.
B66 cannot by itself resolve 53/128/256-bit total-time
ratios of approximately 1.4x, and 128-bit pow3 MSB-first
was still 2.294x in B65. Verify the full CI and ARM64
baseline before any wider algorithmic promotion.

Run from build:
  ./benchmarks/benchmark_compare_integer_backends -t dark --bits 53,128,256,1024,2048,4096 --iterations 20000 --samples 31
B66 native compile, integrated CI and official timing
are not yet verified at commit authoring time.

### B67 own-result short-operand performance decomposition (20 October 2026)

B66's full official x86-64 benchmark (20k iterations, 31 samples,
768 paired GMP-verified inputs and 2688 correct paired operation
results) reported totals/GMP:
  53=1.356x, 128=1.527x, 256=1.490x,
  1024=1.029x, 2048=0.912x, 4096=0.940x,
  ALL=0.979x; the overall worst is 2.447x
  for pow3 at 128 bits. Worst per-size:
  53 pow3=2.264x; 128 pow3=2.447x;
  256 square=2.117x; 1024 quotient=2.028x;
  2048 quotient=1.939x; 4096 quotient=2.074x.
Compared with B63, the MSB-first optimization did
improve actual pow3: 53 72.001->68.925 ns;
1024 232.114->195.173 ns; 2048 342.477->306.423 ns;
4096 638.445->590.928 ns. But other unchanged
operations, particularly construct_decimal,
increased in absolute time: at 1024 444.540->587.105,
2048 986.184->1250.536 and 4096
2657.408->3226.152 ns. Quotient/remainder also
moved without code changes. Since B66 did not
alter construction or division kernels, do not
infer B66 caused their reported slowdown; host noise,
load, clock, allocator and compilation context
remain possible confounders. A single aggregate
benchmark run is not sufficient to reject a
consistently faster production pow kernel.

B67 adds an independent small-operations benchmark
target using both original 32-input
integer_specifications(bits) and
integer_division_specifications(bits). It times
(1) raw owning GMP, (2) ordinary owning signed
Integer API, and (3) direct owning unsigned
BigUnsigned on the identical absolute input
magnitudes for construct_decimal, copy, compare,
mul, square, quotient, remainder, pow3 and
to_decimal. All unsigned results are checked
against an independent GMP absolute-value oracle
before and after timing. The integer and GMP
columns have the exact same signed semantics
as the official workload; the unsigned column
is a non-additive *magnitude-only control*.
Its arithmetic can take different branches
because sign is absent, so the Integer/core
difference is not a direct measurement of
signed-wrapper overhead. It can nevertheless
help distinguish a slow primitive from an
especially costly signed path on 53/128/256 bits.
For pow3 the core controls direct
square(a)*a and is explicitly not identical
to the whole exponentiation algorithm;
B65 already established this is a useful
lower-bound diagnostic, not a generic cube
production shortcut.

No original official 84-case benchmark or its
input pool, production arithmetic, GMP oracle,
tests, or public headers are modified by B67;
a new diagnostic target is appended in CMake.
Acceptance remains each of six
total_time/GMP <1.000x and worst <2.000x over
all 84 cases on the same host, with stable
repeat runs. B67 build/CI and user-host
measurements remain to be verified. Command
from build, one line with user's preferred
theme:
  ./benchmarks/benchmark_integer_small_operations_profile -t dark --bits 53,128,256,1024,2048,4096 --iterations 20000 --samples 31

### B68 same-signed-semantics A/B for decimal adapters and short products (20 October 2026)

The user's B67 full integrated x86-64 GMP-verification PASS
(20000 iterations, 31 samples; the original 32 official
signed input and division-input corpora) reported these
Integer owning / BigUnsigned owning magnitude ratios:
  bits  construct_decimal   mul   square   pow3   to_decimal
    53      1.395x         1.581x  1.313x  1.471x    1.533x
   128      1.250x         1.535x  1.294x  1.334x    1.391x
   256      1.146x         1.298x  1.182x  1.244x    1.231x
At 53/128/256 bits, additional apparent unsigned-core
versus Integer absolute nanoseconds were respectively:
to_decimal +39.820/+41.863/+43.543,
construct_decimal +26.904/+23.005/+22.211,
mul +12.597/+13.212/+13.563.
These are NON-ADDITIVE, not direct signed-wrapper
attribution: unsigned controls omit sign semantics,
return different types and may follow different paths.
B67 also indicates integer/core mul overhead converges
to 1.092x at 4096 bits; to_decimal converges to 1.012x.
Do not infer a proven single root cause or alter
the official benchmark from this comparison.

B68 extends the unchanged independent B67 benchmark with
an additional A/B section whose controls match the complete
SIGNED output type and mathematical semantics, on the
same exact 32 input pairs and same iteration/sample settings:
- construct_signed: normal Integer(String(decimal)) versus
  independently parsed sign with production from_decimal
  and Integer::from_magnitude, avoiding the intermediate
  Ariadne String for the candidate.
- to_decimal_signed: normal std::string(Integer::literal())
  versus production BigUnsigned::to_decimal on the input's
  magnitude plus an exact conditional leading '-' inserted
  into the returned owning std::string; the sign, including
  canonical zero, is preserved.
- mul_signed: normal Integer operator* versus a signed
  result constructed from production BigUnsigned::mul and
  the XOR of the two input signs.
- square_signed: normal self-alias Integer z*z versus
  signed (nonnegative) construction from production
  BigUnsigned::square.

The new candidate results are independently verified
against signed mpz (including input sign, full multiplication,
square, decimal string and signed -0/+0 normalization),
before AND after their own timed loops. The exact GMP
counterparts retain owning result lifetimes. The original
B67 report with unsigned magnitude controls is still
printed intact. These measurements can distinguish
expensive String/signed adaptation from primitive costs
more reliably than B67's unsigned comparisons, but remain
independently timed and susceptible to scheduler/clock
variance; avoid attributing their differences additively.

B68 does not modify Integer/FloatAP production, public
APIs, original 84-operation official benchmark, its
GMP oracle, 32-case input corpus, or existing test
expectations. Promoting an adapter or signed arithmetic
change requires a clear same-host repeated win AND
unchanged 84-case correctness/per-size total and worst
acceptance gates (each total <1x, every worst <2x).
No new production behavior or benchmark weights are
introduced. B68 integrated CI and user-host timings
are not yet established.

Run from build, single line with dark theme:
  ./benchmarks/benchmark_integer_small_operations_profile -t dark --bits 53,128,256,1024,2048,4096 --iterations 20000 --samples 31

### B69 test ARM64 carry-chain transfer against GCC x86 fused at short widths (20 October 2026)

The user provided a second ARM64 official 20,000-iteration/31-sample
benchmark with GMP correctness PASS (768 paired inputs and 2688
operation outputs): total_time/GMP was 0.604x at 53 bits, 0.985x at
128, 0.986x at 256, 0.830x at 1024, 0.971x at 2048 and 1.008x at
4096. Worst/GMP remained below 2x at all sizes (ALL=1.805x).
Thus the native Integer representation and ownership strategy can
perform competitively with GMP on another architecture at short
widths; the x86-64 1.3-1.5x total ratios must not be treated as a
platform-independent unavoidable signed wrapper overhead.
Cross-machine absolute timing comparisons reflect different CPU
and compiler performance and do not isolate an implementation cause.

The same user's B68 x86-64 signed GMP-verified A/B comparisons
(20,000 iterations, 31 samples on original 32 input pairs)
measured signed direct construction candidate/production ratios
0.768x/0.909x/0.934x and signed direct formatter
0.723x/0.835x/0.888x at 53/128/256 bits. The candidate
signed multiplication was 0.949x/1.070x/0.984x and
signed square 1.005x/1.006x/1.010x. This shows repeatable
same-semantics wins for String adapter avoidance, while
manual signed multiply/square construction offers no consistent
benefit. Large x86 operands often show minimal wrapper
differences, so blaming sign handling alone is unwarranted.

B69 makes the proposed ARM64-inspired kernel test on the SAME
GCC x86-64 host, rather than comparing ARM64 and x86 absolute
times. The existing B68 standalone diagnostic gains a third
section ONLY when the production library reports both
selected_multiply_uses_fused() and
has_x86_64_carry_multiply(). For sizes 53, 128 and 256 bits,
it compares the exact existing production fused 128-bit
kernel against the explicitly available x86-64 carry-chain
kernel in TWO independent, non-additive ways:
(1) no-allocation fixed-size caller-owned eight-limb output
arrays, independently timed for nonaliased multiplication
and self-square; and
(2) whole owning SIGNED Integer objects with matching result
type, sign parity and GMP oracle, using both forced kernels.
The whole-result experiment constructs a canonical
BigUnsigned from the eight-word stack buffer; that may
allocate and copy differently from the current production
owning mul(), so its timings are A/B between *forced
variants*, not a direct claim that merely flipping the
production kernel will yield exactly the same performance.

All B69 operations use exactly the existing 32 official
Integer input pairs at each selected size. Additional
300 random/dense/sparse/saturated/zero signed test pairs
at canonical limb lengths 0..4 exercise both forced
kernels for independent GMP quotient-free full products,
including negative signs, squares and nonempty/empty
operands. Both full owning signed products and raw
fixed-buffer output digits must match GMP before AND
after timed experiments. The candidate uses the existing
BigUnsigned::mul_into_fused and
BigUnsigned::mul_into_x86_64_carry entrypoints: no
new architecture assembly, no production dispatch,
ABI, public API, official benchmark workload/weights,
oracle or other test behavior is changed.

Run from build, as a single command:
  ./benchmarks/benchmark_integer_small_operations_profile -t dark --bits 53,128,256 --iterations 20000 --samples 31
A production change should be considered only on
repeatable native GCC x86 improvements in short
owning results, followed by the unmodified 84-case
GMP-verified acceptance gates (all six per-size
totals strictly <1x and every worst ratio <2x),
with no ARM64 regression. Do not confuse the
ARM64 success with proof x86 carry will win.
B69 integrated build and host measurements
remain unverified at commit authoring time.

### B70 short multiplication: selected fused dispatch, ownership and inlining controls (20 October 2026)

The user's B69 GCC x86-64 results (20k iterations, 31 samples,
independent GMP oracle PASS for official and extra short inputs)
disprove the hypothesis that copying ARM64's carry-chain kernel
directly to x86 improves short products. Carry/fused ratios
for matched preallocated multiplies were 1.019x/1.151x/1.031x
at 53/128/256 bits and 1.054x/1.150x/1.033x for square:
carry was SLOWER at all three sizes. Full owning signed
forced-carry/forced-fused was likewise worse in all cases,
e.g. 128-bit mul 32.279/29.876=1.080x.
Keep GCC x86-64 fused production selected, and keep the
ARM64 native carry kernel and other ISA dispatch unchanged.

Nonetheless B69 full signed owning forced-fused vs production
reported:
  bits   mul production -> forced fused   square production -> fused
    53   34.462 -> 27.348 ns (0.794x)    31.467 -> 26.286 (0.835x)
   128   36.850 -> 29.876 ns (0.811x)    34.107 -> 27.449 (0.805x)
   256   58.266 -> 46.562 ns (0.799x)    55.136 -> 44.120 (0.800x)
The two owning results use equal signed types and GMP-verified
values, but B69's forced variant constructs a BigUnsigned from
a stack array *in the benchmark translation unit*, whereas
production Integer::mul and BigUnsigned::mul live in separate
translation units, with different construction/normalization,
function boundaries and possible compiler inlining. Since
production already selects the SAME fused arithmetic kernel,
these 19-21% reductions CANNOT be attributed solely to
the multiplication instruction sequence.

B70 extends the original B67/B68/B69 standalone profiler,
without changing production code, tests, ABI, original
official benchmark or any 32-case GMP corpus. It adds:
- explicit production-selected BigUnsigned::mul_into alongside
  forced mul_into_fused in identical caller-owned 8-limb
  buffers (multiply and square);
- full owning signed results constructed from each selected/
  forced fixed stack buffer;
- an intentionally noinline wrapper of exactly the same
  fused owning variant, to estimate the call/inlining effect;
- a noinline signed wrapper of production owning
  BigUnsigned::mul, contrasted with Integer::mul.
Both selected and forced output digits and signed results
are checked against the independent GMP oracle on all
official 32-pair signed inputs and the existing 300
random/saturated/sparse/zero 0..4-limb stress cases.
Recheck the original oracle after all timed rows.

The noinline attribute is a compiler optimization hint
for this benchmark and does NOT exactly reconstruct
cross-translation-unit ABI/LTO behavior. All rows are
independently timed and not additive, and the production
results may be affected by scheduler/clock drift.
No performance gains or new production acceptance gates
are claimed until same-host B70 data are returned.

From build (one line, -t dark):
  ./benchmarks/benchmark_integer_small_operations_profile -t dark --bits 53,128,256 --iterations 20000 --samples 31
The binding target remains six official size-specific
total_time/GMP strictly <1x and all 84 worst ratios
strictly <2x, no benchmark fixture modifications.
CI/integration verification pending.

### B71 test short owning stack construction inside production library TU (20 October 2026)

B70 on the user's GCC x86-64 host (20,000 iterations,
31 samples, all independent B67/B68/B69/B70 GMP oracles PASS)
found production-selected and forced-fused kernels
approximately matched on the original 32 short operands.
Ratios of selected/fused (candidate/baseline in B70)
were 0.992x / 1.072x / 1.130x for 53/128/256-bit
preallocated mul and 1.036x / 1.015x / 0.944x for
preallocated square. The near-unity and changing
order do not support changing the GCC dispatch.
The B70 own-result selected/fused ratios were
0.996x / 0.915x / 1.014x for multiplication and
0.966x / 0.976x / 1.023x for squaring. These
mixed figures further suggest compiler/code-layout
or benchmark-run effects rather than a consistently
better arithmetic kernel.

B70's noinline fused-owning variants incurred
~7.8%/7.5%/6.1% overhead for mul and ~9.9%/7.3%/7.1%
for square at 53/128/256 bits. Noisy but repeatable
same-run noinline production core wrappers only saved
about 0.972x/0.972x/0.986x in nonaliased mul and
0.948x/0.950x/0.990x in self-square. They do NOT
reconstruct all code generation or ABI effects across
translation units. B70's simple production-vs-forced
owning comparison has a larger gap, but cannot yet
attribute it to instructions alone.

B71 explicitly moves a new diagnostic, exact
short owning BigUnsigned product control INTO THE
SAME src/big_unsigned.cpp translation unit as the
unchanged production BigUnsigned::mul implementation.
For input magnitude lengths at most four limbs, it
calculates the selected-kernel product directly
in a zeroed eight-limb local std::array, constructs
a canonical BigUnsigned from the resulting limb view
and returns the fully owned result. A second function
does exactly the same but explicitly uses the fused
kernel. Both fall back to normal mul() above the
four-limb diagnostic limit or on other toolchains,
so no production dispatch, public API or ARM64
algorithm is modified.

The existing B67/B68/B69/B70 benchmark gains
side-by-side independent median timings of:
(1) owning BigUnsigned::mul vs selected stack helper,
(2) selected vs forced-fused stack helper,
(3) signed owning Integer wrapper calling mul() vs
    the new stack helper, with equivalent non-inlining,
and the same for actual self-square operands.
All new helpers return the same exact owning type
and are linked as production-library functions,
rather than handwritten benchmark-local mul_into
entrypoints. They may still exhibit different code
placement or inlining across build options; this
diagnostic does not create an exact controlled
linker/inlining experiment.

Both unsigned magnitude results and full signed
wrapper outputs are independently checked with
GMP before and after timing, using the original
32 official signed input pairs at each size and
the existing 300 zero/saturated/sparse/random
short 0..4-limb stress pairs. Neither user-facing
benchmark input pools nor the original 84-operation
official timing and GMP oracle are changed. No
new public API or additional CMake build source
is needed. Return of the user's B71 same-host
results is required before proposing any
short-product production change.

The short-size performance strategy remains:
keep fused selected on GCC x86-64 pending robust
contrary evidence; keep native ARM64 carry-chain;
record direct signed decimal adapter improvements
from B68 as a separate possible production change.
Six per-size official total_time/GMP <1x and all
84 worst <2x remain the only acceptance gates.
No reweighted workload was introduced (per user).
Run from build in one line:
  ./benchmarks/benchmark_integer_small_operations_profile -t dark --bits 53,128,256 --iterations 20000 --samples 31
B71 compilation/CI/host measurements are pending.


### B71.1 fix B71 benchmark wrapper declaration order (20 October 2026)

The user reported a native GCC x86-64 benchmark build failure
in `benchmark_integer_small_operations_profile.cpp`:
`verify_b69` referenced `b71_stack_selected_signed` and
`b71_stack_fused_signed` before either was declared.
Move the two existing diagnostic function definitions intact
above `verify_b69` in the same anonymous namespace.
No algorithm, oracle, test pool, production arithmetic,
dispatch, build flag, official 84-case benchmark, or
benchmark timing methodology changes. The fix does not
claim performance or successful native build before
actual CI/user verification.

### B72 production short owning stack product on GCC x86-64 (20 October 2026)

The user's B71.1 on-host independent GMP-verified benchmark
(20,000 iterations, 31 samples; original 32 signed input
pairs and additional 300 short zero/dense/sparse/random cases)
measured the new SAME-library-TU short owning stack
construction against existing BigUnsigned::mul.
At 53 bits, unsigned multiply/square candidate/baseline
were 0.950x/0.954x and signed owning 0.981x/0.978x.
At 128 bits they were 0.978x/0.981x and
0.969x/1.006x. At 256 bits they were
0.928x/0.884x and 0.925x/0.902x.
Selected/fused stack controls were nearly equivalent,
sometimes reversing rank. Thus changing fused to
carry-chain, which the prior B69 tests showed
consistently slower, is not justified. The notable
short owning gain at 256 bits appears linked to
result construction rather than a faster multiplication
kernel. The B70 inline/noinline controls have mixed
results (including reversed rank at 128 bits in B71);
do not treat these as a precise additive call cost.

B72 promotes the precisely measured B71 diagnostic
short owning construction to **production BigUnsigned::mul**
ONLY under the existing GCC x86-64 fused-build guard
and only for nonzero canonical inputs of at most
BigUnsigned::inline_capacity=4 limbs each. It uses
an eight-limb zeroed std::array, writes the exact full
product with the UNCHANGED selected mul_into fused
kernel, then constructs a fully owned normalized
BigUnsigned from a view of the local buffer. Results
larger than four limbs still use the existing owned
vector representation. If either input is >=5 limbs,
the original owning construction/Karatsuba selection
remains untouched. Zero returns at the original
early check; size arithmetic overflow checks retain
their original order. No global buffers, special
benchmark input assumptions, numerical approximation,
public header/ABI, FloatAP arithmetic, ARM64,
Clang or MSVC selected algorithms are changed.
The B71 stack selected/fused diagnostics remain
available as controls, although the B71 production
baseline at <=4 limbs will now use the promoted
variant; compare with B71.1's pre-promotion commit
for historical A/B reference.

Integrated tests add a deterministic independent
GMP oracle over 6x6 operand limb-length boundaries
0/1/2/3/4/5, 32 modes per pair (1152 combinations),
checking mul(a,b), mul(b,a), mul(a,a) and square(a)
for zero, all-one, sparse and random inputs. The
existing 6000-random multiplication oracle and all
previous tests remain intact. The original official
84-operation benchmark and its 32-case input pools
are unchanged. The user must verify both strict
acceptance gates on the SAME GCC x86-64 machine:
total_time/GMP <1x for each of 53/128/256/1024/
2048/4096 bits, and every operation/size ratio <2x
(equivalent to global worst <2x), before considering
this iteration complete. ARM64 baseline preservation
requires its own build/performance verification.
B72 CI, compilation and measured production performance
are pending at publication time. From build:
  ./benchmarks/benchmark_compare_integer_backends -t dark --bits 53,128,256,1024,2048,4096 --iterations 20000 --samples 31


### B73 remove temporary decimal ownership and one signed-result move (20 October 2026)

The supplied x86-64 official baseline reports total_time/GMP
1.429x/1.472x/1.458x/1.034x/0.904x/0.938x at
53/128/256/1024/2048/4096 bits. The short groups require
approximately 30.0%/32.1%/31.4% reductions in their native sums,
not merely a worst-operation fix. The attached B68-B71 profiles
predate production B72. Keep these separate from post-B72 results.
B69's carry kernel remains slower than fused; B70/B71 provide no
basis for changing the selected arithmetic or assigning additive
costs to sign handling, inlining or allocation.

B73 changes ownership and temporary storage, with no new ISA kernel:

- `Integer::literal()` moves the existing `to_decimal()` buffer into
  the standard-string base of its returned `String`, then inserts
  the sign in that same owning string. The old conditional expression
  and `String(std::string const&)` copied the magnitude buffer. This
  retains the exact existing `String` return type and requires no
  utility submodule change. Negative zero still formats as `0`.
- The private signed-magnitude constructor takes `BigUnsigned&&`.
  The public `from_magnitude(BigUnsigned,bool)` signature and layout
  stay the same; its already-owned parameter is moved into the member
  once, instead of through another by-value parameter first. Lvalue
  callers still copy their source; moved-from magnitudes remain zero.
- `from_decimal()` uses one generic block recurrence for both storage
  strategies. At most 77 input digits are accumulated directly into
  the result's four inline limbs: 10^77 < 2^256 bounds every prefix.
  Longer inputs accumulate in the result's own vector, reserved to
  ceil(input_digits/19), which is safe because 10^19 < 2^64. This
  avoids the final second allocation/copy. Long inputs are not
  zero-initialized across their reserved capacity. Leading-zero
  inputs producing at most four limbs are copied into inline storage
  and release their unused vector allocation. The existing 32-digit
  threshold, 19-digit wide recurrence, portable 9-digit recurrence
  and full ASCII validation before parsing remain unchanged. Appends
  occur only after each limb scan, so vector growth/span extension
  cannot invalidate an active scan. The carry append invariant keeps
  stored limbs canonical, including an empty magnitude for zero.
- The radix-10^9 formatter uses local arrays for up to eight source
  limbs and 24 decimal blocks, eliminating its two temporary vector
  allocations at short widths. Each limb contributes fewer than
  three 9-digit blocks. Larger portable leaves retain checked dynamic
  storage. Wide/recursive decimal algorithms and their dispatches,
  multiplication/division kernels, fused/carry selection and B72's
  product result construction are unchanged.

These are shared C++20 ownership improvements on ARM64, x86-64 and
portable targets, not a replacement arithmetic engine. No public
operation, signed division convention, object layout, benchmark
fixture, official input pool, operation weight or oracle is changed.
The private constructor's signature changes; no binary compatibility
for that private symbol is promised. No GMP is added to production.

Validation on Apple M2 Ultra, macOS 26.6.2, AppleClang 21 Release
with warnings-as-errors: all 20 CTest tests pass. The complete native
`test_integer` and `test_big_unsigned` suites also pass a Debug build
with AddressSanitizer and UndefinedBehaviorSanitizer. The complete
BigUnsigned suite compiles with GCC 16.2, warnings-as-errors, native
ARM64 kernels and passes. Additional tests cover signed owning
results at 0/1/4/5/8/9/64 limbs, source preservation, moved-from zero,
negative-zero normalization, retained decimal strings and stream IO;
GMP parsing and formatting at 1/19/31/32/76/77/78/95 digits; and
invalid ASCII/NUL at every position around those boundaries. Existing
randomized and adversarial GMP tests remain intact. A separate forced-portable decimal harness (no __SIZEOF_INT128__,
no native backend) passes 3923 independent GMP/malformed-input checks
under ASan/UBSan, including up to 10000 input digits and the 8/9-limb
formatter boundary. Its compilation suppresses the pre-existing
unused fused-helper warning in this artificially forced configuration;
the normal native builds retain warnings-as-errors. This is not a
native MSVC build. The unchanged 1000x11 four-size FloatAP benchmark
also verifies all 256 paired inputs and 3456 operation results; no
FloatAP performance comparison to B72 was made. Windows/MSVC and
native GCC x86-64 CI remain pending.

Unchanged official 20,000-iteration, 31-sample Integer benchmark,
all six widths, independent GMP verification PASS (768 paired inputs,
2688 paired results) before and in both repeated final candidate runs:

| bits | ARM64 B72 total/GMP | B73 run 1 total/GMP | run 2 total/GMP | B73 worst run 1 / run 2 |
| --- | --- | --- | --- | --- |
| 53 | 0.551x | 0.404x | 0.441x | 0.991x / 1.154x |
| 128 | 0.937x | 0.687x | 0.685x | 1.518x / 1.496x |
| 256 | 1.032x | 0.817x | 0.816x | 1.706x / 1.731x |
| 1024 | 0.882x | 0.847x | 0.847x | 1.462x / 1.433x |
| 2048 | 1.018x | 1.001x | 1.002x | 1.551x / 1.581x |
| 4096 | 1.034x | 0.995x | 0.994x | 1.553x / 1.597x |

Native short conversion nanoseconds (B72 -> B73 two-run range):
53-bit parse 49.110 -> 22.085..22.090, format
74.612 -> 36.923..37.090; 128-bit parse
97.496 -> 62.460..63.244, format
223.083 -> 117.525..119.360; 256-bit parse
127.144 -> 93.123..93.796, format
285.700 -> 170.333..171.842. These are ARM64 observations only;
they do not project x86-64 timings or isolate each change's causal
contribution. The baseline was a single run, candidate two repeats;
host drift exists. Parsing at 2048 bits was about 7% slower locally
(1010.794 -> 1080.692/1084.123 ns) despite lower per-size totals.
Do not describe this as a per-operation regression-free change.

**The requested target has not been achieved or certified by B73.**
Even the ARM64 2048-bit total remains above 1x. Native GCC x86-64
measurements after B72+B73 are required before accepting performance:
each of the six totals must be strictly <1x and every operation/size
ratio strictly <2x, on stable repeated runs on the same machine.
A sum below 1x across all sizes is insufficient. Further work should
follow those post-B73 measurements: short owning square/pow and
nontrivial division remain candidates, with no kernel/threshold
promotion inferred from the old profiles alone.

From the x86-64 Release build (one line):
  ./benchmarks/benchmark_compare_integer_backends -t dark --bits 53,128,256,1024,2048,4096 --iterations 20000 --samples 31


### B74 short owning arithmetic and generic buffered powers (20 October 2026)

The user's post-B73 x86-64 20000x31 official run verifies all 768
paired inputs and 2688 results. Per-size total/GMP is
1.009x/1.250x/1.256x/0.992x/0.893x/0.915x at
53/128/256/1024/2048/4096 bits; worst ratios are
2.003x (pow3), 2.125x (pow3), 2.025x (quotient),
1.982x (quotient), 1.918x (compare), 1.943x (remainder).
The remaining native-sum excesses are 4.352/122.464/178.269 ns
at the three short sizes. GMP's 53-bit timings differ markedly
from the earlier attachment, so do not attribute the entire ratio
improvement to one change or infer stable acceptance from one run.

B74 targets owning intermediate work, retaining the existing kernels:

- One-limb unsigned sums return a scalar magnitude, with a separate
  two-limb carry result when needed. Subtraction whose significant
  range shrinks to one limb returns that exact scalar difference.
  Canonical size, not stale moved-from inline storage, decides zeros.
- Owning `quot`/`rem` with at most four dividend limbs use a private
  shared template. Shorter dividends return q=0/r=a; scalar inputs
  use exact unsigned division/remainder. Equal-width q=0/1 cases
  use B36's existing high-limb criterion before scratch preparation.
  The other cases call the UNCHANGED selected division template with
  local output of four limbs and scratch of nine limbs (n+m+1<=9).
  Zero-divisor checks still precede every dispatch. Joint divrem,
  low-level validation/scratch contracts, >4-limb owning division,
  normalization, preinverse thresholds and ISA kernels are unchanged.
  These storage/scalar paths are shared C++20 on all architectures.
- A new exact `pow(BigUnsigned const&,uint64_t)` uses generic
  MSB-first exponentiation. If base.bit_length()*exponent<=1024,
  two alternating local buffers replace owning intermediate results.
  The bound is checked by division, avoiding overflow; zero, exponent
  zero (including 0^0=1), and exponent one are handled separately.
  Every intermediate exponent is at most the requested exponent.
  Although each canonical intermediate fits sixteen limbs,
  ceil(x/64)+ceil(y/64) may exceed ceil((x+y)/64) by one. Each
  buffer therefore has SEVENTEEN limbs (272 bytes in total), allowing
  the full kernel output before normalization. The input is borrowed
  immutably, buffers/output are disjoint, and every selected kernel
  fully overwrites its output before it is read. No uninitialized
  limb is read and no stack view escapes: only the final BigUnsigned
  owns/copies storage. The >=8-limb symmetric-square selection and
  selected fused/carry/portable multiplication are preserved. Above
  the conservative budget, ordinary owning MSB-first square/multiply
  remains arbitrary-size. The budget is a bounded-storage policy,
  not a measured optimal crossover or an exponent-three special case.
- `Integer::pow` calls this unsigned engine only under the existing
  GCC x86-64 fused-build guard, keeping the same sign parity and
  canonical zero. Other signed front ends retain their old schedule.
  No public existing signature/layout or division convention changes;
  the unsigned pow overload is new. The existing standalone small
  power profile adds a GMP-verified owning SIGNED buffered-power row
  on the original input pool, for direct controls on other targets.

Validation on Apple M2 Ultra / AppleClang 21 Release warnings-as-errors:
20/20 CTest pass; complete Integer and BigUnsigned Debug ASan/UBSan
suites pass; complete GCC 16.2 native ARM64 BigUnsigned suite passes.
Forced-portable (no native backend, __SIZEOF_INT128__ undefined)
ASan/UBSan runs the actual power, short-division, arithmetic-boundary
and common-prefix subtraction tests successfully. The forced build
suppresses only the pre-existing unused fused-helper warning; normal
native builds keep warnings-as-errors. No MSVC/x86 native build is
claimed. New independent checks include 1497 GMP power cases plus
fixed binary/scalar expectations, zero/unit maximum uint64 exponents,
53/64/65/127/128/129/255/256/257/341/342/511/512/513/1023/1024/1025-bit
base boundaries and both sides of the storage budget; 864 division
operand cases at 0..5 limbs (including zero-divisor exceptions),
source preservation, scalar carry, and moved-from stale inline data.
The 65/341-bit dense powers explicitly exercise the seventeenth limb.
All earlier independent and randomized tests remain intact.

Both repeated unchanged ARM64 official 20000x31 runs verify 768
paired inputs and 2688 results. Signed ARM64 power is intentionally
unchanged; the separate same-signed profile verifies buffered powers
before/after timing. Its MSB-control -> buffered ns are
27.156 -> 27.654 (53), 61.435 -> 56.331 (128),
90.904 -> 56.292 (256): neutral at 53, modest at 128, larger at 256.
These ARM64 controls do not project native GCC x86-64 gains and are
not additive components of the official total.

ARM64 total/GMP reference -> two candidate runs:
53: 0.408 -> 0.343/0.374; 128: 0.676 -> 0.645/0.639;
256: 0.809 -> 0.760/0.787; 1024: 0.846 -> 0.844/0.856;
2048: 1.009 -> 1.001/1.002; 4096: 0.963 -> 1.004/1.000.
2048/4096 still fail strict ARM64 totals in these runs. Native
4096-bit summed time nevertheless falls from 11890.528 to
11664.503/11673.024 ns; GMP's reference sum changes from
12345.094 to 11623.273/11669.685 ns. Record absolute times and
repeat on the target rather than claiming either unconditional
parity or a native slowdown from ratios alone.

The unchanged FloatAP 1000x11 four-size benchmark verifies 256
paired inputs and 3456 results for both the saved B73 library and
B74. Its 53-bit relative total changes 1.460x -> 1.739x while native
operation times generally improve and many MPFR times fall sharply
(e.g. exp 1234.709 -> 693.625 ns). These runs do not establish
FloatAP performance acceptance; its pre-existing short-size limits
remain unmet. No floating-point arithmetic/rounding algorithm changed.

B74 is committed directly on remove-gmp per user instruction.
Native GCC x86-64 and Windows CI/performance remain pending. The
binding acceptance is still all SIX official per-size totals <1x
and EVERY operation/size <2x, with stable same-host repeated runs.
No benchmark corpus, weights or oracle was altered. Run from build:
  ./benchmarks/benchmark_compare_integer_backends -t dark --bits 53,128,256,1024,2048,4096 --iterations 20000 --samples 31


### B75 expose small ownership/adaptation bodies to the optimizer (20 October 2026)

The user's B74 GCC x86-64 official 20000x31 run verifies all 768
paired inputs and 2688 operation results. All 84 operation/size
ratios are now below 2x: per-size worst is
1.750/1.856/1.795/1.989/1.942/1.964x. Total/GMP is
0.750/1.091/1.120/0.942/0.854/0.900x at the six widths.
Only 128 and 256 bits still fail total parity: native sums
532.743/776.424 ns versus GMP 488.214/693.014 ns. The remaining
44.529/83.410 ns require about 8.36%/10.74% native reductions.
The overall 0.901x sum does not substitute for these two gates.

B75 makes existing small bodies visible in the class headers:
BigUnsigned's scalar constructor, move constructor/assignment,
limbs() and limb(); Integer's signed-magnitude constructor,
from_magnitude(), destructor, copy constructor/assignment and move
constructor/assignment. Their arithmetic, field transfers, self-move
guards, moved-from canonical zero and negative-zero normalization
are unchanged. The former out-of-line definitions are removed;
there is one ODR-valid inline definition per method. Integer's
default constructor remains out of line. No arithmetic kernel,
unsigned power/division schedule, storage layout, capacity, API
signature, exception specification, ISA dispatch, benchmark fixture,
keep-alive barrier, operation weight or numerical oracle is changed.
No LTO, compiler flag, always-inline attribute or platform dependency.

The purpose is to let each caller optimize ownership transfers and
eliminate redundant call/temporary work, including destruction of
moved-from temporary magnitudes. It is not proof that a specific
GCC ABI/linker rule caused the prior cost: B70's noisy controls did
not establish that. An inline definition permits optimization; it
does not guarantee a particular machine-code sequence or native
x86-64 improvement. The bodies are shared across all toolchains.

**Rebuild the library AND its clients/benchmarks together.** Moving
formerly out-of-line special members to inline definitions can change
which standalone symbols the shared library emits. This iteration
preserves source-level mathematical/lifetime contracts and object
layout; it does not promise linking old precompiled clients against
the new library. The benchmark must be rebuilt normally, not altered
to bypass construction/destruction. From the build directory:
  cmake --build . --parallel --target benchmark_compare_integer_backends

Local Apple M2 Ultra, macOS 26.6.2, AppleClang 21 Release
warnings-as-errors: the complete C++ `everything` target rebuilds,
and 20/20 CTest tests pass. Complete Integer/BigUnsigned Debug
ASan/UBSan suites pass. The complete BigUnsigned suite also passes
GCC 16.2 with native ARM64 kernels; actual storage/boundary/power/
short-division tests pass forced-portable ASan/UBSan without native
kernels or __int128 (only its pre-existing unused fused-helper warning
is suppressed). No tests or assertions were changed. Existing B73
ownership tests and GMP fixed/random arithmetic oracles are retained.
The Integer translation unit also compiles with GCC 16.2 using the
project warning set (including its normal non-template-friend
suppression) and warnings-as-errors. These are ARM64/portable checks,
not native GCC x86-64 or MSVC CI.

The unchanged 20000x31 ARM64 Integer benchmark verifies 768 paired
inputs and 2688 results before and in both final candidate repeats.
A saved B74 executable AND library form the reference, while all
candidate clients are rebuilt with the candidate headers. Native
128-bit sums: 354.681 -> 323.478/335.364 ns; 256-bit sums:
564.974 -> 527.437/550.217 ns. These local reductions do not
predict the required x86-64 improvement or attribute additive costs.

ARM64 per-size total/GMP B74 -> B75 two repeats:
53: 0.358 -> 0.293/0.319; 128: 0.632 -> 0.598/0.598;
256: 0.757 -> 0.738/0.747; 1024: 0.850 -> 0.837/0.843;
2048: 0.998 -> 0.973/0.969; 4096: 0.980 -> 1.001/1.005.
4096-bit native sums are 11786.017 -> 11659.638/11930.548 ns,
with GMP sums 12026.676 -> 11643.269/11873.572 ns. Thus local
ARM64 4096 total acceptance is not established; neither a blanket
nonregression nor target parity is claimed from these mixed runs.

Both matching B74 and rebuilt B75 FloatAP executables/libraries pass
the unchanged 1000x11 four-size numerical check (256 paired inputs,
3456 results). Candidate total/worst at 53/128/256/1024 bits are
1.513/2.764x, 1.070/2.620x, 0.588/2.047x, 0.246/1.719x;
reference is 1.600/4.317x, 1.021/2.975x, 0.586/2.071x,
0.254/1.865x. Existing short-width FloatAP gates remain unmet;
these one-run measurements do not establish performance acceptance.
No FloatAP arithmetic or rounding algorithm was changed.

B75 is committed directly on remove-gmp. Native GCC x86-64 and
MSVC validation, repeated six-size Integer totals strictly <1x and
all 84 ratios strictly <2x remain pending. After rebuilding, run:
  ./benchmarks/benchmark_compare_integer_backends -t dark --bits 53,128,256,1024,2048,4096 --iterations 20000 --samples 31


### R1 — native finite and extended Rational (20 October 2026)

The user accepted Integer as complete after the B75 x86-64 run:
per-size total/worst ratios at 53/128/256/1024/2048/4096 bits were
0.535/1.296, 0.898/1.555, 0.987/1.597, 0.951/1.928,
0.887/1.881 and 0.914/1.930x. Numerical verification passed
768 paired inputs and 2688 results. R1 starts the explicitly requested
Rational migration on remove-gmp; Integer and BigUnsigned are unchanged.

**Representation and invariants.** Rational now privately owns an Integer
numerator, a BigUnsigned denominator and an explicit finite/positive
infinity/negative infinity/NaN kind. Every finite value is reduced, its
denominator is positive and zero is 0/1. Special accessors retain the
mathematically useful 0/0 or +/-1/0 conventions, but classification uses
the explicit kind, never a GMP field. Public _mpq, get_mpq(), the mpq_t
constructor and the unrelated raw-mpz stream helper are removed.
This changes source interfaces for GMP-specific clients and the binary
layout; rebuild the library and all clients together. Ordinary rational
constructors, arithmetic, numerator/denominator accessors, streams and
number/Python integration retain their mathematical interfaces. Builtin
unsigned construction now uses Integer directly, including UINT64_MAX.

**Finite arithmetic.** Addition/subtraction first computes g=gcd(b,d)
for a/b and c/d, forms n=a*(d/g)+/-c*(b/g), and reduces only by
h=gcd(abs(n),g). The denominator is (b/g)*(d/h); no other common
factor can remain because both inputs were reduced. Multiplication
cancels gcd(abs(a),d) and gcd(abs(c),b) before either product; division
cancels numerator and denominator pairs before multiplying and handles
the resulting sign explicitly. Same-object multiplication and sqr use
squares of the already-coprime components. Halving cancels a factor of
two from the numerator or shifts the denominator. Reciprocal exchanges
coprime magnitudes. A private Reduced constructor is used only where
these operations prove coprimality; it still canonicalizes zero.
Powers use shared unsigned powering without repeated fraction gcds;
negative exponents, including Int minimum, do not negate an unrepresentable
signed minimum. Exact fma is implemented as exact multiply then add.
Comparisons use canonical equality or exact cross-products. Floor, ceil
and half-away-from-zero round share one unsigned quotient/remainder.

**Extended semantics and explicit decisions.** Existing ExtendedOperations
continues to dispatch finite versus special arithmetic, comparisons,
min/max, sign, abs and reciprocal. NaN remains unequal to itself and
named cmp returns INCOMPARABLE. No Rational signed zero is introduced.
Finite division by zero is explicitly defined as signed infinity for
nonzero numerators and NaN for 0/0; reciprocal zero is positive infinity.
These cases no longer call GMP on an invalid denominator. Decimal text
accepts signed integers or n/d, ignoring ASCII whitespace, and raises
runtime_error on malformed input; textual denominator zero has the same
explicit extended meaning as the integer-pair constructor. Equal infinities
in Rational/ExactDouble comparisons now compare equal rather than taking
the old unconditional less/greater branch; NaN remains incomparable.
FloatMP mixed comparisons also classify special values before any GMP call.

**Exact binary64 conversion.** ExactDouble import decomposes a finite
binary64 into an integer significand and power of two and cancels trailing
zero bits. get_d computes the binary exponent with exact comparisons,
then truncates an exact scaled integer quotient to at most 53 bits.
Subnormal spacing is 2^-1074; overflow yields signed infinity. Underflow
returns +0, including negative tiny fractions, matching the old GMP
conversion rather than introducing an observable negative zero. Results
are independent of the active host rounding mode. No GMP or floating
approximation is used to decide rounding. The public rational header
passes an include-only compilation without GMP/MPFR headers, and its
arithmetic object has no undefined GMP/MPFR symbols.

**Remaining integration boundaries.** detail/rational_gmp_bridge.hpp is
an independently owned RAII adapter, analogous to the Integer boundary,
for still-GMP-backed Dyadic/FloatAP and MPFR FloatMP only. It is never a
Rational member, arithmetic backend or cache. Export rejects nonfinite
values rather than passing denominator zero to GMP. Dyadic conversion
is implemented beside Dyadic's existing GMP storage; FloatMP conversion
and finite mixed comparison explicitly use temporary exact exports.
FloatAP's Rational conversion imports its significand through the existing
Integer boundary and uses native shifts and normalization; Rational-to-AP
and AP decimal formatting retain explicit transitional GMP exports.
FloatAP's internal GMP significand, Fraction class and fixed-point arithmetic
remain separate future migration work. The production library as a whole
therefore still links GMP/MPFR: R1 does not claim complete link independence.
There are no remaining _mpq/get_mpq call sites in numeric sources, public
headers, tests, Python or benchmark callers. Direct GMP test/benchmark
oracles and FloatAP's independent internal mpq operations are retained.

**Verification and measurements.** After the R1a test-only corrections
documented below, the complete AppleClang 21 Release warnings-as-errors
everything build and all 20 C++ CTest suites pass on
Apple M2 Ultra, macOS 26.6.2. New TestRational cases cover fixed independent
expected fractions, large cancellation, zero/special cases, parsing,
copy/move/self-aliasing, signed minimum exponents, 2000 deterministic
binary64 round trips, 1200 independent random GMP-oracle pairs up to
17 limbs, canonical result gcd/sign checks, floor/ceil and bit-exact
binary64 differential checks at normal/subnormal/overflow boundaries
under all four host rounding modes. Integer, BigUnsigned, Rational,
Float and FloatRounding Debug ASan/UBSan suites pass (5/5).
GCC 16.2 ARM64 passes the native Rational/Integer/BigUnsigned suites
with normal project warnings-as-errors. A complete default GCC build
encountered -Wdisabled-optimization in the existing very large
TestFloat/TestFloatRounding template functions; full-build validation
passes everything and all 20 CTest suites with the larger local GCC
--param=max-gcse-memory=1000000 budget, without changing project CMake
rules or suppressing warnings. Native x86-64 and MSVC are not
available on this host and are not claimed as tested.

The unchanged six-size Rational benchmark (20000 iterations, 31 samples)
passes 384 paired inputs and 2880 operation results. Initial ARM64
per-size total/worst ratios are 1.198/2.211, 1.203/5.643,
2.562/3.567, 2.005/3.066, 1.351/3.435, 0.991/3.890x.
These measurements overlapped compiler/test activity and are diagnostic
only, not performance acceptance or an x86-64 prediction. Rational normalization/GCD, comparison and division
remain optimization work. No fixture, operation weight, GMP oracle,
timing barrier or compiler/ISA selector was changed.

Both saved B75 and rebuilt R1 FloatAP benchmarks pass the unchanged
1000x11 numerical check (256 inputs, 3456 results). Initial reference
versus candidate total/worst ratios at 53/128/256/1024 bits are
1.759/3.953 -> 1.694/4.063, 1.067/2.859 -> 1.068/2.750,
0.600/2.079 -> 0.581/2.140, 0.245/1.628 -> 0.250/1.735x.
These single runs do not establish performance acceptance or blanket
nonregression; the existing short-width AP gates are still unmet.
No FloatAP arithmetic, precision or rounding algorithm was changed.

Rebuild library and clients, then measure native x86-64 from the build directory:
  cmake --build . --parallel --target benchmark_compare_rational_backends
  ./benchmarks/benchmark_compare_rational_backends -t dark --bits 53,128,256,1024,2048,4096 --iterations 20000 --samples 31


### R1a — compile the complete native Rational regression suite (20 October 2026)

The first CI builds of R1 exposed two test-only Clang errors: a local
FloatMP variable named infinity shadowed Ariadne::infinity under -Wshadow,
and the bare GmpRationalBoundary(Rational::nan()) inside ARIADNE_TEST_FAIL
could be parsed as a declaration instead of an expression. Rename the
local and explicitly cast the temporary boundary expression to void.
No Rational arithmetic, representation, exception expectation or oracle
changes. Earlier incremental CTest output had used the previous executable
after a failed rebuild; it is superseded by the following validation,
which checks a successful build exit before starting CTest.

AppleClang 21 Release warnings-as-errors everything including Python
bindings succeeds; all 21 CTest entries pass (20 C++ suites plus the
Python numeric suite with Python 3.14.7). GCC 16.2 ARM64 everything with
the documented local optimization-memory budget succeeds and all 20 C++
suites pass. Rebuilt Debug ASan/UBSan Rational/Integer/BigUnsigned/Float/
FloatRounding suites pass (5/5), including the corrected new tests.
Native x86-64 and MSVC still require remote CI/runtime confirmation.


### R1b — retain division alias coverage on Clang 16 (20 October 2026)

Clang 16 rejects the intentional destination/=destination test under
-Wself-assign-overloaded. Use the existing reference alias as the right
operand, preserving exactly the same object-aliasing case and expected
value without suppressing any warning. AppleClang 21 rebuild and the
complete Rational suite pass. The previous R1a CI already passed macOS
Release/Coverage, Ubuntu GCC 13 x86-64 Release/Coverage and Windows 2025
Release; Clang 16 Release/ThreadSanitizer requires this test spelling fix.
The optional local Clang 15 syntax experiment could not compile against
the installed newer macOS SDK/libc++ combination and is not validation.
Rational arithmetic and Integer/BigUnsigned remain unchanged.

The separate Codecov status thresholds already failed on B75 (project
92.96%, patch 92.56%, target 99.20%). R1a reports project 92.94% and
patch 93.07%; workflow Coverage success is not Codecov gate acceptance.


### R2 — remainder-only Euclidean GCD (20 October 2026)

The initial user-supplied six-size Rational/mpq benchmarks (20,000
iterations and 31 samples) passed 384 paired inputs and 2,880 paired
operation checks on each architecture. Aggregate native/reference
total-time ratios were 1.265x on ARM64 and 1.329x on x86-64;
the most expensive regressions in overall per-size total time were
at 256 bits (2.631x and 3.327x respectively). These values are the
**pre-R2 baselines**, not post-change measurements.

The shared `BigUnsigned::gcd` Euclidean loop now calls the existing
remainder-only `rem(a,b)`, rather than constructing a quotient/remainder
pair through `divrem(a,b)` and discarding the quotient at every
iteration. It preserves the same Euclidean recurrence and edge cases,
including gcd(0,0)=0. It uses the existing tested short and general
remainder kernels and changes no Rational normalization, API,
arithmetic policy, benchmark operand, GMP oracle, or architecture
dispatch. The single shared-core change may also affect other users of
`gcd`; benchmark and functional gates for those users remain applicable.

**Validation status:** no local build or benchmark was available in
the editing environment. No speedup, CI outcome or performance gate
acceptance is claimed for R2 until the updated branch is rebuilt and
tested on ARM64 and x86-64. Re-run the existing TestBigUnsigned,
TestInteger and TestRational differential/alias tests and the unchanged
`benchmark_compare_rational_backends` workload at
53,128,256,1024,2048,4096 bits with 20,000 iterations and 31 samples;
also check the FloatAP regression/performance gates before accepting
this shared-core optimization.


### R3 — bounded stack products for Rational comparison (20 October 2026)

**Measured R2 baseline, not an R3 result.** User-supplied ARM64 and
x86-64 six-size `benchmark_compare_rational_backends` runs both pass
384 paired inputs and 2,880 paired results. Compared with the pre-R2
aggregate totals, ARM64 decreases from 1.265x to 1.073x and x86-64
from 1.329x to 1.021x. At 256 bits, the ARM64/native-to-GMP total
ratio is 1.870x (was 2.631x), and the x86-64 ratio is 2.367x (was
3.327x). The R2 Unix, Windows and Coverage GitHub Actions complete
successfully. Independent native ARM64 and x86-64 R3 measurements
remain pending.

**Change.** The finite `Rational::cmp` cross-product comparison now
uses caller-owned stack limb arrays through the already-selected,
portable-fallback `BigUnsigned::mul_into` kernel when each of the
four unsigned input magnitudes occupies at most eight 64-bit limbs.
Both exact products are compared from their high words, with a
zero-initialized tail to handle different product lengths and any
zero top limb. Above the eight-limb bound, the previous owning
`mul` plus `compare` path is unchanged. Both paths preserve the
existing sign and same-denominator shortcuts, the exact comparison
order and IEEE/extended-special-value dispatch. No approximate
floating comparison, GMP operation, conditional compilation or
architecture-specific instruction is introduced. This specifically
avoids two temporary heap-owned cross-products for the medium-size
input range, while retaining the existing multiplication kernels.

**Regression.** Add deterministic GMP-oracle comparisons around 53,
128, 256, 448 (eight-limb stack boundary), 512 (fallback), and 1024
bits, testing operand order, both-negative reversals and zero cases.
Retain the previously existing randomized GMP differential, aliasing,
NaN/infinity and rounding-mode coverage unchanged.

**Validation status.** R3 has not been locally built or benchmarked:
no native ARM64 or x86-64 execution environment was available in the
repository-editing session. Do not interpret the R2 measurements as
R3 speedup or accepted performance gates. Wait for CI and re-run the
unaltered six-size benchmark on both hosts, plus the full Rational,
Integer, BigUnsigned and FloatAP gates before evaluating another
candidate.


### R4 — proper-fraction floor/ceil shortcuts and strict Rational performance gates (20 October 2026)

**Acceptance policy (user requested):** for **each** of the six widths
53, 128, 256, 1024, 2048, 4096 bits, on **each** of ARM64 and x86-64,
require both the `total_time` ratio strictly below 1.0x and the
`worst` operation ratio strictly below 2.0x versus the directly
measured `mpq_t` backend. The geomean is diagnostic and cannot
substitute for either gate. Keep the 15 operations, 32-pair workloads,
20,000 iterations, 31 samples, numerical oracle and measurement
methodology unchanged. An architecture passes only if all six sizes
pass both gates; do not average away a failure.

**R3 baseline from user-provided runs (before R4):** ARM64 total/worst
per width: 53 0.763/2.091x, 128 0.930/1.812x,
256 1.840/2.432x, 1024 1.560/2.725x,
2048 1.149/3.547x, 4096 0.881/3.945x.
x86-64: 53 1.396/4.066x, 128 1.380/2.574x,
256 2.372/2.931x, 1024 1.695/3.176x,
2048 1.044/3.063x, 4096 0.730/3.889x.
All 384 paired inputs and 2880 paired operation results pass numerical
verification on each architecture. Unix, Windows and Coverage workflow
runs on commit `012b336` succeeded. Only ARM64 128-bit passes both
user gates at R3. This is not an R4 measurement.

**R4 change:** finite `floor` and `ceil` now check the exact unsigned
magnitude order before allocating quotient/remainder storage. For
`|n|<d`, floor gives -1 for negative nonzero and 0 otherwise; ceil
gives 1 for positive nonzero and 0 otherwise. Canonical zero is 0/1.
The existing Euclidean quotient/remainder path is entirely unchanged
when `|n|>=d`. This optimizes mathematically general proper fractions,
not a benchmark-specific shape; it preserves signed results and the
precondition that special values are not valid inputs.

**R4 tests:** deterministic GMP floor- and ceiling-quotient comparisons
for zero, small and large positive and negative proper fractions,
whole numbers and improper fractions at all six widths, including the
4-to-5-limb inline-storage boundary. Existing random GMP differential
and specialized floating-point rounding tests remain intact.

**Validation status:** R4 has not been executed on ARM64 or x86-64 in
this editing session. Its build/CI and six-size performance outcome are
pending; no post-R4 gate is claimed as passing. Afterwards, revisit
`Rational::get_d` (still >2x at 53 bits on both architectures, >4x on
x86-64) and the medium-width construction and four operations that
still dominate the `total_time` failures.


### R5 — exact bounded conversion to binary64 (20 October 2026)

**Acceptance:** preserve R4's per-architecture/per-size strict requirements:
all 6 widths must have total_time <1x and worst <2x relative to
direct mpq_t, with 20,000 iterations, 31 samples and unchanged operation
weights/fixtures. The independent user-provided R4 measurements (ARM64
then x86-64) pass 384 paired inputs and 2,880 operation results each.
Total/worst: ARM64 53 0.723/2.130, 128 0.927/1.890,
256 1.752/2.384, 1024 1.496/2.340, 2048 1.128/2.547,
4096 0.879/2.576x. x86-64 53 1.268/4.009, 128 1.340/2.598,
256 2.286/2.936, 1024 1.642/2.431, 2048 1.014/2.478,
4096 0.713/2.577x. R4 Unix, Windows and Coverage CI pass.

**R5 finite Rational::get_d paths:** for one-limb numerator and
denominator, compute the exact scaled binary64 significand in an
unsigned 128-bit intermediate on GCC/Clang with __int128 support
(including x86-64 and ARM64), avoiding owning big-integer division.
On other platforms, preserve the established portable fallback.

For same-limb-count operands with at least five 64-bit limbs,
normal binary64 values and negative scale of at most 63, first divide
only the top three limbs, shifted by the 53-bit/significand scale.
With n=N*B+nl and d=D*B+dl, B a shared power of two, and
N*2^s=q*D+r, 0<=nl,dl<B, the bounds r>=q and r+2^s<D imply
floor(n*2^s/d)=q. When q=2^52 on normal numbers, the exact
binary exponent already guarantees the lower bound. If neither
proof applies, or the inputs lie outside the selected range,
take the unchanged full-width exact division. This is a generic
exact-prefix certificate, not an approximation or a fixture-specific
shortcut: all results must remain bit-for-bit equal to mpq_get_d,
including signed zeros and all four host rounding modes.

**R5 tests:** extend the existing four-rounding-mode GMP differential
checks with 53,128,256,448,512,1024,2048,4096-bit values,
near-half and near-one ratios, nonzero low tails, negative signs,
and overflow. Keep all earlier randomized tests and benchmark
workload/oracle unchanged.

**Validation:** commit changes without claiming performance or CI
outcomes for R5. Rebuild and run all Unix/Windows/Coverage workflows,
Rational/Integer/BigUnsigned/FloatAP differential gates and unchanged
ARM64/x86-64 six-size benchmarks before considering R5 accepted.


### R6 — exact equal-leading-limb Euclidean subtraction (20 October 2026)

**Context:** R5 benchmark (20,000 iterations, 31 samples, six widths; no fixture or GMP changes) passed both oracle counts on ARM64 and x86-64, but still missed user gates total_time <1.0x and worst <2.0x at several widths. ARM64 total/worst ratios: 53 .705/1.686, 128 .932/1.860, 256 1.760/2.460, 1024 1.510/2.407, 2048 1.096/2.538, 4096 .859/2.569. x86-64: 53 1.210/2.237, 128 1.338/2.497, 256 2.289/3.072, 1024 1.645/2.321, 2048 1.008/2.122, 4096 .707/2.149. R5 Unix/Windows/Coverage CI succeeded.

**R6:** In the shared `BigUnsigned gcd` Euclidean reduction, when both operands occupy the same number (>1) of canonical limbs and their most significant limbs are equal, order them and replace the general remainder operation with the exact subtraction. For ordered a>=b with the same positive high limb, a<2*b, so a mod b=a-b. Other inputs keep the existing remainder-only algorithm unchanged. The input ordering is semantically irrelevant to gcd; zeros preserve the former behavior. This targets division/scratch overhead in a mathematically general case, including differences collapsing from heap to inline limbs at 256 bits, rather than checking any benchmark operand shape.

**Tests:** add deterministic GMP-oracle gcd checks for both input orders at 128/192/256/257/320/512/1024-bit boundaries and multiple low-limb offsets. Keep randomized GMP gcd, Rational arithmetic, FloatAP and the benchmark fixture unchanged.

**Validation:** R6 not locally compiled/benchmarked. R6 CI, cross-platform differential tests, and fresh ARM64/x86-64 baseline-matched performance runs are pending; no target is claimed passed on R6.


### R7 — certified native 128-bit high-prefix binary64 division (20 October 2026)

**R6 baseline (not R7 results):** the independently supplied ARM64
and x86-64 runs both pass 384 paired inputs and 2,880 paired operation
results; R6 Unix, Windows and Coverage workflows succeed.
Per-size total_time/worst ratios at 53,128,256,1024,2048,4096 bits:
ARM64 0.778/1.677, 1.002/1.817, 1.781/2.490,
1.517/2.233, 1.078/2.431, 0.856/2.459.
x86-64 1.275/2.334, 1.328/2.564, 2.208/3.023,
1.539/2.131, 0.954/1.874, 0.670/1.908.
Three of the 12 machine-width cases pass BOTH user acceptance gates;
the requirements stay strictly total_time <1x AND worst <2x at
EACH width and architecture separately, with no benchmark alteration.

**R7 change:** for normal finite binary64 conversion, with a
negative scale of at most 63 bits and equal limb counts, first compute
the high-prefix quotient using native unsigned 128-bit arithmetic.
Use one high limb, and on an ambiguous quotient use two high limbs.
For n=N*B+nl, d=D*B+dl with 0<=nl,dl<B, if
N*2^shift=q*D+r and r>=q, then the full quotient is >=q.
If r+2^shift<D, its full quotient is <q+1. When the
candidate equals 2^52, the exact known normal exponent independently
certifies that the full quotient is >=2^52. Apply both checks before
accepting the computed 53-bit significand. Check shift overflow,
candidate bit range, and prefix validity; otherwise use the existing
R5 three-limb certificate and ultimately the original full exact
division. Preserve the established active-rounding-mode independence,
truncation toward zero, signed zero, overflow/subnormal semantics and
compiler fallbacks. This is independent of particular rational values
and does not weaken correctness guarantees.

**R7 test extension:** in all four floating rounding modes, compare
bit-exact values against GMP for added positive and negative
near-half, near-third, repeating and near-boundary inputs across
the existing 53..4096-bit differential matrix. Leave the existing
random GMP tests, 15 benchmark operations, input pool and medians
unchanged.

**Validation status:** no local build or R7 native benchmark has been
run. The new commit requires Unix/Windows/Coverage, FloatAP/Integer/
BigUnsigned/Rational regressions, and unaltered six-size ARM64/x86-64
runs. Do not attribute measured R6 ratios to R7.


### R8 — platform-select exact conversion and scalar Euclid termination (20 October 2026)

**Acceptance:** unchanged six sizes x two architectures, with strictly
`total_time <1x` and `worst <2x` for EACH machine/size. The user
supplied R7 runs (20,000 iterations, 31 samples; 384 paired inputs,
2,880 operation results PASS per machine). ARM64 total/worst ratios:
53 .779/1.667; 128 1.034/3.987; 256 1.801/3.263;
1024 1.515/2.162; 2048 1.104/2.417; 4096 .872/2.477.
x86-64: 53 1.275/2.381; 128 1.325/2.204; 256 2.184/2.940;
1024 1.531/2.042; 2048 .945/1.908; 4096 .668/1.900.
R7 Unix, Windows and Coverage all pass. These are R7 data, not R8.

**R8 conversion selection:** R7's new unsigned 128/128 top-limb
certificate costs approximately 117 ns versus 36 ns for R6 at 128
bits on ARM64, and also regresses 256+ ARM64 conversion. On x86-64,
R7 materially reduces conversion cost at 128/256/1024/2048 bits.
Therefore enable only the R7 high-prefix unsigned 128-bit conversion
on x86-64 GCC/Clang; ARM64, MSVC and portable targets retain the
R5 three-limb/exact fallback (and the existing one-limb R5 path).
The numerical semantics and shared four-rounding-mode differential
gate remain identical. This selects algorithms by observed host
performance, not precision or benchmark input values.

**R8 fast one-limb comparison:** Rational::cmp now compares exact
unsigned 128-bit cross-products directly when all four numerator/
denominator magnitudes contain exactly one 64-bit limb, avoiding
stack-product arrays, dynamic multiplication dispatch and their
bookkeeping. Longer operands and toolchains without unsigned 128-bit
support retain R3's exact paths. Dedicated GMP comparisons use
near-64-bit-limit, positive and negative fractions.

**R8 scalar GCD termination:** BigUnsigned::gcd orders initial inputs
once, then maintains a>=b invariant throughout Euclid. When the
current divisor fits one 64-bit limb, perform at most one multiprecision
remainder, then finish Euclid with scalar '%' iterations. Special-case
divisor 1, preserving gcd(0,0), symmetry and all other cases. Keep
R6's exact equal-high-limb subtraction on multiprecision inputs.
Add 640 independent randomized scalar GMP comparisons in both orders,
plus widths around the 4/5-limb boundary and 1024/4096-bit dividends.
No changes to benchmark cases, GMP oracle, build flags or timed loops.

**Status:** R8 native CI and six-size ARM64/x86-64 benchmarks pending;
no R8 gate has been claimed passed. Shared-core changes require the
Integer/BigUnsigned/Rational and downstream FloatAP differential gates.


### R8a — fix scalar GCD limb type lookup (20 October 2026)

**Build failure:** Apple Clang ARM64 rejects four unqualified uses of
`Limb` inside the free function `BigUnsigned gcd(BigUnsigned,BigUnsigned)`:
`Limb` is a nested type of `BigUnsigned`, not visible unqualified in
this non-member scope (reported at big_unsigned.cpp around line 1261).

**Fix:** add `using Limb=BigUnsigned::Limb;` inside `gcd`. No algorithm,
branch condition, arithmetic result, benchmark fixture, test expectation,
or compiler flag is changed. The explicit alias resolves the same lookup
problem for Apple Clang, GCC and MSVC. R8a requires CI and fresh
measurements before claiming R8 performance; the R7 benchmark remains
the last measured baseline.


### R8b — repair mixed-width initializer list in GMP GCD test (20 October 2026)

**Blocking build regression:** R8a Unix/Windows/Coverage CI reports build
failures. Inspecting the R8 tests found a second compile-time error:
`for (std::uint64_t d : {1u,2u,3u,17u,65537u,
std::numeric_limits<std::uint64_t>::max()})` tries to deduce one
`std::initializer_list<auto>` from both `unsigned int` and
`std::uint64_t` elements. This is ill-formed under C++20 and was
reproduced with a minimal GCC compilation. Replace it with an explicit
`std::array<std::uint64_t,6u>` of the SAME six values. No GMP oracle,
scenario, numerical behavior, or production source is altered.
CI must be rechecked; the full failure logs were not available through
the repository connector, so other compile problems are not ruled out.


### R9 — one-source BigUnsigned compilation unit (20 October 2026)

**Motivation / user decision:** native architecture conditionals are not
fully separable into platform-specific source files. Division paths,
wide-integer scalar estimation, decimal formatting, Karatsuba and
multiplication dispatch already contain compiler/ISA branches inside
common implementation functions. Keeping external kernels in separate
`.cpp` files therefore creates a misleading partial modular boundary.

**Structural change only:** merge the reference C++20 multiplication
formerly in `big_unsigned_fallback.cpp`, the optional x86-64 kernel
formerly in `big_unsigned_x86-64.cpp`, and the ARM64 native carry
and symmetric-square kernels formerly in `big_unsigned_arm64.cpp`
into `src/big_unsigned.cpp`. Enclose the platform-specific
implementations in the SAME `ARIADNE_BIG_UNSIGNED_X86_64` and
`ARIADNE_BIG_UNSIGNED_ARM64` compile-time guards used before.
Preserve all internal public/diagnostic entry points, implementation
bodies, kernel selectors, GCC fused option, target compiler definitions
and portable reference implementation. Remove the three now-redundant
sources and their CMake source entries; CMake still selects the same
native backend flags based on compiler and target. No algorithm,
operand set, oracle, threshold, runtime condition, test or benchmark
is altered. Preserve historical source attribution in section comments
and Git history.

**Validation gates:** starting from R8b `aa5c520`, whose Unix, Windows
and Coverage CI all succeeded, confirm the refactored single-source
unit builds on GCC x86-64, Clang x86-64, AppleClang ARM64 and MSVC;
run the full differential BigUnsigned, Integer, Rational, FloatAP and
architecture-specific multiply diagnostics. Repeat the six-size,
15-operation GMP comparison on both ARM64 and x86-64 to check for
accidental performance regressions due to TU/codegen changes.
These validations are pending for R9; no improvement is claimed.


### R10 — borrowed initial Euclidean reduction for Rational normalization (20 October 2026)

**Baseline.** R9's single-source BigUnsigned refactor passes Unix,
Windows and Coverage CI. The last user-provided ARM64 and x86-64
R8a/R8b benchmark runs each pass numerical verification, but
`construct_fraction` remains relatively slow at medium/large sizes,
including 256-4096 bits. Continue enforcing per-size and per-host
`total_time <1x` AND `worst <2x` versus direct `mpq_t`.

**Change.** All Rational arithmetic GCD sites pass existing lvalue
magnitudes. Switch `BigUnsigned gcd` from by-value operands to
`const&`, and compute the first exact Euclidean remainder directly
from those borrowed operands, rather than copying both magnitudes.
Return the surviving borrowed input directly if the first remainder
is zero. For one-limb remainder, finish with the established R8
scalar Euclid without copying large inputs. Only if another
multi-limb reduction is necessary, copy the surviving divisor once
and thereafter reuse move-owned remainder buffers. Retain R6's
exact same-high-limb subtraction shortcut and the R2 remainder-only
division path. No benchmark inputs, GMP oracles, expectations, or
runtime thresholds are changed. This function signature changes its
C++ symbol/ABI; dependent binaries must be rebuilt from sources.

**Tests.** Add GMP differential checks across
53/128/255/256/257/512/1024/2048/4096 bits, with both orders,
zero, exact multiples, same-object aliasing, and verification that
borrowed operands are untouched. Keep randomized multi-limb GCD,
Rational construction and downstream FloatAP regressions unchanged.

**Validation.** No local R10 compilation or performance measurement.
Require Unix, Windows, Coverage CI, full Integer/BigUnsigned/Rational/
FloatAP test suites, and unchanged 20,000-iteration 31-sample
ARM64/x86-64 six-width benchmark runs before accepting the change.


### R11 — exact quotient-one/two Euclidean reduction (20 October 2026)

**Baseline.** The newly supplied post-R10 ARM64 and x86-64 runs
(20,000 iterations, 31 samples) each pass the 384 paired-input and
2,880 operation-result GMP checks. Per-width total_time/worst ratios
for 53/128/256/1024/2048/4096 bits are, on ARM64,
0.405/1.229, 0.821/1.750, 1.418/2.309, 1.282/1.663,
0.971/1.901, 0.796/2.063; on x86-64,
0.961/2.248, 1.188/2.299, 1.880/2.459, 1.365/1.621,
0.868/1.544, 0.626/1.586. Five of twelve width-host pairs
meet both strict gates. At the time of R11 creation, R10 Coverage
CI had passed while Unix/Windows CI were still running. No R11
speedup is claimed.

**Change.** Inside BigUnsigned::gcd, R6's equal-high-limb
subtraction shortcut generalizes to ordered positive operands
a>=b with equal canonical limb counts and highest limbs
high_a-high_b <= high_b (computed without overflow). With
B=2^(64*(limb_count-1)), high_b>=1 and
a < (high_a+1)*B <= (2*high_b+1)*B <=3*high_b*B <=3*b.
Thus the Euclidean quotient is either 1 or 2. Compute a-b
once. If high_a<2*high_b, a<2*b, and a-b is the remainder.
Only on high_a==2*high_b is another comparison needed:
if a-b>=b, subtract b once more. This computes the exact
remainder without Knuth division, scratch buffers, or
approximation. Cases outside the exact range still call the
unchanged remainder-only implementation. This applies equally
to dense arbitrary operands and sparse rational inputs,
not a workload signature. Keep R10's borrowed initial operands,
R8's scalar Euclid and all preexisting arithmetic semantics.

**Regression.** Add independent GMP-oracle GCD checks at
53/128/192/255/256/257/320/512/1024/2048/4096 bits
for quotient one, quotient two, exact 2*b, quotient three
fallback, equal operands, reverse order, and verification that
both lvalue arguments remain unchanged. Existing thousands
of randomized GMP differential tests remain untouched.

**Validation.** R11 has not been locally compiled or benchmarked.
Await Unix/Windows/Coverage CI and then unchanged ARM64 and
x86-64 15-operation, six-width native-vs-mpq_t benchmark runs.
Each host-width must satisfy BOTH strict `total_time <1x`
and `worst <2x`, independently; the geomean is diagnostic.


### R12 — exact scalar x86-64 Rational-to-double conversion (20 October 2026)

**Measured R11 baseline (not R12).** Both ARM64 and x86-64
20,000-iteration/31-sample user-supplied runs pass 384 paired inputs
and 2,880 paired operation results. ARM64 total/worst for
53/128/256/1024/2048/4096: .430/1.248, .788/1.791,
1.341/2.544, 1.185/1.398, .894/1.526, .725/1.637.
x86-64: .915/2.307, 1.111/2.147, 1.737/2.237,
1.249/1.464, .797/1.143, .575/1.118.
Six of twelve host/size cases pass both strict requirements
total_time <1 and worst <2. On x86-64 the worst at 53 bits is
to_double (2.307x); on ARM64 the worst at 256 bits is
to_double (2.544x). Construction and GCD reductions are now
substantially faster; 256/1024-bit arithmetic remains a
separate high-priority total_time bottleneck.

**Change.** On GCC/Clang x86-64 only, when both canonical Rational
magnitudes have exactly one 64-bit limb, compute the exact binary
exponent by integer bit-width comparison and bounded native shifts,
without creating BigUnsigned shifts/temporary quotients. The ratio
of two nonzero 64-bit integers is necessarily normal in binary64.
When exponent<=52, form the exact scaled numerator as unsigned
128-bit and divide its two 64-bit halves with baseline x86 DIVQ.
The independently computed exponent proves the quotient has at
most 53 bits, meeting DIVQ's high-word<divisor precondition; an
explicit check retains the portable full division fallback.
When exponent>52, a native 64-bit division by the exact shifted
denominator suffices. The returned integer significand is exactly
representable in double and ldexp preserves truncation toward
zero independent of the host rounding mode. ARM64, Windows/MSVC
and multi-limb inputs keep their previous algorithms unchanged.
No GMP runtime dependency, benchmark special-casing, workload or
threshold modification is introduced.

**Tests.** Add 512 deterministic scalar ratios in each of the
four host rounding modes, testing both signs and near-half,
full-width random, power-of-two and tiny numerator regimes.
Each comparison uses exact binary64 bits against mpq_get_d.
The preexisting randomized multiprecision GMP and special-value
tests are retained. Independently, an x86-64 GCC 14 local
standalone kernel test already compared 320,000 random and
boundary 64-bit ratios to GMP under four rounding modes,
without a mismatch. This verifies the scalar math but is not a
replacement for repository compilation and tests.

**Validation.** R12 native full CI and matched ARM64/x86-64
benchmarks have not run yet. Do not claim a measured improvement
or success on any remaining performance gate before the new data.


### R12a — fix C++ most-vexing-parse in scalar conversion tests (20 October 2026)

**Build failure:** R12 Unix, Coverage, and Windows jobs all failed at
Build. The R12 test statement `const Rational positive(Integer(numerator),Integer(denominator));`
is parsed as a function declaration, not an initialized Rational;
its subsequent unary negation fails to compile. Reproduced the
parsing diagnostic in a local standalone C++20 compiler check.
Replace with direct-list initialization using braces. The exact
same GMP comparison values, test count, and production code are
unchanged. Re-run all build and test gates before accepting R12a.


### R13 — certified small-multiprecision ARM64 binary64 conversion (20 October 2026)

**Motivation:** the user explicitly highlighted the unresolved ARM64
256-bit `to_double` regression. In the supplied R12 results, the
native median is 106.383 ns vs mpq_t 41.617 ns, a 2.556x slowdown;
meeting worst<2x requires <83.234 ns in this matched run.

**Change:** for normal binary64 quotients with negative or zero scale
and equal canonical 5..8-limb magnitudes on ARM64 GNU/Clang, derive
a tentative significand from the top limbs using ordinary floating
arithmetic. Clamp into the valid [2^52, 2^53-1] range. Independently
build the exact shifted numerator and exact divisor-times-candidate
products in fixed-length stack arrays. Only return the candidate if
B*q <= A*2^shift < B*(q+1) is proven by full multiword comparisons.
The test is independent of floating mode, architecture instruction
rounding, approximation quality, and lower limbs. An incorrect or
ambiguous candidate always falls through to the unchanged full-width
conversion; zero, signs, subnormals, overflows, unequal widths and
other targets remain unchanged. This is an exact certification strategy,
not an assumption about benchmark numerator/denominator patterns.

**Regression:** extend existing per-rounding-mode GMP binary64 tests
with 96 deterministic dense inputs at each 5..8-limb length, including
both signs, in addition to the prior 1200 randomized arithmetic pairs
and extensive 53..4096-bit rounding-edge checks. No benchmark inputs,
weights, median methodology or GMP oracle altered.

**Status:** R12a fixes a reproduced most-vexing-parse build regression
and is pending full CI; R13 remains unbuilt and unbenchmarked in this
session. The strict user acceptance criteria remain both total_time<1x
and worst<2x for every six widths on ARM64 and x86-64 independently.


### R13a — avoid rounding-namespace shadow and record R13 performance (20 October 2026)

**Build warning:** AppleClang with `-Wshadow` reports that the local
`rounded` in the R13 ARM64 `Rational::get_d` certificate masks the
existing `Ariadne::rounded` in `numeric/rounding.hpp:261`.
Rename only this local value to `rounded_significand`; there is no
change to arithmetic, rounding, exact certificates or benchmarks.

**R13 measured, user-supplied (not R13a measurements):** ARM64
`to_double` at 256 bits improves from R12's 106.383 ns / 2.556x to
64.521 ns / 1.543x (mpq_t 41.806 ns). The ARM64 worst ratio at 256
bits is now 1.896x (`add`), while total_time remains 1.324x.
Across 53/128/256/1024/2048/4096 bits, ARM64 total/worst ratios:
.424/1.219, .778/1.782, 1.324/1.896, 1.158/1.421,
.897/1.555, .730/1.656. x86-64 total/worst ratios:
.897/1.462, 1.108/2.249, 1.743/2.240, 1.252/1.468,
.797/1.155, .574/1.123. All 384 paired input and 2,880 result
numerical benchmark checks pass on each machine. ARM64 now passes
worst<2x at all widths; only 53,128,2048,4096 pass BOTH gates on
ARM64. x86-64 passes both at 53,2048,4096. Thus 7 of 12 host-width
cases meet both gates after R13, and the strict acceptance target is
still NOT achieved. R13a introduces no performance algorithm changes.

**Validation:** R13's Unix/Windows/Coverage workflows were still
running when R13a was prepared. R13a CI and matched host benchmarks
must be checked independently. No fresh speedup or CI success is
claimed for R13a.


### R14 — relocate self-alias square dispatch from Integer to BigUnsigned (20 October 2026)

**Purpose:** restore the signed/unsigned responsibility boundary,
without changing benchmarks or selected multiplication kernels.
B33/B37 detected same-object `Integer` multiplication and routed
8+-limb operands to symmetric `BigUnsigned::square` on GCC x86-64
fused and ARM64 native backends. This is a magnitude optimization,
not signed arithmetic, so place the exact same platform/size/identity
check in `BigUnsigned mul(const&,const&)`. `Integer mul` now applies
only the XOR of operand signs to the delegated unsigned product.

**Recursion safety:** the `square` fallback calls `mul(a,a)`, but the
fast-path is reachable only at >=8 limbs on backends where `square`
uses a dedicated symmetric kernel at the same threshold. For <8 limbs
or all unsupported targets, `mul` does not dispatch to `square`.
Zero and oversized product checks remain in their existing order.
The new unsigned dispatch may also benefit other genuine self-alias
callers, but no runtime value-pattern or benchmark shortcut is added.

**Differential tests:** new GMP checks for same BigUnsigned object,
distinct-but-equal object, explicit square, source immutability,
zero, dense/all-ones/sparse/random values at 1,4,7,8,9,16,64 limbs.
The preexisting signed self-product, `x*=x`, sign, saturation and
GMP checks remain unchanged. Rebuild and run full Integer,
BigUnsigned, Rational and FloatAP gates. This changes shared
unsigned dispatch; identical performance has not yet been measured.

**Validation status:** before this commit, R13a Unix, Windows and
Coverage CI all passed. R14 compilation, host-specific benchmark
stability and full regression are pending.


### R15 — centralize platform-specific unsigned exponentiation in BigUnsigned (20 October 2026)

**Motivation:** the sign-independent B32 right-to-left exponentiation
schedule remained in `Integer::pow` on ARM64, Clang, MSVC and portable
targets, whereas GCC x86-64 fused builds already called B74's MSB-first
buffered `BigUnsigned::pow`. That made the signed front end responsible
for unsigned algorithms and architecture dispatch. R14 moved the
related self-multiply selection into `BigUnsigned::mul` in a preceding
separate commit.

**Change:** `Integer::pow` now solely converts the exponent, derives
the result sign from negative-base parity and delegates to
`BigUnsigned::pow`. The unsigned function retains B74's original
GCC x86-64 fused MSB-first implementation verbatim under the existing
compiler/backend definition. Other targets use B32's exact
right-to-left loop, moved from `Integer::pow`, preserving first-factor
copy elision, square/multiply ordering, zero and exponent-one handling.
The public `BigUnsigned::pow` API continues to be exact for all
unsigned inputs but now selects the same policy as the signed caller
on each target. This modifies the algorithmic schedule for *direct*
non-GCC callers of `BigUnsigned::pow` (previously MSB-first buffered),
so such callers need independent performance regression measurement.
The result of `Integer::pow` is unchanged for 0^0, +/-0, +/-1,
negative odd/even powers, aliasing, long values and every Nat exponent.
No GMP dependency, benchmark or timed operation is modified.

**Tests:** preserve existing randomized signed and unsigned GMP
exponentiation checks across 0..17 limbs and 0..1025-bit boundaries,
plus 7/8-limb self-multiply and aliasing regression from R14. Add an
independent unsigned GMP power matrix containing bit-pattern exponents
0/1/2/3/5/6/7/8/9/12/15/16/17/31/32/33 with immutable borrowed
inputs. Test and benchmark performance parity separately on ARM64,
GCC/Clang x86-64, MSVC and forced-portable builds.

**Validation:** R13a CI passed before R14. R14 CI was still running
when R15 was authored. R15's Unix/Windows/Coverage workflows and
unchanged signed/unsigned GMP benchmarks have not yet completed;
no semantic regression or performance improvement is claimed absent
those checks. Rebuild and rerun full Integer, BigUnsigned, Rational,
FloatAP and performance gates before accepting.


### R15a — fix MSVC narrowing warning in unsigned power GMP test (20 October 2026)

**Observed CI failure:** the R15 Windows MSVC 14.51 Release build
under `/W4 /WX` fails at `tests/test_big_unsigned.cpp:146` with
C4244: passing a `std::uint64_t` exponent to GMP's
`mpz_pow_ui(..., unsigned long)` potentially narrows on Windows,
where `unsigned long` is 32-bit. Unix and Coverage workflows for
R15 pass. This is solely a test harness type issue, not a production
arithmetic error.

**Fix:** use `const unsigned exponent` in the new R15 GMP-oracle
loop. Every one of the 16 original exponents (0 through 33 with
irregular bit patterns) is unchanged, and the same exponent widens
safely for `BigUnsigned::pow(..., std::uint64_t)` and fits exactly
in GMP's Windows `unsigned long`. No benchmark workloads, oracle,
algorithm, or actual test scenario changes.

**Validation status:** changes pushed directly onto `remove-gmp`;
Windows, Unix and Coverage CI of R15a must be checked. No local
MSVC compiler is available in the editing session.


### R16 — borrow coprime denominators and move normalized sum (20 October 2026)

**Baseline:** R15a (`b7f79d1`) completed Unix, Windows and Coverage
CI successfully. The last supplied side-by-side Rational benchmark
measurements are R13 (before the R14/R15 refactoring); do not present
those as R16 or R15a performance. The user acceptance gates remain
`total_time <1x` AND `worst <2x` at each of 53,128,256,1024,
2048,4096 bits on both ARM64 and x86-64. R13 was compliant on
7 of 12 host/precision combinations.

**Change:** in `FiniteOperations<Rational>::add_sub`, after
`g=gcd(den_a,den_b)`, if g=1, work directly with both borrowed
canonical denominators instead of copying them through
`exact_quotient(den,1)`. For canonical p/q and r/s with gcd(q,s)=1,
the signed cross-sum numerator `p*s +/- r*q` is coprime to q:
modulo q it is p*s (coprime to q), and similarly modulo s
it is +/-r*q (coprime to s). Hence the product q*s already is
fully reduced, provided the numerator is nonzero. Skip the redundant
`gcd(n,1)`, the quotient-by-1 copies and numerator copy; construct
from an owning numerator moved into the result. Canonical
cancellation to zero still follows the existing zero path. When
g!=1, keep the established reduced-denominator algorithm, but if
the final h=gcd(n,g)=1, again move the already-reduced numerator
and borrow the original denominator for its final product. h>1
continues through the prior exact quotient path. All cases preserve
the prior rational sign, denominator, alias and exceptional semantics.
No benchmark corpus, operation count, weights or independent GMP
oracle was modified.

**Tests:** add GMP-oracle exact arithmetic checks at all six widths
for coprime and shared-factor denominators, with g=1, g>1/h=1,
g>1/h>1, both operand orders and signs, cancellation to zero,
self-aliasing via += and -=, and canonical numerator/denominator.
Existing 1200 randomized GMP test pairs and FloatAP gates remain
unchanged.

**Validation:** R16 has not been locally compiled or benchmarked.
Require Unix, Windows and Coverage CI, the complete downstream
Integer/Rational/FloatAP test suite and fresh unchanged benchmark
runs on both hosts before claiming a speedup. The fast-path is an
algebraic identity for all canonical rationals, not benchmark-only.


### R17 — certified two-leading-limb Rational cross-product comparison (21 October 2026)

**Baseline:** R16 Unix and Coverage workflows passed, while Windows
failed during an external vcpkg `z-applocal` step with `Access is denied`
while copying benchmark DLLs after linking. The user classified this
as a temporary GitHub environment issue and explicitly requested that
it be ignored: no CMake or build-tool changes. R16's separately
supplied ARM64 and x86-64 benchmark reports both pass 384 exact input
checks and 2,880 operation-result checks. Total/worst ratios across
53/128/256/1024/2048/4096 bits: ARM64 .423/1.233,
.767/1.876, 1.232/1.565, 1.150/1.401, .888/1.562,
.697/1.553; x86-64 .882/1.422, 1.095/2.170,
1.691/2.039, 1.221/1.448, .774/1.144, .561/1.111.
Seven of twelve host-width pairs pass both required strict gates.
The principal x86-64 remaining worst-case is compare at 128 bits.

**R17 change:** when all four unsigned Rational cross-product
magnitudes contain at least two limbs and the two product limb sums
are identical, compute the products only of each operand's highest
two 64-bit limbs (4-limb local output). Writing operands as
x=U*B^k+r and y=V*B^m+s, with B=2^64 and 0<=r<B^k,
0<=s<B^m, the exact product obeys
UV*B^(k+m) <= xy < (U+1)(V+1)*B^(k+m).
The exact prefix product UV occupies up to four limbs. Form the
exact upper prefix bound UV+U+V+1 by four-limb additions with full
carry propagation; if the bound overflows four limbs, decline
certification. If the upper left bound is <= the lower right bound,
the left product is strictly smaller; symmetrically for right.
Otherwise fall back to the existing complete 8-limb stack/fallback
multiplication and comparison. Size mismatch and one-limb cases
remain on their prior exact paths. This is not a benchmark-specific
heuristic: it proves the result before returning. No approximations,
compiler-specific instructions beyond the existing 128-bit GCC/Clang
arithmetic support, input restrictions or GMP dependencies are added.

**Tests:** add GMP differential comparisons at 128/256/448/512/
1024/2048/4096 bits across contrasting and nearly equal rational
values, signs, swapped operands, and stack/heap thresholds. Retain
all prior randomized GMP tests and unchanged benchmark workloads.

**Validation:** R17 is not locally compiled or benchmarked. Require
Unix/Windows/Coverage checks plus unchanged ARM64/x86-64 31-sample,
20,000-iteration runs before interpreting performance. The target
remains `total_time<1x` and `worst<2x` for every host-width pair.


### R18 — use CMake architecture selectors consistently (21 October 2026)

Replace 40 architecture preprocessor predicates in `src/big_unsigned.cpp` and three in `src/rational.cpp` with the existing CMake-selected `ARIADNE_BIG_UNSIGNED_X86_64` / `ARIADNE_BIG_UNSIGNED_ARM64` target-private flags. This is a lexical-only replacement of `defined(__x86_64__)`, `defined(__aarch64__)`, and `defined(_M_X64)`; all compiler feature guards, algorithm bodies, thresholds, benchmarks and test oracles are unchanged. Document PRIVATE definition scope in `src/CMakeLists.txt`. Native regression and benchmarks pending; do not claim performance changes. R17 Unix and Coverage CI passed prior to this refactor.


### R19 — borrow uncancelled Rational multiply/divide factors (21 October 2026)

**Baseline.** R18 (`e1a8fb2`) Unix, Windows and Coverage workflows
all passed. The most recent separately supplied 6-width host
benchmarks are from R16; neither R17 nor R18 has paired native
performance measurements yet. R16 total_time/worst per width:
ARM64 .423/1.233, .767/1.876, 1.232/1.565, 1.150/1.401,
.888/1.562, .697/1.553; x86-64 .882/1.422, 1.095/2.170,
1.691/2.039, 1.221/1.448, .774/1.144, .561/1.111.
Both reports passed 384 paired inputs and 2,880 paired results.

**Change.** Canonical Rational multiplication and division use two
cross-GCDs (g and h), then produce two products from four exact
cancelled factors. Prior code evaluated `exact_quotient(x,1)` by
copying `x`; each cross-GCD equal to 1 could therefore cause two
unnecessary BigUnsigned copies. New private
`multiply_cancelled(lhs,lhs_divisor,rhs,rhs_divisor)` borrows each
original magnitude when its exact divisor is 1, independently for
both factors, and computes a real exact quotient only when the
divisor exceeds 1. The result is constructed as an owning numerator
and denominator, with numerator sign determined as before by XOR.
Keep the same cross-GCDs, exact multiplication, canonicality,
special-values and self-alias-square dispatch. Remove now-unused
Integer wrappers for `mul` and `exact_quotient`. These are
optimizations for all coprime and partially coprime factors,
not special cases of benchmark inputs.

**Regression.** Add independent GMP-oracle multiply/divide checks
across 53,128,256,1024,2048,4096 bits for all 4 combinations
of trivial and nontrivial cross factors, all 4 sign combinations,
both operand orders, compound operators, self-alias multiply and
self-division. Every result is checked for GMP numeric equality
and canonical positive denominator. Retain randomized GMP
differential suites, FloatAP coverage, benchmarks and workloads.

**Validation.** No local R19 full build or matched ARM64/x86-64
measurement has yet been performed. Require Unix/Windows/Coverage
and unchanged six-width 20,000-iteration/31-sample paired benchmarks
before interpreting performance. Strict acceptance remains
`total_time<1x` and `worst<2x` at every width on BOTH hosts.


### R19a — restore signed exact quotient used by Rational addition (21 October 2026)

**Build-blocking regression caused by R19:** the cleanup removed the
`exact_quotient(Integer const&,BigUnsigned const&)` wrapper as
supposedly unused. However, the `g>1,h>1` normalization branch in
`FiniteOperations<Rational>::add_sub` still calls
`exact_quotient(n,h)` with an Integer numerator. The remaining
BigUnsigned overload cannot accept Integer, causing AppleClang's
two reported compile errors in rational.cpp (around line 253).
R19 Unix, Windows and Coverage workflows all failed in Build.

**Minimal repair:** restore the exact signed wrapper from R18. It
delegates to the unsigned exact quotient of the magnitude and restores
the numerator's original sign with `Integer::from_magnitude`.
This is the existing, verified behavior for nontrivial shared factors
in add/sub. It is intentionally still needed by this branch; the
statement in R19 about removing that wrapper was incorrect. The
now-unused `integer_mul` wrapper stays removed. All of R19's
copy-avoiding factor borrowing in mul/div is retained, without any
change to operand values, GMP oracle, benchmarks or numeric tests.

**Verification:** ensure Unix, Windows and Coverage CI compile,
the GMP differential tests exercise the g>1,h>1 add/sub branch, and
the unchanged six-width benchmark is rerun on ARM64 and x86-64 before
claiming R19 performance. No successful R19a build or speedup is
claimed in advance.


### R19b — fix GMP test variable shadowing on Clang (21 October 2026)

**Observed failure:** R19a's AppleClang ARM64 `-Wshadow -Werror`
build rejects three GMP oracle declarations in
`TestRational::test_gmp_differential`: inner `ga(lhs)` and
`gb(rhs)` (line 651), and later `ga(a)` (line 669), mask the
`GmpRationalBoundary ga,gb,expected` variables declared at the
start of the same function (line 555). R19a Unix, Windows and
Coverage CI did not pass.

**Narrow correction:** in the additional R19 multiplication/division
oracle cases, rename inner `ga`/`gb` to
`oracle_lhs`/`oracle_rhs`, and the inner self-alias `ga`
to `oracle_self`. Update their exact references only.
All test inputs, 4-way sign/operand cases, exact GMP operations,
canonicality checks, benchmark workloads and production code remain
unchanged. No test coverage is weakened.

**Validation:** verify fresh Unix, Windows and Coverage CI.
Fresh R19 performance data cannot be claimed until native benchmark
runs are available.


### R19c — preserve self-division alias test without Clang self-assign warning (21 October 2026)

**Blocking CI error:** Clang 16 with `-Werror` rejects the R19
GMP-oracle test statement `self/=self` as
`-Wself-assign-overloaded` at `tests/test_rational.cpp:676`.
The failure is a warning-policy/test-expression issue, not an
identified Rational arithmetic defect. R19b Unix CI failed;
Coverage passed at the latest check.

**Minimal fix:** use `Rational& self_alias=self;` and invoke
`self/=self_alias`. Both references still name the exact same
object, so the test continues to exercise in-place self-division
and to compare against `mpq_div(a,a)`, followed by canonicality
verification. Do not disable diagnostics or alter the GMP oracle,
number of tests, operands, benchmarks or production implementation.

**Validation:** fresh Unix, Windows and Coverage builds and
regressions are required; no new performance conclusion until
the existing ARM64/x86-64 six-width benchmark is rerun.


### R20 — bounded stack-based signed cross-product sums (21 October 2026)

**Baseline:** R19c (`3ade341`) passed Unix, Windows and Coverage.
Both freshly supplied R19c benchmarks (20,000 iterations, 31 samples)
pass the 384 paired inputs / 2,880 operation results on ARM64 and
x86-64. Total_time/worst ratios for 53/128/256/1024/2048/4096 bits:
ARM64 .404/1.213, .746/1.771, 1.176/1.547,
1.105/1.438, .819/1.447, .681/1.649;
x86-64 .871/1.409, 1.074/1.991, 1.636/2.046,
1.167/1.360, .738/1.145, .539/1.119.
Seven of 12 host/width pairs satisfy both strict gates
total_time<1x AND worst<2x.

**Change:** in `Rational::add_sub`, when each signed numerator
magnitude and unsigned scale factor occupies at most eight 64-bit
limbs, replace two owning multiplications and a third owning
integer add/sub with exact multiplication into two 16-limb local
arrays, followed by `BigUnsigned::add_into` or `sub_into` into
a 17-limb local result. Use full-width magnitude comparison to
choose the subtraction direction, effective sign and cancellation
to canonical zero. Materialize only one final owning magnitude
with `Integer::from_magnitude`. The same helper serves both
coprime and shared-denominator branches. Four entirely scalar
factors, any unit scale factor and long operands preserve
the previous path to avoid undoing existing specializations.
The public API, GMP oracle, all mathematical operations,
platform flags, thresholds of underlying kernels and benchmark
inputs remain unchanged. This is a general allocation-reduction
technique, not an optimization keyed to benchmark values.

**Tests:** add independent GMP add/sub differential cases at
53/128/255/256/320/447/448/511/512/1024 bits. Exercise both
effective operand signs, dense and closely matched fractions,
zero cancellation, carries, subtract-borrow, compound operators
and the 8/9-limb boundary. Preserve the prior randomized
canonical GMP suites and FloatAP regressions.

**Validation pending:** R20 has not been compiled or measured
locally. Require Unix, Windows and Coverage CI, Integer,
BigUnsigned, Rational and FloatAP checks, and native six-width
ARM64 and x86-64 benchmarks before claiming speedup. Keep the
strict user gates per host and per size.

### R21 — bounded borrowed-view Euclidean GCD (21 October 2026)

**Implementation.** For nonzero inputs of
at most 128 limbs, keep Euclid's exact remainder sequence in three
independent stack buffers. Borrow both original inputs until their
first reduction; rotate the result buffer so neither live input ever
aliases the output. Normalize every remainder before the next step.
Use the existing remainder-only division kernel and its portable/ISA
selection; no new instructions, CPU requirements or floating estimates.
Preserve R11's exact quotient-one/two subtraction and trim identical
leading limbs before subtracting. Finish single-limb divisors with the
existing selected remainder kernel (including B61) and scalar Euclid.
Only the final GCD becomes an owning BigUnsigned. Larger operands
retain the existing arbitrary-precision path. The fixed bound limits
scratch, not supported integer precision. No benchmark inputs or
canonical Rational contracts change. Test zero, equality, input
preservation, common factors, dense random values, quotient boundaries
and the 128/129-limb dispatch boundary against independent values/GMP.
Record native timings and sanitizer/compiler checks after execution.

The 128-limb bound uses 3*128 remainder limbs plus 257 division scratch
limbs (5,128 bytes of arrays). Those arrays need no initial zero fill:
each kernel fully writes its active output and every consumer sees
only its normalized span. Euclid maintains a>=b>0; for equally long
inputs ah-bh<=bh proves a<3*b without overflowing 2*bh. At most two
exact subtractions suffice. Trimming equal high limbs is safe because
the first unequal lower limb proves the remaining unsigned order.
The second subtraction uses independent division scratch, then copies
the normalized result back. Three-way rotation never overwrites either
live input. No architecture-dependent arithmetic was added.

**Regression coverage.** New independent cases use
gcd((2^k-1)*f,(2^k+1)*f)=f, f=2^64+3, with k>=126.
The coprimality follows from the odd operands differing by two.
Dense deterministic GMP cases at 2/4/5/8/16/64/127/128/129/160 limbs
exercise long buffer rotations, both argument orders, the exact stack
limit, fallback, and unchanged sources. Existing randomized GCD,
quotient-one/two, equality, zero, scalar, Rational canonicality and
special-value suites are preserved.

**Online measurement.** A dedicated Rational benchmark push/manual
workflow compares HEAD^ with HEAD on the same Ubuntu 24.04 x86-64
runner, separately with GCC 13 and Clang 16. It builds both Release
revisions with warnings-as-errors, tests the candidate, and runs the
unchanged six-width 20,000x31 benchmark sequentially. Raw logs, CPU
and compiler details and both SHAs are uploaded as artifacts; the
summary is informational because shared hosted CPUs have timing noise.
No speed gate is weakened and no benchmark fixture is changed.

**Local validation.** AppleClang 21 Release warnings-as-errors: all
20 C++ CTest suites pass. Debug ASan/UBSan: Integer, BigUnsigned,
Rational, Float and FloatRounding pass (5/5). GCC 16.2 ARM64 Release
warnings-as-errors: Integer, BigUnsigned and Rational pass (3/3).
The workflow YAML parses and git diff --check is clean. Native MSVC,
x86-64 and remote Unix/Windows/Coverage validation remain pending at
publication; do not infer their results from the local checks.

**Native ARM64 measurements.** Apple M2 Ultra (Mac Studio),
macOS 26.6.2, AppleClang
21.0.0.21000334, Release, unchanged benchmark, six widths,
20,000 iterations and 31 samples. R20 parent and R21 were rebuilt
in separate directories; the recorded comparison ran sequentially
without concurrent builds/tests. Both pass 384 paired inputs and
2,880 exact operation results. Rational/GMP total/worst ratios:

| bits | R20 total | R21 total | R20 worst | R21 worst |
|---:|---:|---:|---:|---:|
| 53 | .356 | .367 | 1.245 | 1.235 |
| 128 | .649 | .565 | 1.803 | 1.714 |
| 256 | 1.120 | .894 | 1.697 | 1.750 |
| 1024 | 1.111 | .898 | 1.402 | 1.413 |
| 2048 | .825 | .720 | 1.499 | 1.345 |
| 4096 | .682 | .627 | 1.630 | 1.516 |
| ALL | .772 | .680 | 1.803 | 1.750 |

All six ARM64 widths meet strict total<1 and worst<2 in this run.
Worst cases are to_double at 53/256/1024, compare at 128 and
construct_fraction at 2048/4096. These kernels were not altered;
small absolute timing changes include measurement/code-layout noise.
The sum of native medians at 256 bits falls from 1,987.330 to
1,587.294 ns, and at 1024 from 4,087.420 to 3,302.547 ns.
Examples of absolute R20 -> R21 medians (ns): 256-bit construction
194.525 -> 138.338, add 315.548 -> 237.144, mul 372.442 -> 283.812;
1024-bit mul 863.206 -> 629.456, div 857.354 -> 636.850.
The 53-bit native sum increases 1.9% (411.401 -> 419.288 ns),
where scalar inputs retain the old GCD path; no universal speedup
is claimed. Recheck strict gates on native x86-64.

The paired Integer 20,000x31 comparison passes 768 inputs/2,688
results: aggregate total .958 -> .951 and worst 1.659 -> 1.646.
Existing marginal 2048/4096 total ratios remain around 1x
(1.007/1.013 -> 1.001/1.003); no new Integer gate guarantee.
FloatAP's unchanged 1,000x11 comparison passes 256 inputs/3,456
results. Total/worst at 53/128/256/1024 bits is
1.397/2.483 -> 1.358/2.471,
1.071/2.696 -> 1.072/2.683,
.596/2.053 -> .591/2.055,
.249/1.642 -> .248/1.615.
Its prior short-width performance deficits remain outside R21's scope.

### R22 — exact scalar division and borrowed cancelled factors (21 October 2026)

**Implemented contract.** Add a shared exact
single-limb quotient kernel for the proven-divisibility cases in Rational
normalization and cross-cancellation. Write a disjoint caller buffer of
at least the source length; preserve the source and zero extra output.
Split d=odd*2^s, require the source's low s bits to vanish, and read the
right-shifted source without an owning shift. For odd d, compute its
inverse modulo B=2^64, then recover quotient digits from low to high:
q_i=(a_i-carry)*inverse mod B; carry'=high(q_i*d)+borrow(a_i,carry).
The invariant q_partial*d = a_partial + carry*B^(i+1) proves the digits;
carry<d guarantees the high product plus borrow cannot overflow.
A nonzero final carry rejects nonexact input. The inverse seed
(3*d) XOR 2 is correct modulo 32 for odd d; four Newton iterations
x=x*(2-d*x) double the valid bits to at least 64. All arithmetic is
unsigned modulo B. Power-of-two and scalar source cases retain direct
shifts/native scalar division. A multiprecision divisor keeps the
existing general quotient path under the documented exact-divisibility
precondition. No precision limit or GMP dependency is introduced.

Wide products use compiler-available unsigned 128-bit arithmetic,
MSVC x64 _umul128, or the existing pure C++20 32-bit-half multiply.
No new ISA requirement or dispatch option is added. The inverse fixes
the exact low digit; the wide product supplies the next high carry; inputs, output and all temporary
factor/product buffers are disjoint. Rational's <=8-limb cancelled
product can borrow unit-divisor factors and store other scalar-divisor
quotients in independent eight-limb arrays, then materialize only the
final product. Longer operands/multiprecision divisors retain owning
fallbacks. All gcd and canonical Rational proofs remain unchanged.

Validate independent exact products, all divisor bit lengths, odd/even
and power-of-two divisors, carries, zeros, UINT64_MAX, inline/heap
boundaries, large operands, source preservation and invalid scalar
inputs against GMP; retain full Rational/Integer/FloatAP regressions.
No performance result or threshold acceptance is claimed in advance.

**R21 completed validation and x86-64 baseline.** Commit 97e9ae3 passed
Unix, Windows and Coverage. The user's x86-64 20,000x31 result passes
384 inputs/2,880 results: per-width total/worst at 53/128/256/1024/
2048/4096 bits is .872/1.391, .892/1.995, 1.294/1.730,
.980/1.382, .654/1.064, .500/.975. The 256-bit total remains
above 1x; compare at 128 bits has almost no margin under 2x.
The new hosted R21 workflow passed both compiler jobs, but its actual
performance limits are compiler/host dependent: GCC 13 on AMD EPYC
9V74 reports total/worst .729/1.672, .652/2.210, 1.066/1.667,
.828/1.006, .561/1.065, .462/.953. Clang 16 on AMD EPYC 7763
reports .501/1.163, .668/2.018, 1.100/2.029, 1.552/2.943,
1.881/3.223, 2.473/3.970. Do not treat successful benchmark jobs
as performance-gate acceptance. Hosted Clang's large-width square
path is a separate existing deficit. The artifacts were downloaded
read-only and matched the SHA-256 digests returned by GitHub.

**Validation.** Apple M2 Ultra, macOS 26.6.2, AppleClang 21 Release
warnings-as-errors: all 20 C++ suites pass. Debug ASan/UBSan passes
Integer, BigUnsigned, Rational, Float and FloatRounding (5/5).
GCC 16.2 ARM64 Release warnings-as-errors passes Integer, BigUnsigned
and Rational (3/3). The new direct limb/owning tests cover 6,200 exact
cases constructed with GMP, all divisor bit lengths, every odd seed
residue modulo 32, mixed even/odd divisors, all-ones/random/sparse/zero
quotients and up to 130 dividend limbs. Independent fixed values and
invalid scalar division/zero/undersized buffers are checked separately.
Rational tests add controlled scalar cross-cancellations (including
14 and 15), all sign pairs, compound operators, fixed exact expected
fractions and the eight/nine-limb boundary. The preserved full suite
still covers multiprecision-divisor fallbacks and special values.
A standalone 2,048-case GMP differential build with ARM64 selection
but __SIZEOF_INT128__ undefined exercises the new kernel's existing
C++20 32-bit-half wide-product fallback, including mixed even divisors.
No native x86-64/MSVC performance claim precedes the remote run.

**Native ARM64 performance.** Apple M2 Ultra, AppleClang 21, macOS
26.6.2, Release, unchanged six-width 20,000x31 comparison. Both R21
and final R22 pass 384 inputs and 2,880 paired exact results. Recorded
runs were sequential and had no concurrent builds or tests. The
earlier R22 trials are diagnostic; the table is the final source version.

| bits | R21 total | R22 total | R21 worst | R22 worst |
|---:|---:|---:|---:|---:|
| 53 bits | 0.360 | 0.369 | 1.235 | 1.236 |
| 128 bits | 0.573 | 0.538 | 1.837 | 1.738 |
| 256 bits | 0.902 | 0.882 | 1.701 | 1.687 |
| 1024 bits | 0.905 | 0.884 | 1.420 | 1.443 |
| 2048 bits | 0.729 | 0.709 | 1.298 | 1.153 |
| 4096 bits | 0.621 | 0.594 | 1.472 | 1.276 |
| ALL | 0.680 | 0.656 | 1.837 | 1.738 |

All six ARM64 Rational widths retain total<1x and worst<2x. This is
an incremental improvement; no native x86-64 gate acceptance is inferred.
Absolute 256-bit native medians (R21 -> R22, ns):
- construct_fraction: 138.954 -> 137.442.
- add: 242.127 -> 241.981.
- sub: 241.113 -> 241.244.
- mul: 288.898 -> 274.852.
- div: 283.467 -> 276.931.
Sum of native medians at 53 bits: 415.819 -> 419.818 ns.
Sum of native medians at 128 bits: 1106.676 -> 1044.257 ns.
Sum of native medians at 256 bits: 1604.316 -> 1602.129 ns.
Sum of native medians at 1024 bits: 3325.132 -> 3321.063 ns.
Sum of native medians at 2048 bits: 5470.241 -> 5485.076 ns.
Sum of native medians at 4096 bits: 10695.727 -> 10704.573 ns.

Integer 20,000x31 passes 768 inputs/2,688 results. Aggregate total/
worst is .951/1.614 -> .947/1.619; individual marginal long-width
ratios continue to vary around 1x. No Integer gate guarantee is added.
FloatAP 1,000x11 passes 256 inputs/3,456 exact results for both
revisions. Final R22 total/worst at 53/128/256/1024 is
1.654/3.431, 1.079/2.697, .593/2.069, .247/1.611.
The unchanged R21 baseline itself varied at 53 bits from total 1.337
in the first run to 1.679 in a repeat, largely through the MPFR timing
denominator. Do not interpret the first-run ratio change as a measured
FloatAP arithmetic regression or claim universal nonregression. The
existing short-width FloatAP gates remain unmet; no FloatAP arithmetic
was changed. Raw first/repeated runs are retained for comparison.
Native x86-64 GCC/Clang measurements and Unix/Windows/Coverage CI for
R22 remain pending at publication.
