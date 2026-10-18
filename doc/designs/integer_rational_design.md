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
