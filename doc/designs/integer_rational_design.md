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

