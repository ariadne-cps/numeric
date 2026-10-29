# FloatDP — fixed-precision binary64 floating point

> **Current design perspective — 29 October 2026.** This is the living
> architecture, correctness and optimization journal for the `FloatDP`
> backend. FloatDP is implemented independently of `FloatMP`; the certified
> near-half transcendental kernel in `src/floatdp.cpp` has **no
> `BigUnsigned` dependency** as of D2. This is a documentation reorganization
> of existing history, not an additional numeric optimization.

Keep this file synchronized with subsequent **FloatDP-specific** changes.
For the native arbitrary-precision `FloatMP` implementation, D-series
optimizations, M148 compiler fix and shared dual-backend benchmarking,
see [`floatmp_design.md`](floatmp_design.md).
The shared benchmark is defined in
`benchmarks/benchmark_backend_workload.hpp` and its outputs require
same-machine/same-configuration comparisons before causal attribution.

## Journal scope and chronology

The following historical sections were moved here from
`floatmp_design.md` without changing their technical content:
the 17 October FloatDP-related benchmark maintenance,
the D123 FloatDP exact-binary64 test-construction fix,
the 25 October independent certified transcendental kernel,
and all optimization interventions **D1–D6** with their measured
follow-ups. Their original identifiers, dates, verification caveats
and per-operation benchmark measurements are retained.
Historical descriptions of shared limb multiplication or
`BigUnsigned` in the FloatDP implementation record its *former*
state; D1/D2 supersede that runtime dependency.

The `FloatMP` M148 and M148.2 changes are deliberately retained in
[`floatmp_design.md`](floatmp_design.md). Any results comparing
the two backends are preserved here when they measure a FloatDP
optimization; shared benchmark infrastructure and FloatMP-only
algorithm changes remain in the FloatMP journal.

## FloatDP benchmark: omit redundant ALL row (17 October 2026)

The FloatMP/FloatDP benchmark is limited to binary64/53-bit precision, so its `53 bits` summary
already represents every one of the 27 operation ratios. Remove the redundant `ALL` summary row
**only** from `benchmark_compare_double_precision_backends.cpp`. Retain the aggregated slowdown
accumulator for `--max-slowdown` validation and keep all 27 timing rows, CLI flags and the
1000-iteration/11-sample defaults unchanged. This is an output-only change; no numerical backend,
operation or measurement method is modified.

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

### D123 — use public FloatDP exact-binary64 constructors in decimal tests (23 October 2026)

Fix four CI compilation errors in `tests/test_float.cpp` introduced
by M122: FloatDP(double) is a private constructor, deliberately not a
public import path. Decimal rounding regression tests now construct
FloatDP via its existing public FloatDP(ExactDouble const&,DoublePrecision)
constructor, using the `dp` precision tag. The original exact
binary64 inputs 0.125, 0.375 and -0.0 remain unchanged, including
negative zero's sign bit. Expected oracle strings, rounding modes,
production code and benchmarks are untouched. The compiler previously
failed before running any new decimal checks, so CI still needs to
confirm compilation and execution of these tests.

## FloatDP certified transcendental decoupling (25 October 2026)

**Architecture correction.** The binary64 `FloatDP` transcendental fast path must not construct
`FloatMP` values or initialize its coefficient cache through the `FloatMP` transcendental
engine. The two public floating-point specializations may share exact low-level integer primitives,
but one floating-point type must not be used as the implementation substrate of the other.

`FloatDP` now owns a binary64-specific near-half fixed-point kernel in `src/floatdp.cpp`.
Its centers are outward dyadic enclosures at scale `2^-85`, derived from the repository's
independently generated exact transcendental oracle data. Runtime evaluation uses only binary64 bit
decomposition, two-limb fixed-point coefficients, the shared `BigUnsigned::mul_into` primitive,
explicit remainder/error margins, and final binary64 rounding. No `FloatMP` object, precision tag,
adaptive `FloatMP` series, or `FloatMPShortTranscendentals` cache is involved.

The centered derivatives use exact identities where possible: `tan' = 1+tan^2`,
`atan' = 4/5`, `asin' = -acos' = 2/sqrt(3)`, and `tanh' = 1-tanh^2`.
At binary64 precision this reproduces the original `53+32=85` certification scale and its `|delta| <= 2^-43` centered window. The fast path is accepted only when the lower and upper certified fixed-point endpoints round to
the same binary64 value in the requested rounding mode; otherwise evaluation falls back to the
existing general `FloatDP` implementation.

**Verification before integration.** A direct simulation of the candidate fixed-point kernel against
all versioned 53-bit independent oracle rows for `sin`, `cos`, `tan`, `atan`, `asin`,
`acos`, `tanh`, centered `exp`, and centered `log` produced 289 fast-path decisions with
zero mismatches. Rows outside the certified near-half window, or whose enclosure did not decide a
unique rounded result, remained on the fallback path. Repository CI, especially MSVC Debug, is the
integration gate; no performance result is inferred until those jobs complete.


## D1 — FloatDP native fixed-point correction without BigUnsigned multiplication (29 October 2026)

**Naming policy and incremental plan.** Optimization interventions
after M148 are numbered D1, D2, D3, etc. The user requests that
three proposed FloatDP changes be **implemented and benchmarked one
at a time** and that FloatDP ultimately use **no BigUnsigned**.
D1 changes only multiplication for the near-half fixed-point
correction; D2 will replace the remaining BigUnsigned calls for
applying an endpoint displacement; D3 will replace the final
std::ldexp-based conversion with direct IEEE-754 binary64 bit
assembly. No D2 or D3 implementation is included in this commit,
so latency changes can be attributed to a single local operation.

**D1 scope and invariant.** In src/floatdp.cpp,
floatdp_correction previously called BigUnsigned::mul_into on
a 2-limb derivative bound and a 1-limb exact offset. The existing
compile-time certificate proves the slope is <2^87, the offset
is <=2^11, and the input denominator shift is 53 or 54. Replace
the generic multiplication call with native unsigned 64-bit
arithmetic on two 32-bit partial products of the lower word
plus the derivative's upper word and carry:
bottom=(slope.low&0xffffffff)*offset;
top=(slope.low>>32)*offset+(bottom>>32);
low=(top<<32)|(bottom&0xffffffff);
high=slope.high*offset+(top>>32).
Here low/high are the exact low 128 product bits; by the
certificate the product has fewer than 98 significant bits and
there is no missing third limb. The same exact right shift,
discarded-bit test, and optional outward increment are preserved.
This is strictly integer arithmetic, needs neither a compiler
128-bit extension nor an additional library, and is portable
across GCC, Clang/AppleClang and MSVC.

**Independent arithmetic verification.** An isolated C++20
differential program using boost::multiprecision::cpp_int as
an independent exact full-product oracle completed **614,700**
floor/ceil correction comparisons with GCC and separately
**614,700** with Clang, with zero mismatches. Cases cover
shift=53 and 54, offsets 0..2048, randomized bounded slopes,
zero/highest coefficient words, all-one lower words, exact
products and nonzero discarded remainders. These checks
were run independently of the numeric repository; they do
not constitute a complete build, CI pass, or timing result.

**Regression extension.** The existing frozen
FloatDP/FloatMP shared transcendental golden-oracle checks
now additionally compare IEEE-754 encoded bits of the
actual FloatDP result with the exact dyadic oracle's binary64
encoding, not only Rational value equality. Current floating
rounding mode preservation remains checked. No changes to
oracle source data, certified windows, coefficient tables,
or comparison benchmark.

**Temporary BigUnsigned usage.** D1 does *not* yet remove the
BigUnsigned include or its use inside floatdp_endpoint and
floatdp_round_fixed; those belong to D2. D1 removes only the
hot-path multiplication dependency in order to preserve
one-change-at-a-time performance attribution. The requirement
that FloatDP have **no BigUnsigned** is an explicit acceptance
criterion for completion of D2; D3 is independent of that
cleanup. Until then, the file still uses BigUnsigned.

**Performance gate and future work.** Baseline for D1 is the
user's post-M148 53-bit paired run: FloatMP/FloatDP total_time
0.947x, geomean 1.273x, worst 8.098x. Independently sum
FloatDP transcendental medians and compare to 479.018 ns.
Rerun identical 53-bit inputs, 20,000 iterations and 31 samples
on the same machine; do not infer a speedup from compilation or
the mathematical proof. Run the full FloatDP source-keyed
oracles, Unix/Windows/MSVC, Sanitizers and Coverage before
regarding D1 as validated. D2 and D3 must only be implemented
after observing the D1 benchmark independently.


### D1.1 — user-measured p=53 benchmark after D1 (29 October 2026)

The user supplied a complete same-configuration 53-bit benchmark
(20,000 iterations, 31 samples, 56 paired operations). After D1:
**FloatMP/FloatDP total_time=1.106x**, **geomean=1.651x**,
**worst=8.105x** (fma_nearest). Previous post-M148, pre-D1
measurement: total_time=0.947x, geomean=1.273x, worst=8.098x.
The FloatMP measurements themselves stayed broadly similar, but
the FloatDP certified transcendental operations became much faster.
Summing the 27 separately measured FloatDP transcendental medians
(exp/log/sin/cos/tan/tanh/atan/asin/acos in nearest/down/up) gives
**305.266 ns after D1**, down from **479.018 ns pre-D1**
(36.27% lower). Summing the 18 directed (down/up) entries gives
**137.645 ns after D1** versus **257.387 ns pre-D1**
(46.52% lower). Representative paired FloatDP timings:
exp_down 13.994 -> 7.765 ns; log_down 13.140 -> 6.860 ns;
sin_down 14.558 -> 8.035 ns; cos_down 14.281 -> 7.223 ns;
tan_down 14.562 -> 8.012 ns; atan_down 14.850 -> 7.996 ns.
These are per-operation medians supplied by the user, not
internally executed or independently controlled runs. The
measured acceleration strongly supports continuing incremental
optimization; no claims of broad input-domain speedups beyond
the certified near-half kernel follow from the benchmark.

## D2 — FloatDP native certified endpoint arithmetic, no BigUnsigned (29 October 2026)

**Purpose and one-change isolation.** Following successful D1
measurements, replace only the remaining native-word endpoint
arithmetic in src/floatdp.cpp; retain the D1 slope multiplication,
D1 exact-binary64 oracle checks, the 85-bit fixed-point scale,
the 53/54-bit denominator restrictions, and all compile-time
center-certification conditions. Leave the final
std::ldexp-based double construction UNCHANGED for the separate
D3 experiment. No unrelated FloatMP or benchmark modifications.

**Before D2.** floatdp_endpoint took a FloatDPFixed85 center
(two 64-bit words), created a two-element BigUnsigned::Limb
array, applied one signed displacement with BigUnsigned::add_into
or sub_into, and returned another two-word array. The rounder
floatdp_round_fixed consumed that array. This occurred once
per certified transcendental call. It did not require arbitrary
precision, because D1 and the existing exhaustive consteval
certificates already bound the center and correction to <128 bits.

**After D2.** floatdp_endpoint now returns the existing plain
FloatDPFixed85 struct and applies a signed one-word correction
via portable uint64_t addition/subtraction, comparing the
original/new low words to propagate a 0/1 carry or borrow to
the high word. Absolute displacement is computed with unsigned
modular subtraction to handle even INT64_MIN without signed
overflow. The existing exact compile-time 85-bit certificate
guarantees that actual supported arguments cannot overflow or
underflow the full two-word positive endpoint. The rounder
accesses magnitude.low and magnitude.high instead of array
indices, with its bit extraction, mode-dependent rounding,
ldexp scaling and sign unchanged. Remove the now-unused
numeric/big_unsigned.hpp and <array> includes. The text
BigUnsigned does not appear anywhere in src/floatdp.cpp as
of this change: this completes the stated no-BigUnsigned
requirement for the concrete FloatDP implementation file.
The shared independent frozen oracles, 56-operation benchmark
and compiler-time center proofs are unchanged.

**Correctness checks.** Add five compile-time static_asserts
in src/floatdp.cpp for carries, borrows and a zero resulting
low word. A standalone C++20 differential model compared
the portable two-word operation against the full signed
cpp_int result for **900,000** randomized and adversarial
in-range endpoint/displacement cases with GCC 14.2 and
**900,000** with Clang; no mismatches. Samples included
zero/all-one/high-bit lower words, positive and negative
displacements, boundary carries/borrows and INT64_MIN
subtraction. This verifies the mathematical kernel in
isolation, not the repository's compilation or CI. The
frozen FloatDP transcendental bitwise-oracle checks added
in D1 remain active. CI, the full tests, MSVC and a user
same-machine benchmark are still required.

**D2 performance baseline.** Use D1's independently measured
FloatDP p=53 transcendental sums: 305.266 ns for all 27
medians, and 137.645 ns for the 18 directed operations.
The paired FloatMP/FloatDP total_time ratio is 1.106x,
geomean 1.651x. Assess absolute FloatDP operation latency
and relative paired results after D2, including changes
outside the nine near-half kernels. Do not attribute
small global changes to this patch without controlled
repeated measurements.

**Next: D3, intentionally NOT in D2.** Once D2 performance
is measured, construct the normal IEEE-754 double directly
using std::bit_cast of the exponent/sign/fraction fields,
under the existing strict compile-time no-carry/no-tie
certificate. Compare output bit patterns against frozen
binary64 source-keyed oracles; retain all other floating-point
operations in FloatDP unchanged.


### D2.1 — user-measured 53-bit performance after D2 (29 October 2026)

A new user-supplied paired run at p=53, 20,000 iterations and
31 samples for 56 operations gives **FloatMP/FloatDP
total_time=1.183x**, **geomean=1.949x**, **worst=8.105x**
(fma_nearest); the previous D1 run reported 1.106x,
1.651x and 8.105x. The FloatMP medians were essentially
unchanged on this pair of runs. The 27 FloatDP transcendental
nearest/down/up medians sum to **246.292 ns after D2**,
compared with **305.266 ns after D1** (19.32% lower). Their
18 directed/down-up medians sum to **94.852 ns after D2**,
versus **137.645 ns after D1** (31.09% lower). Representative
FloatDP downward values in ns, D1 -> D2:
exp 7.765 -> 6.579; log 6.860 -> 4.281;
sin 8.035 -> 5.254; cos 7.223 -> 4.771;
tan 8.012 -> 5.473; atan 7.996 -> 5.438.
This is a user measurement on the actual benchmark workload,
not a local CI/build result. The stable FloatMP timings and
the consistent FloatDP improvements across modes support
D2's intended optimization, without proving improvements for
arbitrary transcendental argument ranges.

## D3 — IEEE-754 bitwise FloatDP result construction (29 October 2026)

**Scope and sequencing.** After D1's fixed-point multiplication
and D2's exact two-word endpoint arithmetic, the only requested
remaining intervention is constructing the certified near-half
FloatDP transcendental result directly as its binary64 encoding.
This patch does NOT change D1 or D2, the center tables, compile-
time certificates, input decoder, exact endpoint displacements,
retained-significand calculation, any rounding decisions, or
the public FloatDP APIs. It changes only the final conversion
in floatdp_round_fixed, formerly a double conversion followed
by std::ldexp and application of the sign.

**Correctness argument.** The exhaustive consteval proofs
floatdp_verify_certificates cover *every binary64 input*
inside the certified near-half 53/54 denominator window for
the nine supported transcendental operations and both signs.
For all such inputs, the rounded magnitude has a normalized
53-bit significand, never reaches an exact halfway case,
never carries beyond 53 bits, and has a finite normal
binary64 exponent. In this kernel the fixed-point scale is
B=85; width is within [84,87], giving the unbiased output
exponent 'highest=width-1-B' within [-2,+1]. Conversion
of any 53-bit unsigned significand to double is exact,
and scaling by 2^(highest-52) is exact throughout this
normal exponent range. Constructing the binary64 bit
pattern with sign bit 63, biased exponent (highest+1023)
in bits 62..52 and the low 52 significand bits is therefore
mathematically and bitwise identical to the previous
std::ldexp(static_cast<double>(retained),highest-52).
This holds in every ambient rounding mode and preserves
negative results. Any input outside the certified window
still uses the existing non-certified fallback, unchanged.

**Implementation.** Introduce one internal constexpr
floatdp_normal_bits(bool,int,uint64_t) mapper. Ensure the
platform has IEEE-754 binary64 via static_assert of the
standard double properties. Add constexpr asserts for
the exact encodings of +1.0, -0.5, +0.375 and +2.5,
covering sign, fractional significand and exponent limits.
Call std::bit_cast<double>(floatdp_normal_bits(...)) in
the last step of floatdp_round_fixed. Do NOT replace
the separate FloatDP(TwoExp) constructor's use of ldexp:
that API supports a different argument domain and is
outside this single-change D3 optimization. D3 introduces
no BigUnsigned; D2's no-BigUnsigned invariant stays true.

**Independent validation.** An isolated C++20 test using
-O2 -Wall -Wextra -Werror -frounding-math compared the
former ldexp output (including negative sign) against
the new integer encoding bit-for-bit for 4,000,000
random, normalized, non-overflowing 53-bit significands
covering width=84,85,86,87 and both output signs in
each of the four ambient hardware rounding modes.
All **4,000,000 cases passed on GCC 14.2** and all
**4,000,000 on Clang 17**, with four constexpr encoded
reference-value assertions. This is a standalone
bit-identity check of the conversion, NOT the entire
Numeric repository build, platform CI or a timing result.
The existing D1 source-keyed frozen FloatDP oracles now
also verify exact output bit patterns and remain unchanged.

**Release gate.** After committing, benchmark D3 alone
with the user's p=53, 20,000-iteration/31-sample paired
setup. D2 baselines: FloatDP sum of 27 transcendental
medians 246.292 ns; 18 directed transcendental medians
94.852 ns; paired total_time 1.183x, geomean 1.949x.
Compare individual operation absolute latencies rather
than inferring a gain from the ratio alone. Confirm
bitwise frozen-oracle checks, full source tests,
Windows/MSVC, GCC/Clang/AppleClang, Sanitizers and
Coverage. Do not claim any D3 speedup or integrated CI
status until those measurements and checks are available.

**Next work is separate from D3.** If more improvements
are desired, scope new D4 work only after observing
D3's isolated benchmark. In particular avoid adding
multiprecision facilities to the FloatDP centered
kernel or changing the strong consteval certification
without corresponding exhaustive independent tests.


### D3.1 — user-measured p=53 performance after D3 (29 October 2026)

The user supplied the paired 56-operation benchmark at 53-bit precision
with 20,000 iterations and 31 samples after D3, the IEEE-754
bit-assembly optimization. FloatMP/FloatDP **total_time=1.188x**,
**geomean=2.033x**, and **worst=8.056x** (fma_nearest). The
preceding D2 run measured 1.183x, 1.949x and 8.105x. The very small
change in total_time is not strong evidence of a change in aggregate
performance; per-operation measurements show a more informative
pattern. The FloatMP timings remain broadly stable.

For the nine transcendental functions in the comparison benchmark
(exp, log, sin, cos, tan, tanh, atan, asin, acos), the unweighted
sum of all 27 FloatDP per-operation medians changed from
**246.292 ns (D2)** to **234.770 ns (D3)**, a reduction of
**11.522 ns / 4.68%**. For the 18 directed (downward/upward)
entries, the sum changed from **94.852 ns (D2)** to
**85.397 ns (D3)**, a reduction of **9.455 ns / 9.97%**.
Nearest-mode entries collectively changed from 151.440 ns to
149.373 ns (approximately 1.36% lower). These are sums of
independent median measurements, not measured sequential elapsed time.

Representative individual FloatDP downward measurements (ns),
D2 -> D3:
* exp: 6.579 -> 6.419 (2.4% lower);
* log: 4.281 -> 3.265 (23.7% lower);
* sin: 5.254 -> 4.379 (16.7% lower);
* cos: 4.771 -> 4.258 (10.8% lower);
* tan: 5.473 -> 4.871 (11.0% lower);
* tanh: 5.331 -> 4.804 (9.9% lower);
* atan: 5.438 -> 4.708 (13.4% lower);
* asin: 5.252 -> 4.460 (15.1% lower);
* acos: 5.450 -> 5.169 (5.2% lower).

The results are consistent with removing std::ldexp and integer-to-double
conversion from the certified short path, but do not constitute an
isolated instruction-level causal proof. They reflect one
user-supplied run per version, without cross-run variance or
interleaved A/B trials. D1, D2 and D3 together reduced the
27-operation FloatDP transcendental median sum from
**479.018 ns after M148 (before D1)** to **234.770 ns after D3**
(about 51.0% lower), and the 18 directed median sum from
**257.387 ns to 85.397 ns** (about 66.8% lower). Over that sequence
the paired FloatMP/FloatDP total_time ratio moved from 0.947x
to 1.188x. These comparisons use the actual fixed
near-half-heavy benchmark pool, not arbitrary input distributions.

**Final D1–D3 scope statement.** The three agreed incremental
optimizations are complete and independently benchmarked by the user.
src/floatdp.cpp contains no BigUnsigned dependency. Its certified
near-half path now uses portable native-word multiplication (D1),
portable native-word endpoint carry/borrow (D2), and IEEE-754
binary64 bit assembly (D3), preserving all static consteval
certificates. No D4 implementation is authorized or included by
this benchmark entry. Full repository builds, independent frozen
oracle execution, Windows/MSVC checks, sanitizers, coverage and
additional controlled host repetitions remain outstanding as
verification tasks; nothing here asserts those tests were executed.


## D4 — correctly rounded native binary64 square root (29 October 2026)

**User observation and optimization target.** The post-D3 user benchmark
at 53 bits (20,000 iterations, 31 samples, 56 paired operations) reports
FloatDP sqrt_down=95.900 ns and sqrt_up=96.781 ns, whereas FloatMP
takes 25.873 ns and 37.644 ns, respectively. The corresponding
FloatDP/FloatMP latency ratios are about 3.71x (down) and 2.57x (up).
FloatDP sqrt_nearest=27.158 ns; FloatMP sqrt_nearest=34.560 ns.
The entire post-D3 FloatMP/FloatDP benchmark reported total_time
1.188x, geomean 2.033x and worst 8.056x (fma_nearest).
The aim of D4 is to remove a demonstrably expensive, obsolete
binary64 square-root algorithm; improvements are to be established
in a new user benchmark, not assumed from this analysis.

**Root-cause analysis.** The explicit FloatDP sqrt(mode,x) adapter in
include/numeric/floatdp.hpp calls _apply_rnd(&sqrt_rnd,mode,x), which
uses RoundingSession<FloatDP> to set and restore the requested hardware
rounding mode. src/double.cpp::sqrt_rnd previously implemented square
root via frexp normalization, an unbounded Newton iteration involving
division and averaging until two consecutive approximations compared
equal, and ldexp reconstruction. The repeated division operations
make it costly, especially under directed rounding; the equality
termination also does not itself provide a proof of a correctly
rounded IEEE-754 result. This is not the certified transcendental
path improved in D1--D3; it is an independent primitive.

**D4 implementation.** Keep the existing NaN forwarding and the
nonnegative-domain precondition. Replace only the Newton/frexp/ldexp
body of src/double.cpp::sqrt_rnd with std::sqrt(x). Under the
repository's strict floating-point compiler policy, CMake applies
-frounding-math to GNU/Clang/AppleClang and /fp:strict to MSVC.
On IEC-60559-conforming implementations, the hardware/native
sqrt operation rounds to the ambient rounding direction, which
_apply_rnd has already selected. The existing special-value path
remains, and native sqrt preserves the sign of negative zero
(previously the x==0 branch returned +0 for both signs).
Explicit FloatDP mode selection and restoration remain unchanged.
No BigUnsigned dependency is added to FloatDP; neither D1--D3
nor any 53--64-bit FloatMP kernel is modified.

**Correctness tests added to Numeric.** Extend the existing
test_float_rounding.cpp arithmetic regression with
DoublePrecision-only exact Rational checks over sqrt(2), sqrt(3),
exact squares, powers and adjacent binade endpoints, smallest
subnormal, smallest normal and largest finite binary64. For
each test input, verify directed lower/upper output squared
against the exact positive input; prove lower and upper outputs
are equal or adjacent binary64 values; verify toward-zero equals
downward; derive nearest output from the *exact rational squared
midpoint* and tie-to-even encoding. Verify preservation of the
ambient rounding mode after the four explicitly rounded calls,
negative-zero sign, positive infinity and NaN.

**Independent external arithmetic validation.** A standalone
C++20 differential test compiled and run with GNU GCC 14.2 and
Clang 17, -O2 -frounding-math -fno-fast-math, exercised 39,991
finite positive nonzero binary64 inputs per compiler, including
2,364 subnormals, exact squares, exponent extremes and values
near binade boundaries, for all four runtime rounding modes.
Exact dyadic decode and boost::multiprecision::cpp_int square
comparisons checked directed lower^2<=x<=upper^2 and adjacency.
Exact rational midpoint-squared comparisons tested nearest-even,
including the tie rule, independently of the candidate's
approximation. There were no mismatches on either compiler.
The test also checked std::sqrt(-0) sign and infinities/NaNs.
This is an independent host arithmetic check, NOT a compiled
build of the Numeric repository, not a macOS/Windows check,
and not a user timing measurement.

**Integration and release checks.** Run Numeric's complete
test suite, the new exact reference assertions, and tests
of CurrentRoundingMode and explicit FloatDP sqrt adapters
on Unix/AppleClang, GCC, Windows/MSVC and Sanitizers.
The result depends on conforming runtime rounding semantics:
validate compiler, libm, build flags and fenv interaction
on the target architecture. Reject any compiler/target pair
that uses an erroneously round-to-nearest-only sqrt despite
the explicit rounding mode. Preserve source-keyed golden
oracles and 100%-coverage behavior: do not adapt expected
answers merely to make tests pass. Compare p=53 sqrt_down,
sqrt_up and sqrt_nearest absolute timings to the D3 baselines
95.900, 96.781 and 27.158 ns; compare overall paired
total_time=1.188x and geomean=2.033x without attributing
other-operation noise to D4.

**Scope boundary.** This is D4 alone. Do not conflate the
native binary64 square root with FloatMP's integer sqrtrem
algorithm or D1--D3's certified centered transcendental
interventions. No D5 changes or new dependencies are included.


### D4.1 — user-provided p=53 results after D4 (29 October 2026)

The user supplied the 56-operation comparison benchmark after D4,
at p=53, 20,000 iterations and 31 samples. The measured
FloatMP/FloatDP ratio is **total_time=1.430x**,
**geomean=2.207x**, with worst=**11.385x**
for sqrt_nearest. The pre-D4 (post-D3) run reported
total_time=1.188x, geomean=2.033x and worst=8.056x
(fma_nearest). These metrics are sums and geometric means of
per-operation median latencies, not whole-program elapsed time.

The primary D4 target improved substantially, in ns:

| binary64 operation | pre-D4 | D4 | reduction |
| --- | ---: | ---: | ---: |
| sqrt_nearest | 27.158 | 3.140 | 88.4% |
| sqrt_down | 95.900 | 29.967 | 68.8% |
| sqrt_up | 96.781 | 30.369 | 68.6% |

The sum of the three FloatDP square-root medians decreased from
219.839 ns to 63.476 ns (approximately 71.1% lower). The
corresponding FloatMP measurements in this latest benchmark are
35.746 ns nearest, 25.627 ns downward, and 37.721 ns upward.
FloatDP's nearest square root is therefore about 11.385x faster
than FloatMP's; the worst reported FloatMP/FloatDP ratio now
occurs at this operation. The FloatMP square-root values and
the other FloatMP operations remained broadly stable.

**Interpretation and next bottleneck.** The native IEEE-754 sqrt
primitive has removed the Newton-loop overhead. Its
nearest-mode result now costs only about 3.140 ns, whereas
the directed modes still cost around 30 ns. Inspection of
include/numeric/floatdp.hpp confirms that sqrt(mode,x) goes
through _apply_rnd(&sqrt_rnd,mode,x), which constructs a
RoundingSession<FloatDP>, reads the ambient rounding mode,
sets the requested mode if different, and restores the
original mode on destruction. This also explains the
~30-ns directed arithmetic latencies visible for unrelated
FloatDP add/sub/mul/div/fma operations. Thus rounding-mode
save/set/restore is the leading **hypothesis** for the
remaining ~27-ns directed sqrt difference; this was not
confirmed by instruction-level profiling, and the hardware
rounding mode may already equal the requested mode in
other call patterns.

**Measurement caution.** Relative total_time changed
substantially but still reflects both backends and unrelated
run-to-run noise. The strong and targeted absolute square-root
improvement is the meaningful D4 result. To attribute the
remaining difference to fenv overhead, separately benchmark
native sqrt with CurrentRoundingMode in a single batched
RoundingSession and with per-operation save/set/restore,
using the same inputs. Any future D5 attempt that avoids
changing the hardware rounding mode must independently
prove correctly rounded directed results for all finite
normal/subnormal binary64 values, exact squares, tie cases,
signed zero and exceptional inputs. Do not assume a
round-to-nearest sqrt followed by nextafter is sufficient
without an exact comparison and underflow analysis.

**Status.** D4 was compiled and timed by the user to obtain
these values; this entry does not imply that GCC/Clang/MSVC
CI, the full independent source-keyed numeric oracle suite,
sanitizers, or coverage are passing. The preceding M148
GCC-13 -Werror=shadow repair remains a separate, unverified
CI fix. No D5 source modification is included in this
documentation-only update.


## D5 — enable certified nearest rounding for FloatDP log, tan and atan (29 October 2026)

**Observation and root cause.** The user-measured post-D4 FloatDP
p=53 benchmark (20,000 iterations, 31 samples, 56 operations)
reported:
- log_nearest=24.633 ns vs log_down=3.269 ns, log_up=3.398 ns;
- tan_nearest=40.606 ns vs tan_down=4.871 ns, tan_up=4.981 ns;
- atan_nearest=36.348 ns vs atan_down=4.704 ns, atan_up=4.754 ns.

These three nearest operations were *explicitly excluded* by the
initial guard in floatdp_certified_transcendental in
src/floatdp.cpp: even when their binary64 argument was in the
certified near-half window, nearest returned false and invoked
the much slower general transcendental implementation in
src/double.cpp. By contrast, the already certified directed
modes used the D1/D2/D3 all-native fast path. This is a path
selection discrepancy, not a consequence of an inherently
harder nearest-even rounding operation.

**Proof supporting the change.** Existing compile-time
floatdp_verify_certificates<Function>() static_asserts already
cover *all* near-half binary64 inputs admitted by the 53/54
denominator decoder, both input signs where defined, and all
nine supported transcendental centers. The exhaustive
certificate explicitly checks equality of truncated endpoints,
nonzero discarded remainder, no half-ULP ties for either
endpoint, and identical nearest-increment decisions. The
runtime floatdp_round_fixed already implements nearest
rounding by testing the exact half-bit, relying on the
compile-time proof to exclude ties. Thus the same
mathematical and binary64 enclosure already used for
downward/upward also certifies nearest for log, tangent and
arctangent. No new Taylor constants, tolerances, heuristic
thresholds or approximation algorithms are introduced.

**D5 source change.** Replace the initial five-function
allowlist for nearest with a deliberate exclusion of *only*
FloatDPTranscendental::exponential in nearest mode.
Consequently log, tan and atan now share the verified
certified near-half path with their directed modes and
with sin, cos, asin, acos and tanh nearest. Preserve the
original exponential-nearest exception as an independent
future candidate (no benchmark attribution mixing).
The log input sign check remains: negative x cannot use
log's positive-only near-half center. Out-of-window
arguments retain the unchanged fallback transcendental
implementations. The D1 fixed-point multiplication, D2
native endpoint arithmetic, D3 direct IEEE encoding,
D4 hardware sqrt, and all source-keyed oracle fixtures
remain untouched. No BigUnsigned dependency returns
to src/floatdp.cpp.

**Independent regression coverage.** In the existing
DoublePrecision golden-oracle verification block of
tests/test_float_rounding.cpp, switch the log, tan and
atan records from check_nearest=false to
check_nearest=true. The tests use the existing independently
generated FloatMP transcendental golden-oracle data and
verify both exact dyadic Rational value and complete
binary64 encoding, together with preservation of the
caller's active rounding mode. Leave the exp-nearest
oracle flag false to isolate D5. The existing exhaustive
consteval certificate tests remain active and unchanged.
Do not remove or weaken any independent-oracle checks.

**Benchmark gate.** The anticipated improvement is
specific to 53-bit inputs near one-half. Re-run the user's
unchanged 56-operation paired benchmark with 20,000
iterations and 31 samples. Compare FloatDP log_nearest,
tan_nearest and atan_nearest against 24.633, 40.606
and 36.348 ns, respectively; compare their directed
counterparts to detect any regressions and exp_nearest
to confirm it remains unchanged. Report FloatMP/FloatDP
total_time and geomean, but compare absolute per-operation
times as primary evidence of the intervention. No speedup,
full repository compile or CI success is asserted at
commit preparation. In particular, confirm the frozen
nearest golden oracles, GCC/Clang/AppleClang, MSVC,
sanitizers and coverage before considering this stable.

**Scope boundary.** D5 is only a removal of an unnecessary
nearest-mode fast-path exclusion for three operations.
Do not conflate it with the separate previously suggested
rounding-session optimization; that work remains
unimplemented and would require additional correctness
proofs for directed arithmetic.


## D6 — enable certified nearest for FloatDP exponential (29 October 2026)

**Scope and remaining exclusion.** After D5 allowed log, tan and atan
nearest into the short certified FloatDP transcendental kernel, the
initial guard in src/floatdp.cpp still excluded only
FloatDPTranscendental::exponential in ROUND_TO_NEAREST mode. The
same near-half fast path already served exp_down and exp_up. D6 removes
that last operation-specific nearest exclusion and leaves the
precisely bounded input decoder and general fallback untouched.
Consequently, *all nine existing certified operations* are now
eligible in nearest mode when their input lies inside the already
certified near-half window. This is not an extension to inputs
outside that window and does not alter unsupported domains.

**Correctness basis.** The existing consteval
floatdp_verify_certificates<FloatDPTranscendental::exponential>()
static_assert rigorously verifies every admitted binary64 input
in the two 53/54-denominator windows, both signs and the
rounding behavior of both error-bounded interval endpoints.
The proof includes nonzero discarded fractions, absence of
halfway cases, equality of retained significands and agreement
of the two endpoints' nearest increment decisions. The
D1/D2/D3 native-word fixed-point evaluator already supports
ROUND_TO_NEAREST via floatdp_round_fixed; only the dispatch
guard prevented its use. D6 changes no coefficients, derivative
bounds, rounding formula, binary64 encoding, input windows,
or special-case fallbacks. It introduces no BigUnsigned to
FloatDP. Outside the short certified domain, exp_nearest uses
the original general path unchanged.

**Independent regression coverage.** Switch the existing
independent source-keyed exp golden oracle invocation in
tests/test_float_rounding.cpp
(ARIADNE_FLOATMP_CENTERED_EXP_ORACLE_PATH, 568 rows)
from check_nearest=false to true. For the admitted
p=53 near-half arguments, the test already compares exact
dyadic Rational results, exact IEEE-754 binary64 encodings,
and the preservation of ambient rounding mode. This adds
coverage that was intentionally skipped for exp nearest
without regenerating, weakening or reinterpreting any
oracle fixture. The already-active consteval exhaustive
certificates remain unchanged.

**Baseline and measurement isolation.** In the latest
user-supplied D4 benchmark at 53 bits (20,000 iterations,
31 samples, 56 operations), FloatDP exp_nearest took
**22.935 ns** while exp_down and exp_up took respectively
**6.421 ns** and **6.544 ns**. The last reported paired
FloatMP/FloatDP total_time was **1.430x** and geomean
**2.207x**, but D5 was committed subsequently and its
own benchmark has not yet been supplied in this session.
Therefore those aggregate values are *pre-D5/pre-D6*
references, not a controlled baseline uniquely for D6.
A new run should compare absolute exp_nearest latency
against ~22.9 ns and the directed exp times against
~6.4–6.5 ns; additionally verify that log_nearest,
tan_nearest and atan_nearest exhibit the expected
D5 behavior. As D5 and D6 were not benchmarked in between,
a new global aggregate ratio measures the cumulative D5+D6
effect, not D6 alone.

**Verification status and release gates.** Checked the
source conditions and the presence of the unmodified
exponential static_assert, the original near-half
window restriction, the independent golden-oracle
harness, and absence of BigUnsigned from FloatDP.
This does *not* constitute an integrated Numeric build
or successful CI. Re-run GCC 13 strict -Wshadow/-Werror,
Clang/AppleClang, MSVC with strict floating point,
the complete independent binary64 oracle tests, all
modes and exceptional argument fallbacks, sanitizers,
and coverage checks. In particular avoid changing
the source-keyed fixtures or replacing the strong
consteval proof with benchmark-specific assumptions.

**Next work.** No F7 optimization is included. This
concludes the operation-specific nearest exclusions
of the existing FloatDP short transcendental certifier.


### D6.1 — user-measured combined D5/D6 nearest improvements (29 October 2026)

The user supplied the complete 56-operation FloatDP/FloatMP paired
benchmark at p=53 after both D5 and D6, using 20,000 iterations
and 31 samples. Relative FloatMP/FloatDP results:
**geomean=2.670x, total_time=1.761x, worst=11.054x**
(sqrt_nearest). The last user benchmark before D5 (after D4)
reported geomean=2.207x, total_time=1.430x and worst=11.385x
(sqrt_nearest). There is no intermediate measurement after
D5 but before D6, so the new global result includes both
changes and any concurrent benchmark-host variation.

Absolute FloatDP nearest latencies in nanoseconds, previous
post-D4 -> new post-D6, with reduction relative to the
earlier measured value:
- exp_nearest: 22.935 -> 6.815 ns (70.3% lower, D6 target).
- log_nearest: 24.633 -> 3.154 ns (87.2% lower, D5 target).
- tan_nearest: 40.606 -> 5.065 ns (87.5% lower, D5 target).
- atan_nearest: 36.348 -> 4.429 ns (87.8% lower, D5 target).

Their sum of four separate medians changed from **124.522 ns**
to **19.463 ns** (84.37% reduction); log/tan/atan alone
changed from **101.587 ns** to **12.648 ns** (87.55% reduction).
Those sums are not a contiguous combined benchmark duration,
but they provide a transparent absolute-latency comparison.
The corresponding directed operations in the new run remain
fast (exp_down/up 6.294/6.465 ns,
log_down/up 3.421/3.321 ns,
tan_down/up 4.754/4.717 ns,
atan_down/up 4.369/4.485 ns), broadly matching their
nearest-mode performance. This is strong workload-specific
evidence that D5/D6 removed the previous dispatch overhead.

**Important cross-run confounder.** Some operations not changed
by D5/D6 also shifted materially: for example FloatDP
construct_rational 37.175 -> 20.546 ns and copy 0.631 ->
0.417 ns, while FloatMP construct_rational moved 28.800 ->
30.065 ns. Changes of this scale can reflect a different
runtime machine, compiler, build, OS, thermal/frequency,
measurement setup, or other uncontrolled factors. Therefore
do NOT claim that the entire global total_time ratio change
from 1.430x to 1.761x is causally explained by D5/D6.
The four specifically targeted operations improved by much
larger, consistent factors and now have nearest latencies
comparable to their own directed modes in the *same run*.

**Verification and next step.** The user-supplied benchmark
is a timing result, not evidence of completed full test
execution, CI, independent source-keyed oracle success,
GCC/Clang/MSVC compatibility, or sanitizers/coverage.
Confirm all the nearest source-keyed oracles activated in
D5/D6, retain their existing rigorous compile-time fixed85
certificates and re-run the paired benchmark on a controlled
same-host/same-binary setup for final attribution. The
FloatDP near-half code still has no BigUnsigned usage.
No F7 optimization or changes to numeric sources are
included in this journal update.

## D4.2 — AppleClang strict test-build diagnostics (29 October 2026)

A user-supplied macOS AppleClang Debug build log failed in the
D4 source-keyed FloatDP square-root checks in
`tests/test_float_rounding.cpp` with `-Werror -Wshadow` and
`-Wvexing-parse`. Three independent test-source issues caused
the diagnostics:

- The local `DoublePrecision dp` declaration shadows the existing
  `Ariadne::dp` constant. Rename it to `sqrt_precision` and update
  the D4 test's constructors/special-value calls.
- The local `FloatDP toward_zero` declaration shadows the existing
  `Ariadne::toward_zero` rounding-mode token. Rename it to
  `toward_zero_result` and update its single bitwise check.
- `const Rational x(ExactDouble(exact_input));` is parsed as a
  function declaration (most vexing parse). Use brace initialization
  `const Rational x{ExactDouble(exact_input)};`. The subsequent
  arithmetic comparison errors are downstream consequences of this
  incorrect declaration.

The exact test input set, all four rounding modes, mathematical
inequalities, nearest-even midpoint check, infinity/NaN and signed-zero
checks, and production algorithms remain unchanged. This repairs
**compilation diagnostics in the tests**, not a failing numeric oracle.
The full AppleClang build, CTest and other CI jobs have not been rerun
or asserted successful in this documentation update. This record is
included in the new backend-specific FloatDP journal.
