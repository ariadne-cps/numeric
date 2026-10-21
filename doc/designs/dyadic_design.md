# NewDyadic implementation and optimization journal

**Branch:** `remove-gmp`.
**Companion core journal:** [`integer_rational_design.md`](integer_rational_design.md).
**Experimental baseline:** [`0f7f06c16`](https://github.com/ariadne-cps/numeric/commit/0f7f06c16), R33a (21 October 2026).

This journal tracks the native GMP-free `NewDyadic` experiment, exact arithmetic compatibility with the legacy `Dyadic`, performance on ARM64 and x86-64, and associated shared `BigUnsigned` changes. Preserve existing exact-value differential coverage and the unchanged paired benchmark workload. A candidate optimization is not promoted on benchmark intuition alone: compile, test, compare against exact oracles, and report actual same-host timings. Shared substrate design and Integer/Rational historical changes remain in the companion journal.

### R33a — introduce experimental GMP-free NewDyadic (21 October 2026)

Commit `0f7f06c16` introduces the initial native `NewDyadic` implementation. It is experimental and separate from GMP-backed `Dyadic`; migration is not implicit. Subsequent entries record API parity, dedicated tests, benchmark design, platform issues, and performance changes.

### R33b — paired Dyadic/NewDyadic experiment (21 October 2026)

Keep existing Dyadic/mpf untouched. The proposed NewDyadic candidate owns only an Integer coefficient, nonnegative power-of-two denominator exponent and finite/infinity/NaN kind. It normalizes powers of two after every operation; add, sub, mul, compare, half, round, floor and ceil are exact native arithmetic. A provisional Rational conversion permits independent exact value comparison. The explicit CMake benchmark_compare_dyadic_backends target reports 18 matched operations per width, default 53,128,256,1024,2048,4096 bits with 20000 iterations and 31 samples, paired input pool of 32 deterministic signed dyadics and exact off-timer verification against legacy Dyadic converted to Rational. Report old ns, new ns, New/Old ratio; this is NOT a GMP comparison. No benchmark CI jobs. NewDyadic currently does not replace the full Dyadic public API or its legacy MPFR bridge; promote only after strict correctness and equivalent behavior on edge cases.


### R33d — align NewDyadic public arithmetic API with Dyadic (21 October 2026)

User requested signature/semantic parity and that NewDyadic tests derive from the original test_dyadic.cpp, not an unrelated minimal test. Restore Boolean comparison returns, Comparison cmp, Sign sgn, Boolean eq/lt, Bool special predicates, Positive<NewDyadic> square and abs results, named arithmetic, compound assignments, signed/unsigned pow distinction, numerator/exponent queries, exact String/ExactDouble/TwoExp construction, and canonical native storage. The legacy Dyadic implementation remains unchanged. The paired benchmark now times `cmp` on both sides and `sqr` on both sides. Retain exact Rational oracle outside measurements. Port test_dyadic's literal/conversion/arithmetic/comparison/special sections into the independent NewDyadic test; DyadicBounds/Ball and raw mpf interop still require separate compatibility work before migration/promotion.


### R33e — port original Dyadic test scenarios to NewDyadic (21 October 2026)

**Test provenance.** Rewrite tests/test_new_dyadic.cpp as a parallel
section-by-section adaptation of tests/test_dyadic.cpp: test_literal,
test_conversions, test_arithmetic, test_comparisons, test_infinity.
Carry over exact literal values 3.25, -11.375, 0.375, 2^3/2^-3;
invalid/non-dyadic decimal rejection; all 11 legacy rounding
expectations; negative floor/ceil and infinity exceptions; arithmetic
3/4 +/- (-5/8), signed positive/negative powers, compound updates;
negative comparisons; infinities, NaN, signs and invalid operations.
Add exact cross-checks with the legacy Dyadic and native Rational at
53,128,256,512,1024,2048,4096 bits. Verify public Boolean and
Comparison return signatures at compile time. Do not weaken a
failing numeric expectation or silently substitute a different
numeric rounding rule. Keep the previous test_dyadic.cpp untouched.

**Deliberate gaps before API parity.** The source test sections
test_bounds and test_ball require NewDyadic-specialized bounds,
balls, positive and abstraction wrappers, and are NOT yet ported.
The existing legacy mpf constructor/accessor, default Writer API,
two-exp literal suffixes and full stream interface likewise
remain pending separate compatibility work. This is a test-first
gap inventory, not a claim that NewDyadic already replaces Dyadic.
No benchmark or Rational CI job is re-enabled.


### R33f - fix AppleClang native NewDyadic construction (21 October 2026)

The ARM64 paired benchmark failed to build with incomplete ExactDouble in NewDyadic(double) and ambiguous overload selection for unsigned exponents: a deleted, non-template signed Int overload was an equally ranked implicit conversion as size_t. Include numeric/builtin.hpp (the full ExactDouble definition) and constrain the deleted signature using the BuiltinSignedIntegral concept. Signed arguments remain deliberately rejected; unsigned calls select the existing size_t constructor. No change to reference Dyadic, benchmark fixtures, arithmetic or correctness assertions. Previous Unix and Coverage R33e builds failed during compilation, so validate the next revision independently.


### R33g - Dyadic CLI/report alignment (21 October 2026)

The NewDyadic/Dyadic comparison adopts benchmark_command_line.hpp, Ariadne CommandLineInterface, logging_command_line, and the same parser validation/help behavior as the Integer and Rational comparison executables. Named options: -i/--iterations, -n/--samples, -b/--bits, --summary-only, --details, --max-slowdown, plus logging options. Positional ITERATIONS SAMPLES stay supported for backwards compatibility. Defaults remain 20000, 31 and 53,128,256,1024,2048,4096 bits. The benchmark reports matched absolute medians, optional operation ratios, grouped and overall geomean, weighted total_time and worst NewDyadic/Dyadic ratio; max-slowdown checks the worst observed operation ratio. Do not interpret the legacy Dyadic backend as raw GMP. The fixture, measured operation callbacks and off-clock exact checks are unchanged. No benchmark CI workflow is restored; compilation and runtime still need verification after this change.


### R34 - direct native binary64 and virtual-shift compare (21 October 2026)

The user-supplied ARM64 Dyadic/NewDyadic 20000x31 comparison has 18 operations at each width and the unchanged paired fixture check passes. NewDyadic is substantially quicker at construction, multiplication, squaring and integer rounding. The dominant regressions are to_double (215.242 ns versus 2.856 ns at 4096 bits) and cmp (96.508 ns versus 3.577 ns). The attachment described as x86-64 actually contains a Rational/GMP comparison, so it is NOT a Dyadic result. Replace NewDyadic get_d's owning Rational construction with direct exact 53-bit significand extraction reading at most two limbs; preserve truncation, exponent overflow, subnormal underflow and the legacy mpf_get_d positive-zero convention on negative underflow. Replace cmp's owning shifts with virtual radix-2^64 mantissa alignment and a constant-time high-exponent/sign comparison. Add bitwise GMP-backed legacy Dyadic/get_d checks across normal/subnormal/overflow offsets and exact old/new cmp checks across widths and signs. A standalone 120000-case GMP random test verified nonzero conversion results and discovered positive zero on legacy underflow, which the production candidate follows. Benchmark corpus and public Dyadic/NewDyadic interfaces are unchanged, Rational stays frozen and Rational benchmark CI stays removed. CI and independent full performance repeats remain necessary before accepting R34.


### R35 — virtual-shift, single-output native Dyadic add/sub (21 October 2026)

**Paired observations.** User supplied two independent Dyadic/NewDyadic
runs at 20,000 iterations, 31 samples and six widths, with exact
verification PASS. At 4096 bits legacy versus NewDyadic add/sub
times are 71.027/192.369 and 70.302/231.856 ns on the first
reported system, and 98.457/250.992 and 97.097/271.727 ns
on the second. Their labels "amd64" and "x86-64" are architecture
synonyms, not proof of distinct ISAs. The benchmark's comparison
and binary64 conversion are now substantially improved by R34.

**Implementation.** Preserve the original Dyadic and public
NewDyadic API. On finite add/sub operands with 5..64 maximum
magnitude limbs, each individual size <=64 and binary exponent
gap <=63, replace the construction of two shifted owning
Integer temporaries, the signed sum, and subtraction's negated
intermediate with a single 66-limb stack accumulator. Each
input's shifted radix-2^64 digit is read from borrowed canonical
BigUnsigned limbs (at most two input reads); no source is mutated.
For opposite-sign contributions, determine the exact larger
virtual magnitude from bit lengths and then high-to-low limbs,
choose the result sign, and subtract with borrow. For matching
signs, add with carry and retain the final carry limb.
Canonicalize one final owning magnitude and signed Integer.
All nonfinite, smaller 1..4-limb, >64-limb, or exponent-gap
>=64 cases retain their exact previous algorithm. This
retains exact cancellation, zero, infinity, aliasing and
all special-value semantics.

**Verification.** A separate C++20 Boost.Multiprecision oracle
compiled with -Wall -Wextra -Werror verified the underlying
radix-2^64 aligned addition/subtraction on 200,000 varied random
up-to-64-limb signed operand pairs including cancellation;
this is not an integrated project test. Add ported
test_new_dyadic.cpp legacy-Dyadic/Rational cross-checks at
53,64,65,127,128,255,256,257,320,511,512,1024,
2048,4096,4097 bits and gaps 0,1,7,31,63,64,65,
all four sign combinations, carry, borrow, reverse
subtraction, += and -=. Existing test_dyadic is unchanged,
as are all benchmark operations, sample parameters and
the removal of benchmark CI. Independent Unix/Windows/
Coverage and same-host performance comparisons are required;
no benefit is claimed from unmeasured code.


### R36 — architecture-aware add/sub kernel crossover (21 October 2026)

**Host identity correction.** The user explicitly clarified that
the first paired Dyadic/NewDyadic benchmark host is ARM64 and the
second is x86-64 (the earlier wording "amd64" for the first run
was a mislabel). R35 numerical verification PASS on both hosts,
with unchanged 384 input pairs, 18 operations, 20,000 iterations
and 31 samples.

**Measured tradeoff.** ARM64 R34 -> R35 at 4096 bits:
add 192.369 -> 178.894 ns, sub 231.856 -> 177.590 ns.
At 1024 bits add 90.152 -> 66.335 ns and sub
112.017 -> 65.581 ns. On x86-64 R34 -> R35 at
1024 bits add 135.814 -> 123.249 ns and
sub 153.144 -> 123.163 ns; at 2048 bits add
172.223 -> 191.680 ns and sub 188.964 -> 194.720 ns;
at 4096 bits add 250.992 -> 316.049 ns and sub
271.727 -> 321.672 ns. Thus the exact R35 fused
accumulator is beneficial to 64 limbs on ARM64, but
the x86-64 crossover lies between 16 and 32 limbs.

**R36 policy.** Keep R35's exact shared virtual-shift
algorithm unchanged. Use it for finite operands with
at least one magnitude of five or more limbs, a maximum
exponent gap of 63 bits, and *both* operand magnitudes
up to 64 limbs on ARM64 (__aarch64__/_M_ARM64)
or up to 16 limbs on x86-64 and other architectures.
All other cases take the existing portable owning-shift
fallback. This is a compile-time machine architecture
policy, not an input-fixture lookup, compiler-warning
workaround or semantic change. The cutoff is a measured
heuristic and must be revisited if other x86 compilers
or machines contradict it.

**Coverage.** Extend test_new_dyadic's test port with
legacy-Dyadic/Rational add/sub/reverse-sub checks at
959/960/961, 1023/1024/1025, 1087/1088/1089 and
4031/4032/4033, 4095/4096/4097, 4159/4160/4161
bit boundaries, both signs and exponent gaps 0/1/63/64.
Retain previous randomized oracle and all R35 tests.
Keep Dyadic, Rational, the 18-operation paired benchmark,
its CLI, CI policy and vcpkg/Windows linking unchanged.

**Acceptance.** R36 is a provisional dispatch selection.
After Unix/Windows/Coverage correctness checks, compare
on the *same physical ARM64 and x86-64 machines* to R35.
Expected directions are ARM64 parity with R35 and x86-64
recovering the previous 2048/4096 add/sub timings.
Do not infer measured speedup from the cutoff itself.


### R37 — constant-time canonical fractional check in NewDyadic floor/ceil (21 October 2026)

**Measured context.** R36 ARM64 and x86-64 paired measurements,
20,000x31, keep NewDyadic significantly faster than original
Dyadic in the total weighted 18-operation workload
(ARM64 total_time .402x, x86-64 .404x). ARM64 R36 4096-bit
add/sub = 174.365/173.502 ns vs R35 178.894/177.590 ns;
x86-64 R36 4096-bit add/sub = 271.830/290.305 ns
vs R35 316.049/321.672 ns. The strict tests Unix,
Windows and Coverage are all green on R36. The user
confirmed the first host is ARM64 and the second x86-64.

**Canonicality proof.** NewDyadic's finite invariant is
coefficient/2^e with coefficient odd whenever e>0:
normalize() divides out all common powers of two and
sets e=0 for zero. Therefore for every canonical
nonzero finite NewDyadic, e>0 if and only if it has a
nonzero fractional part. The original has_fraction()
called shift_right, shift_left and BigUnsigned::compare,
temporarily materializing full-width integers solely to
rediscover this invariant. In floor/ceil, replace
has_fraction(magnitude,e) with e!=0, and delete the
redundant helper. Keep the existing truncated sign
correction, exact Integer output, special-value errors
and public overload signatures unchanged.

**Test coverage.** Extend the original Dyadic-derived
NewDyadic test with legacy Dyadic oracle floor/ceil
at 16 widths from 1 to 4097 bits, positive/negative,
six low-word patterns (including even and odd inputs),
nine binary exponent offsets, plus zero and exact
powers of two up to 4096 bits to validate canonical
factor removal. Verify floor<=ceil and
(floor==ceil iff normalized denominator exponent==0).
No changes to tests/test_dyadic.cpp, the paired corpus,
Rational, build/link config, vcpkg or benchmark CI.
R37 is a general semantic simplification, not yet a
measured performance gain; require independent CI and
native same-host benchmark results before promotion.


### R37.1 — compile fix for mixed Boolean/bool test comparison (21 October 2026)

All R37 Unix/Windows/Coverage workflows failed at the compilation step. The user-supplied AppleClang diagnostic pinpoints tests/test_new_dyadic.cpp:187: a comparison of Ariadne Boolean (floor(native)==ceil(native)) to built-in bool integer causes ambiguous operator== overload selection. Explicitly convert the Ariadne Boolean to bool before comparing to the canonicality expectation. Keep the assertion, original Dyadic oracle tests, all R37 arithmetic, benchmark fixture, vcpkg, build workflows and API unchanged. Both R37 user-supplied native benchmark runs passed their 384-input fixture numerical verification and showed floor/ceil improvements; these are independent of integrated CI. Revalidate Unix/Windows/Coverage before considering R37 stabilized.


### R38 — one-allocation shifted signed add/sub for NewDyadic (21 October 2026)

**User-defined goal.** Focus exclusively on add and sub and lower the
maximum measured NewDyadic/Dyadic slowdown below 2 on ARM64 and x86-64.
The R37.1 user benchmarks (20,000x31, unchanged 18 operations and 384
signed pairs) returned PASS and an overall total_time of 0.374x
ARM64 and 0.368x x86-64. Their 4096-bit add/sub ratios were
175.535/67.027=2.619x and 174.450/66.727=2.614x on ARM64,
versus 254.798/95.906=2.657x and 274.429/96.617=2.840x
on x86-64. All Unix/Windows/Coverage R37.1 checks succeeded.

**Proposed algorithm.** Old fallback forms two full shifted owned
Integer operands and one signed result; R35's accelerated kernel
instead uses a 66-limb stack result and then allocates and copies
it into a canonical owned BigUnsigned. Both incur avoidable
large-width work. Add two *new*, isolated BigUnsigned APIs,
add_shifted and sub_shifted: for shifts 0..63, read two borrowed
input views with virtual left alignment, calculate carry/borrow
through one radix-2^64 loop, and write directly into exactly
one result-owned BigUnsigned allocation; subtract rejects a
negative result. The preexisting BigUnsigned arithmetic kernels
and benchmark corpus are unchanged. In NewDyadic use signed
virtual-magnitude compare_scaled to select direction and sign;
then construct exactly one canonical native Integer/2^exponent
result. Route large finite operands with difference<=63 bits
to this kernel: largest operand >=32 limbs on ARM64 or >=17
limbs on x86-64 and other architectures. Shorter inputs retain
R35 stack or original paths, preserving their measured behavior.
No changes to public Dyadic, Rational, Positive wrappers,
multiplication, benchmark CLI, or Windows/vcpkg.

**Exact oracles.** Expand test_big_unsigned.cpp's preexisting
GMP oracle for both shifted operations across lengths
0/1/4/16/17/31/32/33/64/65/129 and shifts
0/1/7/31/63 on both operands, including all-ones,
signed direction and underflow. Extend the original
Dyadic-derived NewDyadic tests with direct old-Dyadic/
Rational add/sub/reverse-sub/+=/-= checks at the
17- and 32-limb boundaries, 4096 and 8192 bits,
both signs and exponent gaps including the 63/64
crossover. Ensure the new kernel never mutates
borrowed operands or weakens correctness assertions.

**Status.** R38 is a candidate, not a claimed success.
Run Unix/Windows/Coverage and remeasure both machines with
the same 20,000x31 workload. The strict goal is worst<2x
for every operation and precision on both ARM64/x86-64;
do not alter benchmark cases or thresholds to satisfy it.

### R38.1 — fix AppleClang macro parsing in BigUnsigned shifted-subtraction test (21 October 2026)

Commit `8f7ea1992` wraps the lambda expression passed to `ARIADNE_TEST_ASSERT` in an extra pair of parentheses. The comma separating the two `BigUnsigned` local declarations inside the lambda had been interpreted by the preprocessor as an additional function-like macro argument under AppleClang. The added parentheses group the expression as one macro argument; the underflow test semantics and checked expectations are unchanged. CI rerun required before marking the fix validated.

### R39 — avoid redundant temporary copy in `abs` and `sqr` (21 October 2026)

Add a checked rvalue overload `Positive<NewDyadic>(NewDyadic&&)` next to the existing const-reference constructor. This lets the negative `abs` path consume the already materialized `-x` temporary without copying its potentially large native Integer coefficient again; the positive `abs` path still copies the caller-owned lvalue, preserving nonmutation, and `sqr` can likewise move its multiplication result. Validation of the positivity invariant is performed on the constructed object. Keep the public signature of `abs`, exact sign/special semantics, original Dyadic and benchmark workload unchanged. This is a structural allocation/copy reduction candidate, not a measured ARM64 gain. Run compilation and full regression tests plus native same-host ARM64 benchmark before accepting a performance improvement.

### R40 — zero-shift direct add/sub path for large NewDyadic on x86-64 (21 October 2026)

When both finite operand denominator exponents coincide, R38's selected owning add/sub path has two zero shifts. Use the existing exact `BigUnsigned` `add`/`sub` in that case, in either subtraction direction selected by `compare_scaled`, instead of the more general virtual-shift combiner. Retain all previous signed-order and zero/cancellation logic. When the denominator exponents differ, keep R38's shifted kernel and all original crossover thresholds. This is an exact fast-path candidate particularly relevant to the x86-64 large-operand path, not a proven speedup; no ISA-specific instruction set is required. Run full correctness and the unchanged paired 20,000x31 benchmark on both hosts, including nonzero shifts and boundary cases, before making performance claims. Do not change benchmark corpus, public Dyadic or Rational.

### R40.1 — fix GMP oracle accessor calls in R38 tests (21 October 2026)

The R38 shifted add/sub tests in `tests/test_big_unsigned.cpp` incorrectly used `shifted_a.value()` and `shifted_b.value()` in two GMP left shifts and two GMP inverse right shifts. Both are `mpz_t` data members of `GmpInteger`, not callable functions; passing `.value` directly corrects the compile error. The test matrix, assertions, reference arithmetic, production kernels, and benchmark workload are unchanged. Confirm Unix, Windows, and Coverage CI after this fix; the implementation alone is not evidence of a green build.

### D41 — x86-64 streamed first-operand shift in owned add/sub (21 October 2026)

**Scope and motivation.** The official unchanged NewDyadic/Dyadic
20,000 x 31 paired workload constructs operands with exponents e and
e+1. Canonicalization of even coefficients can change that relationship,
so the R40 equal-exponent specialization is used for part, not all,
of the signed add/sub samples. The remaining unequal-exponent large
operands often use R38 with the first coefficient virtually shifted.
The R38 loop rereads an earlier limb and tests both shift counts in
every result iteration.

**Candidate.** On x86-64 only, when shift_a is in [1,63] and
shift_b==0, accumulate the shifted first input with a forward spill
register. Each iteration reads the current first input limb at most
once, combines its left-shifted word with the carried high bits, and
uses the second operand unshifted; the original exact carry/borrow
combination remains unchanged. This preserves one owned result, input
nonmutation and the previous underflow behavior. All other shift
arrangements and ARM64 preserve the R38 path. It is a direct shared
BigUnsigned kernel optimization, not an x86 inline-assembly feature.

**Independent preliminary diagnostic (NOT integrated performance).**
A standalone equivalent C++20 radix-2^64 kernel probe compared the
original and streamed loops on GCC 14.2 and Clang 17 x86-64, with
same inputs and a generic-versus-streamed differential PASS across
widths 0..129, shifts 0,1,7,31,63, and both addition and subtraction,
including the underflow exception. In synthetic 1-bit *first-operand*
shift addition, the streamed/generic timing ratios were approximately
0.79..0.97 on Clang and 0.72..0.88 on GCC at 16/32/64/128 limbs.
The second-operand-shift case regressed on some GCC widths, so that
direction was not selected. These diagnostics do not measure the actual
NewDyadic constructor, allocation, normalization or full benchmark;
there is no validated user-visible speedup claim.

**Validation gate.** R38's retained independent GMP shifted add/sub
oracle already exercises 0/1/4/16/17/31/32/33/64/65/129-limb
operands, shifts 0/1/7/31/63 on either operand, all-ones, zero,
carry, borrow and negative-difference rejection. Run Unix, Windows
and Coverage, then compare the unchanged 18-operation NewDyadic/Dyadic
20,000 x 31 benchmark on the same physical x86-64 host to R40.1;
recheck ARM64 parity. Do not use a synthetic isolated ratio as
promotion evidence. No public Dyadic/Rational change, no new benchmark
CI job, and no alteration to the benchmark fixture.

### D42 — direct canonical absolute value without redundant materialization (21 October 2026)

**ARM64 target.** User-supplied unchanged 20,000-iteration, 31-sample
paired NewDyadic/Dyadic benchmark reports abs at 4096 bits of 76.419 ns
versus 36.298 ns, a 2.105x slowdown; previously it was 74.656 ns
versus 36.767 ns, 2.031x. The overall ARM64 worst<2x gate therefore
remains unsatisfied, despite the exact 384-pair benchmark verification
passing. Do not infer D41 caused the small change without matched runs.

**Change.** The old abs conditionally generated -x (copying the signed
Integer coefficient), then constructed Positive<NewDyadic> from a
materialized NewDyadic. Instead, as the already declared friend of
NewDyadic, abs copies the canonical unsigned magnitude once into a
nonnegative Integer and assigns the source denominator exponent and
kind directly, mapping only negative infinity to positive infinity.
Because the original finite representation is canonical, this avoids
a second normalization/copy. The already available checked rvalue
Positive constructor consumes the result without copying its owned
magnitude. The source is never mutated. Preserve positive and negative
finite values, zero, NaN, both infinities, and all original APIs.

**Correctness.** Add test_wide_abs: old Dyadic reference and exact
Rational comparison across 1/53/128/256/1024/2048/4096/8192-bit
magnitudes, different odd low words, five exponents and both signs;
also check canonical exponent, result sign, source immutability, zero,
infinities, and NaN. The preexisting benchmark fixture and test suites
are unchanged.

**Promotion gate.** Require Unix/Windows/Coverage CI and fresh native
ARM64 20,000x31 comparison on the same workload; worst slowdown must
be strictly below 2x. A C++ source change or successful CI alone does
not establish that performance condition. Do not declare NewDyadic
finished until the gate is measured.

### D43 — standalone native dyadic timing benchmark (21 October 2026)

The user has closed the isolated NewDyadic performance campaign and authorized
promoting it into Dyadic, retiring the legacy GMP-backed implementation and
its comparison test. Promotion requires porting downstream bindings and
rounded-float/ball/bounds and writer APIs before deleting the old class:
the experimental NewDyadic implementation alone does not provide these.

As the first independent, compilable step, replace
`benchmark_compare_dyadic_backends.cpp` with
`benchmark_dyadic.cpp` and update the benchmark CMake target.
The new executable measures the same 18 arithmetic and construction
operations at the same six default widths with the same 20,000
iterations and 31-sample medians; it reports absolute ns values and
per-width aggregate sums only. Remove Dyadic/mpf reference timings,
ratios, geomeans, slowdown gates and comparative CLI flags. The
native Rational exact arithmetic checks remain outside the timing
regions; later migration will rename the temporary `NewDyadic`
type identifier to `Dyadic`. CI and benchmark compilation remain
required; this step is not the complete class migration.

### D44 — remove direct Dyadic mpf access on FloatDP/FloatMP input (21 October 2026)

As part of the authorized retirement of the GMP-backed Dyadic,
FloatDP assignment now uses the exact type's `get_d()` conversion
instead of `mpf_get_d(w.get_mpf())`. FloatMP construction and
assignment use the exact mathematical `Rational(w)` through the
existing temporary `GmpRationalBoundary` and MPFR's `mpfr_set_q`,
with explicit handling of NaN and infinities on assignment.
This removes all direct reads of Dyadic's GMP representation from
these input paths without changing MPFR's implementation or
using `double` for multiprecision conversion. The temporary
boundary is GMP-backed because FloatMP/MPFR remains a separate
legacy floating-point backend, not because Dyadic owns GMP data.
Exact output conversions FloatMP->Dyadic and FloatAP->Dyadic still
need migration before old Dyadic files can safely be removed.
Do not conflate this intermediate dependency elimination with
completed promotion or green CI.

### D45 — narrow Dyadic migration to FloatAP, not FloatMP (21 October 2026)

**Scope correction.** FloatMP is scheduled for removal, with FloatAP
remaining as the arbitrary-precision float backend. Consequently, do
not spend migration effort rewriting FloatMP <-> Dyadic conversions.
Revert only the FloatMP input-conversion changes introduced in D44,
restoring `src/floatmp.cpp` to its pre-D44 content pending removal
of that source and its build, test, and API references. The native
`FloatDP` input adjustment from D44 is retained for now.

**Promotion requirement.** Before retiring the old Dyadic class, move
the supported external APIs and conversions to the native replacement:
FloatAP -> Dyadic must preserve the *exact* finite binary rational
(significand and binary exponent), and Dyadic -> FloatAP must respect
its precision and rounding contracts. FloatMP-specific overloads,
storage access, and tests should be removed with FloatMP, not ported.
No Dyadic `get_mpf()`, `_mpf`, or GMP `mpf_t` shim should survive
the final replacement. Keep independent correctness checks for
FloatAP conversion; removing tests does not establish correctness.

**Status.** This is a scope correction, not completion of FloatMP
removal or Dyadic promotion. Existing legacy FloatMP code is
unchanged relative to D43 and still depends on the legacy Dyadic
representation until its separate deletion.

### D46 — exact FloatAP-to-Dyadic export without mpf intermediates (21 October 2026)

FloatAP finite values already carry an integer significand and a signed
binary exponent. Replace the former `mpf_init2 / mpf_set_z /
mpf_mul_2exp / mpf_div_2exp / Dyadic(mpf_t)` conversion with
`integer_from_gmp(_significand)` to a native BigUnsigned, followed
by either an exact left shift for nonnegative exponents or an exact
Dyadic(integer, denominator-exponent) construction for negative
exponents. Preserve the independent sign bit, canonical zero,
NaN and signed infinities. Avoid an implicit narrowing for negative
exponents beyond the current legacy Dyadic Nat constructor limit:
report overflow explicitly; the final native Dyadic API will use
size_t exponents, removing that legacy-only restriction.

Add a D46 exact round-trip and Rational oracle in `test_float.cpp`
for 1/53/127/128/191-bit significands, five exponent choices, both
signs, and special values. The existing FloatAP->Rational conversion
remains a separate exact oracle. This changes no FloatMP adapter and
does not bring back an `mpf_t` Dyadic compatibility API. Because
the final native rename, downstream API transfer and FloatMP removal
have not occurred, this is an intermediate migration commit, not a
claim that the old Dyadic backend has already been removed.

### D47 — parallel native FloatAP/NewDyadic interface before final rename (21 October 2026)

Introduce direct conversions in both directions between FloatAP and NewDyadic
without using legacy Dyadic or mpf_t. The input side converts the exact
Rational value and distinguishes exact construction (rejecting inexact
precision) from explicitly rounded construction. Assignment preserves the
exact-conversion contract. The output side exports FloatAP's integer
significand and signed binary exponent as an exact native dyadic, including
NaN and signed infinity. Add a test matrix of 1/53/127/128/191 bit operands,
five exponents, both signs, exact Rational verification, assignment,
explicit rounding and nonfinite values. These overloads are transitional:
when NewDyadic becomes Dyadic, eliminate the old FloatAP/Dyadic overloads
and retain just these native implementations. Do not migrate FloatMP.
Validation through Unix/Windows/Coverage CI is required.

### D48 — stage native interval and enclosure wrappers (21 October 2026)

Extract the existing Dyadic bounds/ball/lower-bound/upper-bound/
approximation type declarations into an independent
`numeric/new_dyadic_extensions.hpp` for the experimental native type,
preserving value/raw access, lattice operations and exact enclosure
semantics. This is intentionally a parallel API, not a rename of the
production Dyadic type. The inherited FloatMP-specific constructors
and `get(MultiplePrecision)` declarations are omitted: FloatMP will
be retired, not adapted. The `operator ValidatedNumber()` conversion
is also deferred until the final name migration because its
NumberWrapper/ExactNumber instantiation requires replacing the
production class.

Provide exact four-vertex interval multiplication and parity-aware
integer powers directly, rather than relying on legacy GMP-bound
source definitions. Add focused tests in `test_new_dyadic.cpp` for
arithmetic enclosures, crossing-zero square and powers, subset
refinement, ball conversion, and lower/upper/approximation wrappers.
No storage or arithmetic function references GMP. This stages the
interfaces; it is NOT the final removal of `dyadic.hpp`,
`dyadic.cpp`, legacy test files, or FloatMP, and CI must validate
the changes before declaring this layer ready.

### D49 — revert D48's separate extensions header (21 October 2026)

D48 introduced `include/numeric/new_dyadic_extensions.hpp`, a parallel wrapper header with no counterpart in the existing Dyadic organization. This was an architectural deviation and is reverted in full: delete the new header and restore `tests/test_new_dyadic.cpp` to the D47 version. Preserve the D48 entry above for historical traceability, but treat its implementation as **retracted**.

The replacement must use the existing file layout exclusively: `include/numeric/dyadic.hpp` for the native Dyadic class, Positive specialization, writers, and bounds/ball/lower/upper/approximation declarations; `src/dyadic.cpp` for required out-of-line definitions; `tests/test_dyadic.cpp` for the retained and updated tests. Delete `include/numeric/new_dyadic.hpp` and `tests/test_new_dyadic.cpp` at final promotion. Do not create `dyadic_extensions.hpp` or similar files; do not port FloatMP adapters, which are scheduled for retirement. Existing mathematical coverage must be retained through the final renamed test.

### D50 — integrate native enclosure wrappers into the existing Dyadic header layout (21 October 2026)

D48 was reverted in D49 because a separate `new_dyadic_extensions.hpp`
file was not part of the project's Dyadic layout. Restore its already
CI-tested native wrapper implementation in the **same native class
header**, currently `include/numeric/new_dyadic.hpp`; this is to be
renamed in-place to `include/numeric/dyadic.hpp` during final
promotion. No additional headers or new source files are created.
Retain the D48 interval/ball/lower/upper/approximation tests in the
existing native test file, but remove its include of the rejected
extensions header. The existing wrapper specializations contain
no accesses to GMP Dyadic storage and no FloatMP-specific overloads.

Production `dyadic.hpp`, `dyadic.cpp` and `test_dyadic.cpp`
remain unchanged until final type replacement; removal of these
three legacy bodies and final rename remain explicit follow-up work.
A passing D48 CI validates only the former split-header variant;
Unix/Windows/Coverage must validate the new in-header arrangement.

### D51 — eliminate all FloatMP access to Dyadic's GMP representation (21 October 2026)

FloatMP is scheduled for deletion, so no new Dyadic/GMP compatibility
API is added. Delete the `mpfr_set_f(w.get_mpf())` paths from
`FloatMP(Dyadic,...)` and its assignment; the existing Rational/MPFR
boundary is reused only for the residual legacy API until FloatMP is
removed. Delete `mpfr_get_f(res._mpf)` and raw `_mpf[0]` special-value
tags in the reverse conversion. Its remaining legacy operator routes
the *exact* binary value through Rational and reconstructs
the dyadic coefficient/denominator, with explicit NaN and infinity
handling. Crucially FloatMP itself no longer knows the storage of
Dyadic or needs `get_mpf()`.

The source-of-truth FloatMP->Rational conversion is now a direct exact
`mpfr_get_z_2exp` extraction of the MPFR significand and signed binary
exponent, handled with native Integer/Rational; Rational(FloatMP) uses
that conversion, not the former FloatMP->Dyadic->Rational chain.
`cast_integer(FloatMP)` also works directly on Rational, without
round-tripping through Dyadic. Preserve the old public FloatMP
signatures solely because the library and generic tests still
instantiate them; they should be deleted wholesale as part of the
forthcoming FloatMP removal, not migrated into native Dyadic.
No changes to FloatMP representation, test fixtures, or benchmark
workloads. CI must confirm exactness, including finite/special cases.

### D52 — parity of Rational/Integer/Dyadic interoperability with FloatAP (21 October 2026)

Before replacing the public Dyadic backend, maintain parity between the
soon-to-be-removed FloatMP interfaces and the retained FloatAP backend.
`Rational(FloatAP const&)` now mirrors the existing explicit
`Rational(FloatMP const&)` constructor, delegating to FloatAP's exact
`operator Rational()` rather than using double or any Dyadic/MPFR
storage conversion. FloatAP already provides `cast_integer(FloatAP)`,
`FloatAP(Integer,...)` and both exact and explicitly rounded Dyadic
conversion APIs; no duplicate Integer overload is added. Add
independent tests of finite/large dyadic rational values, both
infinities, NaN, FloatAP->Dyadic roundtrip, and exact/non-integral
FloatAP->Integer casts. The production Dyadic enclosure getter APIs
accept MultiplePrecision and will be matched with ArbitraryPrecision
during the final native promotion, without preserving FloatMP-specific
adapters. This is an API parity step, not final Dyadic promotion.

### D53 — build fix for FloatMP infinity sign extraction (21 October 2026)

D51 introduced two invalid calls to `sgn(*this)` inside FloatMP's
conversion operators; unlike Dyadic, FloatMP does not expose that
free function in this context. Restore the existing comparison-based
sign selection `(*this > 0) ? Sign::POSITIVE : Sign::NEGATIVE` in both
FloatMP->Dyadic and FloatMP->Rational special-value branches.
No arithmetic or rounding path changes. Build and CI validation
remain mandatory before treating D51/D52 as tested.

### D54 — native Dyadic promoted to canonical filenames (21 October 2026)

Replace canonical dyadic.hpp/dyadic.cpp with Integer-magnitude and denominator-exponent native backend, fold wrappers directly into dyadic.hpp, preserve writers, remove old mpf_t implementation. Remove new_dyadic.hpp and test_new_dyadic.cpp and update FloatAP and benchmark to canonical type and header. Preserve legacy test_dyadic.cpp plus native wrapper regression checks. FloatMP remains scheduled for removal; build and tests must validate compatibility. This is a large integration candidate: do not equate commit publication with green CI.

### D55 — restore Dyadic public integration after canonical migration (21 October 2026)

Address compile errors from D54: checked Natural exponent extraction uses the inherited Integer::get<std::size_t>() API, not nonexistent get_ui(); remove three obsolete DyadicBounds definitions from rational.cpp now implemented in dyadic.hpp; include complete float enclosure types where out-of-line getters are defined; retain the already used FloatMP getter declarations until its separate deletion without reintroducing mpf_t into Dyadic; convert legacy FloatMP lower/upper inputs via their exact Rational value. Restore missing rational-valued reciprocal/division contracts for Dyadic and explicit mixed RationalBounds/DyadicBounds operations and comparisons required by NumberWrapper. The exact interval arithmetic is unchanged. No test fixtures suppressed. These are build fixes and still require native compilation/CI validation.
