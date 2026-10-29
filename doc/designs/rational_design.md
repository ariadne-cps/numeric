# Rational — native exact rational design journal

**Repository:** `ariadne-cps/numeric` · **branch:** `main`
**Scope:** `Rational` finite and extended arithmetic, numerator/denominator canonicalization, GCD and cancellation use, comparisons, conversions, exact oracles, and R-series performance and correctness history.
**Companion:** [Integer](integer_design.md) · [BigUnsigned](big_unsigned_design.md) · [Dyadic](dyadic_design.md) · [FloatMP](floatmp_design.md)

Each historical milestone is preserved verbatim below with its original date and identifier; older entries document the state at the time, not necessarily the final implementation. Full validation and CI must be verified independently.

## Current architecture and correctness contract

A finite `Rational` has an `Integer` numerator and a positive `BigUnsigned` denominator; zero is canonical and nonzero finite fractions are reduced. Positive/negative infinity and NaN are explicit states. Preserve signed arithmetic, cancellation before large products, exact comparisons, floor/ceil/round, IEEE-754 conversion, and independently certified test fixtures. The portable and native unsigned primitives, division and GCD substrate are documented in [`big_unsigned_design.md`](big_unsigned_design.md).

## Implementation and optimization timeline

### R1 — native finite and extended Rational (20 October 2026)


Numerical verification passed 768 paired inputs and 2688 results.

**Representation and invariants.** Rational now privately owns an Integer numerator, a BigUnsigned denominator and an explicit finite/positive infinity/negative infinity/NaN kind. Every finite value is reduced, its denominator is positive and zero is 0/1. Ordinary rational constructors, arithmetic, numerator/denominator accessors, streams and number/Python integration retain their mathematical interfaces. Builtin unsigned construction now uses Integer directly, including UINT64_MAX.

**Finite arithmetic.** Addition/subtraction first computes g=gcd(b,d) for a/b and c/d, forms n=a*(d/g)+/-c*(b/g), and reduces only by h=gcd(abs(n),g). The denominator is (b/g)*(d/h); no other common factor can remain because both inputs were reduced. Multiplication cancels gcd(abs(a),d) and gcd(abs(c),b) before either product; division cancels numerator and denominator pairs before multiplying and handles the resulting sign explicitly. Same-object multiplication and sqr use squares of the already-coprime components. Halving cancels a factor of two from the numerator or shifts the denominator. Reciprocal exchanges coprime magnitudes. A private Reduced constructor is used only where these operations prove coprimality; it still canonicalizes zero. Powers use shared unsigned powering without repeated fraction gcds; negative exponents, including Int minimum, do not negate an unrepresentable signed minimum. Exact fma is implemented as exact multiply then add. Comparisons use canonical equality or exact cross-products. Floor, ceil and half-away-from-zero round share one unsigned quotient/remainder.

**Extended semantics and explicit decisions.** Existing ExtendedOperations continues to dispatch finite versus special arithmetic, comparisons, min/max, sign, abs and reciprocal. NaN remains unequal to itself and named cmp returns INCOMPARABLE. No Rational signed zero is introduced. Finite division by zero is explicitly defined as signed infinity for nonzero numerators and NaN for 0/0; reciprocal zero is positive infinity. Decimal text accepts signed integers or n/d, ignoring ASCII whitespace, and raises runtime_error on malformed input; textual denominator zero has the same explicit extended meaning as the integer-pair constructor. Equal infinities in Rational/ExactDouble comparisons now compare equal rather than taking the old unconditional less/greater branch; NaN remains incomparable.

**Exact binary64 conversion.** ExactDouble import decomposes a finite binary64 into an integer significand and power of two and cancels trailing zero bits. get_d computes the binary exponent with exact comparisons, then truncates an exact scaled integer quotient to at most 53 bits. Subnormal spacing is 2^-1074; overflow yields signed infinity. Results are independent of the active host rounding mode.

**Verification and measurements.** After the R1a test-only corrections documented below, the complete AppleClang 21 Release warnings-as-errors everything build and all 20 C++ CTest suites pass on Apple M2 Ultra, macOS 26.6.2. Integer, BigUnsigned, Rational, Float and FloatRounding Debug ASan/UBSan suites pass (5/5). GCC 16.2 ARM64 passes the native Rational/Integer/BigUnsigned suites with normal project warnings-as-errors. A complete default GCC build encountered -Wdisabled-optimization in the existing very large TestFloat/TestFloatRounding template functions; full-build validation passes everything and all 20 CTest suites with the larger local GCC --param=max-gcse-memory=1000000 budget, without changing project CMake rules or suppressing warnings. Native x86-64 and MSVC are not available on this host and are not claimed as tested.

The unchanged six-size Rational benchmark (20000 iterations, 31 samples) passes 384 paired inputs and 2880 operation results. These measurements overlapped compiler/test activity and are diagnostic only, not performance acceptance or an x86-64 prediction. Rational normalization/GCD, comparison and division remain optimization work.

These single runs do not establish performance acceptance or blanket nonregression; the existing short-width MP gates are still unmet. No FloatMP arithmetic, precision or rounding algorithm was changed.

### R1b — retain division alias coverage on Clang 16 (20 October 2026)


Clang 16 rejects the intentional destination/=destination test under -Wself-assign-overloaded. Use the existing reference alias as the right operand, preserving exactly the same object-aliasing case and expected value without suppressing any warning. AppleClang 21 rebuild and the complete Rational suite pass. The previous R1a CI already passed macOS Release/Coverage, Ubuntu GCC 13 x86-64 Release/Coverage and Windows 2025 Release; Clang 16 Release/ThreadSanitizer requires this test spelling fix. The optional local Clang 15 syntax experiment could not compile against the installed newer macOS SDK/libc++ combination and is not validation. Rational arithmetic and Integer/BigUnsigned remain unchanged.

The separate Codecov status thresholds already failed on I75 (project 92.96%, patch 92.56%, target 99.20%). R1a reports project 92.94% and patch 93.07%; workflow Coverage success is not Codecov gate acceptance.

### R2 — remainder-only Euclidean GCD (20 October 2026)


The shared `BigUnsigned::gcd` Euclidean loop now calls the existing remainder-only `rem(a,b)`, rather than constructing a quotient/remainder pair through `divrem(a,b)` and discarding the quotient at every iteration. It preserves the same Euclidean recurrence and edge cases, including gcd(0,0)=0. The single shared-core change may also affect other users of `gcd`; benchmark and functional gates for those users remain applicable.

**Validation status:** no local build or benchmark was available in the editing environment. No speedup, CI outcome or performance gate acceptance is claimed for R2 until the updated branch is rebuilt and tested on ARM64 and x86-64.

### R3 — bounded stack products for Rational comparison (20 October 2026)


Compared with the pre-R2 aggregate totals, ARM64 decreases from 1.265x to 1.073x and x86-64 from 1.329x to 1.021x. The R2 Unix, Windows and Coverage GitHub Actions complete successfully. Independent native ARM64 and x86-64 R3 measurements remain pending.

**Change.** The finite `Rational::cmp` cross-product comparison now uses caller-owned stack limb arrays through the already-selected, portable-fallback `BigUnsigned::mul_into` kernel when each of the four unsigned input magnitudes occupies at most eight 64-bit limbs. Both exact products are compared from their high words, with a zero-initialized tail to handle different product lengths and any zero top limb. Above the eight-limb bound, the previous owning `mul` plus `compare` path is unchanged. Both paths preserve the existing sign and same-denominator shortcuts, the exact comparison order and IEEE/extended-special-value dispatch. This specifically avoids two temporary heap-owned cross-products for the medium-size input range, while retaining the existing multiplication kernels.

**Validation status.** R3 has not been locally built or benchmarked: no native ARM64 or x86-64 execution environment was available in the repository-editing session. Do not interpret the R2 measurements as R3 speedup or accepted performance gates. Wait for CI and re-run the unaltered six-size benchmark on both hosts, plus the full Rational, Integer, BigUnsigned and FloatMP gates before evaluating another candidate.

### R4 — proper-fraction floor/ceil shortcuts and strict Rational performance gates (20 October 2026)


Keep the 15 operations, 32-pair workloads, 20,000 iterations, 31 samples, numerical oracle and measurement methodology unchanged. An architecture passes only if all six sizes pass both gates; do not average away a failure.

All 384 paired inputs and 2880 paired operation results pass numerical verification on each architecture. Unix, Windows and Coverage workflow runs on commit `012b336` succeeded. Only ARM64 128-bit passes both user gates at R3. This is not an R4 measurement.

**R4 change:** finite `floor` and `ceil` now check the exact unsigned magnitude order before allocating quotient/remainder storage. For `|n|<d`, floor gives -1 for negative nonzero and 0 otherwise; ceil gives 1 for positive nonzero and 0 otherwise. Canonical zero is 0/1. The existing Euclidean quotient/remainder path is entirely unchanged when `|n|>=d`. This optimizes mathematically general proper fractions, not a benchmark-specific shape; it preserves signed results and the precondition that special values are not valid inputs.

**Validation status:** R4 has not been executed on ARM64 or x86-64 in this editing session. Its build/CI and six-size performance outcome are pending; no post-R4 gate is claimed as passing.

### R5 — exact bounded conversion to binary64 (20 October 2026)

 The independent user-provided R4 measurements (ARM64 then x86-64) pass 384 paired inputs and 2,880 operation results each. R4 Unix, Windows and Coverage CI pass.

**R5 finite Rational::get_d paths:** for one-limb numerator and denominator, compute the exact scaled binary64 significand in an unsigned 128-bit intermediate on GCC/Clang with __int128 support (including x86-64 and ARM64), avoiding owning big-integer division. On other platforms, preserve the established portable fallback.

For same-limb-count operands with at least five 64-bit limbs, normal binary64 values and negative scale of at most 63, first divide only the top three limbs, shifted by the 53-bit/significand scale. With n=N*B+nl and d=D*B+dl, B a shared power of two, and N*2^s=q*D+r, 0<=nl,dl<B, the bounds r>=q and r+2^s<D imply floor(n*2^s/d)=q. When q=2^52 on normal numbers, the exact binary exponent already guarantees the lower bound. If neither proof applies, or the inputs lie outside the selected range, take the unchanged full-width exact division.

Keep all earlier randomized tests and benchmark workload/oracle unchanged.

**Validation:** commit changes without claiming performance or CI outcomes for R5. Rebuild and run all Unix/Windows/Coverage workflows, Rational/Integer/BigUnsigned/FloatMP differential gates and unchanged ARM64/x86-64 six-size benchmarks before considering R5 accepted.

### R6 — exact equal-leading-limb Euclidean subtraction (20 October 2026)


R5 Unix/Windows/Coverage CI succeeded.

**R6:** In the shared `BigUnsigned gcd` Euclidean reduction, when both operands occupy the same number (>1) of canonical limbs and their most significant limbs are equal, order them and replace the general remainder operation with the exact subtraction. For ordered a>=b with the same positive high limb, a<2*b, so a mod b=a-b. Other inputs keep the existing remainder-only algorithm unchanged. The input ordering is semantically irrelevant to gcd; zeros preserve the former behavior. This targets division/scratch overhead in a mathematically general case, including differences collapsing from heap to inline limbs at 256 bits, rather than checking any benchmark operand shape.

**Validation:** R6 not locally compiled/benchmarked. R6 CI, cross-platform differential tests, and fresh ARM64/x86-64 baseline-matched performance runs are pending; no target is claimed passed on R6.

### R7 — certified native 128-bit high-prefix binary64 division (20 October 2026)

 **R6 baseline (not R7 results):** the independently supplied ARM64 and x86-64 runs both pass 384 paired inputs and 2,880 paired operation results; R6 Unix, Windows and Coverage workflows succeed.

**R7 change:** for normal finite binary64 conversion, with a negative scale of at most 63 bits and equal limb counts, first compute the high-prefix quotient using native unsigned 128-bit arithmetic. Use one high limb, and on an ambiguous quotient use two high limbs. For n=N*B+nl, d=D*B+dl with 0<=nl,dl<B, if N*2^shift=q*D+r and r>=q, then the full quotient is >=q. If r+2^shift<D, its full quotient is <q+1. When the candidate equals 2^52, the exact known normal exponent independently certifies that the full quotient is >=2^52. Apply both checks before accepting the computed 53-bit significand. Check shift overflow, candidate bit range, and prefix validity; otherwise use the existing R5 three-limb certificate and ultimately the original full exact division. Preserve the established active-rounding-mode independence, truncation toward zero, signed zero, overflow/subnormal semantics and compiler fallbacks. This is independent of particular rational values and does not weaken correctness guarantees.

**Validation status:** no local build or R7 native benchmark has been run. The new commit requires Unix/Windows/Coverage, FloatMP/Integer/ BigUnsigned/Rational regressions, and unaltered six-size ARM64/x86-64 runs.

### R8 — platform-select exact conversion and scalar Euclid termination (20 October 2026)

The user supplied R7 runs (20,000 iterations, 31 samples; 384 paired inputs, 2,880 operation results PASS per machine). R7 Unix, Windows and Coverage all pass. These are R7 data, not R8.

**R8 conversion selection:** R7's new unsigned 128/128 top-limb certificate costs approximately 117 ns versus 36 ns for R6 at 128 bits on ARM64, and also regresses 256+ ARM64 conversion. On x86-64, R7 materially reduces conversion cost at 128/256/1024/2048 bits. Therefore enable only the R7 high-prefix unsigned 128-bit conversion on x86-64 GCC/Clang; ARM64, MSVC and portable targets retain the R5 three-limb/exact fallback (and the existing one-limb R5 path). The numerical semantics and shared four-rounding-mode differential gate remain identical. This selects algorithms by observed host performance, not precision or benchmark input values.

**R8 fast one-limb comparison:** Rational::cmp now compares exact unsigned 128-bit cross-products directly when all four numerator/ denominator magnitudes contain exactly one 64-bit limb, avoiding stack-product arrays, dynamic multiplication dispatch and their bookkeeping. Longer operands and toolchains without unsigned 128-bit support retain R3's exact paths.

**R8 scalar GCD termination:** BigUnsigned::gcd orders initial inputs once, then maintains a>=b invariant throughout Euclid. When the current divisor fits one 64-bit limb, perform at most one multiprecision remainder, then finish Euclid with scalar '%' iterations. Special-case divisor 1, preserving gcd(0,0), symmetry and all other cases. Keep R6's exact equal-high-limb subtraction on multiprecision inputs.

**Status:** R8 native CI and six-size ARM64/x86-64 benchmarks pending; no R8 gate has been claimed passed. Shared-core changes require the Integer/BigUnsigned/Rational and downstream FloatMP differential gates.

### R9 — one-source BigUnsigned compilation unit (20 October 2026)


**Motivation / user decision:** native architecture conditionals are not fully separable into platform-specific source files. Division paths, wide-integer scalar estimation, decimal formatting, Karatsuba and multiplication dispatch already contain compiler/ISA branches inside common implementation functions. Keeping external kernels in separate `.cpp` files therefore creates a misleading partial modular boundary.

**Structural change only:** merge the reference C++20 multiplication formerly in `big_unsigned_fallback.cpp`, the optional x86-64 kernel formerly in `big_unsigned_x86-64.cpp`, and the ARM64 native carry and symmetric-square kernels formerly in `big_unsigned_arm64.cpp` into `src/big_unsigned.cpp`. Enclose the platform-specific implementations in the SAME `ARIADNE_BIG_UNSIGNED_X86_64` and `ARIADNE_BIG_UNSIGNED_ARM64` compile-time guards used before. Preserve all internal public/diagnostic entry points, implementation bodies, kernel selectors, GCC fused option, target compiler definitions and portable reference implementation. Remove the three now-redundant sources and their CMake source entries; CMake still selects the same native backend flags based on compiler and target. No algorithm, operand set, oracle, threshold, runtime condition, test or benchmark is altered. Preserve historical source attribution in section comments and Git history.

**Validation gates:** starting from R8b `aa5c520`, whose Unix, Windows and Coverage CI all succeeded, confirm the refactored single-source unit builds on GCC x86-64, Clang x86-64, AppleClang ARM64 and MSVC; run the full differential BigUnsigned, Integer, Rational, FloatMP and architecture-specific multiply diagnostics. These validations are pending for R9; no improvement is claimed.

### R10 — borrowed initial Euclidean reduction for Rational normalization (20 October 2026)


**Baseline.** R9's single-source BigUnsigned refactor passes Unix, Windows and Coverage CI. The last user-provided ARM64 and x86-64 R8a/R8b benchmark runs each pass numerical verification, but `construct_fraction` remains relatively slow at medium/large sizes, including 256-4096 bits.

**Change.** All Rational arithmetic GCD sites pass existing lvalue magnitudes. Switch `BigUnsigned gcd` from by-value operands to `const&`, and compute the first exact Euclidean remainder directly from those borrowed operands, rather than copying both magnitudes. Return the surviving borrowed input directly if the first remainder is zero. For one-limb remainder, finish with the established R8 scalar Euclid without copying large inputs. Only if another multi-limb reduction is necessary, copy the surviving divisor once and thereafter reuse move-owned remainder buffers. Retain R6's exact same-high-limb subtraction shortcut and the R2 remainder-only division path. This function signature changes its C++ symbol/ABI; dependent binaries must be rebuilt from sources.

Keep randomized multi-limb GCD, Rational construction and downstream FloatMP regressions unchanged.

**Validation.** No local R10 compilation or performance measurement. Require Unix, Windows, Coverage CI, full Integer/BigUnsigned/Rational/ FloatMP test suites, and unchanged 20,000-iteration 31-sample ARM64/x86-64 six-width benchmark runs before accepting the change.

### R11 — exact quotient-one/two Euclidean reduction (20 October 2026)


Five of twelve width-host pairs meet both strict gates. At the time of R11 creation, R10 Coverage CI had passed while Unix/Windows CI were still running. No R11 speedup is claimed.

**Change.** Inside BigUnsigned::gcd, R6's equal-high-limb subtraction shortcut generalizes to ordered positive operands a>=b with equal canonical limb counts and highest limbs high_a-high_b <= high_b (computed without overflow). With B=2^(64*(limb_count-1)), high_b>=1 and a < (high_a+1)*B <= (2*high_b+1)*B <=3*high_b*B <=3*b. Thus the Euclidean quotient is either 1 or 2. Compute a-b once. If high_a<2*high_b, a<2*b, and a-b is the remainder. Only on high_a==2*high_b is another comparison needed: if a-b>=b, subtract b once more. This computes the exact remainder without Knuth division, scratch buffers, or approximation. Cases outside the exact range still call the unchanged remainder-only implementation. This applies equally to dense arbitrary operands and sparse rational inputs, not a workload signature. Keep R10's borrowed initial operands, R8's scalar Euclid and all preexisting arithmetic semantics.

**Validation.** R11 has not been locally compiled or benchmarked.

### R12 — exact scalar x86-64 Rational-to-double conversion (20 October 2026)

 **Measured R11 baseline (not R12).** Both ARM64 and x86-64 20,000-iteration/31-sample user-supplied runs pass 384 paired inputs and 2,880 paired operation results. On x86-64 the worst at 53 bits is to_double (2.307x); on ARM64 the worst at 256 bits is to_double (2.544x).

**Change.** On GCC/Clang x86-64 only, when both canonical Rational magnitudes have exactly one 64-bit limb, compute the exact binary exponent by integer bit-width comparison and bounded native shifts, without creating BigUnsigned shifts/temporary quotients. When exponent<=52, form the exact scaled numerator as unsigned 128-bit and divide its two 64-bit halves with baseline x86 DIVQ. The independently computed exponent proves the quotient has at most 53 bits, meeting DIVQ's high-word<divisor precondition; an explicit check retains the portable full division fallback. When exponent>52, a native 64-bit division by the exact shifted denominator suffices. The returned integer significand is exactly representable in double and ldexp preserves truncation toward zero independent of the host rounding mode. ARM64, Windows/MSVC and multi-limb inputs keep their previous algorithms unchanged.

**Validation.** R12 native full CI and matched ARM64/x86-64 benchmarks have not run yet. Do not claim a measured improvement or success on any remaining performance gate before the new data.

### R13 — certified small-multiprecision ARM64 binary64 conversion (20 October 2026)


**Motivation:** the user explicitly highlighted the unresolved ARM64 256-bit `to_double` regression.

**Change:** for normal binary64 quotients with negative or zero scale and equal canonical 5..8-limb magnitudes on ARM64 GNU/Clang, derive a tentative significand from the top limbs using ordinary floating arithmetic. Clamp into the valid [2^52, 2^53-1] range. Independently build the exact shifted numerator and exact divisor-times-candidate products in fixed-length stack arrays. Only return the candidate if B*q <= A*2^shift < B*(q+1) is proven by full multiword comparisons. The test is independent of floating mode, architecture instruction rounding, approximation quality, and lower limbs. An incorrect or ambiguous candidate always falls through to the unchanged full-width conversion; zero, signs, subnormals, overflows, unequal widths and other targets remain unchanged. This is an exact certification strategy, not an assumption about benchmark numerator/denominator patterns.

**Status:** R12a fixes a reproduced most-vexing-parse build regression and is pending full CI; R13 remains unbuilt and unbenchmarked in this session.

### R13a — avoid rounding-namespace shadow and record R13 performance (20 October 2026)


**Build warning:** AppleClang with `-Wshadow` reports that the local `rounded` in the R13 ARM64 `Rational::get_d` certificate masks the existing `Ariadne::rounded` in `numeric/rounding.hpp:261`. Rename only this local value to `rounded_significand`; there is no change to arithmetic, rounding, exact certificates or benchmarks.

All 384 paired input and 2,880 result numerical benchmark checks pass on each machine. Thus 7 of 12 host-width cases meet both gates after R13, and the strict acceptance target is still NOT achieved. R13a introduces no performance algorithm changes.

**Validation:** R13's Unix/Windows/Coverage workflows were still running when R13a was prepared. R13a CI and matched host benchmarks must be checked independently. No fresh speedup or CI success is claimed for R13a.

### R14 — relocate self-alias square dispatch from Integer to BigUnsigned (20 October 2026)


**Purpose:** restore the signed/unsigned responsibility boundary, without changing benchmarks or selected multiplication kernels. I33/B37 detected same-object `Integer` multiplication and routed 8+-limb operands to symmetric `BigUnsigned::square` on GCC x86-64 fused and ARM64 native backends. This is a magnitude optimization, not signed arithmetic, so place the exact same platform/size/identity check in `BigUnsigned mul(const&,const&)`. `Integer mul` now applies only the XOR of operand signs to the delegated unsigned product.

**Recursion safety:** the `square` fallback calls `mul(a,a)`, but the fast-path is reachable only at >=8 limbs on backends where `square` uses a dedicated symmetric kernel at the same threshold. For <8 limbs or all unsupported targets, `mul` does not dispatch to `square`. Zero and oversized product checks remain in their existing order. The new unsigned dispatch may also benefit other genuine self-alias callers, but no runtime value-pattern or benchmark shortcut is added.

Rebuild and run full Integer, BigUnsigned, Rational and FloatMP gates. This changes shared unsigned dispatch; identical performance has not yet been measured.

**Validation status:** before this commit, R13a Unix, Windows and Coverage CI all passed. R14 compilation, host-specific benchmark stability and full regression are pending.

### R15 — centralize platform-specific unsigned exponentiation in BigUnsigned (20 October 2026)


**Motivation:** the sign-independent B32 right-to-left exponentiation schedule remained in `Integer::pow` on ARM64, Clang, MSVC and portable targets, whereas GCC x86-64 fused builds already called B74's MSB-first buffered `BigUnsigned::pow`. That made the signed front end responsible for unsigned algorithms and architecture dispatch. R14 moved the related self-multiply selection into `BigUnsigned::mul` in a preceding separate commit.

**Change:** `Integer::pow` now solely converts the exponent, derives the result sign from negative-base parity and delegates to `BigUnsigned::pow`. The unsigned function retains B74's original GCC x86-64 fused MSB-first implementation verbatim under the existing compiler/backend definition. Other targets use B32's exact right-to-left loop, moved from `Integer::pow`, preserving first-factor copy elision, square/multiply ordering, zero and exponent-one handling. The public `BigUnsigned::pow` API continues to be exact for all unsigned inputs but now selects the same policy as the signed caller on each target. This modifies the algorithmic schedule for *direct* non-GCC callers of `BigUnsigned::pow` (previously MSB-first buffered), so such callers need independent performance regression measurement. The result of `Integer::pow` is unchanged for 0^0, +/-0, +/-1, negative odd/even powers, aliasing, long values and every Nat exponent.

Test and benchmark performance parity separately on ARM64, GCC/Clang x86-64, MSVC and forced-portable builds.

**Validation:** R13a CI passed before R14. R14 CI was still running when R15 was authored. Rebuild and rerun full Integer, BigUnsigned, Rational, FloatMP and performance gates before accepting.

### R16 — borrow coprime denominators and move normalized sum (20 October 2026)

 **Baseline:** R15a (`b7f79d1`) completed Unix, Windows and Coverage CI successfully. The last supplied side-by-side Rational benchmark measurements are R13 (before the R14/R15 refactoring); do not present those as R16 or R15a performance. R13 was compliant on 7 of 12 host/precision combinations.

**Change:** in `FiniteOperations<Rational>::add_sub`, after `g=gcd(den_a,den_b)`, if g=1, work directly with both borrowed canonical denominators instead of copying them through `exact_quotient(den,1)`. For canonical p/q and r/s with gcd(q,s)=1, the signed cross-sum numerator `p*s +/- r*q` is coprime to q: modulo q it is p*s (coprime to q), and similarly modulo s it is +/-r*q (coprime to s). Hence the product q*s already is fully reduced, provided the numerator is nonzero. Skip the redundant `gcd(n,1)`, the quotient-by-1 copies and numerator copy; construct from an owning numerator moved into the result. Canonical cancellation to zero still follows the existing zero path. When g!=1, keep the established reduced-denominator algorithm, but if the final h=gcd(n,g)=1, again move the already-reduced numerator and borrow the original denominator for its final product. h>1 continues through the prior exact quotient path. All cases preserve the prior rational sign, denominator, alias and exceptional semantics.

**Validation:** R16 has not been locally compiled or benchmarked. Require Unix, Windows and Coverage CI, the complete downstream Integer/Rational/FloatMP test suite and fresh unchanged benchmark runs on both hosts before claiming a speedup. The fast-path is an algebraic identity for all canonical rationals, not benchmark-only.

### R17 — certified two-leading-limb Rational cross-product comparison (21 October 2026)


**Baseline:** R16 Unix and Coverage workflows passed, while Windows failed during an external vcpkg `z-applocal` step with `Access is denied` while copying benchmark DLLs after linking. The user classified this as a temporary GitHub environment issue and explicitly requested that it be ignored: no CMake or build-tool changes. R16's separately supplied ARM64 and x86-64 benchmark reports both pass 384 exact input checks and 2,880 operation-result checks. Seven of twelve host-width pairs pass both required strict gates. The principal x86-64 remaining worst-case is compare at 128 bits.

**R17 change:** when all four unsigned Rational cross-product magnitudes contain at least two limbs and the two product limb sums are identical, compute the products only of each operand's highest two 64-bit limbs (4-limb local output). Writing operands as x=U*B^k+r and y=V*B^m+s, with B=2^64 and 0<=r<B^k, 0<=s<B^m, the exact product obeys UV*B^(k+m) <= xy < (U+1)(V+1)*B^(k+m). The exact prefix product UV occupies up to four limbs. Form the exact upper prefix bound UV+U+V+1 by four-limb additions with full carry propagation; if the bound overflows four limbs, decline certification. If the upper left bound is <= the lower right bound, the left product is strictly smaller; symmetrically for right. Otherwise fall back to the existing complete 8-limb stack/fallback multiplication and comparison. Size mismatch and one-limb cases remain on their prior exact paths. This is not a benchmark-specific heuristic: it proves the result before returning.

**Validation:** R17 is not locally compiled or benchmarked. Require Unix/Windows/Coverage checks plus unchanged ARM64/x86-64 31-sample, 20,000-iteration runs before interpreting performance.

### R18 — use CMake architecture selectors consistently (21 October 2026)


Replace 40 architecture preprocessor predicates in `src/big_unsigned.cpp` and three in `src/rational.cpp` with the existing CMake-selected `ARIADNE_BIG_UNSIGNED_X86_64` / `ARIADNE_BIG_UNSIGNED_ARM64` target-private flags. This is a lexical-only replacement of `defined(__x86_64__)`, `defined(__aarch64__)`, and `defined(_M_X64)`; all compiler feature guards, algorithm bodies, thresholds, benchmarks and test oracles are unchanged. Document PRIVATE definition scope in `src/CMakeLists.txt`. Native regression and benchmarks pending; do not claim performance changes. R17 Unix and Coverage CI passed prior to this refactor.

### R19 — borrow uncancelled Rational multiply/divide factors (21 October 2026)

 **Baseline.** R18 (`e1a8fb2`) Unix, Windows and Coverage workflows all passed. The most recent separately supplied 6-width host benchmarks are from R16; neither R17 nor R18 has paired native performance measurements yet. Both reports passed 384 paired inputs and 2,880 paired results.

**Change.** Canonical Rational multiplication and division use two cross-GCDs (g and h), then produce two products from four exact cancelled factors. Prior code evaluated `exact_quotient(x,1)` by copying `x`; each cross-GCD equal to 1 could therefore cause two unnecessary BigUnsigned copies. New private `multiply_cancelled(lhs,lhs_divisor,rhs,rhs_divisor)` borrows each original magnitude when its exact divisor is 1, independently for both factors, and computes a real exact quotient only when the divisor exceeds 1. The result is constructed as an owning numerator and denominator, with numerator sign determined as before by XOR. Keep the same cross-GCDs, exact multiplication, canonicality, special-values and self-alias-square dispatch. Remove now-unused Integer wrappers for `mul` and `exact_quotient`. These are optimizations for all coprime and partially coprime factors, not special cases of benchmark inputs.

**Validation.** No local R19 full build or matched ARM64/x86-64 measurement has yet been performed. Require Unix/Windows/Coverage and unchanged six-width 20,000-iteration/31-sample paired benchmarks before interpreting performance.

### R19a — restore signed exact quotient used by Rational addition (21 October 2026)


**Build-blocking regression caused by R19:** the cleanup removed the `exact_quotient(Integer const&,BigUnsigned const&)` wrapper as supposedly unused. However, the `g>1,h>1` normalization branch in `FiniteOperations<Rational>::add_sub` still calls `exact_quotient(n,h)` with an Integer numerator. The remaining BigUnsigned overload cannot accept Integer, causing AppleClang's two reported compile errors in rational.cpp (around line 253). R19 Unix, Windows and Coverage workflows all failed in Build.

**Minimal repair:** restore the exact signed wrapper from R18. It delegates to the unsigned exact quotient of the magnitude and restores the numerator's original sign with `Integer::from_magnitude`. This is the existing, verified behavior for nontrivial shared factors in add/sub. It is intentionally still needed by this branch; the statement in R19 about removing that wrapper was incorrect. The now-unused `integer_mul` wrapper stays removed.

No successful R19a build or speedup is claimed in advance.

### R20 — bounded stack-based signed cross-product sums (21 October 2026)

 **Baseline:** R19c (`3ade341`) passed Unix, Windows and Coverage. Both freshly supplied R19c benchmarks (20,000 iterations, 31 samples) pass the 384 paired inputs / 2,880 operation results on ARM64 and x86-64.

**Change:** in `Rational::add_sub`, when each signed numerator magnitude and unsigned scale factor occupies at most eight 64-bit limbs, replace two owning multiplications and a third owning integer add/sub with exact multiplication into two 16-limb local arrays, followed by `BigUnsigned::add_into` or `sub_into` into a 17-limb local result. Use full-width magnitude comparison to choose the subtraction direction, effective sign and cancellation to canonical zero. Materialize only one final owning magnitude with `Integer::from_magnitude`. The same helper serves both coprime and shared-denominator branches. Four entirely scalar factors, any unit scale factor and long operands preserve the previous path to avoid undoing existing specializations. This is a general allocation-reduction technique, not an optimization keyed to benchmark values.

Exercise both effective operand signs, dense and closely matched fractions, zero cancellation, carries, subtract-borrow, compound operators and the 8/9-limb boundary.

**Validation pending:** R20 has not been compiled or measured locally. Require Unix, Windows and Coverage CI, Integer, BigUnsigned, Rational and FloatMP checks, and native six-width ARM64 and x86-64 benchmarks before claiming speedup. Keep the strict user gates per host and per size.

### R21 — bounded borrowed-view Euclidean GCD (21 October 2026)


**Implementation.** For nonzero inputs of at most 128 limbs, keep Euclid's exact remainder sequence in three independent stack buffers. Borrow both original inputs until their first reduction; rotate the result buffer so neither live input ever aliases the output. Normalize every remainder before the next step. Use the existing remainder-only division kernel and its portable/ISA selection; no new instructions, CPU requirements or floating estimates. Preserve R11's exact quotient-one/two subtraction and trim identical leading limbs before subtracting. Finish single-limb divisors with the existing selected remainder kernel (including B61) and scalar Euclid. Only the final GCD becomes an owning BigUnsigned. Larger operands retain the existing multiple-precision path. The fixed bound limits scratch, not supported integer precision. No benchmark inputs or canonical Rational contracts change. Record native timings and sanitizer/compiler checks after execution.

The 128-limb bound uses 3*128 remainder limbs plus 257 division scratch limbs (5,128 bytes of arrays). Those arrays need no initial zero fill: each kernel fully writes its active output and every consumer sees only its normalized span. Euclid maintains a>=b>0; for equally long inputs ah-bh<=bh proves a<3*b without overflowing 2*bh. At most two exact subtractions suffice. Trimming equal high limbs is safe because the first unequal lower limb proves the remaining unsigned order. The second subtraction uses independent division scratch, then copies the normalized result back. Three-way rotation never overwrites either live input. No architecture-dependent arithmetic was added.

**Regression coverage.** New independent cases use gcd((2^k-1)*f,(2^k+1)*f)=f, f=2^64+3, with k>=126. The coprimality follows from the odd operands differing by two. Existing randomized GCD, quotient-one/two, equality, zero, scalar, Rational canonicality and special-value suites are preserved.

**Online measurement.** A dedicated Rational benchmark push/manual workflow compares HEAD^ with HEAD on the same Ubuntu 24.04 x86-64 runner, separately with GCC 13 and Clang 16. It builds both Release revisions with warnings-as-errors, tests the candidate, and runs the unchanged six-width 20,000x31 benchmark sequentially. Raw logs, CPU and compiler details and both SHAs are uploaded as artifacts; the summary is informational because shared hosted CPUs have timing noise. No speed gate is weakened and no benchmark fixture is changed.

**Local validation.** AppleClang 21 Release warnings-as-errors: all 20 C++ CTest suites pass. Debug ASan/UBSan: Integer, BigUnsigned, Rational, Float and FloatRounding pass (5/5). GCC 16.2 ARM64 Release warnings-as-errors: Integer, BigUnsigned and Rational pass (3/3). The workflow YAML parses and git diff --check is clean. Native MSVC, x86-64 and remote Unix/Windows/Coverage validation remain pending at publication; do not infer their results from the local checks.

**Native ARM64 measurements.** Apple M2 Ultra (Mac Studio), macOS 26.6.2, AppleClang 21.0.0.21000334, Release, unchanged benchmark, six widths, 20,000 iterations and 31 samples. R20 parent and R21 were rebuilt in separate directories; the recorded comparison ran sequentially without concurrent builds/tests. Both pass 384 paired inputs and 2,880 exact operation results.

| bits | R20 total | R21 total | R20 worst | R21 worst |
|---:|---:|---:|---:|---:|
| 53 | .356 | .367 | 1.245 | 1.235 |
| 128 | .649 | .565 | 1.803 | 1.714 |
| 256 | 1.120 | .894 | 1.697 | 1.750 |
| 1024 | 1.111 | .898 | 1.402 | 1.413 |
| 2048 | .825 | .720 | 1.499 | 1.345 |
| 4096 | .682 | .627 | 1.630 | 1.516 |
| ALL | .772 | .680 | 1.803 | 1.750 |

Worst cases are to_double at 53/256/1024, compare at 128 and construct_fraction at 2048/4096. These kernels were not altered; small absolute timing changes include measurement/code-layout noise. The sum of native medians at 256 bits falls from 1,987.330 to 1,587.294 ns, and at 1024 from 4,087.420 to 3,302.547 ns. Examples of absolute R20 -> R21 medians (ns): 256-bit construction 194.525 -> 138.338, add 315.548 -> 237.144, mul 372.442 -> 283.812; 1024-bit mul 863.206 -> 629.456, div 857.354 -> 636.850. The 53-bit native sum increases 1.9% (411.401 -> 419.288 ns), where scalar inputs retain the old GCD path; no universal speedup is claimed. Recheck strict gates on native x86-64.

The paired Integer 20,000x31 comparison passes 768 inputs/2,688 results: aggregate total .958 -> .951 and worst 1.659 -> 1.646. FloatMP's unchanged 1,000x11 comparison passes 256 inputs/3,456 results. Its prior short-width performance deficits remain outside R21's scope.

### R22 — exact scalar division and borrowed cancelled factors (21 October 2026)


**Implemented contract.** Add a shared exact single-limb quotient kernel for the proven-divisibility cases in Rational normalization and cross-cancellation. Write a disjoint caller buffer of at least the source length; preserve the source and zero extra output. Split d=odd*2^s, require the source's low s bits to vanish, and read the right-shifted source without an owning shift. For odd d, compute its inverse modulo B=2^64, then recover quotient digits from low to high: q_i=(a_i-carry)*inverse mod B; carry'=high(q_i*d)+borrow(a_i,carry). The invariant q_partial*d = a_partial + carry*B^(i+1) proves the digits; carry<d guarantees the high product plus borrow cannot overflow. A nonzero final carry rejects nonexact input. The inverse seed (3*d) XOR 2 is correct modulo 32 for odd d; four Newton iterations x=x*(2-d*x) double the valid bits to at least 64. All arithmetic is unsigned modulo B. Power-of-two and scalar source cases retain direct shifts/native scalar division. A multiprecision divisor keeps the existing general quotient path under the documented exact-divisibility precondition.

Wide products use compiler-available unsigned 128-bit arithmetic, MSVC x64 _umul128, or the existing pure C++20 32-bit-half multiply. No new ISA requirement or dispatch option is added. The inverse fixes the exact low digit; the wide product supplies the next high carry; inputs, output and all temporary factor/product buffers are disjoint. Rational's <=8-limb cancelled product can borrow unit-divisor factors and store other scalar-divisor quotients in independent eight-limb arrays, then materialize only the final product. Longer operands/multiprecision divisors retain owning fallbacks. All gcd and canonical Rational proofs remain unchanged.

No performance result or threshold acceptance is claimed in advance.

**R21 completed validation and x86-64 baseline.** Commit 97e9ae3 passed Unix, Windows and Coverage. The 256-bit total remains above 1x; compare at 128 bits has almost no margin under 2x. Clang 16 on AMD EPYC 7763 reports .501/1.163, .668/2.018, 1.100/2.029, 1.552/2.943, 1.881/3.223, 2.473/3.970. Do not treat successful benchmark jobs as performance-gate acceptance. Hosted Clang's large-width square path is a separate existing deficit. The artifacts were downloaded read-only and matched the SHA-256 digests returned by GitHub.

**Validation.** Apple M2 Ultra, macOS 26.6.2, AppleClang 21 Release warnings-as-errors: all 20 C++ suites pass. Debug ASan/UBSan passes Integer, BigUnsigned, Rational, Float and FloatRounding (5/5). GCC 16.2 ARM64 Release warnings-as-errors passes Integer, BigUnsigned and Rational (3/3). Independent fixed values and invalid scalar division/zero/undersized buffers are checked separately. Rational tests add controlled scalar cross-cancellations (including 14 and 15), all sign pairs, compound operators, fixed exact expected fractions and the eight/nine-limb boundary. The preserved full suite still covers multiprecision-divisor fallbacks and special values. No native x86-64/MSVC performance claim precedes the remote run.

**Native ARM64 performance.** Apple M2 Ultra, AppleClang 21, macOS 26.6.2, Release, unchanged six-width 20,000x31 comparison. Both R21 and final R22 pass 384 inputs and 2,880 paired exact results. Recorded runs were sequential and had no concurrent builds or tests. The earlier R22 trials are diagnostic; the table is the final source version.

| bits | R21 total | R22 total | R21 worst | R22 worst |
|---:|---:|---:|---:|---:|
| 53 bits | 0.360 | 0.369 | 1.235 | 1.236 |
| 128 bits | 0.573 | 0.538 | 1.837 | 1.738 |
| 256 bits | 0.902 | 0.882 | 1.701 | 1.687 |
| 1024 bits | 0.905 | 0.884 | 1.420 | 1.443 |
| 2048 bits | 0.729 | 0.709 | 1.298 | 1.153 |
| 4096 bits | 0.621 | 0.594 | 1.472 | 1.276 |
| ALL | 0.680 | 0.656 | 1.837 | 1.738 |

Absolute 256-bit native medians (R21 -> R22, ns):
- construct_fraction: 138.954 -> 137.442.
- add: 242.127 -> 241.981.
- sub: 241.113 -> 241.244.
- mul: 288.898 -> 274.852.
- div: 283.467 -> 276.931. Sum of native medians at 53 bits: 415.819 -> 419.818 ns. Sum of native medians at 128 bits: 1106.676 -> 1044.257 ns. Sum of native medians at 256 bits: 1604.316 -> 1602.129 ns. Sum of native medians at 1024 bits: 3325.132 -> 3321.063 ns. Sum of native medians at 2048 bits: 5470.241 -> 5485.076 ns. Sum of native medians at 4096 bits: 10695.727 -> 10704.573 ns.

Integer 20,000x31 passes 768 inputs/2,688 results. No Integer gate guarantee is added. FloatMP 1,000x11 passes 256 inputs/3,456 exact results for both revisions. The existing short-width FloatMP gates remain unmet; no FloatMP arithmetic was changed. Raw first/repeated runs are retained for comparison. Native x86-64 GCC/Clang measurements and Unix/Windows/Coverage CI for R22 remain pending at publication.

### R23 — fused unsigned product accumulation for Rational sums (21 October 2026)

 **Baseline.** R22 (2f9e710) passed Unix, Windows, Coverage and both hosted x86-64 benchmark jobs. User-supplied 20,000x31 reports, arm64 then x86-64, both pass 384 inputs/2,880 exact results. On x86-64: .896/1.398, .875/1.966, 1.234/1.762, .959/1.413, .647/1.107, .499/.985. Only the 256-bit total fails.

**Shared unsigned primitive.** accumulate_product_into(a,b,out,subtract) updates an existing disjoint accumulator by +/-a*b modulo B^out.size(), B=2^64, and returns the one-bit overflow or borrow. For nonempty factors require out.size()>=a.size()+b.size(); check size arithmetic and storage before mutation. A zero factor leaves the accumulator unchanged. Both read-only factors may alias each other but not output. No allocation, scratch storage, finite precision limit or new ISA requirement is introduced; no existing kernel selector changes.

Each row forms a_i*b_j+carry, applies its low limb to the accumulator, and carries its high limb plus addition overflow/subtraction borrow. The combined intermediate is at most (B-1)^2+2*(B-1)=B^2-1, so the next carry fits a limb. Propagate the last row carry/borrow through all higher accumulator words. Since a*b<B^out.size(), adding it to an unsigned accumulator can overflow at most once; subtraction can borrow at most once. The OR of escaping row bits is therefore the exact final overflow/borrow. Skip exact-zero multiplier rows; other rows retain the same full-width recurrence. Optional compiler unsigned 128-bit multiplication, MSVC x64 _umul128, and the existing portable 32-bit-half multiply provide the same low/high products. There is no CPU feature probe or assembly addition in this change.

**Rational integration.** For bounded cross sums, form the first product in width+1 limbs, where width covers both full products, and accumulate the second directly. One guard limb proves the unsigned sum exact. For opposite effective signs, a returned borrow identifies a negative difference; convert its two's-complement buffer back to an unsigned magnitude in place, then combine its sign with that of the first numerator. Signed interpretation remains in Rational. Zero is canonicalized by Integer and the existing Rational zero branch. This replaces two product arrays plus a separate add/sub output by one 17-limb array and removes the separate product comparison. All scalar, unit-factor and >8-limb operand fallbacks are retained.

Completed nonzero canonical add/sub components also move directly into the result, followed by setting kind=finite. Both owning values are fully evaluated before member assignment, preserving aliases. Member move assignment is noexcept; public Rational move construction, move assignment and their existing source-value semantics are unchanged. The direct assignment alone did not show an ARM64 speedup in the short trial, so it is not claimed as an independently measured optimization. No GCD, reduction proof, benchmark operand, timing weight or oracle changes.

**Functional validation.** Apple M2 Ultra, macOS 26.6.2, AppleClang 21 Release warnings-as-errors: all 20 C++ CTest suites pass. Debug ASan/UBSan: Integer, BigUnsigned, Rational, Float and FloatRounding pass (5/5). GCC 16.2 ARM64 Release warnings-as-errors: Integer, BigUnsigned and Rational pass (3/3). New public replacement/self-alias tests cover both source orders, zero cancellation, every denominator-GCD reduction branch and replacing an existing special destination. No x86-64/MSVC performance or remote CI success is claimed before those jobs run.

**Native measurements.** Sequential R22/R23 Release runs on the same M2 Ultra, without concurrent builds/tests: Rational 20,000x31, all six widths, both pass 384 paired inputs/2,880 exact results.

| Bits | R22 total | R23 total | R22 worst | R23 worst |
| --- | ---: | ---: | ---: | ---: |
| 53 | .368 | .362 | 1.257 | 1.238 |
| 128 | .537 | .523 | 1.762 | 1.746 |
| 256 | .869 | .851 | 1.680 | 1.723 |
| 1024 | .886 | .884 | 1.417 | 1.552 |
| 2048 | .703 | .705 | 1.121 | 1.244 |
| 4096 | .610 | .611 | 1.273 | 1.282 |
| ALL | .664 | .662 | 1.762 | 1.746 |

All Rational ARM64 gates pass. Absolute add/sub medians (R22 -> R23, ns): 53 bits 90.767/89.033 -> 81.604/80.842; 128 bits 188.167/186.790 -> 168.244/171.060; 256 bits 242.121/242.142 -> 217.710/220.954; 1024 bits 615.410/600.740 -> 577.708/578.765. The 256-bit add/sub improvement is about 10%/9%; it does not establish that the user's x86-64 total gate now passes. Native x86-64 reports and remote Unix/Windows/Coverage results remain pending at publication.

Integer 20,000x31 passes 768 inputs/2,688 results in both revisions. FloatMP 1,000x11 passes 256 inputs/3,456 exact results in both revisions. These single sequential runs include timing variation; they do not prove universal nonregression. Existing FloatMP short-width deficits remain. No FloatMP arithmetic or existing Integer kernel selection changed.

### R24 — diagnostic of common short GCD costs (21 October 2026)

The 256-bit total still violates the strict per-width <1 gate; compare at 128 bits is close to the <2 worst-case gate. R23 revision 904d41fe passed Unix, Windows, Coverage and hosted Rational benchmark workflows. A successful hosted benchmark job does not, by itself, constitute a performance-gate pass.

**Diagnostic only.** Add \`benchmark_rational_gcd_profile\` to the normal benchmarks build, never to CTest. Recreate the unchanged official 32-pair Rational operand fixture at all six official widths. Measure separately the two construction reductions, common denominator GCD, the final gcd(n,g) for addition and subtraction, the two multiplication cross-GCDs, and numerator/denominator division cross-GCDs.

A diagnostic stack-backed binary GCD (Stein) runs only for magnitudes of at most 8 limbs, using independent local buffers for shift and subtraction and a final exact common-power-of-two restoration. For larger inputs, the binary candidate falls back to the production algorithm for independent correctness boundary checks and is reported as inapplicable in timing data. No input set, oracle, operation weight, threshold, or reporting field of the original official Rational benchmark is changed.

**Hosted x86-64 measurement.** The existing two-compiler benchmark workflow now builds and runs the R24 diagnostic on its candidate only, after the unchanged R23-parent/R24 official comparison, and uploads the diagnostic with host information and raw official logs. These are **not** projected speedups. Consider production promotion only after inspecting measured x86-64 GCC/Clang costs, separately for both the GCD and complete Rational operations, with no ARM64 regression and all Unix/Windows/Coverage checks.

**Status at authoring:** no local native x86-64 execution or benchmark result has been claimed, and no GCD optimization is promoted.

### R24.1 — Production-equivalent experimental GCD measurements (21 October 2026)

The 3000-iteration, 11-sample R24 diagnostic reports show that the
experimental eight-limb Stein GCD loses to the existing Euclidean
implementation in every measured category: at 256 bits, binary/native
ratios range from 11.3–23.9x on ARM64 and 8.0–17.1x on x86-64.
The ARM64 256-bit construct_b case is close to parity
(111.458 vs 112.917 ns), but larger widths remain slower.

The diagnostic only measures secondary reduction GCD calls when
`g=gcd(den_a,den_b)>1` and the cross-sum is nonzero, matching
the actual `Rational::add_sub` path. In the 256-bit fixture, 8 of
32 operand pairs have `g>1` (all `g=2`), while 24 have `g=1`.
Timings that include the skipped calls are not representative
of the production workload. **Decision:** do not select the
experimental Stein GCD over the established Euclidean kernel.
### R25 — nonproduction bounded Lehmer experiment (21 October 2026)


**Observed R24.1 input.** The user supplied corrected ARM64 and x86-64 R24 diagnostics (3,000 iterations, 11 samples, all GCD checks pass). Binary Stein loses in every measured applicable category on both architectures and is not a promotion candidate. At 128 bits, the conditional add_reduce GCD costs 79.068 vs 34.492 ns but runs only for 10 of 32 benchmark pairs. At 256 bits, add_reduce runs for 8 of 32. Consequently, optimizing only a low-limb GCD tail cannot by itself close the complete 256-bit Rational total-time gap. The earlier R24.1 CI succeeded for Unix and Coverage; Windows and the hosted benchmark were not yet complete when the user reported these measurements.

**Experiment.** Add a self-contained header under benchmarks/ only, with a compiler-128-bit bounded Lehmer GCD candidate (at most eight 64-bit limbs). For equal-length multi-limb inputs, extract the same 127 leading bits from both, compute up to 20 agreed quotient steps using a matrix with small signed coefficients, and apply the resulting unimodular matrix to *entire* limb arrays using signed base-2^64 carry arithmetic. Reject nonrepresentable, negative, or nondecreasing transforms; otherwise their determinant +/-1 proves exact gcd preservation, even though the prefix approximations are truncated. The fallback performs an exact Euclidean step through existing public BigUnsigned::divrem_into, deliberately charging its *both-output* cost rather than altering any private production remainder selector. A scalar tail uses exact 128-bit remainder and scalar std::gcd. Coefficient bounds of 2^19 and at most 20 steps bound all signed 128-bit intermediates. Inputs and all work arrays are disjoint. Non-GCC/Clang platforms without __int128 and operands >8 limbs retain production gcd, with no misleading candidate timing reported.

**Acceptance rule.** A candidate must be meaningfully faster on native x86-64 in the costly 256-bit GCD categories, without harm to ARM64, before considering a selective production kernel. Do not alter the official benchmark or claim R25 performance improvement before these runs.

### R26 — prove and benchmark quotient-three shortcut in bounded GCD (21 October 2026)


On x86-64, Lehmer/native at 256 bits: construct_a=1.848, construct_b=1.193, add_den=1.789, mul_cross_a=1.609, mul_cross_b=1.294, div_num=1.157, div_den=1.715. On ARM64 the corresponding range is 1.414-2.321. Neither binary GCD nor this conservative Lehmer implementation is selected for production.

**Candidate change.** Expand the existing R11/R21 short-quotient subtraction in *bounded* BigUnsigned::gcd from at most two to at most three subtracts for ordered, equal-size canonical input views. Writing a=high_a*B^k+low_a and b=high_b*B^k+low_b, a>=b and high_a<=3*high_b imply a<(high_a+1)*B^k<= (3*high_b+1)*B^k<4*b (the last strict inequality follows from b>=high_b*B^k and the integrality of the quotient; for high_b=1, a<4*B^k<=4*b). Consequently floor(a/b)<=3. This is a general arithmetic bound, not a corpus-specific recognition test. Test the condition with unsigned differences high_delta=high_a-high_b, high_delta<=high_b || high_delta-high_b<=high_b using short-circuit evaluation, without overflowing 3*high_b. The first subtraction stores in the next disjoint GCD buffer; each additional subtract uses the independent division workspace and copies the trimmed output back. A lexicographic comparison ensures no subtraction occurs when the remainder is already <b. The normal general remainder-only division remains the fallback for larger quotients and unequal sizes. No change to quotient or remainder public kernels or Rational invariants.

Existing randomized GCD, decimal, multiplication, Rational and FloatMP tests remain enabled.

**Status.** R26 is a measured candidate, **not** an asserted speedup: the x86-64 256-bit full Rational total and 128-bit worst operation must improve or at least stay within the project gates without regressing ARM64 and high widths; use the hosted same-run parent/candidate x86-64 comparison and independent ARM64 native reports. If the extra comparisons/subtractions cost more than the avoided general divisions, revert production dispatch while retaining useful diagnostic and design findings.

### R26.1 — reject quotient-three production dispatch after dual-host measurement (21 October 2026)

Earlier unchanged-benchmark R23/R25 baselines from the same respective hosts were x86-64 .856/.819/1.176/.938/.624/.483 and ARM64 .402/.557/.852/.876/.695/.604. The apparent 256-bit x86-64 reduction 1.176x -> 1.170x is approximately 0.5% relative, too small to establish a dependable performance improvement from one set of hardware measurements.

**Decision.** Revert *only* the R26 quotient-three modification to src/big_unsigned.cpp and restore its precise R25 production contents (unaltered R21/R11 quotient-one/two bounded GCD). Preserve all prior R24/R25 diagnostics and source documentation, clearly labelling their candidates as rejected/unselected. No benchmark fixture, sample count, threshold, operation, oracle or architecture gate changed. The four GitHub checks on R26 (Unix, Windows, Coverage and x86-64 Rational benchmark) completed successfully; they establish functional correctness, not sufficient performance. The R26.1 revert must pass the normal checks independently. Future efforts should target the actual x86-64 five-limb / 256-bit Rational bottleneck rather than widening repeated-subtraction dispatch without independent repeatable end-to-end benefit.

### R27 — short GCD scalar-tail and equal-prefix candidate (21 October 2026)


**Basis.** The R25 binary/Lehmer diagnostics reject both expensive short-GCD candidates. Investigate the frequent transition to one-limb gcd tails instead: for a canonical multi-limb dividend A and scalar divisor d, the first high word A_top needs no 128/64 division when A_top<d. More generally, for *all* nonzero d, an exact first reduction A_top%d, followed by the existing per-limb recurrence, preserves A mod d. For powers of two, the complete remainder is A_low&(d-1). For two ordered equal-length operands with identical limbs 1..n-1, gcd(A,B)=gcd(A,A_low-B_low), and equal operands return A unchanged. These laws hold for arbitrary inputs, not just the workload's 2^k+small offsets.

**Candidate implementation.** Introduce a nonallocating local scalar Euclid helper in big_unsigned.cpp for dividend lengths <=8, using the exact DIVQ-backed divide_wide on remaining words and ordinary native scalar Euclid; skip the top-word division when it is already smaller than d. Use it for the bounded GCD's scalar tail and the large-fallback scalar-tail cases of up to eight limbs. For >8 limbs, retain B61's preinverted remainder dispatch and existing scalar algorithm. Before entering bounded GCD for 2..8 limbs, check whether the complete high prefix matches. If so, return the full original magnitude on equality, otherwise compute gcd(A,A_low-B_low) using the same scalar helper without allocating the three rotating remainder buffers. Other multi-limb paths and quotient-one/two shortcuts are unchanged.

**Measurement rationale and limitations.** A separate local GCC 14 x86-64 microbenchmark used the exact algebraic Rational 256-bit fixture to construct the scalar-tail inputs outside timed loops. On that *different* Xeon 8573C host, n=4-5 scalar remainders took about 13-20 ns for regular DIVQ recurrence, versus about 13-16 ns for top-divide-skipping and about 14-18 ns for preinverted pi1, depending on GCD site. The gains are small and must NOT be extrapolated to the user's x86-64 CPU or to complete Rational operations. More elaborate reciprocal-32 variants were slower on the same local host.

Do not modify workload operands, the numerical oracle, default build flags, 20,000x31 sample schedule or Rational correctness. The new production change is **provisional** pending Unix, Windows, Coverage and same-host x86-64 parent/candidate CI as well as independent ARM64 comparison. A tiny win within ordinary measurement noise is not enough to justify a broadly slower path.

### R28 — Rational add/sub and to_double candidates (21 October 2026)


ARM64 total at 256 bits is 0.876x. Prior binary, Lehmer and quotient-three GCD attempts did not close this gap, so prioritize these three operations without modifying mul/div dispatch.

**R28a: get_d exponent ordering.** The previous exact comparison of n with denominator*2^exponent, or numerator*2^(-exponent) with denominator, materialized a shifted owning BigUnsigned. Replace it by exact borrowed-limb comparison. Compare the bit lengths first; if equal, reconstruct each virtual shifted limb including the 64-bit offset and cross-word carry, from most to least significant. Protect size_t bit-length addition from overflow. No integer values are approximated, no rounding mode depends on the helper, and the pre-existing binary64 truncation, prefix certificates, specials, underflow and overflow stay untouched. This can remove a heap-backed intermediate at the 4/5-limb boundary.

**R28b: bounded positive product accumulation.** For GCC x86-64 with unsigned __int128, when both factor views have at most eight limbs and R23 accumulate_product is adding (not subtracting), fuse a_i*b_j+output[i+j]+carry into one exact 128-bit operation and split its low/high limbs. This sum is less than B^2 for B=2^64. Keep zero-row skipping, carry propagation and the escaped-carry bit as in R23. All other targets, negative accumulations, large factors, APIs and all Rational dispatch paths remain unchanged. Any speedup remains to be measured independently and cannot be inferred from compiler or operation counts.

The six-width 15-operation official 20,000x31 benchmark and its 384 inputs/2,880 expected results are unchanged. Run Unix, Windows, Coverage and paired x86-64 GCC/Clang official benchmarks; request an independent ARM64 run. R28 is provisional until both the whole Rational workload and these isolated operations show reproducible benefit with no architecture regression.

### R29 — retain R28 double conversion, retire unproven add accumulator, certify sparse prefixes (21 October 2026)


**R28 measurements.** The user provided native ARM64 and x86-64 official Rational 20,000x31 runs. Both pass 384 paired inputs and 2,880 results. R27 -> R28 at 256 bits: ARM64 total .876x -> .835x, to_double 64.390 -> 41.237 ns, add 203.394 -> 198.042 ns, sub 190.767 -> 192.492 ns; x86-64 total 1.164x -> 1.154x, to_double 104.773 -> 85.729 ns, add 355.474 -> 357.068 ns, sub 348.538 -> 349.180 ns. These tiny add/sub differences do not establish a benefit from R28b GCC accumulation, which added an extra bounded kernel path. All four R28 CI workflows passed. The 128-bit x86-64 worst-case compare was 2.019x (strictly above the <2x gate), while 256-bit total was still 1.154x. Compare code did not change; do not attribute this small timing variation without paired evidence.

**Exact binary64 R29 candidate.** Inside the R7 two-limb high-prefix certificate, if both leading two-limb prefixes have exactly zero low word, write N=n_top*B, D=d_top*B, s=-scale (s<=63). If n_top*2^s fits in 64 bits, compute exact q=floor(n_top*2^s/d_top) and r=(n_top*2^s)%d_top using native 64-bit arithmetic. The full-prefix remainder is R=r*B. Accept only if 2^52<=q<2^53 and (r>0 or q==2^52). This is the same existing R7 exact certificate: r>0 implies R>=B>q, while q==2^52 has a lower bound from the exact exponent. The strict upper bound holds because D-R=(d_top-r)*B>=B>2^s. Ambiguous or out-of-bound cases fall back to the existing 128-bit prefix division. The predicate recognizes a general sparse leading-prefix property, not benchmark operand values. All normal, subnormal, special, negative and other precision cases keep the existing semantics, with no floating estimate or rounding mode dependency. ARM64, MSVC and Windows remain unchanged.

The unchanged Rational 15-operation official 20,000x31 benchmark and all strict cross-platform gates still apply. R29 remains a provisional performance candidate pending CI and independently measured ARM64/x86-64 end-to-end results.

### R30 — sparse-radix fraction constructor reduction (21 October 2026)


The x86-64 full Rational 256-bit total remains 1.143x. Both exact official benchmarks pass. R29 Unix, Windows, Coverage, Rational benchmark workflows passed.

**Provable constructor-only fast path.** Preserve the public numerator/denominator constructor, special and zero semantics, sign normalization, Reduced constructor and exact quotients. Before GCD normalization, inspect both nonzero magnitudes for the same 2..128 limb width with all intermediate limb positions 1..width-2 equal to zero. Write N=p*B^k+n0, D=q*B^k+d0 with B=2^64 and k=width-1. Choose the operand with the larger high coefficient as H, the other as L. Only when high(H)=t*high(L) exactly and t*low(L) fits in uint64, delta=abs(low(H)-t*low(L)) also fits in uint64. Euclid gives gcd(N,D)=gcd(L,delta). If delta=0, return L; if delta=1, return 1. For delta>65535 use existing gcd(L,one-limb delta). For 1<delta<=65535, compute B^k modulo delta with binary exponentiation, and the scalar remainder from the two nonzero limbs of L, then return std::gcd(delta,remainder). Residues <=65534 keep all products and additions safely in uint64; no 128-bit integer or ISA assumptions.

All predicates describe general integers rather than fixture values. Dense inputs bail out on their lowest intermediate word; unmatched limb counts, unmatched high coefficients, low-product overflow and inputs beyond 128 limbs use the unchanged full Euclidean GCD. No BigUnsigned or other Rational operator implementation is changed.

Keep the original six-width Rational benchmark, 20,000x31 samples, 384 paired inputs and 2,880 expected results unchanged.

**Status and acceptance.** R30 is a provisional constructor optimization, not yet an integrated or production performance pass.

### R30.1 — Native construction benchmark confirmation (21 October 2026)

Against R29 on the same respective hosts, R30
`construct_fraction` improved at 256, 2048 and 4096 bits:
ARM64 **114.196 -> 94.010**, **386.390 -> 155.429**
and **692.450 -> 219.613 ns**; x86-64
**211.755 -> 136.206**, **336.809 -> 214.528**
and **537.964 -> 278.859 ns**.
ARM64 aggregate 256-bit runtime ratio changed
from 0.851x to 0.827x. These are independent same-host
full-workload measurements, not synthetic proxy estimates.
### R31 — eliminate final Rational mul temporary; extend sparse GCD to add/sub and mul (21 October 2026)


**Motivation.** R19 already borrows unit-cancelled multiply factors; R20/R23 already fuse signed addition's two cross products into one bounded scratch array and one owned final magnitude. R22 already divides scalar cancellation factors into caller-owned buffers. Thus additional blind removal of product temporaries is unlikely to explain the x86-64 versus ARM64 gap. The R24 GCD breakdown showed that multi-limb GCD is disproportionately expensive on x86-64.

**Change 1: direct mul result commitment.** Multiply all four cancelled factors into two fully owned BigUnsigned outputs, then create the signed Integer numerator and pass both to FiniteOperations::set_sum. This avoids construction of a temporary Rational via reduced() and its move/swap assignment. The caller destination is only mutated after the products are calculated, retaining multiply self-alias correctness. Sqr and div paths are intentionally unchanged.

**Change 2: generic sparse-radix GCD reuse.** Rename the R30 constructor-only helper to gcd_sparse_radix and use it also for add/sub's gcd(den_a,den_b) and mul's two cross-GCD operations. For arbitrary high coefficients n=p*B^k+n0, d=q*B^k+d0, compute the exact determinant delta=q*n0-p*d0 if both scalar multiplications fit 64 bits. Because q*n-p*d=delta, gcd(n,d) divides delta. If delta!=0 and gcd(p,delta)=1, then gcd(n,d)=gcd(n,delta): a divisor of n and delta divides p*d and is coprime to p. Symmetrically if gcd(q,delta)=1, use gcd(d,delta). Otherwise fall back unchanged to the full BigUnsigned::gcd. For a gap<=65535, the specialized sparse-scalar step computes B^k mod gap in logarithmic time by exact modular exponentiation; larger gaps use the existing general scalar GCD kernel. All inputs and the result remain canonical, with no new ISA assumptions, precision limits or dependence on specific benchmark operand values. The original zero and special Rational cases remain unchanged. This general determinant route particularly handles denominators with high coefficients 2 and 3 without assuming a high-coefficient multiple.

No automated Rational benchmark job is restored: the user-requested workflow deletion remains intact.

**Status.** R31 is a candidate. Verify Unix/Windows/Coverage and request independent native ARM64/x86-64 20,000x31 measurements before attributing end-to-end gains. Never treat a workflow success as a performance-gate acceptance.

### R32 — reuse exact sparse GCD and direct result commitment in Rational div (21 October 2026)


**R31 native measurements.** On the unchanged official Rational 20,000x31 workload the user supplied independent ARM64 and x86-64 reports. Both pass numerical verification (384 paired input cases; 2,880 exact operation results). R30 -> R31 at 256 bits: ARM64 add 199.883 -> 180.285 ns, sub 192.454 -> 176.277 ns, mul 244.844 -> 217.831 ns, total 0.827x -> 0.777x. x86-64 add 353.754 -> 311.886 ns, sub 346.689 -> 304.128 ns, mul 475.300 -> 412.156 ns, total 1.103x -> 1.030x. R31 passed Unix, Windows and Coverage GitHub Actions; the retired Rational benchmark workflow remains absent. The x86-64 128-bit compare worst is 1.978x and 256-bit total still fails strict <1.

**Candidate change.** R31 introduced a proved generic sparse-radix GCD used by constructor/add/sub/mul, but division still performed its two cross-cancellation GCD calls through the generic multi-limb gcd. Replace both calls with the exact same gcd_sparse_radix, preserving its full fallback for dense, incompatible, overflowing or non-certifiable operands. Once the cancelled numerator and denominator are owned, commit them directly using FiniteOperations::set_sum, as in R31 multiplication. This removes the temporary reduced Rational and its move/swap. All inputs remain borrowed until both outputs are computed, keeping aliasing safe. No change to the arithmetic identity, numerical representation, zero/special cases, algorithms in BigUnsigned, or benchmark inputs/cadence.

Run Unix, Windows and Coverage independently. Measure ARM64 and x86-64 with unchanged 20,000x31 official benchmark for promotion; a successful CI alone does not show performance gain. R32 remains provisional until evaluated.

### R32.1 — dual-host native benchmark confirms nominal Rational gates (21 October 2026)

The x86-64 256-bit total 0.990x has only a ~1% margin and the x86-64 128-bit worst 1.980x has only ~1% margin. Neither is a repeatability claim or proof of robust platform superiority. R32 is therefore a successful *measured candidate*, but should be remeasured under matched-run stable host conditions before declaring a reliable performance gate.

**CI distinction.** The original R32 Unix workflow's Ubuntu GCC 13 Release, Ubuntu Clang 16 Release, and Ubuntu Clang 16 ThreadSanitizer jobs all succeeded. The macOS AppleClang 21 Release job initially failed during GitHub Actions 'Set up job' (before compilation); this is a runner setup issue rather than evidence of a C++ source failure. Coverage succeeded. Windows and any retried macOS checks must be tracked separately until completed. The user-requested stand-alone Rational benchmark CI remains deleted; no performance thresholds, sample counts, workload or test cases have been changed by this journal entry.

**Next action.** Freeze the R32 arithmetic implementation for now. Prefer repeated independent ARM64 and x86-64 native 20,000x31 measurements and confirmation of all CI platform tests over risky further production dispatch alterations when the strict gates are only narrowly met.

## NewDyadic continuation

The experimental `NewDyadic` work beginning with commit [`0f7f06c16`](https://github.com/ariadne-cps/numeric/commit/0f7f06c16) (Y33a), including Y33b–Y38, is maintained separately in [`dyadic_design.md`](dyadic_design.md). The shared BigUnsigned foundation is documented in [the BigUnsigned journal](big_unsigned_design.md), and signed Integer history in [the Integer journal](integer_design.md); Rational historical records remain in this journal. New shared-kernel changes that affect both projects should be cross-referenced in each journal.

### R33 — oracle migration and removal of reverse Rational external reference bridge methods (22 October 2026)


Do not remove or hide that live dependency until those FloatMP paths are made native. This is test-only and bridge surface cleanup; no arithmetic implementation, architecture or performance gates were modified. CI validation remains required.

### R33.1 — remaining reverse-conversion callers (22 October 2026)


Additional exact-conversion callers in Rational required fixture updates. Rational R4 moves all 54 canonical floor/ceil cases to an independently computed TSV. These changes do not change numerical kernels or performance requirements. Build/CI validation required before claiming closure.

### R33.2 — R16 reference fixture migration (22 October 2026)


R16 is now fully source-keyed in rational_r16.tsv: all 84 original configurations are retained, not reduced because the seven GCD classes and both operand orders are semantically distinct. Preserve R23 assignment/alias and special-destination checks. No production arithmetic changes; CI validation pending.

### R33.3 — R0 randomized reference migration (22 October 2026)


The original 1,200-case Rational random differential loop now tests 120 source-keyed rows from rational_random.tsv, selected adaptively one from each ten visits with full 17-width coverage, raw sign combinations, zero operands and comparison classes. Retain the full mt19937_64 advancement without discarded arithmetic to keep downstream R29/R13 input identity. No production kernels touched. Verify the new commit on Unix, Windows, Coverage and Sanitizers.

### R33.4 — R20 exact cross-add/sub oracle (22 October 2026)


Verify direct and compound assignment arithmetic against canonical TSV numerators and denominators; preserve cancellation and carry pattern classes. No production changes. CI validation pending.

### R33.5 — R19 factor cancellation reference migration (22 October 2026)


Source-keyed TSV contains canonical operands and independent product/division outcomes for both orientations; test checks direct and compound operations and self-aliasing. CI validation pending.

### R33.6 — R22 scalar cancellation oracle migration (22 October 2026)


Freeze all 160 R22 cases in tests/data/rational_r22.tsv, retaining ten original 53/127/128/255/256/447/448/511/512/1024-bit thresholds, factors 1, 1/15, 14, 14/15, all four sign masks, and the two exact cancellation identities. Oracle uses independent signed big integers with gcd normalization. Compare native direct/in-place multiplication and division to exact normalized numerator/denominator from TSV. Build/CI still required.

### R33.7 — freeze R30, R17 and R3 rational independent exact-reference fixtures paths (23 October 2026)


Preserve ALL R30 10 widths, 7 high-coefficient pairs, 6 low-tail pairs, 4 sign masks and 9 dense cases: 1689 rows (including the 8192-bit and 128/129-limb cases), with both canonical numerator and denominator frozen separately. Preserve all 252 R17 7 x 6 x 6 comparison signs and 6 R3 boundary comparisons, including both operand orders, negative reversals and zero ordering.

Readers reject missing, misordered and extra rows. Do not weaken canonical numerator/denominator assertions to only Rational equivalence: the constructor's normalized representation remains an explicit tested contract.

Cross-platform CI and source-keyed historical test suites must pass before accepting this as production test migration.

### R33.8 — source-keyed R31/R32 canonical arithmetic oracle (23 October 2026)


Retain all ten widths 128/192/255/256/320/448/512/1024/2048/4096 bits, all 24 tail configurations and all four operand-sign combinations at each width (960 original operand pairs). Keep direct add/sub/mul/div and in-place +=, *=, /= checks: 6720 operation outputs, plus canonical source left/right operand checks (1920 source fractions). Compare canonical signed numerator and positive denominator components, not only rational equality, so reduction errors cannot be disguised by expected-value normalization.

Files tests/data/rational_r31_<bits>_<chunk>.tsv cover all 960 cases without sampling. The 2048-bit rows are divided into two 12-index chunks and 4096-bit rows into four 6-index chunks because the GitHub blob connector limits very large single-file transfers. All other precisions use one chunk. Runtime readers enforce the original iteration/index/sign sequence, missing rows and no trailing rows in every chunk; the file split does not change numerical test semantics. No change to benchmark, CMake, production Rational or FloatMP. Treat this commit as unvalidated until fresh cross-platform CI confirms it.

### R33.9 — independently frozen R7 binary64 conversion oracle (23 October 2026)


Preserve all 8 original widths (53/128/256/448/512/1024/2048/4096), all 12 original near-boundary/repeating/negative/overflow rational probes and all four FE_* host rounding modes: 384 expected values, no sampling. Inputs still come from the original explicit Rational construction expressions; expected bits were computed independently from exact signed multiple-precision integer numerator/denominator pairs by finding the exact binary exponent and truncating the exact significand quotient toward zero, with full IEEE overflow and subnormal handling. Unlike a numeric float comparison, raw-bit equality catches sign, last-bit, infinity and rounding-mode errors.

The reader checks mode/width/probe source keys and rejects missing, out-of-order and extra rows. Test-independent CMake and benchmark sources are unchanged.

### R33.10 — remove R12 runtime external reference conversion oracle (23 October 2026)


All 512 original samples are retained for each of the four FE rounding modes, and both signs per sample are checked: 2048 source-keyed rows and 4096 bit-exact binary64 outcomes (no thinning or skipped cases). The existing std::mt19937_64 engine, seed, draws and all four input classes remain unchanged. The TSV explicitly freezes the original numerator and denominator next to each expected positive/negative uint64 bit pattern. Test readers verify original source values, mode/index order, missing rows and excess rows.

The std::mt19937_64 sequence was validated against a separate compiled C++20 generator for both seeds used in the conversion group. All original R5/R7/R13/R29 conversion checks and subsequent random-stream consumption remain unchanged.

Production numeric code, benchmarks and CMake are untouched. Cross-platform CI is required before declaring this migration verified.

### R33.11 — migrate R29 and R13 binary64 conversion oracles (23 October 2026)


Migrate the complete remaining R29 sparse-prefix and R13 dense random Rational::get_d() differential blocks to independently computed source-keyed TSV expectations. R29 preserves 4 host modes, 6 widths, 3 numerator and 3 denominator high-word coefficients, and both signs: 216 source rows, 432 exact binary64 outcomes. R13 preserves 4 host modes, 4 widths (5..8 limbs), all 96 original cases per width: 1536 source rows and exact binary64 outputs. No random cases are skipped. R13 TSVs are split by mode and limb count for manageable GitHub transfers, without changing the original 1536 cases. Each reader rejects missing, misordered and extra rows. Cross-platform CI is required before treating this migration as verified.

### R33.12 — migrate remaining R5 get_d independent exact-reference fixtures (23 October 2026)


Preserve the original 12 shifts (0/1/52/53/1022/1023/1024/1073/1074/ 1075/1076/2048), four FE rounding modes and four signed numerator/ denominator scale patterns per visit: 48 exact source pairs and 192 raw binary64 checks. No sampling, operand substitution or changes to the historical std::mt19937_64 stream. Offline source generation reproduces the previously advanced R0 and R29/R13 generator consumption, and was checked against a separately compiled C++20 std::mt19937_64 source producer including five R5 samples. Compare all four raw uint64 patterns and check exact generated source inputs, row order and no extra or missing rows. Maintain explicit mismatch diagnostics for mode, shift, operation, actual and expected patterns.

Production arithmetic, benchmark reference, and all CMake files are unchanged. Require Unix/Windows/Sanitizers/Coverage CI before declaring the test migration validated.

### R33.14 — frozen GCD oracles for BigUnsigned R6/R8/R10/R11 (23 October 2026)


Preserve R6 7 widths x 4 offsets x 2 orders = 56 cases; R8 all 640 original mt19937_64 source pairs plus 6 x 6 wide/small checks = 676 cases; R10 nine widths x five borrowed-source pairs and one self-alias case = 54; R11 eleven widths x five Euclidean quotient-boundary cases = 55. Total 841 records, no sampling. R8 source stream was checked against a separate compiled C++20 generator for its first twelve output words; skipped conditional draws remain skipped.

TSV rows include original source-key parameters, exact canonical hexadecimal left/right operands, and independently computed exact unsigned Euclidean gcd values. Each test verifies canonical inputs, GCD(a,b) and GCD(b,a) against the same frozen oracle, then verifies source values again to preserve the original borrowed operand immutability assertions. The TSV readers reject missing, reordered, malformed and excess rows. The new bidirectional GCD assertions strengthen some original checks without eliminating any original source scenario, storage threshold or alias case.

No, CMake, production source or benchmark edits. Require four-workflow CI validation.

### R33.15 — freeze BigUnsigned GCD R26/R27/R21 runtime independent exact-reference fixtures (23 October 2026)


R26 retains 9 limb widths (2/3/4/5/8/16/ 127/128/129), 40 trials per width and all 9 quotient/variant combinations per trial: 360 frozen base sources and 3240 result checks. The TSV records source low/b values and all 9 expected gcd results; the source mt19937_64 still consumes exactly the original draws. R27 retains 7 x 96 shared-prefix and 8 x 96 scalar-tail cases = 1440. R21 retains all 10 x 4 dense samples = 40, plus existing exact fixed-factor identities at every width. Total 4720 original GCD source/result checks, without sampling or case substitution.

The three mt19937_64 seeds were independently crosschecked against compiled C++20 engine outputs. Tests enforce exact source hexadecimal, order/keys, missing/excess rows, canonical limbs and bit length, both argument orders and immutability. The R26 cases are split by width for manageable data files only, not for sampling. No edits to CMake, benchmarks or production numeric code. Require four-workflow CI before declaring R33.15 validated.

### R80 — independent rational R28/R30/R22 fixture generators (23 October 2026)


Add standalone standard-C++20 `tests/data/generator_rational_exact.cpp` to reproduce 36 R28 exact signed add/sub rows, all 1,689 R30 exact canonical normalization rows (including nine 128..8192-bit dense fallbacks), and all 160 R22 scalar signed multiplication/division rows. R28's raw operand columns remain unreduced as in the original frozen fixture, while resulting sum and difference columns are canonical. Seeds, source parameters, sign classes and row ordering are unchanged; no sampling or new test cases.

The C++20 generator was compiled and --verify executed with GCC and Clang using CI-style strict warnings including -Wsign-conversion -Wshadow -Werror. Regenerated outputs matched all previously committed Git object hashes: R28 7292632988a8b8d7e02f9e2a866704f46fce73d2; R30 3381e112eac6cbc56c30dcaa18a858053894014e; R22 fe5f8415bab7db3f4980c34bf41360e13d4b67b0. Read-only verification also rejected a deliberately corrupted fixture. Update tests/CMakeLists.txt only to compile the generator as a dependency of the aggregate `tests` target; it is not automatically run. Update README-generators.md. Coverage rises from 36 to 39 versioned source-reproducible TSVs. Full Unix, Windows, Sanitizers and Coverage CI remains required.

Existing source-tree TSVs remain mandatory in the meantime.

### R81 — version exact independent R16/R19/R20 arithmetic sources (23 October 2026)


Extend `tests/data/generator_rational_exact.cpp`, reusing its independent unsigned-limb signed rational, binary GCD, aligned-subtraction division and shift/add multiplication support. Historical pre-TSV C++ test inputs were recovered and represented directly, not read from existing fixture inputs. Reproduce R16 all 6 x 7 x 2 = 84 signed add/sub, doubling and alias rows; R19 all 6 x (4 x 4 + 1) = 102 signed cross-factor multiply/divide and self-alias rows; R20 all 10 x 4 x 4 = 160 signed cross-add/sub rows. In R20 pair_index=3 deliberately repeats the same Rational on both sides; this preserves all original exact cancellation cases.

Independent JavaScript BigInt exact arithmetic regenerated all three complete byte streams, independently confirming their Git blob SHA-1s: R16 ed55680b1d2ced1b911c5581276f62c9009b0902; R19 09fafdc815bd0bc58c0d41fd0385743f4218dd17; R20 f7895eb61aef5b04e4a23fa2d66bd8a26b2fe720. Existing original TSVs are unchanged; verification is read-only by default and explicit `--write` remains required for regeneration. `tests/CMakeLists.txt` already registers this generator; no new target or `generate-data` rule is added. Source-reproducible coverage reaches 42 of 115 TSV files. Check full Unix/Windows/Sanitizers/Coverage CI after commit. The preceding R80 macOS Coverage failure was in the Upload Code Coverage step, while its Ubuntu Coverage job and Windows workflow succeeded; this should not be misclassified as a generator compilation failure.

### R88 — independent Rational R7/R12 binary64 truncation fixtures (23 October 2026)


Add standalone `tests/data/generator_rational_binary64.cpp` for the original R7 4 x 8 x 12 = 384 source-keyed rational double-conversion cases and R12 4 x 512 = 2048 one-limb random conversion cases (both positive and negative outputs). Reproduce source numerators and denominators by exact powers of two and the original mt19937_64 (0x1200530012812026) replay, respectively. The expected raw IEEE-754 binary64 patterns are derived by exact bit lengths, bit-aligned integer quotient, normal/subnormal exponent/significand construction and explicit truncation toward zero.

GCC and Clang compiled and verified all output with C++20 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror. Both generated TSVs match their original committed blob SHA-1 values exactly: R7 9d8e75190d09b49bf91b183cd68e2c6afaa7f4b4; R12 0ead6f9d1aac895ac4712dfb8a73041f2ea9002f. Read-only verification detects deliberate fixture corruption. Register the eleventh executable in explicit generators and generate-data targets, not the tests target. Coverage rises to 61/115 TSVs. Full platform CI remains to be checked.

### R89 — shared RNG exact R29/R13/R5 generator (23 October 2026)


Add standalone `tests/data/generator_rational_shared_stream.cpp` for `rational_r29.tsv` (216 source-keyed sparse-prefix binary64 cases), all 16 `rational_r13_mode{0..3}_limbs{5..8}.tsv` files (96 exact signed dense-limb conversion cases each), and `rational_r5.tsv` (48 shifts, four signed conversion expectations per row). Reconstruct every source input using the original mt19937_64(0x6672616374696f6e) sequence: first consume every original 1,200 R0 visit schedule (four operands, each with all limbs plus one sign draw). For each of four subsequent host rounding modes, generate R29 followed by R13 limbs 5,6,7,8 followed by R5's 12 shift probes. This original interleaving is essential; advancing all R29 modes before R13 changes the random-source values.

GCC and Clang C++20 compilations succeeded with -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror. An additional Clang build used the strict extended CI warning set. All 18 generated fixtures matched the existing Git blob SHA-1 IDs, as checked against the current main tree. --verify rejects deliberately corrupted data without rewriting it. Register the twelfth executable with the explicit generators and generate-data targets; tests and default all still exclude generators. Coverage rises from 61 to 79 of the 115 unchanged, versioned TSV files. Full platform CI remains to be checked.

### R90 — complete deterministic Rational R31/R32 exact fixture regeneration (23 October 2026)


Add standalone standard-C++20 `tests/data/generator_rational_r31.cpp` using only the independent radix-2^64 `Nat` helper. Generate all 10 original widths (128, 192, 255, 256, 320, 448, 512, 1024, 2048, 4096 bits) x 24 original input index configurations x all four sign masks = 960 source-keyed rational add/sub/mul/div test rows. Numerators and denominators are reconstructed from the original test's formulas; all outputs are independently canonicalized using binary limb GCD and exact division, independent shift/add multiplication and signed fraction arithmetic. Retain all 14 original TSV file partitions: eight single chunks, two 2048-bit chunks of 12 indices, and four 4096-bit chunks of six indices, including exact original index-range header comments.

The generator compiled with GCC and Clang under strict C++20 warnings, including -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror; Clang passed the extended CI-style warning set. All 14 regenerated TSVs match the original Git object SHA-1 values byte for byte; read-only verification correctly rejects a tampered file. Source replay is purely deterministic and reads no original TSV as input. The complete generator ran locally in about 2.3 seconds excluding compilation on Linux; this is a diagnostic, not a bound for other machines.

Register the thirteenth standalone generator in `ORACLE_GENERATORS` for the explicit `generators` and `generate-data` targets, leaving `tests` decoupled. Update README. Source-reproducible coverage rises from 79 to 93 of 115 versioned TSVs. Full platform CI is still required.

### R91 — independent source-keyed Rational R0 random fixture (23 October 2026)


Add standalone `tests/data/generator_rational_random.cpp` to reproduce the complete preexisting `rational_random.tsv` from its original mt19937_64(0x6672616374696f6e) stream. Replay all 1,200 original visits (47,860 RNG draws), including four operand limb vectors and four sign draws at each iteration, keeping the 120 exact adaptive one-per-ten selected indices. Source integers, signed canonical fractions, add/subtract/multiply, square, half, optional division and reciprocal, floor/ceil, comparison and raw hexadecimal IEEE-754 binary64 bits are all computed with independent standard-C++20 base-2^64 limb arithmetic.

An independent exact BigInt replay reproduced all 120 original 19-field rows. GCC 14 and Clang 17 compiled with strict C++20 warnings including -Wsign-conversion -Wshadow -Werror; an extended CI-style Clang strict-warning compilation also passed. Both executables regenerated 579,270 bytes, exactly the original Git blob SHA-1 `1f1db8b7fe27af6c30dcaa18c7feb6c3ee6ce7de`. Read-only --verify detected deliberate corruption. Register the fourteenth generator in explicit generators and generate-data targets without coupling tests; update README. Coverage rises from 93 to 94/115 versioned TSV files. Full platform CI pending.

