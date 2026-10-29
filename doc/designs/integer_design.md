# Integer — native signed-magnitude design journal

**Repository:** `ariadne-cps/numeric` · **branch:** `main`
**Scope:** `Integer` sign semantics, signed-magnitude representation, canonical zero, quotient/remainder rules, signed parsing/formatting, self-aliasing, exponentiation, ownership and Integer-specific benchmarks.
**Companion journals:** [BigUnsigned](big_unsigned_design.md) · [Rational](rational_design.md) · [Dyadic](dyadic_design.md) · [FloatMP](floatmp_design.md)

This is the living implementation journal for the **signed `Integer` interface**. `Integer` builds on the independently maintained [`BigUnsigned` arithmetic substrate](big_unsigned_design.md). B-series identifiers denote historical milestones, not ownership: an unsigned kernel may have originated in an Integer benchmark and vice versa. Every milestone retained below appears with its original heading, date and technical text; shared low-level kernel milestones are retained in the BigUnsigned journal, not duplicated here.

## Current architecture and correctness contract

`Integer` is a signed-magnitude arbitrary-precision exact type with a canonical unsigned magnitude and a sign. Signed zero is normalized to zero, signed division truncates toward zero, the remainder has the dividend's sign, and arithmetic preserves alias and ownership semantics. The underlying unsigned multiplication, division, GCD, shifts and decimal kernels are designed and validated in [`big_unsigned_design.md`](big_unsigned_design.md). Signed tests and independent exact fixtures must remain warning-clean on supported C++20/ISA toolchains.

## Implementation and optimization timeline

### G0 — Branch creation and design baseline (18 October 2026)


Chose a **common unsigned arithmetic core with specialized higher-level implementations**, not multiple incompatible integer engines and not a general Rational allocation path for FloatMP. The first milestone is the living design only. No new numerical implementation or benchmark has yet been added. The FloatMP journal remains the performance reference and is to be updated in parallel when core changes affect FloatMP.

### G0 clarification — wrapper compatibility is not a goal (18 October 2026)


The only explicitly emphasized compatibility requirement is the handling of `Rational` NaN and ±infinity, in addition to exact mathematical correctness. Follow a staged caller migration rather than maintaining a permanent compatibility facade. Decide quotient/remainder semantics explicitly and update consumers as needed. Begin with isolated `BigUnsigned` and independent tests before replacing either type; FloatMP performance/correct rounding remains the key regression gate.

### G0b clarification — keep benchmarks out of CTest (18 October 2026)


**Current design rule (supersedes previous G0b smoke-test integration):** CTest and the custom `tests` target contain only functional/unit tests and their build prerequisites. Benchmarks are separate executables managed solely by `benchmarks/CMakeLists.txt` and the `benchmarks` target. Specifically, remove all four `smoke_benchmark_*` CTest entries, the special `benchmark_smoke_executables ALL` target, and the benchmark dependency of the `tests` target. Retain the four named benchmark executables and their comparison/absolute data paths.

The earlier G0b journal entries below record the historical smoke-test experiment and its CI workarounds; they are no longer the intended configuration. A benchmark can still fail by throwing on mismatched post-timing outputs **when invoked explicitly**. Such verification is part of the benchmarking program, not a unit test automatically registered with CTest. Performance runs use the same-host standard `1000x11` setup unless the user explicitly specifies otherwise. No numerical code, input pool, expected result or measured operation is modified by this separation. Revalidate the ordinary Unix, Windows and Coverage builds after removal.

### Implementation record template


Never claim CI/performance results before observing them.

### G0b implementation log — benchmarks introduced (18 October 2026)


Added four independent executables and shared workload/report helpers while leaving current Integer, Rational and FloatMP implementation untouched. Aggregate and absolute timings are separated; numerical agreement is verified after the timing loops. CMake registers all targets and four lightweight CTest smoke tests. The initial commit only introduces benchmark infrastructure: green CI, successful smoke tests and same-host 1000x11 numbers must be checked and recorded separately.

### G0b follow-up — ensure CTest benchmark executables are built by coverage (18 October 2026)


The first new benchmark commit `f0a60f78` compiled and passed the Release pipelines on GCC 13, Clang 16, AppleClang 21 and Windows. Both Coverage jobs, however, failed in their combined build/test phase, and the Unix ThreadSanitizer job failed at Test. The connector exposes failure *steps* but not full job logs, so the specific failure message is **not yet established**.

The Coverage target builds a project-defined tests target, which originally depended only on the unit-test executables. Four new CTest smoke tests invoke the new benchmark executables, so those executables must be built whenever the `tests` target is built. Add four explicit benchmark target dependencies to `tests`. This does not change any benchmark code or operation times, and avoids relying on an all-target build in Coverage. Re-run Unix/Windows/Coverage and distinguish any remaining sanitizer failures; no all-green claim is warranted until observed.

### G0b ThreadSanitizer: build smoke targets with default ALL (18 October 2026)


The supplied ThreadSanitizer CTest output from build `build-tsan` shows **all 19 numerical tests passing**, followed by four `smoke_benchmark_*` tests marked `Not Run` because their benchmark executables do not exist. This is definitively a missing-build-target failure, not a numerical comparison failure or detected data race. The workflow calls `cmake --build . --parallel 2` (default `ALL` target) before `ctest --output-on-failure`, whereas the previous commit only attached these executables to the *custom* `tests` target; the two build roots are different.

The solution is a single `benchmark_smoke_executables ALL` custom CMake target, depending on exactly the four registered benchmark binaries. The `tests` target also explicitly depends on it. Hence both TSAN's default build and Coverage's `tests`-based build produce all executables that CTest needs. This does **not** remove the smoke tests, skip TSAN, touch numerical code, change benchmark inputs or affect the Release performance workload. The preceding `24e4a290` commit is confirmed Coverage and Windows green, but Unix remains red from TSAN due to missing binaries. The current change must pass Unix/GCC/Clang/macOS, Windows and Coverage before claiming completion.

### G1 priority change — FloatMP before Integer/Rational (18 October 2026)


The previously recommended shared unsigned foundation and G1a `BigUnsigned` facade are therefore **not an immediate milestone**. The `FloatMPNative::Unsigned` type in `src/floatmp_native_unsigned.hpp` is deliberately private to FloatMP; its correctness and FloatMP performance are more important than premature reuse. The new `test_floatmp_native_unsigned` is a real functional CTest test, distinct from the four benchmark executables which must never be registered with CTest.

### Migration priority refinement — one FloatMP source pair (18 October 2026)


The user explicitly requires only `floatmp.hpp` and `floatmp.cpp` for FloatMP, not a separate header per arithmetic algorithm. All unsigned kernels and quotient/remainder rounding have been moved into `floatmp.cpp`, with no changes to the public header. Do not build a temporary compatibility layer or factor out code for Integer/Rational before measuring an actual need. Private implementation tests have been replaced with tests of the FloatMP public mathematical behavior, under CTest. Benchmark targets remain separate.

### B17 native signed-magnitude foundation, before Integer facade migration (19 October 2026)


B16 CI `3ea38f4` is green on Unix, Windows and Coverage. A one-commit replacement could break unrelated numeric semantics and leave test failures hard to localize.

B17 therefore introduces a separate `SignedInteger` native signed-magnitude class in `include/numeric/signed_integer.hpp` and `src/signed_integer.cpp`. It contains a canonical `BigUnsigned` magnitude and a sign flag, normalizing negative zero. Supported operations: signed/unsigned 64-bit construction (including INT64_MIN), strict decimal parse/output, sign, comparison, negate, add/sub/mul, square, exponentiation by squaring, unsigned-magnitude signed divrem, shifts and checked 64-bit extraction. Division is **truncation toward zero**, with remainder having the dividend's sign: `a=q*b+r` and `abs(r)<abs(b)`, including all sign combinations. Right shift truncates toward zero; `pow(0,0)=1`.

Full B17 cross-platform project CI remains pending.

### B17.1 decision: SignedInteger is temporary verification scaffolding (19 October 2026)


There must be **only one signed multiple-precision type in the final architecture: `Integer`**.

Check that no code or build reference to `SignedInteger` remains. This deletion is a required completion gate for the migration, not optional cleanup.

B17 remains an intermediate development and test stage; it does not establish a second production API. This documentation-only clarification does not alter B17 code or its CI results.

### B18 direct native Integer storage and arithmetic (19 October 2026)

**Implementation.** Migrate `Integer` itself to a canonical `BigUnsigned` magnitude plus sign bit. Implement signed construction, comparison, addition/subtraction, multiplication/square, quotient and remainder, powering, decimal parsing/printing, increment/decrement, absolute value, bit-length and checked numeric conversions using the shared exact core. Moving an Integer resets the moved-from sign and cannot introduce negative zero.

**Division contract.** Signed division truncates toward zero; the remainder has the dividend's sign. For every nonzero divisor, `a = b*q + r` and `abs(r) < abs(b)`. Checked extraction reports overflow instead of silently narrowing. Keep these observable rules explicit in code and tests.

**Tests and staging.** Retain fixed signed boundaries (including the 64-bit extrema), negative/zero normalization, invalid decimal inputs, all sign combinations, powers and self-aliasing, plus independently computed expected results. The earlier `SignedInteger` class remained temporary test scaffolding at B18; its eventual retirement is recorded separately. B17 Unix/Windows/Coverage succeeded at `31c81e1`; B18 itself required new cross-platform validation.

### B21 complete retirement of transient signed class (19 October 2026)


After successful B20 native ARM64 and x86-64 benchmark verification and green Unix/Windows/Coverage CI for `c767b45`, retire the independent B17 migration scaffold. The production `Integer` class continues to own canonical `BigUnsigned` magnitude and sign. Delete exactly `include/numeric/signed_integer.hpp`, `src/signed_integer.cpp`, and `tests/test_signed_integer.cpp`; remove its source from `src/CMakeLists.txt` and its test executable from `tests/CMakeLists.txt`. Do not replace the scaffold with another signed numeric wrapper.

Preserve B18's existing distinct 1700-case random differential corpus. Check the fixed sign/negative-zero, signed extrema, parsing failures, checked conversions, zero-power, negative square, signed 20-case division matrix, division-by-zero, and 130-shift pattern from B17. Do not add a new permanent signed class or public API.

The deletion is an architecture completion gate, not merely a build-time convenience. All production code continues to reside under `namespace Ariadne`. Full Unix/Windows/Coverage and CTest results for this migration are pending; do not assert complete validation until green. The benchmark remains unchanged (all six sizes, fourteen operations) and can be rerun without migration specific modifications.

### B33 alias-aware Integer self-multiplication dispatch to B32 square (19 October 2026)


Source-level inspection explains the latter: the benchmark's unchanged `IntegerOperations::square(a)` is defined as `a*a`; the inline `Integer::operator*` forwards both references directly to `Integer mul(Integer const&,Integer const&)`, which until B33 always called generic `BigUnsigned::mul`. The B32 symmetric `BigUnsigned::square` was not reached for that ordinary self-multiplication use case. This is a genuine production API opportunity, not a request to alter the benchmark workload.

B33 detects **object identity only** (`&z1==&z2`) in `Integer::mul`, and uses the already verified B32 `BigUnsigned::square` when the magnitude has >=8 limbs, only in GCC x86-64 builds with the existing `ARIADNE_BIG_UNSIGNED_GCC_FUSED` definition. Distinct-but-equal objects continue through generic `BigUnsigned::mul`; different values are unaffected. All smaller GCC inputs, Clang, MSVC, ARM64 Apple M2/M4 and portable configurations retain their previous binary-multiplication path. The square of a negative Integer is positive, so the alias case constructs an unsigned magnitude with sign=false; zero stays canonical. This also correctly implements `x*=x`: the original input is read completely before the assignment of the resulting Integer. The check is O(1), based on semantic aliasing of operand references, and does not recognize benchmark-generated values or inspect their numeric pattern. B32's existing full-product size guard and kernel dispatch remain authoritative.

A fresh stand-alone x86-64 GCC 14.2 isolated fixed-pool test of the exact B32 limb recurrences measured schoolbook -> symmetric at 8,16,32,64 limbs around 66.4->49.5 / 230.8->147.3 / 903.6->515.8 / 3625.6->1922.1 ns (approximately 25/36/43/47% less kernel time). Clang's isolated prototype also favored symmetry, but this does **not** change its production selector, which remains unmodified. These new production tests and an integrated B33 binary have **not yet been run** by the assistant. Thus no integrated x86-64 `square` gain may be claimed before repeating the user's 84-operation benchmark on the same x86-64 host, ideally with a GCC fused build and optionally its preexisting option-OFF control.

ARM64 B32 measurements were recorded on Apple **M2**, whereas preceding B31.1 values came from Apple **M4**; absolute cross-host differences must not be interpreted as B32 code regressions. No ARM64 code path is changed by B33. Neither CI nor Codecov was consulted at the user's request.

### B34 alias-dispatch investigation: targeted real-production benchmark (19 October 2026)


The B33 full x86-64 benchmark verified 1024/2048/4096-bit native square improvements of approximately 32/39/40% over B32. Its generic `mul` results at 53/128/256 bits, however, increased by approximately 19/13/10% even though B33 changed only a pointer identity guard in `Integer::mul` and retained the exact pre-B33 multiplication body for distinct operands.

Tested B33's current identity-first guard against size-first reordering, an explicit unlikely branch hint and both combined, using alternating-order GCC 14.2 x86-64 stand-alone probes with 64-bit inline capacity four, the exact fused schoolbook arithmetic and the B32 symmetric square recurrence. Measured native 1/2/4-limb distinct-operation times clustered around 20/23/44-47 ns; comparing a generic call without an alias guard to a call with the B33 guard generally yielded differences around 0-1 ns with variable sign. At 8-64 limbs no alternative had consistent performance wins across replicates. These isolated probes are *not* the complete production Integer implementation, so cannot conclusively explain the user's 5-6 ns small-`mul` shifts. They do show that blindly reordering the guard or attaching a branch hint lacks a credible performance justification. Binary layout or inter-run variability remain hypotheses, not established causes.

**Decision:** Preserve B33's current production `Integer::mul` and symmetric kernel rather than force an unproven B34 dispatch. Add a separate, opt-in `benchmark_integer_alias_dispatch` executable, built against the **same actual NUMERIC_SRC** as the official Integer benchmark. The new diagnostic records median ns for public `a*b` on distinct operands, distinct objects with identical values, and same-object `a*a`. Its references call `BigUnsigned::mul` via `Integer::from_magnitude`, and use the actual B32 `BigUnsigned::square` as a third self-square control. Report the library's actual `selected_multiply_uses_fused()` selection.

To run from the same build on x86-64: `./benchmarks/benchmark_integer_alias_dispatch --bits 53,128,256,512,1024,2048,4096 --iterations 20000 --samples 31`. Interpret public/reference differences cautiously because they include call structure as well as the alias guard; run B32 and B33 builds on the **same CPU**, alternating order if possible, for a causal performance conclusion. The C++ source and new CMake target were added but a full production build and this actual diagnostic have not yet been executed by the assistant. All Integer, BigUnsigned, Rational, FloatMP implementation paths, existing tests and official benchmarks remain untouched. B32 ARM64 numbers are Apple M2, B31.1 ARM64 numbers were Apple M4; never use cross-machine absolute times as evidence of a code regression. CI and Codecov remain unconsulted by user request.

### B34.1 diagnostic warning fix and ARM64/x86-64 results (19 October 2026)


The user built and executed B34's production-linked diagnostic on both Apple ARM64 (reported `fused=false`) and GCC x86-64 (reported `fused=true`), with 20,000 iterations, 31 samples. Clang on ARM64 reported `-Wswitch-default` because the diagnostic's enum switch did not contain an explicit `default` case. B34.1 moves the preexisting `std::logic_error("Unknown diagnostic mode")` throw into the `default` label of that switch, removing only the warning.

The x86-64 GCC results demonstrate B33's genuine benefit on aliased self-products by comparing its production `a*a` with the diagnostic's explicit generic multiply: 512 bits: 94.453 vs 102.012 ns (7.4% faster); 1024 bits: 181.105 vs 253.586 ns (28.6%); 2048 bits: 503.167 vs 809.378 ns (37.8%); 4096 bits: 1834.987 vs 3047.787 ns (39.8%). The direct `BigUnsigned::square` control at 4096 bits was 1827.186 ns and confirms that the alias path is reaching the specialized square rather than generic mul.

For distinct operands at 53,128,256 bits, the same x86-64 run yielded native/generic pairs 51.296/49.827, 42.475/40.296, and 69.793/67.834 ns, respectively, corresponding to gaps of about 1.5, 2.2 and 2.0 ns. Those differences include compiler code generation and call/dispatch structure, so **must not be treated as the pure alias-check cost**. The 5-6 ns B33/B32 small-multiplication regression has not been causally established. For 1024,2048,4096 bits, distinct native/generic were 257.940/257.924, 821.206/819.870 and 3063.295/3063.707 ns, indicating near parity for large inputs.

ARM64 native/generic measurements were very close, as expected for a platform that does not enable B33's GCC fused alias specialization (4096-bit distinct 2359.202/2355.177 ns, aliased 2362.954/2361.898 ns). These diagnostic numbers identify an ARM64 host only, not whether it was M2 or M4: never infer its specific CPU without explicit user confirmation. The earlier official ARM64 B32/B33 runs were Apple M2, whereas B31.1 was measured on Apple M4.

**Decision:** keep B33 unchanged, neither reorder its alias guard nor add a branch hint. The user has executed the B34 diagnostic, but the B34.1 warning-only commit has not undergone a new full build by the assistant. GitHub CI and Codecov were not consulted.

### B41 avoid temporary negation and trim common high limbs in subtraction (19 October 2026)


2.13x/2.35x/2.77x. The signed benchmark operand corpus has `a=+/- (2^(bits-1)+13+2*i)` and `b=+/- (2^(bits-1)+31+5*i)`, with 32 deterministic sign combinations. Half of the subtractions have opposite signs and therefore require a true n-limb **addition** of unsigned magnitudes. The other half subtract nearly identical n-limb magnitudes and produce a small signed result. Treating every timed `sub` as a long unsigned borrow chain would be technically incorrect.

Before B41, `Integer::sub(z1,z2)` invoked `add(z1,neg(z2))`. `neg` produces a new owning `Integer` and **copies the entire `BigUnsigned` magnitude**, even though the signed subtraction only requires reading the sign and magnitude of z2. B41 performs direct signed-magnitude dispatch without materializing this negated copy. If signs differ, add magnitudes and keep z1's sign; if signs match, compare magnitudes, return canonical zero when equal, otherwise subtract larger-minus- smaller and select z1's sign for larger z1, its inverse for larger z2. All four sign combinations, zero, x-x and x-=x preserve the exact original result and canonical sign.

Unsigned `BigUnsigned::sub(a,b)` formerly pre-scanned the full numbers with `compare`, allocated `a.size()` output limbs, subtracted all limbs using `sub_into`, then normalized the entire result. B41 first rejects a shorter minuend without allocation. For equal-width canonical inputs, it scans from the high limb to find the highest **different** position k. If no position differs, return zero directly; if a[k]<b[k], throw `underflow_error` without modifying an operand. Otherwise the high equal prefix cancels, and the entire nonnegative difference is represented by the low k+1 limbs. Allocate precisely that result width and call the existing exact `sub_into` on only those limbs, then normalize. For a strictly longer minuend there is no high equal-width prefix, and the full-width subtraction remains exact. This avoids duplicate whole-number comparison within `BigUnsigned::sub` and reduces work/allocations for all equal-prefix values, not only the benchmark pool. Full low-range borrow propagation, underflow behavior, canonical-zero and output-normalization invariants remain preserved. The lower-level public `sub_into` signature, bounds checks, output-tail clearing and borrow-out behavior are not modified. No new ISA, platform-specific selection, CMake flags, benchmark operands, library dependencies, scratch changes or API changes are introduced.

Note: the production tests and full CI are **pending** on the B41 commit. The same-M2 integrated benchmark must measure B41 versus B40, especially `sub` at 1024/2048/4096 bits, while using `add` and `negate` as controls. Avoid claiming a speedup before same-host timing data; inspect Windows/macOS/Linux CI separately from external Codecov.

### B41.1 Clang 16 compatibility of same-object subtraction regression (19 October 2026)


The native Apple M2 B41 full 31-sample official benchmark confirmed large speedups without changing the official workload. B40 -> B41 `sub` at 1024/2048/4096 bits improved 71.156 -> 27.744 ns (61.0%), 97.917 -> 44.462 ns (54.6%) and 159.254 -> 73.335 ns (54.0%).

B41's Unix CI failed during **compilation** on Ubuntu Clang 16 Release and ThreadSanitizer at `tests/test_integer.cpp:463`, not during an arithmetic test: `-Werror,-Wself-assign-overloaded` diagnosed the intentional `aliased-=aliased` regression statement. Both macOS AppleClang, Ubuntu GCC, GCC/AppleClang Coverage, and Windows workflows passed B41. This compiler-specific diagnostic does not by itself establish an arithmetic or sanitizer defect.

B41.1 changes **only** the intentional same-object test's spelling to: `Integer const& same_object=aliased;` followed by `aliased-=same_object;`. Both names still reference the exact same `Integer` object; the test retains its aliasing, canonical-zero and nonnegative-sign assertions. No source arithmetic, API, compiler warning settings, test scenarios, workload, oracle, scratch contracts, or other platform paths change. A local Clang 17 `-Wself-assign-overloaded -Werror` compilation of the same reference-alias idiom passed; full Clang 16 release/TSan and all CI remain pending until the B41.1 workflows finish.

### B64 short-operand pow/square ownership diagnostic (20 October 2026)


The overall worst is 2.476x (pow3, 128 bits). B63 reduced the 1024/2048 quotient to 176.383/306.194 ns, from B61 292.537/554.274 ns. The 53/128/256 short-size total-time gates are still far from <1x; the 4096-bit worst compare margin is small and could be sensitive to repeat-run noise.

The direct cube is NOT to be promoted solely because the official workload includes exponent 3. It also measures production alias z*z versus an owning BigUnsigned::square-based signed Integer construction, retaining the exact output semantics.

The two timing groups measure whole owned result objects, not non-additive kernel stages. Different implementations may allocate internally; this is precisely one of the effects under investigation and must not be claimed as a proven bottleneck until measured. The explicit generic lazy algorithm is a diagnostic candidate, not a B64 production change. No change is made to any production path, public API, existing tests or FloatMP; this new target is appended to benchmarks/CMakeLists.txt.

Same-host command, from build and with user's preferred theme and single-line formatting: ./benchmarks/benchmark_integer_small_power_profile -t dark --bits 53,128,256,1024,2048,4096 --iterations 20000 --samples 31 The result will determine whether deferred-copy generic exponentiation can reduce the 53/128-bit pow3 worst without the unsupported conclusion that one micro optimization will bring all three small-size totals under 1x. B64 integrated compilation/CI and actual GCC user-host timings are not yet confirmed.

### B65 assess ARM64 B37-B40 transfers at short widths and MSB-first power (20 October 2026)


Following the user's question about lessons from the earlier ARM64 square/power campaign, inspect both B37-B40 and x86 B32/B53/B54:
- ARM64 B37's triangular/symmetric square was already present on GCC x86-64 since B32, with the SAME >=8-limb production gate.
- B38's skip-zero outer rows and B39/B40's adaptive dense/sparse protection were already transferred to x86 by B53 and B54.
- The existing balanced/dense Karatsuba and zero-row dispatch serves large inputs; it does not address 1/2/4-limb work. There is no missing ARM64 assembly kernel that can be copied verbatim to x86, and the x86 production square under eight limbs still uses generic multiply.

The lazy 3-10% gain is much smaller than the 12-43% direct control's gain. A cube-only production shortcut would be a fragile and workload-specific response, not a general solution.

For exponent three, the generic recurrence naturally computes square(base) then mul(square,base), matching the direct cube control, but works for any unsigned exponent (including exponent zero and signed inputs). (2) A benchmark-only forced invocation of the existing GCC x86-64 symmetric-square kernel at 1..7 limbs, with a bounded 14-limb local output before constructing the same owning BigUnsigned. Production still dispatches to generic multiplication under eight limbs; ARM64, Clang, MSVC and other architectures remain unchanged. The forced comparison is reported at 53/128/256 bits; 1024+ retains its original selected square path and has no new short-kernel row.

MSB-first exponentiation changes the intermediate multiplication order at general exponents, so it must not be promoted from a cube-only measurement without wider correctness and same-host performance coverage and ARM64 regression checks. Likewise the symmetric short square should only be promoted if it improves the *whole owning* square result at all relevant small sizes, not merely its preallocated arithmetic. B65 native timings and integrated CI are pending. Command from build (one line; user-preferred -t dark): ./benchmarks/benchmark_integer_small_power_profile -t dark --bits 53,128,256,1024,2048,4096 --iterations 20000 --samples 31

### B66 promote independently verified MSB-first signed power on GCC x86-64 (20 October 2026)


The diagnostic direct-cube sequence was faster, but it is not a general justification to special-case exponent 3. The B65 short symmetric square control improved 53 bits by only 5.5%, changed 128 bits by -0.9% and regressed 256 bits by 0.2%; DO NOT lower the >=8-limb square production gate based on these data.

It handles exponents 0 and 1 explicitly, preserving 0^0=1, negative-sign parity and zero normalization. The exponent is the general Nat argument, not benchmark-only exponent three. The original B32 implementation remains selected on ARM64, AppleClang/Clang, MSVC and other builds.

Verify the full CI and ARM64 baseline before any wider algorithmic promotion.

### B67 own-result short-operand performance decomposition (20 October 2026)


Worst per-size: 53 pow3=2.264x; 128 pow3=2.447x; 256 square=2.117x; 1024 quotient=2.028x; 2048 quotient=1.939x; 4096 quotient=2.074x. Compared with B63, the MSB-first optimization did improve actual pow3: 53 72.001->68.925 ns; 1024 232.114->195.173 ns; 2048 342.477->306.423 ns; 4096 638.445->590.928 ns. But other unchanged operations, particularly construct_decimal, increased in absolute time: at 1024 444.540->587.105, 2048 986.184->1250.536 and 4096 2657.408->3226.152 ns. Quotient/remainder also moved without code changes. A single aggregate benchmark run is not sufficient to reject a consistently faster production pow kernel.

B67 adds an independent small-operations benchmark target using both original 32-input integer_specifications(bits) and integer_division_specifications(bits). Its arithmetic can take different branches because sign is absent, so the Integer/core difference is not a direct measurement of signed-wrapper overhead. It can nevertheless help distinguish a slow primitive from an especially costly signed path on 53/128/256 bits. For pow3 the core controls direct square(a)*a and is explicitly not identical to the whole exponentiation algorithm; B65 already established this is a useful lower-bound diagnostic, not a generic cube production shortcut.

B67 build/CI and user-host measurements remain to be verified. Command from build, one line with user's preferred theme: ./benchmarks/benchmark_integer_small_operations_profile -t dark --bits 53,128,256,1024,2048,4096 --iterations 20000 --samples 31

### B68 same-signed-semantics A/B for decimal adapters and short products (20 October 2026)


B67 also indicates integer/core mul overhead converges to 1.092x at 4096 bits; to_decimal converges to 1.012x. Do not infer a proven single root cause or alter the official benchmark from this comparison.

B68 extends the unchanged independent B67 benchmark with an additional A/B section whose controls match the complete SIGNED output type and mathematical semantics, on the same exact 32 input pairs and same iteration/sample settings:
- construct_signed: normal Integer(String(decimal)) versus independently parsed sign with production from_decimal and Integer::from_magnitude, avoiding the intermediate Ariadne String for the candidate.
- to_decimal_signed: normal std::string(Integer::literal()) versus production BigUnsigned::to_decimal on the input's magnitude plus an exact conditional leading '-' inserted into the returned owning std::string; the sign, including canonical zero, is preserved.
- mul_signed: normal Integer operator* versus a signed result constructed from production BigUnsigned::mul and the XOR of the two input signs.
- square_signed: normal self-alias Integer z*z versus signed (nonnegative) construction from production BigUnsigned::square.

The original B67 report with unsigned magnitude controls is still printed intact. These measurements can distinguish expensive String/signed adaptation from primitive costs more reliably than B67's unsigned comparisons, but remain independently timed and susceptible to scheduler/clock variance; avoid attributing their differences additively.

No new production behavior or benchmark weights are introduced. B68 integrated CI and user-host timings are not yet established.

Run from build, single line with dark theme: ./benchmarks/benchmark_integer_small_operations_profile -t dark --bits 53,128,256,1024,2048,4096 --iterations 20000 --samples 31

### B70 short multiplication: selected fused dispatch, ownership and inlining controls (20 October 2026)


Full owning signed forced-carry/forced-fused was likewise worse in all cases, e.g. 128-bit mul 32.279/29.876=1.080x. Keep GCC x86-64 fused production selected, and keep the ARM64 native carry kernel and other ISA dispatch unchanged.

Since production already selects the SAME fused arithmetic kernel, these 19-21% reductions CANNOT be attributed solely to the multiplication instruction sequence.

Recheck the original oracle after all timed rows.

The noinline attribute is a compiler optimization hint for this benchmark and does NOT exactly reconstruct cross-translation-unit ABI/LTO behavior. All rows are independently timed and not additive, and the production results may be affected by scheduler/clock drift. No performance gains or new production acceptance gates are claimed until same-host B70 data are returned.

CI/integration verification pending.

### B73 remove temporary decimal ownership and one signed-result move (20 October 2026)


The short groups require approximately 30.0%/32.1%/31.4% reductions in their native sums, not merely a worst-operation fix. The attached B68-B71 profiles predate production B72. Keep these separate from post-B72 results. B69's carry kernel remains slower than fused; B70/B71 provide no basis for changing the selected arithmetic or assigning additive costs to sign handling, inlining or allocation.

B73 changes ownership and temporary storage, with no new ISA kernel:

- `Integer::literal()` moves the existing `to_decimal()` buffer into the standard-string base of its returned `String`, then inserts the sign in that same owning string. The old conditional expression and `String(std::string const&)` copied the magnitude buffer. This retains the exact existing `String` return type and requires no utility submodule change. Negative zero still formats as `0`.
- The private signed-magnitude constructor takes `BigUnsigned&&`. The public `from_magnitude(BigUnsigned,bool)` signature and layout stay the same; its already-owned parameter is moved into the member once, instead of through another by-value parameter first. Lvalue callers still copy their source; moved-from magnitudes remain zero.
- `from_decimal()` uses one generic block recurrence for both storage strategies. At most 77 input digits are accumulated directly into the result's four inline limbs: 10^77 < 2^256 bounds every prefix. Longer inputs accumulate in the result's own vector, reserved to ceil(input_digits/19), which is safe because 10^19 < 2^64. This avoids the final second allocation/copy. Long inputs are not zero-initialized across their reserved capacity. Leading-zero inputs producing at most four limbs are copied into inline storage and release their unused vector allocation. The existing 32-digit threshold, 19-digit wide recurrence, portable 9-digit recurrence and full ASCII validation before parsing remain unchanged. Appends occur only after each limb scan, so vector growth/span extension cannot invalidate an active scan. The carry append invariant keeps stored limbs canonical, including an empty magnitude for zero.
- The radix-10^9 formatter uses local arrays for up to eight source limbs and 24 decimal blocks, eliminating its two temporary vector allocations at short widths. Each limb contributes fewer than three 9-digit blocks. Larger portable leaves retain checked dynamic storage. Wide/recursive decimal algorithms and their dispatches, multiplication/division kernels, fused/carry selection and B72's product result construction are unchanged.

No public operation, signed division convention, object layout, benchmark fixture, official input pool, operation weight or oracle is changed. The private constructor's signature changes; no binary compatibility for that private symbol is promised.

Validation on Apple M2 Ultra, macOS 26.6.2, AppleClang 21 Release with warnings-as-errors: all 20 CTest tests pass. The complete native `test_integer` and `test_big_unsigned` suites also pass a Debug build with AddressSanitizer and UndefinedBehaviorSanitizer. The complete BigUnsigned suite compiles with GCC 16.2, warnings-as-errors, native ARM64 kernels and passes. Its compilation suppresses the pre-existing unused fused-helper warning in this artificially forced configuration; the normal native builds retain warnings-as-errors. This is not a native MSVC build. The unchanged 1000x11 four-size FloatMP benchmark also verifies all 256 paired inputs and 3456 operation results; no FloatMP performance comparison to B72 was made. Windows/MSVC and native GCC x86-64 CI remain pending.

Native short conversion nanoseconds (B72 -> B73 two-run range): 53-bit parse 49.110 -> 22.085..22.090, format 74.612 -> 36.923..37.090; 128-bit parse 97.496 -> 62.460..63.244, format 223.083 -> 117.525..119.360; 256-bit parse 127.144 -> 93.123..93.796, format 285.700 -> 170.333..171.842. These are ARM64 observations only; they do not project x86-64 timings or isolate each change's causal contribution. The baseline was a single run, candidate two repeats; host drift exists. Parsing at 2048 bits was about 7% slower locally (1010.794 -> 1080.692/1084.123 ns) despite lower per-size totals. Do not describe this as a per-operation regression-free change.

**The requested target has not been achieved or certified by B73.** Even the ARM64 2048-bit total remains above 1x. A sum below 1x across all sizes is insufficient. Further work should follow those post-B73 measurements: short owning square/pow and nontrivial division remain candidates, with no kernel/threshold promotion inferred from the old profiles alone.

### B75 expose small ownership/adaptation bodies to the optimizer (20 October 2026)


The user's B74 GCC x86-64 official 20000x31 run verifies all 768 paired inputs and 2688 operation results. The remaining 44.529/83.410 ns require about 8.36%/10.74% native reductions. The overall 0.901x sum does not substitute for these two gates.

B75 makes existing small bodies visible in the class headers: BigUnsigned's scalar constructor, move constructor/assignment, limbs() and limb(); Integer's signed-magnitude constructor, from_magnitude(), destructor, copy constructor/assignment and move constructor/assignment. Their arithmetic, field transfers, self-move guards, moved-from canonical zero and negative-zero normalization are unchanged. The former out-of-line definitions are removed; there is one ODR-valid inline definition per method. Integer's default constructor remains out of line. No arithmetic kernel, unsigned power/division schedule, storage layout, capacity, API signature, exception specification, ISA dispatch, benchmark fixture, keep-alive barrier, operation weight or numerical oracle is changed. No LTO, compiler flag, always-inline attribute or platform dependency.

The purpose is to let each caller optimize ownership transfers and eliminate redundant call/temporary work, including destruction of moved-from temporary magnitudes. It is not proof that a specific GCC ABI/linker rule caused the prior cost: B70's noisy controls did not establish that. An inline definition permits optimization; it does not guarantee a particular machine-code sequence or native x86-64 improvement. The bodies are shared across all toolchains.

**Rebuild the library AND its clients/benchmarks together.** Moving formerly out-of-line special members to inline definitions can change which standalone symbols the shared library emits. This iteration preserves source-level mathematical/lifetime contracts and object layout; it does not promise linking old precompiled clients against the new library. The benchmark must be rebuilt normally, not altered to bypass construction/destruction.

Local Apple M2 Ultra, macOS 26.6.2, AppleClang 21 Release warnings-as-errors: the complete C++ `everything` target rebuilds, and 20/20 CTest tests pass. Complete Integer/BigUnsigned Debug ASan/UBSan suites pass. The complete BigUnsigned suite also passes GCC 16.2 with native ARM64 kernels; actual storage/boundary/power/ short-division tests pass forced-portable ASan/UBSan without native kernels or __int128 (only its pre-existing unused fused-helper warning is suppressed). No tests or assertions were changed. The Integer translation unit also compiles with GCC 16.2 using the project warning set (including its normal non-template-friend suppression) and warnings-as-errors. These are ARM64/portable checks, not native GCC x86-64 or MSVC CI.

The unchanged 20000x31 ARM64 Integer benchmark verifies 768 paired inputs and 2688 results before and in both final candidate repeats. A saved B74 executable AND library form the reference, while all candidate clients are rebuilt with the candidate headers. Native 128-bit sums: 354.681 -> 323.478/335.364 ns; 256-bit sums: 564.974 -> 527.437/550.217 ns. These local reductions do not predict the required x86-64 improvement or attribute additive costs.

Thus local ARM64 4096 total acceptance is not established; neither a blanket nonregression nor target parity is claimed from these mixed runs.

Existing short-width FloatMP gates remain unmet; these one-run measurements do not establish performance acceptance. No FloatMP arithmetic or rounding algorithm was changed.

### B84 — version exact structural Integer/BigUnsigned generators (23 October 2026)


Add the self-contained C++20 `tests/data/generator_integer_structural.cpp` as a seventh standalone generator. It independently derives the original source operands and exact expected outputs for `integer_b66.tsv` (360 rows), `integer_sub.tsv` (192 rows), `big_unsigned_sub.tsv` (100 rows), `integer_b33.tsv` (72 rows), `integer_b37.tsv` (72 rows) and `integer_b38.tsv` (30 rows). Replay the entire original mt19937_64 seeded streams, including the visits discarded by the historical structural sampling. Preserve exponent boundaries, all signed masks, alias cases, underflow markers, and canonical zero.

GCC and Clang compiled the same generator with C++20 `-O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror`; both wrote identical byte streams. Git blob object hashes of the regenerated TSVs match the previously committed originals exactly: integer_b66 9c0a5bfc956b78aead0e7a1b13d93a7e57f896e2; integer_sub b1ba7f3454dbc190755e8c7b66bc05a0918a6fdf; big_unsigned_sub f57078e92fd83ba82fef41aceb9bd148604da08d; integer_b33 ff92f80126f7e409904d933ff2eee065626ba557; integer_b37 c1a2861fc851a2c048f73dd1e1817f2ecdcdc9e3; integer_b38 cf359723cc0f5d70984fab635ea03a49ec32442b.

Register the program in the `ORACLE_GENERATORS` list so explicit `generators` compiles it and explicit `generate-data` runs it with `--write`; unit-test builds remain independent. Update generator README. 50 of the 115 versioned TSVs are now source-reproducible. No fixture contents changed. The complete platform CI remains to be checked.

### B87 — restore source-keyed signed randomized B17/B18 generators (23 October 2026)


Add `tests/data/generator_integer_random.cpp`, standalone standard-C++20 source replaying all 2500 B17 and 1700 B18 original mt19937_64 visits, concatenated decimal string operands, sign rules and zero cases. The historical one-per-ten sampling positions are preserved exactly: `rep=10*k+(7*k+3)%10`, yielding 250 B17 and 170 B18 cases. The independent radix-2^64 limb implementation computes signed add/subtract, shift/add product, bit-aligned long quotient/remainder with truncation toward zero and sign-of-dividend remainder, and B17 signed right/left shifts including all 0..129 bit boundaries.

A separate exact JavaScript BigInt/MT19937-64 implementation reproduced both original TSV contents independently, followed by direct C++20 GCC and Clang execution and byte-for-byte comparison. C++ builds used `-Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror`. Original Git blob object IDs matched: B17 `542d28c6079f163c8d21f36adc35f5f57ae2d91f`; B18 `28f074d21007fe6412004a8d3a0432cb6478d20d`. Register the tenth executable in ORACLE_GENERATORS for `generators` and `generate-data`, not `tests`; 59 of 115 TSVs now reproduced. Full platform CI remains to be checked.

