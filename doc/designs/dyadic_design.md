# Dyadic — native implementation journal

**Repository:** `ariadne-cps/numeric` · **branch:** `main`
**Scope:** exact binary-rational `Dyadic`, its shared integer substrate, floating-point conversions, and validated enclosure wrappers.
**Companion core journal:** [Integer and BigUnsigned](integer_design.md) · [Rational](rational_design.md)

This journal records the evolution from the early `NewDyadic` prototype to the canonical `Dyadic` implementation. Historical revision labels are preserved, but each record focuses on algorithms, semantic behavior, independent tests, benchmarks of the native implementation, and outstanding CI requirements. Earlier experimental class names do not imply that those types remain active on `main`.

## Current design and invariants

- Finite values are exact signed integer coefficients multiplied by powers of two, with canonical normalization. Special cases include NaN and signed infinities.
- Addition and subtraction align powers of two with virtual shifts and owned-output kernels; comparisons avoid unnecessary large temporaries. Multiplication, squaring, halving, powers, floor and ceil must remain exact.
- `get_d()` follows the documented binary64 conversion contract, including normal and subnormal boundaries and overflow/underflow. `FloatMP` conversions preserve exact dyadic values where specified and distinguish exact from explicitly rounded construction.
- Bounds, balls, lower/upper bounds and approximations follow the existing public header layout and preserve enclosure semantics. Interval tests cover signs, zero crossing, powers and exceptional values.
- Verify integer/exponent invariants, independent exact values, API integration, platform builds and sanitizers. For performance, use native absolute timings at fixed widths and repeated runs, not unrelated baselines.

## Implementation and optimization timeline

### R33a — introduce experimental native NewDyadic (21 October 2026)


Commit `0f7f06c16` introduces the initial native `NewDyadic` implementation. Subsequent entries record API parity, dedicated tests, benchmark design, platform issues, and performance changes.

### R33d — align the native arithmetic API (21 October 2026)

Complete the prototype's public signatures: Boolean comparison results, `cmp`, `sgn`, `eq`, `lt`, special-value predicates, `Positive<NewDyadic>` results for `sqr` and `abs`, named and compound arithmetic, signed/unsigned power overloads, coefficient/exponent accessors, and exact construction from String, ExactDouble and TwoExp. Port semantic tests by functionality and verify the exact Rational representation outside any timed operation.

### R33e — extend the exact Dyadic regression suite (21 October 2026)

Organize the prototype tests into literals, conversions, arithmetic, comparisons and exceptional values. Exercise exact dyadic literals 3.25, -11.375, 0.375 and powers of two, rejection of invalid/non-dyadic decimals, rounding cases, negative floor/ceil, `3/4 ± (-5/8)`, powers and compound assignment, NaN, infinities and invalid operations. Test widths 53, 128, 256, 512, 1024, 2048 and 4096 bits with independent exact Rational checks and compile-time verification of the Boolean and Comparison result types.

Bounds, balls, writer operations and complete stream coverage were explicitly deferred at this stage; the later D-series entries trace their integration.

### R33f — AppleClang constructor compilation repair (21 October 2026)

Fix the incomplete `ExactDouble` type in `NewDyadic(double)` by including `numeric/builtin.hpp`. Constrain the deleted signed-exponent overload with `BuiltinSignedIntegral` so unsigned `size_t` arguments select the correct constructor. Signed exponent arguments remain intentionally rejected. Prior Unix/Coverage compilation failed; compilation must be checked again after this change.

### R33g — benchmark CLI and report setup (21 October 2026)

Use the shared `benchmark_command_line.hpp` parser, `Ariadne::CommandLineInterface` and logging conventions. Support iterations, samples, widths, summary/detail modes and positional iteration/sample compatibility. Use the deterministic finite/special corpus and check results with exact arithmetic outside the timed region. The experimental default widths were 53, 128, 256, 1024, 2048 and 4096 bits, with 20,000 iterations and 31 samples. Benchmark execution was not registered as an automatic CI test.

### R34 — direct binary64 conversion and virtual-shift comparison (21 October 2026)

Replace an allocating `Rational` conversion in `get_d()` with direct exact extraction of the 53-bit significand from at most two magnitude limbs. Account explicitly for finite signs, exponent overflow, normal/subnormal spacing and underflow; the chosen tiny-negative conversion convention yields positive zero. Replace the allocating aligned-integer comparison with radix-`2^64` virtual-limb alignment and early high-exponent/sign decisions.

Add boundary tests for normal/subnormal/overflow conversion and exact signed comparisons across multiple widths. Cross-platform build and repeated native latency validation remain separate requirements.

### R35 — stack-based virtual-shift add/sub (21 October 2026)

For finite operands of 5–64 limbs with exponent difference at most 63 bits, combine shifted radix-`2^64` limbs directly in one 66-limb stack accumulator. Compare virtual magnitudes before opposite-sign subtraction, propagate carry/borrow exactly, and materialize only the final canonical signed integer. Preserve zero/cancellation, aliasing and special-value semantics; outside the bounded regime retain the exact general implementation.

An independent Boost.Multiprecision C++20 check covered 200,000 signed operand pairs, including cancellation. It was a standalone diagnostic rather than an integrated project CI result.

### R36 — architecture-specific add/sub crossover (21 October 2026)

The R35 stack accumulator reduced native ARM64 add/sub latencies at larger tested widths, but regressed at 2048 and 4096 bits on the measured x86-64 system. Restrict the bounded stack algorithm to magnitudes of at most 64 limbs on ARM64 and 16 limbs on x86-64/other targets, with exponent gap at most 63; otherwise use the owning-shift path. The architecture cutoff is a measured heuristic, not part of the arithmetic semantics.

R36's machine-specific timing observations require confirmation after every change to compiler or dispatch. Keep all finite/special correctness cases unchanged.

### R37 — constant-time fractional detection in floor/ceil (21 October 2026)

For a normalized nonzero finite coefficient `c / 2^e`, the coefficient is odd whenever `e>0`. Therefore the number has a fractional part exactly when `e != 0`. Replace full-width shift-and-compare detection with this constant-time invariant in `floor` and `ceil`, retaining sign corrections and exceptional-value errors. Verify `floor(x) <= ceil(x)` and equality exactly for integral normalized values. The R36 Unix, Windows and Coverage workflows were reported green; a new performance benefit from R37 was not yet established.

### R37.1 — compile fix for mixed Boolean/bool test comparison (21 October 2026)

 All R37 Unix/Windows/Coverage workflows failed at the compilation step. Explicitly convert the Ariadne Boolean to bool before comparing to the canonicality expectation. Keep the assertion, original Dyadic oracle tests, all R37 arithmetic, benchmark fixture, build workflows and API unchanged. Both R37 user-supplied native benchmark runs passed their 384-input fixture numerical verification and showed floor/ceil improvements; these are independent of integrated CI. Revalidate Unix/Windows/Coverage before considering R37 stabilized.

### R38 — single-owned-allocation shifted add/sub (21 October 2026)

Introduce `BigUnsigned::add_shifted` and `sub_shifted` for virtual left shifts from 0 to 63 bits. Read borrowed operand views without mutation, propagate carry/borrow through a single radix-`2^64` loop, and write into one owned result. Use exact signed scaled-magnitude comparison to select direction and sign. For larger finite magnitudes, select the owned kernel at 32 limbs on ARM64 or 17 limbs on x86-64/other targets; shorter inputs retain the bounded stack path.

Test widths and transitions near 16/17 and 31/32 limbs, 4096/8192-bit inputs, signed cancellation and exponent gaps 63/64. No success is inferred until compilation, exact tests and native measurements complete.

### R38.1 — fix AppleClang macro parsing in BigUnsigned shifted-subtraction test (21 October 2026)


Commit `8f7ea1992` wraps the lambda expression passed to `ARIADNE_TEST_ASSERT` in an extra pair of parentheses. The comma separating the two `BigUnsigned` local declarations inside the lambda had been interpreted by the preprocessor as an additional function-like macro argument under AppleClang. The added parentheses group the expression as one macro argument; the underflow test semantics and checked expectations are unchanged. CI rerun required before marking the fix validated.

### R39 — avoid redundant temporary copy in `abs` and `sqr` (21 October 2026)


Add a checked rvalue overload `Positive<NewDyadic>(NewDyadic&&)` next to the existing const-reference constructor. This lets the negative `abs` path consume the already materialized `-x` temporary without copying its potentially large native Integer coefficient again; the positive `abs` path still copies the caller-owned lvalue, preserving nonmutation, and `sqr` can likewise move its multiplication result. Validation of the positivity invariant is performed on the constructed object. Keep the public signature of `abs`, exact sign/special semantics, original Dyadic and benchmark workload unchanged. This is a structural allocation/copy reduction candidate, not a measured ARM64 gain. Run compilation and full regression tests plus native same-host ARM64 benchmark before accepting a performance improvement.

### R40 — zero-shift direct add/sub path for large NewDyadic on x86-64 (21 October 2026)


When both finite operand denominator exponents coincide, R38's selected owning add/sub path has two zero shifts. Use the existing exact `BigUnsigned` `add`/`sub` in that case, in either subtraction direction selected by `compare_scaled`, instead of the more general virtual-shift combiner. Retain all previous signed-order and zero/cancellation logic. When the denominator exponents differ, keep R38's shifted kernel and all original crossover thresholds. This is an exact fast-path candidate particularly relevant to the x86-64 large-operand path, not a proven speedup; no ISA-specific instruction set is required. Run full correctness and the unchanged paired 20,000x31 benchmark on both hosts, including nonzero shifts and boundary cases, before making performance claims. Do not change benchmark corpus, public Dyadic or Rational.

### R40.1 — repair test oracle value access (21 October 2026)

Fix four incorrect accessor expressions in the shifted-arithmetic tests: the fixture values are data members rather than callable functions. This is a test-compilation repair only; the numerical test cases, arithmetic kernels and benchmark inputs are unchanged. Confirm Unix, Windows and Coverage after the edit.

### D41 — streamed first-operand shift on x86-64 (21 October 2026)

For `shift_a` in 1–63 and `shift_b == 0`, keep a forward spill register for the high bits of the first operand instead of rereading earlier limbs. Each loop iteration reads that operand's current limb at most once; the second operand is unshifted. Preserve exact carry/borrow, one owned destination and input immutability. Other shift arrangements and ARM64 retain the established kernel.

A standalone C++20 diagnostic on GCC/Clang x86-64 checked widths 0–129 and shifts 0, 1, 7, 31 and 63, including underflow. The second-operand-shift candidate was not selected because some widths regressed. This diagnostic does not establish integrated Dyadic performance; project CI and full native timing still need verification.

### D42 — direct canonical absolute value (21 October 2026)

Construct the nonnegative `Integer` coefficient from one copy of the already canonical unsigned magnitude and transfer the original denominator exponent/kind directly into the result. Map negative infinity to positive infinity, leave sources unmodified and use the checked rvalue `Positive` constructor to avoid a second owning copy. Cover zero, NaN, infinities, both signs and coefficient widths up to 8192 bits, checking the exact Rational value and canonical exponent. Compilation and timed native measurements remain distinct acceptance checks.

### D43 — standalone native Dyadic benchmark (21 October 2026)

Introduce `benchmarks/benchmark_dyadic.cpp` and its CMake target. Measure 18 construction and arithmetic operations at the default widths 53, 128, 256, 1024, 2048 and 4096 bits using 20,000 iterations and 31-sample medians. Report **absolute nanoseconds** and per-width native aggregate times. Exact Rational checks run outside timed callbacks.

This benchmark becomes the independent native performance log. The public class rename and downstream API integration are tracked separately in D44–D56.

### D44 — FloatDP input uses exact Dyadic conversion (21 October 2026)

Route `FloatDP` assignment from a dyadic value through the type's public `get_d()` conversion. Keep the previous conversion semantics and validate compilation; this is an integration step, not evidence of complete Dyadic migration.

### D45 — define precise FloatMP interoperability (21 October 2026)

The `FloatDP` input adjustment is retained. `FloatMP → Dyadic` must preserve an exact finite binary fraction (integer significand and signed binary exponent). `Dyadic → FloatMP` distinguishes exactly representable construction from explicitly rounded construction and respects the requested precision. Keep independent tests for finite and exceptional values.

### D46 — exact FloatMP-to-Dyadic export (21 October 2026)

Read the multiprecision float's exact integer significand and binary exponent. For nonnegative exponents, shift the magnitude left exactly; for negative exponents, form a Dyadic coefficient with a nonnegative denominator exponent. Preserve sign, canonical zero, NaN and signed infinities, and check that exponent conversions cannot silently narrow. Add exact round-trip and Rational tests for 1, 53, 127, 128 and 191-bit significands, five exponent choices and both signs.

### D47 — native bidirectional FloatMP/Dyadic conversion (21 October 2026)

Add direct conversions between `FloatMP` and the native prototype. Import uses an exact Rational representation and distinguishes exact construction (rejecting inexact target precision) from explicitly rounded conversion. Export uses exact integer significand and binary exponent, including NaN and infinities. Test assignment, exact/rounded conversion, finite values of 1–191 bits and special values. The temporary prototype overloads were later consolidated into the canonical public type.

### D48 — stage native interval and enclosure wrappers (21 October 2026)


Extract the existing Dyadic bounds/ball/lower-bound/upper-bound/ approximation type declarations into an independent `numeric/new_dyadic_extensions.hpp` for the experimental native type, preserving value/raw access, lattice operations and exact enclosure semantics. This is intentionally a parallel API, not a rename of the production Dyadic type. The `operator ValidatedNumber()` conversion is also deferred until the final name migration because its NumberWrapper/ExactNumber instantiation requires replacing the production class.

Add focused tests in `test_new_dyadic.cpp` for arithmetic enclosures, crossing-zero square and powers, subset refinement, ball conversion, and lower/upper/approximation wrappers.

### D49 — revert D48's separate extensions header (21 October 2026)


D48 introduced `include/numeric/new_dyadic_extensions.hpp`, a parallel wrapper header with no counterpart in the existing Dyadic organization. This was an architectural deviation and is reverted in full: delete the new header and restore `tests/test_new_dyadic.cpp` to the D47 version. Preserve the D48 entry above for historical traceability, but treat its implementation as **retracted**.

The replacement must use the existing file layout exclusively: `include/numeric/dyadic.hpp` for the native Dyadic class, Positive specialization, writers, and bounds/ball/lower/upper/approximation declarations; `src/dyadic.cpp` for required out-of-line definitions; `tests/test_dyadic.cpp` for the retained and updated tests. Delete `include/numeric/new_dyadic.hpp` and `tests/test_new_dyadic.cpp` at final promotion. Existing mathematical coverage must be retained through the final renamed test.

### D50 — integrate native enclosure wrappers into the existing Dyadic header layout (21 October 2026)


D48 was reverted in D49 because a separate `new_dyadic_extensions.hpp` file was not part of the project's Dyadic layout. Restore its already CI-tested native wrapper implementation in the **same native class header**, currently `include/numeric/new_dyadic.hpp`; this is to be renamed in-place to `include/numeric/dyadic.hpp` during final promotion. No additional headers or new source files are created. Retain the D48 interval/ball/lower/upper/approximation tests in the existing native test file, but remove its include of the rejected extensions header.

A passing D48 CI validates only the former split-header variant; Unix/Windows/Coverage must validate the new in-header arrangement.

### D51 — exact conversion adapter consolidation (21 October 2026)

Consolidate the exact Dyadic conversion route through `Rational`: reconstruct the binary coefficient and denominator exponent, with explicit NaN and signed-infinity handling. This integration step requires finite and exceptional round-trip validation; it does not alter the Dyadic arithmetic kernels.

### D52 — parity of Rational/Integer/Dyadic interoperability with FloatMP (21 October 2026)


FloatMP already provides `cast_integer(FloatMP)`, `FloatMP(Integer,...)` and both exact and explicitly rounded Dyadic conversion APIs; no duplicate Integer overload is added. Add independent tests of finite/large dyadic rational values, both infinities, NaN, FloatMP->Dyadic roundtrip, and exact/non-integral FloatMP->Integer casts. This is an API parity step, not final Dyadic promotion.

### D53 — infinity sign handling in conversion paths (21 October 2026)

Restore sign classification through an exact numeric comparison against zero when producing infinity, preserving positive/negative infinity and NaN semantics. No arithmetic or rounding algorithm is changed. Validate the conversion path with the integration tests.

### D54 — promote native Dyadic to canonical files (21 October 2026)

Replace the public `include/numeric/dyadic.hpp` and `src/dyadic.cpp` implementation with the canonical Integer-magnitude / denominator-exponent representation. Keep bounds, balls, approximation wrappers and writers in the project's established Dyadic header/source organization. Retire the temporary prototype header and test files; update `FloatMP` and the standalone benchmark to use the canonical `Dyadic` type. Preserve the existing regression scenarios and the new wrapper tests. Publication of this integration change does not imply green CI.

### D55 — restore Dyadic public integration after canonical migration (21 October 2026)


Restore missing rational-valued reciprocal/division contracts for Dyadic and explicit mixed RationalBounds/DyadicBounds operations and comparisons required by NumberWrapper. The exact interval arithmetic is unchanged. No test fixtures suppressed. These are build fixes and still require native compilation/CI validation.

### D56 — public API and regression-test integration repair (21 October 2026)

Restore `Dyadic::inf()` for Python alongside `inf(int)` and `inf(Sign)`. Preserve exact Rational-valued `Dyadic / Dyadic`, `RationalBounds` division of dyadic bounds, and mixed exact comparisons between Rational, Dyadic and Integer in either order; NaN comparisons return `Comparison::INCOMPARABLE`. Convert 13 wrapper checks to `ARIADNE_TEST_ASSERT` and verify canonical integer/exponent invariants. Keep `dyadic.hpp`, `dyadic.cpp` and `test_dyadic.cpp` as the canonical file organization. Full compilation and regression remain the acceptance gate.
