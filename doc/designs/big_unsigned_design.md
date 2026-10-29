# BigUnsigned — native unsigned multiprecision design journal

**Repository:** `ariadne-cps/numeric` · **branch:** `main`
**Scope:** `BigUnsigned` unsigned limb representation, native add/sub/multiply/square/divide/GCD/power, scratch and caller-owned buffers, architecture-specific kernels, decimal I/O, source-keyed oracle generation, and historical A/B-series optimizations.
**Companion journals:** [Integer](integer_design.md) · [Rational](rational_design.md) · [Dyadic](dyadic_design.md) · [FloatMP](floatmp_design.md)

This is the living design and correctness journal for the **shared unsigned substrate** used by `Integer`, `Rational` and other numeric backends. Original milestone titles, dates, proof obligations, performance evidence and verification caveats are retained. Work on signed wrappers or signed benchmark semantics is cross-referenced in [`integer_design.md`](integer_design.md) rather than duplicated. Many historical B-series entries originated in the joint development of both components; the distinction here follows implementation ownership, not the numeric milestone prefix.

## Current architecture and correctness contract

Unsigned limb storage and normalized views use exact C++20 arithmetic with canonical zero, consistent carry/borrow and aliasing, and portable fallbacks for selected ARM64 and x86-64 specializations. `BigUnsigned` provides low-level multiplication, squaring, Euclidean division, GCD, power, decimal conversion, scratch storage and caller-provided output. The algorithm-specific benchmarks and independently generated golden fixtures remain validation gates. Signed quotient and sign rules belong to [`Integer`](integer_design.md); finite/extended numerator-denominator and binary64 conversion rules belong to [`Rational`](rational_design.md).

## Implementation and optimization timeline

### BigUnsigned phase B0 — portable storage and foundation (18 October 2026)


`FloatMP`, `Integer` and `Rational` retain their current production representations; no operations have been rewired to use BigUnsigned.

**Storage contract.** Nonnegative exact values are stored as a little-endian sequence of 64-bit unsigned limbs; each represents a coefficient of `(2^64)^i`. Zero is the empty sequence. The highest stored limb is always nonzero for nonzero values. Four limbs (256 bits) fit in the value's inline buffer; larger values use dynamically sized contiguous storage. `size()` returns the number of significant limbs, `limb(i)` returns zero beyond that length, `limbs()` exposes a read-only, invalidated-on-mutation span. Constructor normalization strips upper zero limbs. All normalizing mutations switch back to inline storage if at most four limbs remain. Copy owns independent storage; move leaves the source as zero, and `swap` swaps all representation fields.

**Portable operations now implemented:** `compare` returns -1/0/+1; `add` and `sub` compute exact nonnegative magnitudes with carry and borrow limited to 0 or 1; `sub` throws `std::underflow_error` if the result would be negative. `shift_left` computes exact multiplication by `2^k` and checks representable result-size overflow; `shift_right` is exact floor division by `2^k`. All operations return owning values and leave input values unchanged. No CPU-specific intrinsics, inline assembly or runtime ISA dispatch are present yet.

**Cost and deliberate limits.** Allocations are avoided for the common result range of at most four limbs, including the **normalization of carry-free add/shift results** that provisionally require an extra limb. Data for larger values uses standard contiguous vector storage. Still to be implemented: writable limb views and external scratch, in-place arithmetic, multiplication and squaring, exact division/remainder, small-divisor APIs and GCD; only then test optional ARM64 and x86-64 variants. This phase makes no performance claims.

The actual repository CTest must also pass in CI before calling the change fully integrated across platforms. These standalone runs used a minimal local macro shim so they are not a substitute for repository CMake/CTest.

**Architecture specific code** remains strictly subject to the mandatory per-kernel journal and ISA requirements in `floatmp_design.md`; there is currently no specialization to document. Neither benchmarks nor the frozen ARM64 20,000x31 FloatMP workload was changed.

### Authorized conditional ARM64/x86-64 source split (18 October 2026)


The user explicitly approved the **exact future filenames** `src/big_unsigned_arm64.cpp` and `src/big_unsigned_x86-64.cpp` **if and when** distinct architecture-specific implementations are justified. These files are **not part of phase B0**, do not exist yet and must not be created as empty placeholders.

`include/numeric/big_unsigned.hpp` remains the single architecture-neutral interface. `src/big_unsigned.cpp` owns the portable C++20 reference implementation, storage, normalization and any architecture-neutral selection logic. Genuinely ARM64/AArch64-specific low-level kernels may eventually live in `src/big_unsigned_arm64.cpp`; genuinely x86-64-specific low-level kernels may eventually live in `src/big_unsigned_x86-64.cpp`. The hyphen in `x86-64` is intentional in this approved filename. No third-party arithmetic engine and no duplicated floating-point semantics.

**Split criterion.** Create a dedicated source file only when a first kernel using native ISA instructions demonstrates a meaningful advantage under a reproducible, full-operation workload and its portable equivalent has passing bit-exact differential tests. Document separately: contract and carry bounds, instruction/CPU-feature requirements, compiler flags, runtime/compile-time dispatch and portable fallback, aliasing and scratch lifetime, test forcing of each path, ISA-appropriate CI and sanitizers, benchmark host/compiler, absolute timings, and operation-level crossover. Architecture files must be compiled only on supported targets; generic builds including Windows/MSVC must continue to build through portable fallbacks. No unguarded instruction assumptions. The per-kernel documentation belongs in this journal and in `floatmp_design.md` with the same source change.

**Current validation note.** The original B0 commit `be34905` passed Unix and macOS CI jobs but its Windows Release job failed during **Build**, before tests. The precise MSVC diagnostic has not yet been verified; this documentation-only approval does not resolve the Windows failure. Do not call phase B0 cross-platform validated until Windows passes.

### BigUnsigned B0 Windows fix and style (18 October 2026)


The initial `be34905` Windows 2025 Release job failed in MSVC 14.51.36231 when parsing `BigUnsigned::_inline` in the class declaration. The resulting C2226/C2334 and follow-on C3409/C2143/C2059 diagnostics point to an identifier collision: `_inline` is interpreted as a compiler-specific inline keyword/extension. The member is renamed everywhere to `_inline_limbs` with no change to the storage representation, arithmetic algorithms, tests, or ABI-visible public methods.

The same maintenance commit expands whitespace, blank-line separation and control-flow formatting in `src/big_unsigned.cpp` and the dedicated test, following neighboring Ariadne numeric sources. Existing differential and adversarial tests are retained unaltered in their mathematical coverage. MSVC/Windows CI must pass before claiming the fix is verified; GCC/Clang regressions must also be checked. No architecture-specific kernels or FloatMP arithmetic changes are introduced.

### BigUnsigned B0 Windows independent exact-reference fixtures bit-count fix (18 October 2026)


Their shift parameter is `mp_bitcnt_t`, which can be 32-bit on Windows, whereas the randomized native shift was held in 64-bit `std::size_t`. No warning is disabled or demoted, no test distribution/case coverage changes, and no production `BigUnsigned` code is modified. Revalidate the Windows Build and CTest alongside the existing Linux and macOS checks.

### BigUnsigned B1: portable full and truncated multiplication (18 October 2026)


Implemented in the common numeric substrate, not FloatMP: `mul(a,b)` computes the exact nonnegative product; `square(a)` is currently `mul(a,a)`; `mul_low(a,b,k)` computes the product modulo `2^(64*k)` without calculating higher limbs. Inputs are immutable, output canonical; overlarge sizes raise `std::length_error`.

**Arithmetic proof:** break 64-bit factors into 32-bit halves, form four 32x32 products, then combine their high/low halves with a bounded middle sum to obtain exact 128-bit `high:low`. In row multiplication with radix `B=2^64`, the next-limb carry arises from the high product plus two overflow bits from accumulating the existing destination limb, low product and prior carry. The full combined sum is bounded by `(B-1)^2 + 2(B-1) = B^2-1`; therefore the next carry fits 64 bits. `mul_low` ignores contributions with `i+j>=k`, so its modular result is exact. The destination is zero-initialized and disjoint from source spans.

The integrated CTest adds 6,000 more deterministic randomized pairs, including all-one and zero limbs, varying sizes up to 34 limbs, low truncation to 0/1/4/5 or random limb lengths, full product, square and commutativity. The project CI must still validate MSVC and complete CTest. No performance improvement is claimed; a specialized fast square, in-place views, scratch reuse and ARM64/x86-64 kernels remain future work. The existing Integer and FloatMP 20,000x31 benchmarks are unchanged.

### BigUnsigned B2: limb views, reusable scratch, direct-output kernels (18 October 2026)


Portable shared substrate adds `LimbView` and `MutableLimbView` (non-owning `std::span` aliases), `BigUnsigned::Scratch` with sixteen inline 64-bit limbs and vector-backed reuse for larger requests, plus `BigUnsigned::assign(LimbView)`. `assign` normalizes its result and supports self-referential subspans. `Scratch::acquire(count)` zeroes the returned range, invalidates previous views, and reuses vector capacity when sufficient; scratch belongs to the caller, without globals or concurrency sharing.

The new `add_into` and `sub_into` accept two read-only limb views and one disjoint writable destination whose length must cover the larger input. They write the low output limbs and zero higher available destination limbs, returning the final 0/1 carry or borrow. `sub_into` computes a modular subtraction even when the mathematical difference is negative; the public `sub` keeps its checked nonnegative contract. `mul_low_into` computes the exact product modulo `2^(64*output.size())`, including clearing the destination on zero input. All direct-output kernels forbid input/output overlap and allocate no internal memory; caller-provided scratch can be reused. Existing owning `add`, `sub` and `mul_low` now delegate to these same kernels. The multiplication carry proof from B1 applies unchanged.

GCC and Clang strict C++20 builds as well as Clang ASan/UBSan passed locally in a shim-assisted standalone build. No speedup is claimed.

### BigUnsigned B3: selectable x86-64 multiplication (18 October 2026)


**Scope and identity.** B3 introduces `mul_into` (full exact product into a disjoint caller-owned buffer of at least `a.size()+b.size()` limbs), `mul_into_portable`, and `mul_low_into_portable`, alongside the existing `mul_low_into` selected kernel. Empty factors zero all output limbs; truncated output computes the exact residue modulo `2^(64*output.size())`; all kernels zero the full output and allocate no scratch. Inputs are immutable, little-endian 64-bit limbs. `mul` uses the direct full-product output kernel. No FloatMP, Integer, Rational, or existing benchmark workload is modified.

**ISA and dispatch.** `src/big_unsigned_x86-64.cpp` is compiled only on x86-64/AMD64 with GCC, Clang, AppleClang, or MSVC. Its 64x64->128 unsigned product uses GCC/Clang unsigned `__int128` in this conditional source or MSVC x64 `_umul128`; neither is required by the C++20 reference in `src/big_unsigned.cpp`. x86-64 baseline integer multiplication is sufficient: no BMI2, ADX, CPUID, runtime dispatch, assembly, or optional instruction flags. The selected path calls the architecture kernel only for two or more limbs in *both* inputs, otherwise falls back to the portable kernel. Other targets always use the reference; `has_native_multiply()` reports whether the optional implementation was built. Explicit portable entry points are available to tests and microbenchmarks. The CMake CPU and compiler conditions are a build-time feature gate.

**Proof and memory contracts.** Represent radix `B=2^64`. For one inner product, `a[i]*b[j]+out[i+j]+carry <= (B-1)^2+2(B-1)=B^2-1`; thus the carry is a 64-bit limb. Adding low-product, old destination and prior carry generates at most two carry-out flags, which are added to the high product. Full and truncated products use identical recurrence and exact modular truncation. Caller owns output, which must not overlap either source (including partial overlaps); no internal allocation, global scratch, or mutable dispatch state. Contract violations for too-small nonempty full output raise `std::length_error` before touching the result; full/low outputs for empty factors are zero. No alignment beyond naturally valid 64-bit spans is required. Each call is independent and thread-safe when buffers are disjoint.

The selected and explicitly portable kernels are compared word-for-word for full and low products; output overwrite, scratch reuse and insufficient full-buffer exception are exercised. A matching isolated x86-64 algorithm probe passed 45,736 edge/random cases with GCC 14.2, Clang 17, and Clang ASan/UBSan. These local checks are **not** a substitute for the project's full CMake/CTest and Windows/MSVC/Apple ARM64 validation.

**Performance and acceptance.** Previously reported x86-64 preallocated-output 20,000x31 microbenchmark medians (GCC 14.2, -O3) were portable/optimized: 53 bit 10.29/10.30 ns, 128 bit 17.97/16.14 ns, 256 bit 57.78/32.27 ns, 1024 bit 877.04/386.57 ns, 2048 bit 3576.64/1558.36 ns; these are **user-provided prior-session measurements**, not new measurements of this integration. In particular, 1024-bit speedup was 2.27x in a kernel-only comparison. No adaptive performance crossover is claimed from these results. The source paths are `src/big_unsigned.cpp` and `src/big_unsigned_x86-64.cpp`; these reference kernels must continue to agree bit-for-bit. ISA reference: x86-64 baseline unsigned `MUL` (Intel SDM); MSVC `_umul128` compiler intrinsic contract. No third-party arithmetic source was copied.

### BigUnsigned B3: first integrated Ubuntu x86-64 benchmark and CLI alignment (18 October 2026)


**User-reported result, B3 commit `f8e56dc`.** On Ubuntu x86-64, `benchmark_big_unsigned_multiply` ran with 20,000 iterations and 31 median samples; the selected x86-64 path was active. Timings are nanoseconds per operation, with the originally frozen deterministic operand construction:

Defaults remain 20,000 / 31 and operand sizes 53/128/256/1024/2048 bits. Only argument handling and presentation differ. Tests and benchmarks remain separate.

### BigUnsigned B3: Ubuntu x86-64 benchmark after CLI/report integration (18 October 2026)


**User-reported same-host rerun** of the updated benchmark at `f87dcd9`, with the unchanged five deterministic operand sizes, 20,000 iterations and 31 median samples. The selected native x86-64 path was active. These are not performance measurements obtained by the assistant.

**Interpretation constraints.** The 2048-bit owning operation appears about 3.2% faster than its nominally cheaper selected direct-output kernel despite owning a result; the measurement order, runtime state and compiler may contribute. Do not assume this reflects an allocation benefit. This benchmark uses one deterministic operand pair per size repeatedly: broaden operand patterns in a **separate diagnostic**, without changing or quietly rebasing this existing workload. Profile actual instruction-level inner loops before attempting an adaptive algorithm or new ISA optimization. Do not alter the mathematical oracle or report microbenchmark gains as full public-operation gains.

### BigUnsigned multiplication benchmark: explicit allocation/dispatch labels (18 October 2026)


Following user review of the B3 report, the benchmark labels its three BigUnsigned cases as `prealloc_portable`, `prealloc_selected`, and `new_selected`. The first calls `BigUnsigned::mul_into_portable` with caller-provided output, the second calls `BigUnsigned::mul_into` with caller-provided output, and the third calls the public `mul(a,b)`, which also dispatches to the selected architecture kernel but owns/normalizes a newly created result. Output column widths are adjusted for the longer labels. Row labels, corresponding diagnostic labels and aggregate names change, **not** operand generation, timing order, tested algorithms, mathematical oracle, warm-up, iterations, sample counts, allocation or measurement semantics.

### BigUnsigned B4 experiment: x86-64 column-wise Comba candidate (18 October 2026)


**Experimental only:** `BigUnsigned::mul_low_into_comba` and `mul_into_comba` provide a forced, allocation-free comparison kernel for the same disjoint read-only limb views and caller-owned output views as B3. The selected production `mul_into`, `mul_low_into`, owning `mul` and `square` continue to use exactly the B3 dispatch. The Comba candidate is implemented only in `src/big_unsigned_x86-64.cpp`, compiled behind the existing x86-64 build gate; on other platforms its explicit entry point delegates to the portable reference. No runtime feature probing, optional CPU extensions or architecture-specific floating-point behavior is introduced.

**Arithmetic and carry proof.** Column `k` sums all partial products `a[i]*b[k-i]` where `0<=i<a.size()` and `0<=k-i<b.size()`, plus the carry from the prior column. The sum is held in three unsigned 64-bit words `(top,hi,lo)` (radix `B=2^64`). Each exact product `(product_hi,product_lo)` adds to `(hi,lo)` using 64-bit wrapping additions, with each carry-out propagated to the next word. For at most `min(a.size(),b.size())` products in any column, the three-word sum is sufficient for any representable input spans: the worst case is `min(n,m)*(B-1)^2` plus carry of at most `min(n,m)*(B-1)`, giving a high-word count bounded by the number of products and representable span lengths. After emitting `lo` to output[k], shift `hi` and `top` down to serve as the following column's initial carry, and zero `top`. The algorithm computes the exact full result or its low `output.size()` limbs; output is completely zeroed beforehand and nonoverlapping input/output is required. No extra scratch, heap allocation, or mutable global state; `mul_into_comba` checks full-buffer length before writing and handles empty operands like B3. x86-64 wide products use the same baseline `MUL` lowering via conditional GNU/Clang unsigned `__int128` or MSVC `_umul128`; no additional ISA feature is needed.

On non-x86-64 this also exercises the documented reference fallback.

**Performance and report controls.** The benchmark retains default 20,000 iterations, 31 samples and the existing 53/128/256/1024/2048-bit operand set. The default benchmark must remain comparable with the B3 user-reported Ubuntu timings. An isolated probe with GCC/Clang showed inconsistent crossovers between Comba and existing schoolbook depending on compiler and system state; do not claim Comba faster from those exploratory times. Require same-host integrated measurements and green cross-platform CI before considering any production dispatch change. User-facing mode labels follow B3 naming.

### BigUnsigned B4 measurements and B5 fused-accumulation experiment (18 October 2026)


**B4 user-reported Ubuntu x86-64 measurements (20,000 iterations, 31 samples).** All timings in ns per multiplication, using the frozen deterministic inputs. Comba is faster only at 53 bits and slower at the other four sizes, so **B3 remains the production selector**. Windows, Unix and Coverage CI all passed on B4 commit `47776c9`.

**B5 hypothesis.** A native x86-64 schoolbook inner step can be expressed as a single unsigned 128-bit mathematical sum `total=(Wide)a[i]*b[j]+out[i+j]+carry`; write its low 64 bits to the destination and retain the high 64 bits as the next carry. The bound `(B-1)^2+2(B-1)=B^2-1` for radix `B=2^64` proves that the intermediate is exactly representable in 128 bits (no overflow). This differs from B3's explicitly separated 64-bit low/high product and two overflow flags only in how the compiler sees the accumulation; it does not change the arithmetic algorithm or operation count. Conditional x86-64 GCC/Clang use `unsigned __int128`, while MSVC x64 delegates this optional entry point to its existing, verified native B3 multiply, and non-x86-64 platforms delegate to the portable C++20 reference. No BMI2/ADX, runtime CPU dispatch, global scratch or new allocation is introduced. Both input/output aliasing and full/truncated destination constraints remain those of B3.

Compiler-dependent results motivate a gated benchmark comparison before any production path change.

**B5 implementation and benchmark protocol.** Expose optional `mul_low_into_fused` and `mul_into_fused`, with no changes to `mul_low_into`/public owning `mul`. Keep 20,000x31 defaults and the order/content of the four B3 baseline measurements intact; experimental Comba is likewise optional. A local 12,000-case oracle does not establish full project CI correctness, nor does performance on one compiler establish a dispatch crossover. Gate any future selector change on successful Unix/Windows/Coverage CI and repeated same-host GCC/Clang measurements. Both shared design journals must remain synchronized.

### BigUnsigned B6: GCC x86-64 selected fused multiplication (18 October 2026)


**Decision and measurements.** A user-reported integrated Ubuntu x86-64 20,000-iteration/31-sample B5 benchmark measured `prealloc_selected` vs forced `prealloc_fused` medians (ns): 53 bits 15.442/13.076, 128 bits 18.834/13.169, 256 bits 29.297/20.148, 1024 bits 323.861/205.060, 2048 bits 1263.749/807.183. These are user-reported integrated figures on Ubuntu GCC x86-64, not an independent cross-compiler acceptance claim. Clang probe performance in the B5 journal was inconsistent; no Clang/MSVC/ARM64 speedup may be inferred.

**B6 build and dispatch.** When `CMAKE_SYSTEM_PROCESSOR` identifies x86-64/AMD64 and `CMAKE_CXX_COMPILER_ID` is `GNU`, the CMake option `NUMERIC_BIG_UNSIGNED_GCC_FUSED_MULTIPLY` defaults to `ON` and defines the private compile-time `ARIADNE_BIG_UNSIGNED_GCC_FUSED` switch. `BigUnsigned::mul_low_into` then selects the proven B5 unsigned `__int128` fused inner loop, including one-limb cases. The public owning `mul` and `mul_low` inherit that selection through the existing `mul_into` and `mul_low_into` entry points. GCC x86-64 with the option `OFF`, Clang/AppleClang x86-64 and MSVC x64 retain the B3 selected kernel and its one-limb portable fallback. Other architectures retain the C++20 portable reference. There is no runtime CPUID, nonbaseline instruction, ABI change, shared mutable state, input/output aliasing or scratch allocation. B5's carry bound `(B-1)^2+2(B-1)=B^2-1` applies without modification. The `selected_multiply_uses_fused()` query reports the build-time selection.

**Controlled comparison and regression.** `mul_low_into_schoolbook` and `mul_into_schoolbook` explicitly force the former B3 selected schoolbook (or portable fallback elsewhere), in addition to pre-existing portable, Comba and fused explicit kernels. `--include-schoolbook` adds an opt-in `prealloc_schoolbook` row after the four standard report entries and any other experimental comparisons; normal benchmark defaults, operands, measure order and 20,000x31 counts stay unchanged. `prealloc_selected` and `new_selected` **now legitimately use fused on GCC x86-64 default builds**; the older B3 timing is available via explicit schoolbook and an option-OFF control build. Note this changes the meaning of `selected` relative to previous commits but not what its name promises: selected is always architecture/compiler-dependent. Acceptance of B6 requires complete Unix/Windows/Coverage CI, ideally a matched GCC option-OFF CTest and repeat GCC/Clang/Apple ARM64 benchmark checks before any further promotion. B6 does not migrate Integer/Rational or FloatMP.

### BigUnsigned B6: selected fused Ubuntu x86-64 integrated results (18 October 2026)


**User-provided measurements** for B6 commit `63ce810` (GCC x86-64 with `selected kernel=fused`, 20,000 iterations and 31 median samples). Units: ns per multiplication, fixed operand pool unchanged since B3.

B6 selected and forced fused agree closely, validating routing/performance in the integrated Ubuntu benchmark.

**Validation.** GitHub Actions Unix, Windows and Coverage succeeded for `63ce810`. Independent B6 user results, not assistant-side execution.

### BigUnsigned B7: Ubuntu x86-64 result-lifetime benchmark (18 October 2026)


Values below are medians in ns per multiplication. This is not an independent assistant-side timing or a measurement of ARM64/Clang/MSVC.

The additional reported modes (53/128/256/1024/2048 bits) were portable 15.370/23.819/45.232/591.363/2304.662 ns, forced fused 13.081/13.682/19.615/203.619/756.209 ns, and forced B3 schoolbook 17.314/20.246/30.337/320.331/1221.233 ns. All three GitHub Actions workflows (Unix/Windows/Coverage) passed on B7 commit `6722b74`.

Thus at 256 bits lifecycle overhead is comparatively important, but at 1024 and particularly 2048 bits most of the absolute gap is already present with reusable output buffers. These differences of medians are descriptive, not an isolated attribution of allocator vs canonicalization, and paired timing order and compiler optimizations remain potential confounders.

**B8 direction (not yet implemented).** Preserve the fully validated GCC fused B6 selector. A dedicated operand-distribution diagnostic should span balanced and unbalanced lengths and avoid quietly replacing the five-case frozen benchmark. Do not widen `BigUnsigned` inline storage globally or change public representation solely to optimize the 256-bit benchmark: that affects every `Integer`/`Rational`/future `FloatMP` object.

### BigUnsigned B8: opt-in allocation-free Karatsuba experiment (18 October 2026)


B8 therefore tests an independently implemented Karatsuba full-product kernel, while **preserving production B6 selected fused GCC x86-64, existing B3/B4/B5 forced kernels and all four default benchmark modes**. The experimental `BigUnsigned::mul_into_karatsuba(a,b,out,scratch)` computes an exact full product into an output view of at least `a.size()+b.size()` limbs. Inputs, output and scratch must be mutually disjoint; all result words including excess output are overwritten; neither this kernel nor recursive calls allocates. `karatsuba_scratch_limbs(a.size(),b.size())` computes a checked upper bound on required reusable 64-bit scratch words. Insufficient output or scratch raises `std::length_error` before modifying the destination (nonempty inputs). An empty factor zeroes the destination.

**Algorithm and arithmetic invariant.** Split at `k=floor(max(n,m)/2)`; the fast path applies only when both operands have nonempty upper halves and `min(n,m)>8` limbs. Let `z0=a0*b0`, `z2=a1*b1` and `t=(a0+a1)(b0+b1)-z0-z2=a0*b1+a1*b0`. Recursively write z0 and z2 directly into nonoverlapping output slices. Use two scratch-resident half sums, each with one extra high limb, and one scratch-resident cross-product buffer. Subtractions are exact nonnegative operations and cross-term is added to the shifted destination with carry propagation. The recursion scratch requirement is `2*(sa+sb)+max(work(z0),work(z2),work(sumA*sumB))`, where `sa=max(|a0|,|a1|)+1`, `sb=max(|b0|,|b1|)+1`; recursive calls reuse scratch sequentially. Every sum and carry is computed from unsigned 64-bit limbs with overflow checks; the full output is exact at radix 2^64. The recursion's base case uses the current selected (not necessarily portable) B6 schoolbook multiply. This is a full-product implementation only: no general truncated Karatsuba path or production crossover is claimed. The file `src/big_unsigned_karatsuba.cpp` is portable C++20 and requires no `__int128`, third-party source or CPU-specific instructions.

Clang ASan/UBSan passed **4,616** focused differential cases. Those integrated modifications still require Unix/Windows/Coverage CI acceptance. No mutation of the production selector or test oracle is permitted to rescue an unsuccessful experiment.

**Isolated performance (not project results).** A standalone 20,000-iteration/31-sample GCC x86-64 microbenchmark recorded fused/Karatsuba (ns): 1024 bits 316.92/410.54, 2048 bits 1185.21/1480.69, 4096 bits 4703.88/5086.85, 8192 bits 19268.6/16420.2. Clang on the same probe recorded 1024 bits 290.36/398.26, 2048 bits 1078.98/1535.44, 4096 bits 4306.24/5110.39, 8192 bits 17939/16803.1. These are **exploratory isolated measurements**, not the unchanged integrated Numeric benchmark, and show that this initial implementation is slower through at least 4096 bits; only at 8192 bits did a modest speedup appear. The crossover will depend on compiler, input distribution and runtime environment; there is insufficient evidence to make Karatsuba production-default, particularly for the project's 53/128/256/1024/2048-bit performance gate.

For assessing large-input crossover, pass `--bits 1024,2048,4096,8192` explicitly and check `prealloc_karatsuba` versus `prealloc_selected` on the same Ubuntu host. Do not conflate a high-bit diagnostic with replacement of frozen baselines or a whole FloatMP performance result. B8's purpose is to establish evidence, not enable subquadratic dispatch prematurely.

### BigUnsigned B8 build fix: testing macro shadow (18 October 2026)


Ubuntu GCC found an integrated CTest **test-source compilation error**, not evidence of an arithmetic failure, in the newly added B8 wide-limb loop: `ARIADNE_TEST_ASSERT(...)` internally declares `auto result=(expression)`, so a local scratch view named `result` caused name shadowing and use-before-deduction when passed inside the assertion. The build issue was reported with GCC `-Wshadow`, and the full Unix/Windows/Coverage CI for the corrected commit remains the acceptance gate.

### B8.1: Karatsuba cutoff experiment (18 October 2026)


The user-reported Ubuntu 20,000-iteration, 31-sample timings show Karatsuba cutoff 8 slower than fused selected at 1024 (346.542 vs 206.249 ns), 2048 (1256.520 vs 797.097 ns), 4096 (4203.549 vs 2975.416 ns), and 8192 bits (13583.344 vs 12589.667 ns). This contradicts the isolated-prototype crossover suggestion; no production promotion is justified. Same threshold must be used for sizing and multiplication. No runtime workspace allocations, production changes, or FloatMP changes. The CI still needs to be checked for the new commit.

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


An architecture-specific AArch64 multiplication kernel is added in `src/big_unsigned_arm64.cpp`, compiled for aarch64/arm64 with GCC, Clang or AppleClang. Its explicit `mul_into_fused` path uses unsigned 128-bit multiplication and accumulation, lowering to baseline ARM64 `MUL`, `UMULH` and carry operations; compilation to those instructions was checked with Clang 17 targeting AArch64. This confirms code generation, **not performance**. The arithmetic bound `(2^64-1)^2+2(2^64-1)=2^128-1` proves each wide intermediate fits in 128 bits. Truncated and full multiplication, disjoint inputs and outputs, fully overwritten result, and no internal allocations preserve the B5 semantics. `has_native_multiply` now means x86-64 OR ARM64 native support, while `has_x86_64_multiply` specifically identifies availability of Comba/B3 schoolbook. The benchmark permits `--include-fused` on ARM64 but correctly restricts x86-specific modes.

**Dispatch discipline:** ARM64 production selection remains the portable kernel by default. The CMake option `NUMERIC_BIG_UNSIGNED_ARM64_FUSED_MULTIPLY=ON` enables fused ARM64 dispatch explicitly for controlled experiments; the default is OFF, as there are no integrated user macOS ARM64 speed measurements. GCC x86-64 fused production selection and Windows/non-ARM platforms remain unchanged. Accept production promotion only after correct macOS CI/local CTest and repeated same-host measurements demonstrate stable benefit. No FloatMP code change or ARM64 performance claim.

### B9 native Apple Silicon ARM64 benchmark, user-reported (18 October 2026)


The benchmark used 20,000 iterations, 31 samples, the five frozen operand sizes, and unchanged deterministic inputs. All times in ns per multiplication:

Same-run portable-selected to forced-fused speedups: **1.208x**, **1.239x**, **1.626x**, **1.935x**, **1.719x** at the respective sizes, equivalently -17.2%, -19.3%, -38.5%, -48.3%, -41.8% execution time. Thus the explicit ARM64 `__int128` kernel has an important integrated speedup over the baseline portable path on the user's Apple Silicon machine, especially at 1024/2048 bits. Treat these cases as possibly affected by compiler optimization, timing noise or benchmark design; do not assume isolated cause or universal advantage. The commit `3b34604` passed Unix/Windows/Coverage CI, but those workflow results alone do not establish the new kernel was exercised on macOS ARM64.

**Next gate, not yet a default-dispatch change.** Repeat the same user macOS ARM64 benchmark using a separate build configured with `-DNUMERIC_BIG_UNSIGNED_ARM64_FUSED_MULTIPLY=ON` and run project correctness tests on that build. Prefer multiple repeated paired runs before enabling ARM64 fused by default. The present OFF build and all non-ARM64 production code remain unchanged.

### B9 Apple Silicon ARM64 fused selector ON: user-reported owning results (18 October 2026)


**Source:** User reran the integrated benchmark on macOS Apple Silicon, 20,000 iterations, 31 median samples, five original operand sizes, after setting `NUMERIC_BIG_UNSIGNED_ARM64_FUSED_MULTIPLY=ON` and rebuilding. Banner explicitly reported `native available=true | selected kernel=fused`. Units ns; matched rows within this run:

**Comparison against previous OFF run on the same user's Mac** (separate benchmark invocations, not paired trials): `new_selected` OFF=10.842/12.069/44.527/368.960/1396.992 ns versus ON=11.065/11.508/32.946/202.990/836.858 ns at 53/128/256/1024/2048 bits. ON/OFF speedup is approximately 0.98x/1.05x/1.35x/1.82x/1.67x. Small-input 53- and 128-bit owning differences are not reliable signals; the 256-bit and especially the 1024/2048-bit gains are large. Within the ON run, `prealloc_selected` and `prealloc_fused` match closely at 53,256,1024,2048 bits; the ~9.2% difference at 128 bits suggests microbenchmark variability or call-path effects, not a distinct algorithm. These discrepancies underscore the need for repeated measurements and ensuring compiler optimizations do not render fixed-operand benchmarks unrepresentative.

GitHub Actions Unix/Windows/Coverage are green on B9 `3b34604` and documentation commit `e827a57`; macOS ARM64 project CTest in the selected-fused ON build has not been reported and therefore must not be claimed. Keep ARM64's default option OFF pending a clean macOS ARM64 test run and ideally repeated same-host measurements. Candidate for a future default change is now supported by *real* user integrated timings, not by a standalone prediction.

### B9.1: promote native ARM64 fused multiplication without a CMake toggle (18 October 2026)


**Acceptance:** The user reported the macOS ARM64 project tests green with fused dispatch active. Paired selected/forced-fused integrated measurements on Apple Silicon were approximately 16.600/16.592 ns at 256 bits, 180.769/180.544 ns at 1024 bits and 805.240/802.896 ns at 2048 bits; the owning public API also improved materially compared with the earlier portable selection. The earlier B9 notes describe the historical OFF/ON experiments; those settings are now **retired**.

**Final production design:** For supported native ARM64 GNU/Clang/AppleClang builds, `ARIADNE_BIG_UNSIGNED_ARM64` alone enables compilation *and* selection of the fused ARM64 kernel. Remove both the public CMake option `NUMERIC_BIG_UNSIGNED_ARM64_FUSED_MULTIPLY` and the redundant private macro `ARIADNE_BIG_UNSIGNED_ARM64_FUSED`; no transitional flag, runtime toggle, or compatibility branch is retained. The separately exposed `mul_into_portable` and `mul_low_into_portable` remain sufficient for reference comparisons. A direct ARM64 test assertion now checks that the native fused path is selected. x86-64 GCC retains its independent selector/option, while other x86 compilers keep B3 and unsupported platforms remain portable. No changes to FloatMP code or claims about FloatMP performance. New revision requires Unix/Windows/Coverage CI verification; the reported macOS tests covered the prior option-enabled but arithmetically identical implementation.

### B10 ARM64 explicit carry experiment (18 October 2026)


New B10 experimental AArch64 row-schoolbook multiplication `mul_low_into_arm64_carry` / `mul_into_arm64_carry`: baseline instruction sequence MUL, UMULH, ADDS, ADC with two carry propagations per partial product, preserving the 128-bit bound and allocation-free full/truncated semantics. Clang 17 assembled an isolated instruction probe for both AArch64 Linux ELF and Apple arm64 Mach-O; **no native Apple Silicon timing has been obtained**. The production fused ARM64 selector is unchanged. The forced kernel uses portable fallback off ARM64, the integrated `--include-arm64-carry` benchmark mode requires ARM64 and appends an optional comparison. Any promotion requires correctness checks and same-host user benchmark. Compare `prealloc_selected` against `prealloc_arm64_carry`. No performance improvement is asserted.

### B10.1: source layout consolidation (18 October 2026)


Per repository architecture rule, only architecture-specific implementations have dedicated translation units. The experimental allocation-free Karatsuba implementation, its cutoff selection, sizing, recursive helpers and public method definitions have been moved **without algorithmic changes** into `src/big_unsigned.cpp`; `src/big_unsigned_karatsuba.cpp` is deleted and removed from `src/CMakeLists.txt`. GCC/Clang's identical `unsigned __int128` fused accumulation implementation is now maintained **once**, as `multiply_limbs_fused` in the shared `src/big_unsigned.cpp`. The duplicate bodies and entry points in `src/big_unsigned_x86-64.cpp` and `src/big_unsigned_arm64.cpp` are removed. ARM64 source retains only its explicit native `MUL/UMULH/ADDS/ADC` experiment; x86 source retains B3 schoolbook/Comba and the MSVC `_umul128` implementation. The shared fused helper keeps the x86-64 MSVC fallback to B3, the forced fused wrapper's architecture availability, and GCC x86-64/ARM64 default dispatch unchanged. On toolchains without `__int128`, the helper's fallback remains portable; unsupported platforms select the existing portable reference. This is a structural refactor, **not** a claimed speedup: recheck actual compiler output and same-host performance after consolidation, as moving a kernel between translation units can affect optimizer behavior. The user requires no nonarchitecture extra .cpp files for multiplication algorithms; do not recreate a standalone Karatsuba source. CI and macOS ARM64 test/benchmark results for this new layout are not yet verified.

### B10.1 Apple Silicon ARM64 carry-chain user benchmark (18 October 2026)


All figures are median ns:

Same-invocation fused-selected / explicit carry speedup = 1.337x/1.075x/0.998x/1.182x/1.413x for 53/128/256/1024/2048 bits. Larger balanced 1024- and 2048-bit cases show 15.39% and 29.24% lower time, while 256-bit times are effectively identical. The new native assembly offers actual user measured benefit at 1024/2048 bits. The selected fused kernel itself remains near previous 180.769/805.240 ns measurements at 1024/2048 bits, giving a coarse (not conclusive) cross-run indication of no consolidation regression.

**Decision:** keep the current ARM64 production fused selector unchanged until the user runs macOS ARM64 `ctest --output-on-failure` on the newly consolidated and carry-enabled source tree. Unix CI alone cannot certify native ARM64 inline assembly correctness. Rerun representative timings on the same machine, preferably reversed bit-size ordering and including 4096 bits, before selecting the assembly by default. The benchmark uses fixed operands and fixed per-backend measurement ordering, so it is not a universal performance conclusion and notably tiny-operand owning results can be affected by optimization and run-to-run noise. Do not add a CMake experimentation flag, and do not create extra nonarchitecture source files. The two architecture files remain the only exceptional compilation units; the fused and Karatsuba algorithms are kept in the shared file.

### B10.1 second ARM64 carry-chain benchmark, including 4096 bits (18 October 2026)


The user ran `benchmarks/benchmark_big_unsigned_multiply --include-arm64-carry --bits 256,1024,2048,4096 -t dark` from `build-release` on Apple Silicon. Selected remained `fused` (not the carry candidate); 20,000 iterations, 31 samples. Same-run median ns:

Carry vs selected fused same-run speedups: 1.268x (256), 1.166x (1024), 1.395x (2048), 1.623x (4096). Carry execution-time savings: 21.2%, 14.2%, 28.3%, 38.4%.

**Repeatability and qualification:** B10.1's first Apple Silicon benchmark showed carry 16.617/154.744/571.048 ns at 256/1024/2048 bits versus selected fused 16.585/182.898/806.990 ns. The second run reproduces the strong carry benefits at 1024 and 2048 bits, adds a large 4096-bit benefit, but shows variable relative behavior at 256 bits (first run tied, second run faster). One interpretation is timing variability at short durations; do not infer a reliable 256-bit threshold or claim universal speedup. Because carry remains opt-in, user `new_selected` still measures fused and does not yet represent owning carry performance. The user has **not yet reported** native macOS ARM64 `ctest` results for the consolidated B10.1 code, so do not change production dispatch at this stage; retain current fused ARM64 selected kernel. GitHub Actions for the B10.1 source-layout refactor `5d63d7e` passed Unix/Windows/Coverage; those jobs are not proof of ARM64 inline assembly execution. Architecture-only extra translation-unit rule remains in force.

### B10.2: promote ARM64 carry-chain to production (18 October 2026)


Two user-run integrated macOS benchmark invocations measured real carry-chain wins over selected fused at 1024 and 2048 bits; the second also measured 4096 bits: selected fused 3673.671 ns, carry 2263.667 ns (1.62x faster); at 2048 bits fused 796.533 ns, carry 570.777 ns (1.40x faster). At 256 bits results varied across runs and no size-specific threshold is asserted. This is measured integrated evidence, not an isolated prototype prediction.

**Change:** On all supported ARM64 GCC/Clang/AppleClang native builds, `BigUnsigned::mul_low_into` and therefore `mul_into` and owning `mul` select the existing explicit ARMv8-A `MUL/UMULH/ADDS/ADC` carry-chain implementation unconditionally. No CMake options or extra translation units are created. The shared portable `mul_into_portable` and common `mul_into_fused` remain available for comparisons. `selected_multiply_uses_fused()` is now **false on ARM64**, correctly indicating that selected is no longer the shared fused kernel; benchmark output reads `selected kernel=ARM64 carry` based on native architecture availability. The GCC x86-64 fused selection, Clang/MSVC x86-64 B3 fallback, and other architecture portable fallback are unchanged.

Performance is not assumed to carry over mechanically from preallocated forced calls; compare to same-host historical baselines. Full end-to-end FloatMP performance remains unmeasured.

### B10.2 Apple Silicon production carry dispatch benchmark (18 October 2026)


The banner confirmed `native available=true | selected kernel=ARM64 carry`. User-reported median times in nanoseconds:

**Owning speedups verified against the previous fused-selected user run (separate invocations):** at 1024 bits, 202.990 -> 172.594 ns (1.176x speedup, ~15.0% less time); at 2048 bits, 836.858 -> 605.933 ns (1.381x, ~27.6% less time). The tiny widths 53/128/256 bits have no robust owning performance trend. At 256, 1024 and 2048 bits, selected and forced carry are nearly identical, verifying the new dispatch. At 53 and 128 bits, selected appears slower than forced carry by about 46% and 12%; because their arithmetic implementation is the same this is a call-path, optimizer, benchmark-order or timing issue to investigate, **not** evidence of two different algorithms. Forced fused remains available as a performance control.

**Validation state:** The user earlier reported macOS ARM64 `ctest` green with the carry kernel explicitly forced, prior to the B10.2 production selector change. No separate local macOS CTest result has yet been reported for this post-promotion build; benchmark's correctness oracle is not a replacement for full CTest. Unix/Windows/Coverage CI statuses should be inspected for this new commit. Both design journals are updated for the measured owning promotion. No extra translation units or performance CMake toggles. Does not demonstrate an end-to-end FloatMP speedup.

### B10.2 validation closed: selected ARM64 carry, post-promotion (18 October 2026)


**User confirmation:** `ctest --output-on-failure` on the updated macOS ARM64 build was reported **fully green after** the B10.2 production dispatch change. This supersedes the earlier pending native-CTest qualification in the B10.2 benchmark record. **GitHub Actions on functional dispatch commit `d01c505`: Unix, Windows, Coverage all completed successfully.** The subsequent journal-only `4bf1886` CI jobs were still running when checked and should not be confused with the completed functional-gate CI. The native ARM64 carry kernel is now the accepted default. Shared fused and portable kernels remain explicit controls. No new CMake flags or nonarchitecture source files are added; no FloatMP end-to-end speed claim.

### B11 experimental baseline x86-64 carry-chain assembly candidate (18 October 2026)


In the existing architecture-only `src/big_unsigned_x86-64.cpp`, B11 adds `mul_low_into_carry_x86_64`, an **opt-in** row-schoolbook with explicit GNU/Clang x86-64 `MULQ`, `ADDQ`, `ADCQ`: each product computes the exact 128-bit total `a[i]*b[j]+output[i+j]+carry`, using RDX:RAX for the product, then two ADD/ADC pairs to propagate the two carries into the high word. GCC/Clang extended asm ties the low input/output to `+&a` (RAX), the high output to `=&d` (RDX), keeps other operands in nonoverlapping registers, and clobbers condition codes. No BMI2, ADX, CPU feature probes, scratch allocations or new translation units. MSVC retains the existing B3 native schoolbook fallback; non-x86 architectures retain portable fallback. This candidate does not modify any production dispatch, including ARM64 carry and GCC x86-64 fused.

**Correctness approach:** Locally on Linux x86-64, the exact isolated carry instruction sequence was compiled with GCC 14.2 and Clang and matched unsigned `__int128` for 1,500,000 random/edge tuples, including all-ones. This isolated check is **not** a substitute for project tests. The integrated benchmark adds `--include-x86-64-carry`, which prints optional `prealloc_x86_64_carry` after existing modes; it is accepted only on GNU/Clang x86-64 builds, preserving frozen defaults and measurement ordering.

**Performance gate:** No integrated x86-64 measurements of this kernel are yet available and no claim of speedup is made. The x86-64 selected production path is unchanged pending repeated same-host evidence, including public owning results if ever promoted. Maintain architecture-only additional .cpp file policy.

### B11.1 separate portable fallback translation unit (18 October 2026)


Per user-requested three-backend source organization, the **portable C++20 limb multiplication** (`multiply_wide` 32-bit decomposition and `multiply_limbs` row-schoolbook) was moved **without changing the arithmetic** from `src/big_unsigned.cpp` to `src/big_unsigned_fallback.cpp`. `src/big_unsigned_x86-64.cpp` still holds all x86-64-native B3/Comba/experimental `MULQ/ADCQ` kernels and `src/big_unsigned_arm64.cpp` the selected ARM64 `MUL/UMULH/ADDS/ADC` carry kernel. These are the three backend-specific .cpp files. `src/big_unsigned.cpp` stays shared for representation, additions/subtraction, full product size checks, public dispatch, owning helpers, Karatsuba, and the one shared GCC/Clang `__int128` fused implementation (which must **not** be duplicated between x86 and ARM as noted in B10.1). Its no-`__int128` path now calls `mul_low_into_portable` in the fallback TU. CMake always adds `big_unsigned_fallback.cpp`, and uses one `if()/elseif()` branch to include **at most one** native backend and corresponding macro for supported architecture/compiler. Unsupported targets use only the fallback. Neither production selectors nor benchmark/test APIs nor data allocation/lifetime change. Existing x86 GCC fused switch continues unchanged; no new flags.

Do not create any algorithm-specific extra TU, including Karatsuba: only `big_unsigned.cpp`, `big_unsigned_fallback.cpp`, `big_unsigned_x86-64.cpp`, and `big_unsigned_arm64.cpp` are expected for BigUnsigned implementations.

### B11 Ubuntu x86-64 integrated benchmark: MULQ/ADCQ loses to GCC fused (19 October 2026)


The precise same-invocation medians (ns) are:

The explicit x86 carry kernel time relative to GCC fused selected is +10.7%, -0.7%, +6.6%, +7.8%, +23.3%, +6.1% across these widths (positive means **slower**). 128 bits is effectively tied and no general speedup exists. This suggests that higher-level multiplication algorithms/crossover (already-experimental Karatsuba) may be a more promising direction at large sizes than forcing a baseline MULQ/ADCQ instruction sequence, but requires new integrated data.

**Decision:** B11 explicit x86-64 `MULQ/ADDQ/ADCQ` stays **opt-in only**; no promotion to GCC x86-64 production dispatch. ARM64 production carry remains selected and unchanged. Possible causes include GNU/Clang inline-asm register restrictions and compiler optimization/scheduling differences; no causal claim has been proved from instruction-level profiling. CI for B11 kernel commit `7c3f9a7` and subsequent B11.1 fallback-refactor `d86c877` is fully green (Unix, Windows, Coverage each). The user supplied performance numbers, but did not explicitly report a separate local full `ctest` outcome for this new x86-64 build; do not infer full local test status from the benchmark alone. No new source files, ISA flags, CMake switches, or production dispatch changes are justified by these measurements.

### B12 experimental Karatsuba write-volume optimization (19 October 2026)


Existing B8.1 Ubuntu GCC x86-64 integrated cutoff-sweep data show cutoff=32 was the most promising measured candidate (Karatsuba ~2629 ns at 4096 bits and ~8462 ns at 8192 bits), with an estimated balanced-operand crossover between 2048 and 4096 bits; the observation is historical, not an independently reproduced result for this change. This B12 experiment modifies **only** the experimental Karatsuba implementation in the shared `src/big_unsigned.cpp`: (1) `add_shifted` adds the actual cross-term length and propagates any remaining carry only until it clears, rather than scanning every high destination word; it still validates all discarded high cross-term words, preserving overflow checks; (2) recursive internal nodes no longer zero the output before two child products, since `should_split` implies both low halves have exactly `half` limbs and `z0` and `z2` entirely cover disjoint portions of the exact output; (3) the public full-product Karatsuba clears only extra caller destination limbs, leaving all exact product limbs to the recursive writer. Base cases still use the existing architecture-selected low multiplication. The algorithm, recursion splits, minimum cutoff 8, workspace sizing, signatures, and production dispatch do **not** change.

Additional deterministic all-ones/randomized limb patterns test odd and asymmetric sizes including 17x17, 31x33, 63x67, 79x80, 127x129, 255x257, and the non-splitting 1x180 and 8x8 cases; output is prefilled nonzero and oversized, and scratch is bounded by sentinel limbs on either side. This is correctness coverage, not artificial coverage metric inflation. Recheck source build and Unix/Windows/Coverage CI and the full test suite; project CTest has not yet been run on a user's host for this experiment. Neither correctness nor performance improvement is asserted before validations. No new `.cpp` files, CMake flags, or other architecture changes. Do not select Karatsuba automatically from unverified new timings.

### B12 Ubuntu GCC x86-64 integrated Karatsuba cutoff=32 benchmark (19 October 2026)


The user ran the post-B12 `benchmark_big_unsigned_multiply --include-karatsuba --karatsuba-cutoff 32 --bits 2048,2560,3072,3584,4096,8192 -t dark` in their native Ubuntu x86-64 release build. Selected production kernel remains fused GCC. Median ns for the **same invocation**:

Relative to same-run fused preallocated, Karatsuba time is +0.72%, -17.72%, -9.91%, -12.11%, -14.58%, and -33.29% respectively (negative is lower/better). Observed crossover for these **balanced** operands on this one host lies **between 2048 and 2560 bits**, but is not established for asymmetric operands, other targets, or end-to-end owning `mul`. B12's B8.1 historical cutoff=32 measurements for 2048/4096/8192 bits were 815.782/2629.425/8461.684 ns (different invocations). Post-B12 values 764.395/2582.760/8264.658 are somewhat lower, but the differences are small enough compared with cross-run noise that **B12's own optimization effect is not proved** without controlled repeated A/B measurements.

User has not separately reported native `ctest` results for this commit. The current opt-in `prealloc_karatsuba` benchmark preallocates reusable output **and scratch** outside timing and precomputes scratch sizing, while public owning `mul` allocates its own result and does not currently select Karatsuba. Preserve ARM64 carry and x86 fused production dispatch and existing source structure; no additional algorithm-specific .cpp files.

### B13 BigUnsigned multi-operation integrated benchmark (19 October 2026)


Per the agreed sequencing, stop optimizing multiplication and measure already-implemented unsigned operations before designing division. Added `benchmarks/benchmark_big_unsigned_arithmetic.cpp` to the **benchmark-only** executable list in `benchmarks/CMakeLists.txt`, leaving production sources, runtime dispatch, `BigUnsigned` APIs, CMake architecture selection, and all old multiplication benchmark cases unchanged. This is an executable benchmark, **not** an extra algorithm/backend source file. Default sizes 53/128/256/1024/2048 bits, 20,000 iterations, 31 samples; `--bits`, `--iterations`, `--samples`, `--summary-only` and Ariadne logging controls mirror the established integrated benchmark CLI. Deterministic equal-bit-size operands with nonzero high bits; each pair is ordered so unsigned subtraction cannot underflow.

As with the multiplication benchmark, fixed deterministic inputs and benchmark order should be interpreted cautiously; later add asymmetrical/zero/equal/carry-saturated cases for performance modeling and richer workloads as needed. Build/CI and native host measurements remain pending.

**Next implementation step only after baseline measurement:** define and test generic unsigned `divrem` in shared `BigUnsigned`; start with divisor one limb and progress to multiprecision, preserving exact quotient/remainder identity and explicit zero-divisor errors. Then `gcd` and decimal parsing/output, then new signed-magnitude `Integer` integration.

### B14 initial single-limb unsigned division in shared core (19 October 2026)


**Source-location rule clarified by user:** `big_unsigned.cpp` implements generic algorithms for which there is **no native ISA-specialized version**. `big_unsigned_fallback.cpp` contains the generic fallback **only for an algorithm that also has an x86-64 and/or ARM64 specialized implementation** (currently portable schoolbook multiplication). Native `big_unsigned_x86-64.cpp` and `big_unsigned_arm64.cpp` contain their own specialized kernels. The experimental single-limb division was therefore placed exclusively in the shared `src/big_unsigned.cpp` instead of `big_unsigned_fallback.cpp`, and the fallback file was **left unchanged**. This rule supersedes any prior ambiguous statement that all portable primitives belong in fallback; it also applies to prospective generic `divrem`, `gcd`, and conversion functions until genuine architecture specializations exist.

B14 adds exact arbitrary-length unsigned division by a **single nonzero 64-bit limb**: `BigUnsigned::divrem_limb_into(dividend_view,divisor,quotient_view)` returns a remainder and fully overwrites the quotient view (zeroes extra words), with no heap allocation and a disjoint-views contract. The owning free function `divrem_limb(BigUnsigned const&,Limb)` returns a normalized quotient/remainder pair. Invalid divisor 0 raises `std::domain_error`; too-short quotient raises `std::length_error` before mutating output. The implementation uses a pure C++20 radix-2 bitwise division of the 128-bit remainder/digit pair (no `__int128` or compiler-ISA dependencies), maintaining remainder < divisor each bit, including when a logical 65th bit is set. General multilimbed-divisor division is **not** yet implemented.

The B13 add/sub/compare/shift measured baselines remain as supplied by user; there is no justified microoptimization of them in this commit. After CI and native project tests pass, continue with multi-limb `divrem`, not premature Integer migration.

### B14.1 Windows single-limb independent exact-reference fixtures width fix (19 October 2026)


A deterministic case also explicitly checks a nonzero remainder above 32 bits. No production implementation, kernel selection, API, build flags, or benchmarks change. Re-run Unix, Windows and Coverage CI and report their actual outcomes before closing the Windows validation gate.

### B15 common multi-limb unsigned division (19 October 2026)


After the user-confirmed B14.1 Windows green CI, added exact unsigned `divrem` to the shared `src/big_unsigned.cpp`, leaving native multiplication and portable ISA fallback unchanged. Owning `divrem(a,b)` returns canonical nonnegative quotient and remainder satisfying `a=b*q+r`, `r<b`. The no-allocation `divrem_into` accepts pairwise-disjoint canonical limb inputs, output buffers and scratch; it fully overwrites outputs, zeroes extra limbs, validates sizes before mutation, and throws on a zero divisor. For normalized lengths n>=m>=2, `divrem_scratch_limbs` returns n+2*m+2 with checked size arithmetic; other cases need zero scratch.

Algorithm: normalize the multi-limb divisor to set its highest bit, estimate each base-2^64 quotient digit using the existing portable (high:low)/limb division, cap at one limb, compute the full product with the existing selected multiplication kernel, correct the estimate downward while its product exceeds the current dividend window, then subtract and unnormalize the final remainder. For single-limb divisors reuse B14; for smaller dividends return zero quotient and input remainder. No new backend source or ISA assumption is introduced.

**Project CTest and Unix/Windows/Coverage CI are pending** for this commit; the standalone test is not a substitute. Future work: gcd and decimal import/export, then signed-magnitude Integer, only after production validation.

### B15.1 shared subtract helper collision fix (19 October 2026)


B15 failed to compile on AppleClang ARM64 with `-Werror`: the new division-local `subtract_in_place(Write,Read) noexcept` conflicted with the existing later Karatsuba `subtract_in_place(Write,Read)`, which can throw `std::logic_error` if the subtraction's nonnegative precondition is violated. Both helpers were in the same unnamed namespace and had the same argument types; the exception-specification mismatch caused the diagnostic. This is a compile-time regression, **not** evidence that the division algorithm returned an incorrect value.

Remove the redundant division-local implementation, declare the shared helper before the division routine, and reuse the existing checked implementation for both quotient-product subtraction and Karatsuba. This preserves division's exact correction loop and its nonnegative precondition while keeping the existing borrow guard. The previously described standalone prototype checks do **not** verify the integrated B15.1 fix. Re-run Unix/Windows/Coverage CI and host CTest before claiming integrated correctness.

### B16 BigUnsigned GCD and strict decimal conversion (19 October 2026)


B15.1 Unix, Windows, and Coverage workflows passed at `feb15cc`. GCD uses Euclidean remainder with canonical zero rules `gcd(0,0)=0`, `gcd(x,0)=x`. Parsing accepts only a nonempty ASCII unsigned decimal sequence, including leading zeroes; signs, whitespace, prefixes and nondigits raise `std::invalid_argument`. Output always uses canonical digits, with zero as `"0"`.

The decimal parser consumes up to nine digits per block and multiplies existing radix-2^64 limbs by at most 10^9 using checked 32-bit half-product arithmetic, without a required 128-bit compiler extension. Decimal output repeatedly divides a temporary mutable limb copy by 10^9 using the existing portable wide-limb divider, accumulating and zero-padding nine-digit blocks. These are correctness-first baselines, not competitive performance claims.

**Integrated Unix/Windows/Coverage CI remains pending** on this new commit; the isolated harness is not a project-level result. This concludes the planned prerequisite primitive set for the initial signed-magnitude Integer design, but does not constitute Integer migration.

### B19 decimal format speed and representative division benchmark (19 October 2026)


Native B18 benchmark measurements on macOS ARM64 and x86-64 (20,000 iterations, 31 samples, 53/128/256/1024/2048/4096 bits) verified 384 paired input values and 2304 paired operation results. These were measurement issues, not demonstrated arithmetic defects.

B19 makes the **targeted common-only** decimal-format change in `src/big_unsigned.cpp`: while dividing the copied limb vector by `10^9`, split each 64-bit limb into high and low 32-bit halves and perform two ordinary 64-bit quotient/remainder steps, maintaining the invariant `0<=remainder<10^9<2^32`. Each `(remainder<<32)|half` is bounded below `2^64`; this is portable C++20 on Windows LLP64 and AppleClang/GCC/Clang without `__int128` or new ISA/backend sources. The existing arbitrary-divisor `divrem_limb` and `divrem` implementations remain untouched. An isolated Clang warnings-as-errors differential check compared 800,000 high/low inputs and partial remainders against a 128-bit reference and passed; this does not establish an integrated build or performance improvement.

Existing original quotient/remainder rows and the original 32-value pools remain unchanged, while two additional `quotient_nontrivial`/`remainder_nontrivial` rows use a separate deterministic 32-value signed division pool with single-limb, half-width, near-equal and small divisors. Both backend pools are constructed outside timing. No changes to Rational benchmark paths, CI settings, timed sample counts, or data exclusions. **4096 bits remains supported and untouched**; for the next provisional benchmark invocation select `--bits 53,128,256,1024,2048` explicitly for speed, then restore 4096 after measuring `to_decimal` improvements.

Integrated Unix/Windows/Coverage CI and native ARM64/x86-64 benchmark results for B19 are pending. Keep `SignedInteger` temporary and migrate its unique tests before final deletion, independently of benchmark throughput.

### B19.1 scoped type fix (19 October 2026)


AppleClang reported two unqualified `Limb` uses in free function `Ariadne::to_decimal`. The function was already in `namespace Ariadne`; `Limb` is a nested alias of `BigUnsigned`. Qualify both as `BigUnsigned::Limb`. No other source, test, benchmark, algorithm or namespace is changed. Integrated CI remains pending.

### B19.2 equal-prefix compare improvement and accurate benchmark verification count (19 October 2026)


The user supplied B19.1 native ARM64 and x86-64 benchmark runs (20,000 iterations, 31 samples). The shared `to_decimal` optimization reduced 2048-bit conversion from 111251 to 5113 ns on ARM64 (21.76x) and from 240108 to 7589 ns on x86-64 (31.64x). The x86-64 run also included 4096 bits: conversion dropped from 952785 to 27951 ns (34.09x). These are same-user before/after runs, not simultaneous A/B trials.

`BigUnsigned::compare` repeatedly invoked `limb()` for equal high prefixes, including repeated bounds and inline-vs-heap selection. B19.2 changes only this shared implementation to acquire both `LimbView` spans once and compare their limb values directly, preserving length-first comparison, signed ordering, canonical rules and `noexcept`. No ISA-specific source or kernel selection changes. An isolated C++20 compare-view implementation passed 250000 deterministic reference comparisons when compiled with both Clang and GCC, with warnings-as-errors. That does not replace project CTest or establish an observed speedup.

The new benchmark division-input oracle in B19 checked all 32 additional input pairs inside the pre-existing 32-element input loop, incorrectly reporting 10560 distinct inputs for five widths or 12672 for six. B19.2 moves that validation loop outside the old input loop, preserving exactly the same comparisons, with no weakening and no change in timed measurements. The accurate counts are 640 distinct input values for five widths or 768 for six, with unchanged 2240 or 2688 per-operation result checks respectively. This is not a performance change.

The B19.1 Unix/Windows/Coverage CI was green at `a7b7bd1`, and its native benchmark outputs passed numerical verification. The B19.2 full project and cross-platform CI plus new ARM64/x86-64 performance results are **pending**. Keep 4096 available but use `--bits 53,128,256,1024,2048` for interim repeats.

### B20 selected native 128/64 unsigned division and independent exact-reference fixtures (19 October 2026)


B20 changes **only** the shared 128/64 digit division selector inside `src/big_unsigned.cpp`. The invariant for `divide_wide(high,low,divisor)` is `divisor!=0 && high<divisor`, so both quotient and remainder fit in 64 bits and `0<=remainder<divisor`. If `high==0`, plain C++ 64/64 quotient and remainder suffice. With GNU/Clang `__SIZEOF_INT128__` the nontrivial path computes `(uint128(high)<<64)|low` and divides directly, suppressing only the accepted compiler-extension `-Wpedantic` diagnostic; with genuine MSVC x64 it uses the Visual Studio 2019+ `_udiv128` intrinsic from `<immintrin.h>`. Other compilers and architectures retain the unchanged 64-iteration pure-C++20 `divide_wide_portable` fallback, annotated `[[maybe_unused]]` to avoid spurious warnings when native operations are compiled. There is no new backend translation unit or ISA assembly, no changes to long-division correction / remainder algorithms, no new public API and no functions outside `namespace Ariadne`. Both `divrem_limb_into` and `divrem_into` call this one selected helper. The existing `to_decimal` fast radix-`10^9` path remains unchanged. Quotient and remainder at `Integer` level continue to share `divrem` and the existing signed truncation convention.

Outside the project, the selected GNU/Clang implementation was compiled with GCC and Clang at `-std=c++20 -O3 -Wall -Wextra -Wpedantic -Werror` and compared against the original bitwise fallback over **1,200,252** exact 128/64 cases per compiler (1,200,000 deterministic random cases plus 252 boundary cases), all passing. This is an isolated check, not evidence of project CI, of MSVC intrinsic correctness or of speedups.

The previous 6000+6000 random oracle cases and checked preallocated-buffer invariants remain untouched. No benchmark corpus, timed loops, sizes or special mode switches are changed.

### B22 balanced large-integer decimal formatting (19 October 2026)


B22 optimizes only `BigUnsigned::to_decimal` in the shared `src/big_unsigned.cpp`. Preserve B19's exact radix-`10^9`, two-64/32-division formatter as `to_decimal_blocks` and select it unchanged for all values of at most **16 64-bit limbs** (including zero, signed-magnitude single limb, the inline/heap boundary, and every 1024-bit input). For larger inputs, build a private per-call table `P_i=10^(9*2^i)` by exact `BigUnsigned::square`, caching the exact decimal widths `9*2^i`. Avoid constructing the expensive final unnecessary square when its minimum bit length already exceeds half the original input's bit length. Recursively split a positive value at the largest `P_i` whose bit length is at most half the current value's bit length, using the existing exact `BigUnsigned::divrem(value,P_i)`. Convert the high and low parts recursively, left-padding the low decimal string with exactly `9*2^i` digits before concatenation. Because `0 <= low < P_i`, the suffix never exceeds the fixed width; leading high zeros are absent, and low zeros/long internal zero runs are exact. Leaves of at most 16 limbs use `to_decimal_blocks`. Zero is returned as `"0"`. The helper functions are in an anonymous namespace **inside** `namespace Ariadne`, never in the global namespace. This algorithm reuses native divrem and should be evaluated on both architectures for speed rather than assuming a crossover.

B22 integrated cross-platform CI and the actual ARM64/x86-64 performance benchmarks are pending.

Do not optimize or alter the benchmark corpus to manufacture improvements. Keep 4096 bits included.

### B23 specialized single-output unsigned division (19 October 2026)


Starting with the B22 native decimal formatting baseline, split only the *result ownership* of unsigned division while retaining the exact normalized radix-2^64 quotient estimation and correction algorithm. Inside `BigUnsigned`, a private compile-time-selected `divrem_selected_into<WantQuotient,WantRemainder>` shares the same normalization, quotient-digit correction, and working scratch.

New native magnitude `quot(BigUnsigned,BigUnsigned)` and `rem(BigUnsigned,BigUnsigned)` produce only the requested result. The quotient buffer has only `n-m+1` limbs where `n>=m`, rather than `n`; remainder-only allocates no quotient object and skips all quotient output stores, while quotient-only allocates no remainder object and skips final remainder unnormalization. For a one-limb divisor, remainder-only consumes digit divisions without allocating or storing quotient limbs. A dividend shorter than its divisor yields canonical zero quotient or the untouched canonical dividend as remainder. Zero divisor still throws. All helpers/functions remain inside `namespace Ariadne`. `Integer::quot` and `Integer::rem` now call the corresponding single-result unsigned primitive, preserving truncation toward zero, dividend-signed remainder, and canonical nonnegative zero.

Full Unix, Windows, Coverage CI and ARM64/x86-64 performance measurements remain pending; compare B23 against B22 at all 53..4096-bit widths with 20,000 iterations and 31 samples. Inspect all 4 quotient/remainder workloads.

### B23.1 benchmark operation naming clarification (19 October 2026)


Following the user-requested clarification, rename the four integer division benchmark *labels only*: historical `quotient` to `quotient_trivial`, historical `remainder` to `remainder_trivial`, `quotient_nontrivial` to `quotient`, and `remainder_nontrivial` to `remainder`. The original 32-pair corpus satisfies `|a|<|b|`, hence its quotient is zero and remainder equals the dividend; the additional 32-pair division corpus includes small divisors, half-width divisors and near-equal magnitudes, and produces nonzero quotients. Keep every operand, iteration, sample, operation order, verified result, operation count and timed routine unchanged: only report names differ. Old and new benchmark output rows must be matched by this four-way mapping before comparisons. The current B23.1 naming-only commit requires fresh CI validation, but no numerical performance gains are claimed from its label change.

### B24 wide-radix decimal construction (19 October 2026)


After B23.1's report-label-only update, GitHub Unix, Windows, and Coverage CI all succeeded on `68bc5e9`.

B24 changes only native `BigUnsigned::from_decimal` in `src/big_unsigned.cpp`. Preserve the complete B19 base-`10^9` split-32-bit implementation and the shared ASCII digit validation for all strings of **309 characters or fewer**, and for all compiler targets that do not support `unsigned __int128`. For longer strings on GCC/Clang targets defining `__SIZEOF_INT128__`, accumulate decimal chunks of up to **18 digits** (radix up to `10^18`) using one exact 128-bit intermediate per 64-bit limb: `product=(Wide)digit*radix+carry`, `digit=(Limb)product`, `carry=(Limb)(product>>64)`. Because `radix<=10^18<2^64` and `carry<radix`, no product overflows 128 bits. Both the first possibly short chunk and subsequent full chunks preserve leading zeros and all limb carries. Helpers and all modified production functions are inside `namespace Ariadne`. The 309-digit cutoff leaves the complete <=1024-bit benchmark pool on the old implementation.

Its timings, particularly ~1.8x around 617 digits and ~2.7-3.0x around 1234 digits, are isolated **proxies**, not project-integrated ARM64/x86-64 speedups. Keep the benchmark corpora, 14 operations, sample count and architecture test fixtures unchanged; measure actual improvement on each host.

Add four malformed >309-digit cases spanning first/middle/last invalid character positions. Existing 2800 randomized magnitude oracle cases, 1800 parsing oracle cases, 50 B22 split boundaries and signed Integer tests are untouched. Unix/Windows/Coverage CI and native ARM64/x86-64 performance results for B24 remain pending until the commit is validated. Compare B24 vs B23 using unchanged six sizes, 20,000 iterations, 31 samples; pay special attention to `construct_decimal` at 2048 and 4096 bits and ensure no regressions in `to_decimal` or division.

### B25 wide-radix large-leaf decimal formatting (19 October 2026)


Building on B22 balanced decimal formatting and B24's successful wide-radix *parsing*, B25 optimizes the leaf formatter inside `BigUnsigned::to_decimal`, without changing the B22 split-power recursion, quotient/remainder, signed Integer behavior, or benchmarks. For magnitudes of **9 limbs or more** on GCC/Clang targets with `__SIZEOF_INT128__`, use an unsigned-128-bit 128/64 division by the constant `10^18` to process one 64-bit input limb per decimal block pass. A 64-bit running remainder is always below `10^18` and the combined `(Wide(remainder)<<64) | limb` fits in 128 bits; each quotient limb is below `2^64`. Collect 18-digit base-`10^18` remainders and pad all nonleading blocks to 18 digits. This should reduce scans relative to B19's 9-digit method. The original two-64/32-division base-`10^9` algorithm is unchanged for **up to 8 limbs** and on compilers without native unsigned 128-bit integers, including MSVC. The B22 recursive cutoff remains 16 limbs; large values therefore gain the wide leaf inside the existing split recursion. No public API changes, external dependencies, global caches, floating-point estimates, or new architecture-specific assembly.

The implementation's fast-path selection occurs in the existing private `to_decimal_blocks` helper; `to_decimal_blocks_wide` is an internal helper in the anonymous namespace nested inside `namespace Ariadne`.

They justify trying B25 but do not establish actual acceleration. The output fast path differs from B24's parsing fast path. Compare B25 with B24 separately on native ARM64 and x86-64, particularly `to_decimal` at 1024/2048/4096 bits, other operations, and aggregation. Run Unix/Windows/Coverage CI; all B25 integrated validation and production speedups are pending until verified.

### B25.1 ARM64-only wide-radix output, no CMake flags (19 October 2026)


B25's full native benchmarks verified 768 independent paired inputs and 2688 paired operation results on both ARM64 and x86-64. However, B25's base-10^18 decimal output on x86-64 regressed strongly relative to B24: native `to_decimal` rose from approximately 2239 to 4212 ns at 1024 bits, 6725 to 9126 ns at 2048 bits, and 19731 to 23495 ns at 4096 bits. On ARM64, B25 improved the same three sizes from roughly 1590 to 1346 ns, 4556 to 4186 ns and 13138 to 12596 ns. These measurements justify an architecture-specific compile-time selection; do not use the isolated B25 leaf microbenchmark as a production speedup claim.

B25.1 retains the exact B25 base-10^18 unsigned 128-bit leaf implementation **only when compiler target macro `__aarch64__` is predefined** and the compiler defines `__SIZEOF_INT128__` and is GCC/Clang. It continues to use the B19 base-10^9 leaf for ARM64 magnitudes of <=8 limbs. For **x86-64 regardless of `__int128` support**, MSVC, and other targets, the existing portable radix-10^9 path is always selected. Both compile-time guards in the same `src/big_unsigned.cpp` implementation are updated. No CMakeLists.txt file is modified, no user-visible CMake option, compile definition, new feature flag or build-setting dependency is introduced: selection uses the compiler's standard target macro directly.

Do not remove any existing regression tests.

B25.1 Unix/Windows/Coverage CI results and native host benchmarks are pending. Expect the B25 to_decimal results to be retained on ARM64 and the B24 results to be restored on x86-64; these are expectations, not assertions of measured performance. Compare all six established bit sizes with the same 20,000 iterations, 31 samples and unchanged operation pool. Do not modify CMake for architecture selection.

### B26 x86-64 baseline DIVQ for exact 128/64 division (19 October 2026)


ARM64 at the same size was approximately 2.93x and 2.88x. In contrast, `pow3` was approximately 2.64x on x86-64 versus 2.68x on ARM64, so B26 targets division rather than modifying multiplication or exponentiation.

B26 changes *only* the internal `divide_wide` helper in `src/big_unsigned.cpp` for baseline GNU/Clang x86-64 (`__x86_64__`, not `_MSC_VER`). Under the already-established precondition `divisor!=0 && high<divisor`, place the two-limb dividend in `RDX:RAX` and execute a single `divq` with an explicit register input for the divisor. The quotient returns in RAX and remainder in RDX. Use `+a` and `+d` read/write constraints, plus the `cc` clobber; there is no memory access inside the assembly. The high=0 shortcut remains the normal 64-bit division. The precondition ensures that DIVQ cannot overflow its 64-bit quotient. This baseline instruction requires no BMI2, ADX or runtime CPU feature detection. MSVC x64 retains its existing `_udiv128` intrinsic; non-x86 platforms retain their `unsigned __int128` or bitwise portable fallback.

The change affects the expensive 128/64 quotient estimates and single-limb division, not trivial division cases where `|a|<|b|` and the helper is not invoked. It does not change `Integer::pow`, the multiplication dispatch, decimal input/output, benchmarks, build flags, or CMake in any way. No new CMake option, preprocessor build definition, instruction set requirement or ABI dependency is introduced.

The GCC generated instruction stream for the isolated helper contains `div` rather than the wide-division support call. These are local isolated *correctness and code-generation checks*; they do not establish performance on the repository's actual production benchmark.

The existing 6000 multiprecision differential cases continue to exercise all quotient-estimate correction paths. All tests are compiler-portable; no inline assembly in test sources.

Prioritize actual x86-64 `quotient` and `remainder` speedups, report changes in `quotient_trivial` separately, and confirm ARM64 stability. CI and integrated performance are pending at commit preparation.

### B27 fused multiprecision quotient-step subtraction (19 October 2026)


B27 replaces the normalized multi-limb division's separate `mul_into(v,qhat,product)`, full product-vs-window comparison, product correction and subtraction passes with a **single multiply-subtract pass** updating the normalized dividend in place. For each divisor limb calculate the exact 64x64 product plus the incoming high carry, subtract its low limb from the partial dividend and accumulate the subtraction borrow into the high carry. Subtract the final carry from the extra high dividend limb. If this indicates the initial quotient estimate was too large, decrement it and add the divisor back, propagating the addition carry; repeat only when the high limb still denotes a negative result. The proven 64-bit normalized quotient estimator can overestimate by up to two. No underestimated digit is introduced. This is the exact radix-2^64 long division recurrence, not a new approximation.

The product/borrow operations use `unsigned __int128` on GCC/Clang, MSVC x64 `_umul128` and explicit low-word carry handling, and a C++20 half-limb 32-bit product fallback on other toolchains. B26's x86-64 `divq` quotient estimator and its existing platform fallbacks remain unchanged. The architecture-independent fused division step uses **no inline assembly** and needs no additional compiler flags, CMake edits, ISA extensions or dispatch. Quotient-only, remainder-only and joint-output modes share the same implementation; B23 ownership and signed Integer truncation semantics are unchanged.

The exact public scratch requirement for n dividend limbs and m divisor limbs, n>=m>=2, shrinks from `n+2*m+2` to **`n+m+1`**, retaining only the normalized dividend (n+1) and divisor (m). There is no temporary product slice. Required scratch for single-limb divisors or dividends smaller than the divisor stays zero. All pre-mutation output and scratch checks remain in place. Document the new scratch contract in the public header.

Cases covered magnitudes up to 128 limbs, 64 distinct normalization shifts, saturated limbs, high-zero inputs, near-equal operands, and quotients requiring one or two add-back corrections. The oracle recorded **7930 two-correction events** and no case needing more than two. This validates the independent local core algorithm, not yet the whole repository build.

The existing 6000 randomized multiprecision and 64*10 normalization cases are left intact. Benchmark corpus, 14 operation labels, six widths, 20000 iterations and 31 samples must not change. B27 GitHub Unix/Windows/Coverage CI and native ARM64/x86-64 integer benchmarks remain pending until run, and no integrated performance improvement is claimed yet. Compare B27 with B26 separately on each architecture, especially general `quotient` and `remainder`, and check `quotient_trivial` for run-to-run noise.

### B28 exact equal-leading-limb quotient 0/1 fast path (19 October 2026)


The B27 benchmark's 32 meaningful division operands consist of 16 single-limb small divisors, eight approximately half-width divisors and eight near-equal dividend/divisor pairs. An independent GCC x86-64 probe also found that replacing one baseline `divq` with two 32-bit `divl` instructions for small divisors was slower; that candidate was **not** integrated.

For normalized n-limb magnitudes a,b with n>=2 and the same most significant limb h, write a=h*B^(n-1)+a0 and b=h*B^(n-1)+b0, where B=2^64 and 0<=a0,b0<B^(n-1). Since h>=1, a<2*b. Consequently floor(a/b) is exactly 0 or 1. B28 adds an early return inside `BigUnsigned::divrem_selected_into` after all existing pre-mutation output and scratch checks, and after the single-limb case. Quotient-only scans from the highest limb and returns 0 for a<b or 1 otherwise without subtraction. For remainder-only and joint-output, one unsigned multi-limb subtraction obtains a-b and its final borrow. A zero borrow produces quotient 1 and remainder a-b; a nonzero borrow produces quotient 0 and remainder a. All output spans remain zero-initialized as before. The shared B27 normalization/estimation/fused subtraction and correction path is unchanged for other operands. No ISA feature, external dependency, CMake option or public API is added, and the exact scratch requirement n+m+1 is deliberately unchanged even though this fast path does not touch it.

Verify all output modes and check that public output tails are cleared and scratch guards are unmodified. Retain all previous B27 and earlier regressions; benchmark corpus, operations, iterations and samples must remain unchanged.

**Performance evidence is isolated, not integrated:** an x86-64 AMD EPYC 9V74 environment using GCC -O3, fixed operand pools and preallocated output/scratch views measured the near-equal 4096-bit subset at approximately 113.3 -> 32.3 ns (quotient only) and 133.3 -> 84.1 ns (remainder only). In the same local probe the full 32-case mix was approximately 581.4 -> 562.9 ns (quotient) and 607.0 -> 573.0 ns (remainder). Clang -O3 likewise improved the near-equal subset but full-mix remainder timings were variable (599.5 -> 618.9 ns in one run). The isolated test excludes owning Integer allocations, signed wrappers and the exact repository benchmark harness; its absolute numbers must not be compared against native user ARM64/x86-64 benchmark values. No integrated B28 improvement or ARM64 result is claimed.

### B29 selective zero-limb fused subtraction for sparse divisors (19 October 2026)


Disaggregating the **unchanged** 32-pair general division corpus reveals that eight half-width divisors of the form `2^(bits/2)+1+i` contain long zero-limb runs: their quotient/remainder steps dominate the isolated native cost despite accounting for only one quarter of operands. This workload is mathematically representative of sparse bignums, not a new fixture contrived to score better.

B29 keeps B27's fused normalized long-division recurrence and B28's exact equal-leading-limb quotient-0/1 shortcut unchanged. The shared product/subtraction helper now has two compile-time instantiations. For a zero limb of the normalized divisor, the sparse specialization subtracts only the incoming carry/borrow from the dividend limb, then propagates the resulting borrow: `window[i] -= carry; carry = old < carry`. This is algebraically identical to multiplying that zero limb by the current estimated quotient digit, including the high-product zero and carry rules. The dense specialization retains the exact original multiplication code on GCC/Clang `unsigned __int128`, MSVC `_umul128`, and the portable 32-bit half-limb fallback.

**Selection occurs once before the quotient-digit loop**, after B27 normalization: use the sparse specialization only if the divisor has at least eight limbs, there are at least eight quotient digits, and at least half its normalized limbs are zero. Instantiate the two quotient-digit loops using a C++20 templated lambda; do not add a per-digit mode branch to the dense path. This cutoff avoids replacing the dense kernel with a branch on each multiplication, which performed worse on small or dense divisors in probes. The preexisting output and scratch size checks, correction bound, public `n+m+1` scratch allocation, signing/rounding boundaries, benchmark input pool, six sizes, 14 operation labels, 20000 iterations and 31 samples are all preserved. No build/CMake flags, optional CPU instruction or new public API are introduced.

On the 32-pair corpus at 4096 bits, selecting a sparse multiply-subtract kernel reduced measured quotient/remainder work by roughly 10-30% across GCC 14 and Clang 17 runs, while isolated half-width sparse-divisor cases often improved by approximately 20-40%. Absolute isolated kernel timings are not the official owning `Integer` benchmark and can be sensitive to compilation/code layout; **no integrated B29 speedup is claimed**. The single-limb-preinverse alternative showed benefits only in some 4096-bit GCC runs and was not incorporated.

GCC -O3 and Clang 17 AddressSanitizer/UndefinedBehaviorSanitizer runs both passed. The integrated production test suite and native ARM64/x86-64 six-size performance gates remain outstanding; CI and Codecov are not being consulted at user request.

### B30 x86-64 radix-10^18 decimal leaves via B26 DIVQ (19 October 2026)


B25's attempt to divide an unsigned 128-bit numerator by 10^18 regressed on x86-64 in the integrated benchmark and was disabled by B25.1. B30 does **not** restore that implementation: it reuses the exact B26 hardware 128/64 `divide_wide` helper, which already executes baseline x86-64 `divq` under GNU/Clang. Each decimal radix-10^18 quotient step has numerator `(remainder << 64) | limb`, where `0 <= remainder < 10^18`; hence the hardware `RDX < divisor` precondition holds and both the quotient and next remainder are exact 64-bit values.

The leaf formatter `to_decimal_blocks_wide` now selects this existing helper only for GNU/Clang x86-64 with at least nine 64-bit limbs. Preserve B25.1's original `unsigned __int128` implementation on ARM64, B19's radix-10^9 implementation for up to eight limbs and for MSVC/other architectures, the B22 recursive split powers and cutoff of sixteen limbs, B24 decimal parsing, B27/B29 division, B28 quotient-0/1, all output contracts and all Integer/FloatMP/Rational/Dyadic semantics. No additional assembly, new ISA feature, CMake change, build definition or public API is introduced.

An isolated x86-64 AMD EPYC 9V74 probe measured median leaf times for 16-limb operands of approximately 2075 -> 1081 ns under GCC 14.2 and 2260 -> 1152 ns under Clang 17. These are algorithmic proxies, **not** the repository's native `Integer::to_decimal` production timing or proof of cross-host improvements. They address the prior B25 leaf-only-versus-integrated discrepancy more thoroughly but cannot settle it; require the real six-size 20,000 iterations/31 samples ARM64 and x86-64 benchmark before claiming an integrated benefit. The native comparator and benchmark corpus must remain unchanged.

Integrated cross-platform test execution and ARM64/x86-64 production timings remain pending. GitHub CI and Codecov are not being consulted by user request.

### B30.1 revert x86-64 wide decimal leaf after integrated regression (19 October 2026)


Numerical verification passed on both hosts (768 paired inputs, 2688 paired results per run). It falsifies the promising isolated B30 proxy performance result for GNU/Clang x86-64: native `Integer::to_decimal` times at 1024, 2048 and 4096 bits worsened respectively from B29's **2247.362, 5736.002, 14893.242 ns** to B30's **3715.747, 7599.686, 18242.587 ns** (+65.3%, +32.5%, +22.5%). All sizes up to 256 bits use the old formatter and were not expected to change. On ARM64, where the code path was unchanged, corresponding formatting values were **1407.481 -> 1346.456**, **3565.394 -> 3501.423** and **9505.252 -> 9339.346 ns**: the small variations do not establish a causal speedup.

A key diagnostic is that 1024 bits fits the B22 16-limb nonrecursive cutoff, so the pronounced regression is already present in B30's x86-64 decimal **leaf** and cannot be blamed on recursive division or power construction. The native B29 base-10^9 path divides 32-bit halves by a fixed radix, which a compiler can strength-reduce; B30's `DIVQ` on 10^18 is not necessarily cheaper on the user's x86-64 CPU despite halving the number of radix blocks. CPU differences and generated code may explain the discrepancy with the independent AMD EPYC proxy, but the real project benchmark is decisive.

**Decision:** B30.1 restores `src/big_unsigned.cpp` byte-for-byte from the proven B29 commit `3d425d8`. Thus on x86-64 the original 10^9 block formatter is selected again. ARM64 retains B25.1's 10^18 wide-division formatter, and MSVC/portable targets retain their original B19 path. All B29 sparse division, B28 exact quotient-0/1, B27 fused arithmetic, parsing, recursive decimal splitting, public APIs, build flags and benchmark inputs remain unchanged. No performance special casing of benchmark inputs is introduced.

The B30 integrated functional benchmark passed prior to the revert; B30.1 restores the earlier B29 code, but a new integrated run of the B30.1 tree and its unit tests remains outstanding. Do not describe the B29 benchmark timings as new B30.1 measurements. CI and Codecov remain deliberately unconsulted at user request.

### B31 pairwise decimal string assembly and allocation reduction (19 October 2026)


B30.1 re-established the reliable B29 x86-64 radix-10^9 leaf after B30's official native 1024/2048/4096-bit regression. Profiling the reconstructed **exact B30.1 leaf** on an x86-64 AMD EPYC 9V74 with GCC 14 and Clang 17 showed that GCC strength-reduces the constant 1e9 divisions to integer multiplication and shifts, rather than executing a runtime hardware DIVQ. Replacing the divisor is not the current target.

B31 retains that B19 arithmetic, B22 recursive power/split logic, B25.1 ARM64 1e18 wide-radix path, B27-B29 division and all public APIs and signed semantics. It changes **only the rendering and allocation strategy** of `to_decimal_blocks` with radix 1e9: pre-reserve an upper bound of three base-1e9 blocks per 64-bit limb; render the highest block once via `std::to_chars`; allocate an exactly sized result string; write the lower blocks backward from least to most significant using four base-100 digit pairs and one leading digit. A local constexpr 200-character lookup maps 0..99 to fixed two-digit strings. This preserves exactly nine characters per lower block (including leading zeros) without the previous repeated `std::to_string`, `append`, and temporary-string operations. The result-length calculation is overflow-checked; the constexpr table and C++17 `std::to_chars` require no new architecture support, build flag, or CMake change. The ARM64 wide leaf is untouched; on ARM64 the radix-1e9 renderer is still used through eight limbs.

An exploratory alternating-order microbenchmark on the exact prior 1e9 extraction loop and the new renderer, with both random and the official near-power-of-two `2^(bits-1)+c` operand shapes, measured 16-limb GCC baseline ~2198-2206 ns versus candidate ~1769-1798 ns (roughly 18-20% faster). Clang ranged from ~2054-2131 ns to ~1719-1788 ns in short replicated runs (roughly 13-19% faster). Eight- and nine-limb leaves often improved ~25-34%; 32-limb stand-alone leaves improved ~6-11%, though production B22 uses the recursive path beyond 16 limbs. These figures isolate a leaf implementation, not production `Integer::to_decimal`: **no integrated B31 speedup has yet been shown**. Their limited transferability is especially relevant following B25/B30's misleading standalone wins.

Full repository tests and the ARM64/x86-64 integrated benchmark remain pending; GitHub CI and Codecov are deliberately not consulted at user request.

### B31.1 shadow warning fix (19 October 2026)


Rename the local counter constructing B31's constexpr digit-pair table from `value` to `pair_value` to avoid shadowing the `BigUnsigned const& value` argument in `to_decimal_blocks` under `-Wshadow`. Only an identifier changes: the 100 table entries, output bytes, arithmetic, benchmark workload and APIs remain identical. No new performance or test results are claimed; GitHub CI and Codecov remain unconsulted by user request.

### B32 selected GCC x86-64 symmetric squaring and generic pow identity elision (19 October 2026)


B31.1's repeated native ARM64 benchmarks confirmed its decimal conversion improvement and showed the first 256-bit `pow3` spike was measurement noise (181.7 ns falling to 124.5 ns). The remaining major 4096-bit costs are `pow3`, `mul` and `square`. Inspection of `Integer::pow` found that the first set exponent bit multiplied the initial accumulator `1` by a power of the base even though assignment is sufficient. The exponentiation schedule now copies that first factor and executes all subsequent squarings/multiplications as before. This works for every nonnegative Nat exponent, including exponent zero, zero base, negative inputs and even/odd sign; it does not special-case the benchmark's exponent three.

`BigUnsigned::square` previously called `mul(a,a)`, which multiplies every pair of limbs twice. B32 introduces a **GCC x86-64 fused-build-only**, allocation-free symmetric square at eight or more limbs. The result's existing 2n limbs hold the sum of products with i<j, calculated once for each unordered pair; one left shift doubles all cross terms; adding each a[i]^2 at positions 2i and 2i+1 completes the square with explicit carry propagation. The cross-term value cannot overflow 2n limbs because it is less than a^2. Per-row carry placement at i+n is safe: previous rows reach at most i+n-1. The full 2n-limb result cannot overflow, including saturated-limb inputs. GCC x86-64 chooses this path only when its existing `ARIADNE_BIG_UNSIGNED_GCC_FUSED` selection is enabled; turning off the existing CMake option still exercises the unchanged reference schoolbook kernel. ARM64 retains its verified assembly carry-chain multiplier, Clang its selected B3 schoolbook, and MSVC/portable their original paths. No new public API, scratch, ISA assumption, CMake option, benchmark workload or native dependency is introduced.

An independently built x86-64 GCC 14.2 -O3 probe compared the existing fused schoolbook recurrence against the symmetric implementation with fixed-size operand pools: at 8/16/32/64 limbs the times were approximately 67.7->50.6 / 232.2->148.6 / 929.2->531.5 / 3855.9->2044.8 ns. These correspond to 25/36/43/47% lower isolated kernel time. **Integrated B32 performance has not yet been measured** and the isolated x86 numbers must not be compared directly with the six-size Integer benchmark. This does not establish any ARM64 speedup from the symmetric algorithm, which is not selected there.

The integrated production test suite and user-host ARM64/ x86-64 benchmarks remain outstanding; GitHub CI and Codecov are intentionally not consulted.

### B35 selective existing Karatsuba for large balanced GCC x86-64 multiplication (19 October 2026)


B34.1 preserves the successful B33 alias-to-symmetric-square dispatch. Before introducing any new assembly kernel, B35 measured the existing B8 `BigUnsigned::mul_into_karatsuba` recurrence against the GCC-selected B6 fused schoolbook. The existing exact Karatsuba recursion uses caller-owned scratch, recursively computed sums, and the same architecture-selected schoolbook base cases; no new numeric algorithm is added.

At 32 limbs, cutoff 32 offered no useful speedup; at 64 limbs, cutoff 32 improved a dense-mixed preallocated kernel by approximately 12-15% and a benchmark-shaped near-power-of-two pool by approximately 16%. In repeated **owning** tests that allocate output and Karatsuba workspace each iteration, 64-limb sparse operands improved from approximately 3730 to 3148 ns (15.6%). At 96/128 limbs, sparse owning improvements were approximately 24%/30%. GCC -O3 and Clang ASan/UBSan both passed.

B35 therefore changes only **owning** `BigUnsigned::mul` on GCC x86-64 builds already selecting the existing `ARIADNE_BIG_UNSIGNED_GCC_FUSED` kernel, and only when the *smaller* operand has >=64 normalized limbs and the operand lengths are within 25% of the larger. The comparison uses `smaller >= larger - larger/4` to avoid integer overflow. Compute the exact scratch size via the existing `karatsuba_scratch_limbs(n,m,32)`, allocate one private vector, and call the existing exact `mul_into_karatsuba(...,scratch,32)` into the usual owned result. No global scratch, changes in output ownership, new ISA, CMake setting, public API, or benchmark-specific input tests. The original fused dispatch remains selected for 53/128/256/1024/2048-bit benchmarks and strongly unbalanced operands; 4096 bits corresponds to 64 limbs and can use Karatsuba. All ARM64 Apple M2/M4, Clang, MSVC, and portable architectures keep B33's original multiplication; GCC control builds with the fused option OFF also retain the previous kernel. B32 symmetric `BigUnsigned::square` and B33 aliased Integer self-products are unchanged; preallocated `BigUnsigned::mul_into` stays fused.

The original 6000 randomized product tests, B32 symmetric square cases, Karatsuba scratch/guards and other benchmarks remain unchanged. GitHub CI and Codecov were not consulted by user request.

### B35.1 Clang -Werror regression fix and integrated measurements (19 October 2026)


B35.1 changes **only** that B35 loop to `for (const auto& [n,m]: ...)`. The user's reported Clang diagnostic is the reason for the change; successful CI compilation remains to be verified on the new commit.

Two native x86-64 integrated B35 benchmark runs, with 20,000 iterations and 31 samples, measured 4096-bit `Integer::mul` at **2752.212** and **2750.796 ns**. B33 on the same architecture had measured **3065.983 ns**, so B35's 4096-bit multiplication benefit was approximately **10.2-10.3%**, and the two B35 measurements differ by 0.05%. B35's 4096-bit `square` (1830.722/1832.047 ns) and `pow3` (8142.558/8146.656 ns) remained stable. Both x86-64 runs and the ARM64 B35 benchmark reported numerical verification PASS with 768 paired inputs and 2688 paired operation results. Do not infer causality for timing drift in unrelated operations. Earlier ARM64 B35 host was Apple M2; never compare it directly with the Apple M4 B31.1 baseline.

At the user's explicit request of 19 October 2026, **GitHub Actions CI may be inspected again starting with this B35.1 commit**. Prior exclusions of CI checks applied to older commits; this permission does not itself establish a passing build. Existing numeric test coverage and correctness gates are preserved, and Codecov status should be described separately from compilation and functional CI.

### B36 exact high-limb quotient-0/1 extension (19 October 2026)


After B35.1 cleared Unix and Coverage CI (including Clang 16, AppleClang 21 and ThreadSanitizer), the initial Windows run was still active and later completed successfully. The user explicitly reauthorized checking GitHub Actions CI.

Skipping multiply/subtract for zero quotient digits is algebraically valid, but the 4096-bit half-width sparse-divisor corpus has no zero digits in its quotient. Do not select this candidate to manufacture benchmark speedups.

B36 instead extends the architecture-independent B28 exact 0-or-1 quotient shortcut for any nonzero canonical `n==m>=2` magnitudes with different high limbs: if `a.back()<b.back()` then `a<b` and quotient 0, remainder a. If `a.back()>b.back()` and `a.back()-b.back()<b.back()`, then `b<=a<2*b`, giving quotient exactly 1 and remainder one full n-limb subtraction. The subtraction's final borrow is zero in this branch. The subtraction-free quotient-only case directly writes 1; remainder-only does not create any quotient. Equal high limbs use the existing B28 scan/borrow algorithm with unchanged semantics, including a<b determined from lower limbs. All other cases use the existing B27/B29 normalized digit-estimation, sparse multiply-subtract and correction algorithm. The high-limb comparison avoids overflow from evaluating `2*b.back()`; no new CPU features, allocations, scratch, APIs, build flags or CMake changes are required. Shared logic preserves the exact full public scratch/output prechecks and tail initialization on all code paths.

An independent exact-arithmetic Python stress probe verified the high-limb classification against true multiple-precision integer division on 70,000 random equal-length inputs (2..128 limbs): 34,845 zero-quotient cases, 17,690 strict high-dominance one-quotient cases, and 17,465 fallback cases. This validates the shortcut's arithmetic classification, not compiled C++ integration or any runtime performance improvement.

All B28/B29, B35 and earlier cases are retained. The official 14-operation benchmark and all its operand pools remain unchanged. The half-width divisor inputs and B28 equal-high near-equal cases are already handled by their previous paths, so do **not** predict a substantial benchmark benefit from B36 without same-machine integrated measurements. Verify no regressions on Apple ARM64 M2 vs M2 and GCC/Clang x86-64. GitHub Actions may be checked; Codecov is not required to establish numeric correctness.

### B37 native ARM64 symmetric square for Integer self-products and powers (19 October 2026)


The corresponding `pow3` operation was 590.819/1892.640/7303.648 ns versus 352.962/921.167/2818.185 ns (1.67x/2.05x/2.59x). The official benchmark's `square` implementation is `a*a`, not the named `sqr` API. The power algorithm itself is already a binary method with the B32 first-identity-multiply elision; `pow3` requires one `square(base)` and one general product. Prioritize the underlying square rather than a benchmark-specific `pow3` shortcut.

Until B37, ARM64 `BigUnsigned::square(a)` called generic `mul(a,a)`, performing all n*n products; B32's specialized symmetric square dispatch was restricted to GCC x86-64. B37 adds `BigUnsignedDetail::square_into_symmetric_arm64` in the existing ARM64 TU. It computes only off-diagonal `a[i]*a[j]` for i<j using exactly the selected baseline `MUL/UMULH/ADDS/ADC` carry-chain assembly, doubles the cross-term limb array, and adds all diagonal squares with exact 128-bit unsigned accumulation and carry propagation. The sum is the mathematical square and fits exactly in 2*n limbs. There are n*(n+1)/2 wide products versus n*n for general multiplication. Inputs/output must be disjoint; the own-result wrapper provides the exact output size. No runtime CPUID, non-baseline ARM ISA, global storage, new API, CMake setting, ABI change or additional heap scratch is introduced.

Dispatch ARM64 AppleClang/GCC/Clang with selected native `ARIADNE_BIG_UNSIGNED_ARM64` only for magnitudes >=8 limbs; preserve the previously selected carry-chain schoolbook for smaller operands. On ARM64, extend the existing B33 self-alias `Integer::mul(z,z)` to call `BigUnsigned::square` for >=8-limb same-object operands, with a size-first branch; numerically equal but separate objects retain generic multiplication, including the usual sign semantics. All paths on x86-64 GCC/Clang, Windows/MSVC and portable targets are unchanged. The new kernel also accelerates `sqr` and internal `Integer::pow` squaring automatically.

The six-instruction AArch64 multiplication-accumulation assembly block also cross-compiled with Clang 17's `--target=aarch64-none-elf`, yielding the intended MUL, UMULH, ADDS and ADC instructions. These are local algorithmic and assembly-syntax checks, **not** a native Apple M2 performance measurement or complete project build.

Preserve B32/B33 and all previous tests. The full project test suite, GitHub Actions for B37 and same-M2 31-sample integrated performance comparison against B36 remain pending until run. Specifically gate on `square` at 1024,2048,4096 bits, `pow3` at all widths, `mul` as a control, and confirm no regression at 53/128/256 bits. No integrated speedup is claimed before those measurements. Keep M2 and M4 baseline results separate. The user has reauthorized GitHub CI inspections; report third-party Codecov statuses separately.

### B38 skip exact-zero rows in native ARM64 multiplication and symmetric square (19 October 2026)


The B37 Apple M2 integrated benchmark improved 4096-bit `square` from B36's 2298.379 to 1476.506 ns (35.8%) and `pow3` from 7303.648 to 6380.048 ns (12.6%). At 1024/2048 bits, `pow3` improved to 506.440/1667.752 ns. All B37 Unix, Windows and Coverage GitHub Actions later completed successfully (including AppleClang, Clang 16, GCC and ThreadSanitizer). The generic ARM64 product in `pow(x,3)` is an unbalanced product of a roughly n-limb operand and its roughly 2n-limb square. For values containing long zero-limb runs, both the native generic schoolbook kernel and B37's symmetric square were needlessly executing wide multiplications by exact zero. This is a general arithmetic property, not a `pow3` or benchmark-operand-specific recognition.

B38 preserves B37's six-instruction MUL/UMULH/ADDS/ADC accumulation for each **nonzero** outer-limb row but skips the entire row when `a[i]==0` in `mul_low_into_carry_arm64`. Since the output is zero-filled up front, a skipped row adds exactly zero; its final carry slot was never touched by earlier rows. The same exact-zero row elimination applies to `square_into_symmetric_arm64` for the off-diagonal upper triangle. It is safe because slot `output[i+n]` is not touched before row i. Doubling and diagonal square accumulation still process **all** limbs, including those skipped in the off-diagonal step. There is no division, dynamic input classification, new allocation or new ISA. This also preserves full, truncated, oversized-output and asymmetric `BigUnsigned::mul_low_into` contracts, including the reversed order of the operands; no special cases of `Nat(3)` or of the official operand pool were introduced. The generic multiplication API, `Integer::pow`, x86-64, AppleClang selection, CMake build flags and official benchmarks remain unchanged. A dense input incurs only an outer-row zero comparison, so check actual Apple M2 small and dense cases before accepting a performance claim.

All passed. On the same proxy, 64x128-limb sparse products reduced the *inner schoolbook work* from 64 rows to 2 and measured about 6471 -> 249 ns, while dense 64x128 results varied by about 1%. These are **x86 proxy** results, not actual ARM64 or production integrated `pow3` timings.

All previously tested saturated, no-zero, high/low-boundary, x*=x and power cases remain intact. The complete production build/CI and B38 same-host integrated benchmarks are **pending**; do not attribute the x86 proxy's large speedups to Apple Silicon. Compare B38 only against B37 on the same Apple M2, and use B37 x86-64 as an unchanged-backend control.

### B39 adaptive ARM64 dense/sparse schoolbook and selective owning Karatsuba (19 October 2026)


The same-M2 dense benchmark with cutoff-32 Karatsuba compared B37 and B38 directly. At 1024/2048 bits, B38's selected times were 154.600/576.502 ns versus B37 154.896/573.425 ns. Karatsuba cutoff 32 was slightly slower than selected schoolbook at 1024/2048 bits, but at 4096 bits it improved B38 preallocated time 2548.485 -> 1972.463 ns (22.6%) and B37 2284.963 -> 1964.965 ns (14.0%). The Karatsuba control reuses caller-owned scratch, so those timings do **not** establish an equally large speedup when the owning `mul` allocates scratch.

B39 addresses this with two exact compile-time instantiations of the **existing** native ARM64 MUL/UMULH/ADDS/ADC row-schoolbook recurrence. `mul_low_into_carry_arm64_impl<false>` has no zero-row guard, matching B37's dense inner loop. `...<true>` retains B38's exact-zero-row skip. A cheap linear prescan selects the sparse version only when at least `max(2,floor(active_rows/8))` of the active outer rows are zero; the scan stops once that count is reached. For active row counts <8 or multipliers <8 limbs, retain the B38 sparse implementation without any prescan. In full and truncated cases active rows are `min(a.size(),output.size())`. No reordered operand assumptions, allocations or altered alias/scratch/output contracts are introduced.

The B37 symmetric square gets the same template separation: the cross-product triangle is dense and branch-free unless the zero-row prescan selects B38's sparse version; doubling and the full diagonal update remain identical on both paths. The last input limb is excluded from the cross-row density count because it is never an outer cross row. The branch only chooses an equivalent exact algorithm, not a different mathematical result.

On ARM64 only, **owning** `BigUnsigned::mul` now dispatches to the existing generic B8 Karatsuba recurrence (cutoff 32) if both input magnitudes have >=64 limbs, are balanced within 25% of the longer magnitude, and each contains no more than floor(length/8) zero limbs. The test is early-exit for sparse operands, so the two-nonzero-limb 4096-bit official benchmark still uses sparse B38, and the approximately n x 2n general product inside `pow3` remains on the original ARM64 schoolbook path. Only this owning path allocates the exact private scratch; the preallocated `mul_into` path and the public scratch/output contracts do not change. Non-ARM64 implementations including the existing GCC x86-64 B35 Karatsuba dispatch remain unchanged.

Actual native M2 timings of the adaptive dense path and of owning Karatsuba, along with integrated pow3 results and CI, remain **unmeasured on B39**. Compare B39 dense 4096-bit `prealloc_selected` against B37 2284.963 ns, and new_selected against B37 2318.360 ns; do not conflate preallocated Karatsuba and owning timing. Confirm the sparse 4096-bit `mul`, `square` and `pow3` do not regress relative to B38 158.777/228.442/565.898 ns.

### B40 eliminate medium-size density prescans and Karatsuba leaf overhead (19 October 2026)


All three B39 GitHub Actions workflows completed successfully.

A separate same-M2 dense multiplication benchmark identified *new* B39 dispatch costs at moderate sizes. Compared with B38, dense preallocated selected products rose 154.600 -> **166.696 ns** at 1024 bits (+7.8%) and 576.502 -> **600.090 ns** at 2048 bits (+4.1%). For the 4096-bit dense case B39 recovered much of B38's regression: 2548.485 -> 2343.858 ns, while large owning products used Karatsuba and improved 2575.535 -> 2265.110 ns. Repeated B39 density classification in cutoff-32 Karatsuba's 32/33-limb recursive leaf products is also consistent with the measured preallocated Karatsuba changes from B38 to B39: 158.329 -> 170.012 ns (1024), 580.327 -> 603.987 ns (2048), 1972.463 -> 2065.662 ns (4096). These observations motivate but do not individually prove the source of all timing drift.

B40 conservatively **bypasses the B39 linear density prescan** in the ARM64 carry-chain multiplier whenever the active outer row count `min(a.size(),output.size())<64` **or** the multiplier contains fewer than 64 limbs. Those cases dispatch directly to the already validated B38 exact-zero-row implementation; the dense B37-like branch-free specialization is still available at >=64x64 limbs after B39's exact zero-density scan. This restores the B38 dispatch behavior in 1024/2048-bit dense full products, cutoff-32 Karatsuba's small recursive base cases, and ordinary low-truncated multiplication. At high sizes no density threshold, algorithm, instruction set or output contract changes.

For ARM64 symmetric squaring, there are only `n-1` off-diagonal outer rows; B40 bypasses classification at `n<=64`, invoking the B38 symmetric kernel directly at 4096 bits and below, including cases with nonzero dense diagonals. At larger sizes B39's adaptive choice stays in place. There are still two compile-time native MUL/UMULH/ADDS/ADC variants, no new allocations, no scratch changes, no new features or runtime flags. Dense balanced >=64-limb *owning* multiplication retains B39's cutoff-32 Karatsuba selection and density gates. Small, sparse or imbalanced owning products continue down the selected ARM64 schoolbook. The Karatsuba recursion itself is unchanged; only the base-case multiplier ceases to repeat an O(n) scan for n=32/33. Caution: this should not be described as a measured B40 speedup until native M2 data are collected.

Check owning multiplication in both orders, low products in both orders, the public symmetric square, and for large shapes also explicit cutoff-32 Karatsuba with caller-owned scratch. Cover dense, all-ones, moderately sparse, threshold-zero, and highly sparse limb patterns. Full production C++ CI and same-M2 dense/sparse benchmarks for B40 are pending; do not compare ARM M2 absolute timings to any other CPU.

### B42 single-allocation ARM64 radix-10^18 decimal formatting (19 October 2026)


B41.1 passes all Unix, Windows and Coverage GitHub Actions, including Clang 16 Release and ThreadSanitizer. B42 retains the exact existing recursive splitting and wide 128/64 quotient/remainder algorithm but optimizes formatting of the ARM64 radix-10^18 leaves: the decimal block vector reserves at most two blocks per input limb, then one exactly sized result string is filled using `to_chars` for the leading block and constexpr digit-pair lookup for all nine pairs in each remaining 18-digit block. No per-block `std::to_string` or temporary allocation is used. Other platforms' radix-10^9 path, all arithmetic kernels and the official benchmark are unchanged. B42 production CI and Apple M2 timings are pending.

### B43 append-only recursive decimal output and one root string (19 October 2026)


B42 on Apple M2 reduced native `to_decimal` versus B41.1 from 1514.429 -> **1178.808 ns** at 1024 bits (-22.2%), 3675.169 -> **3085.894 ns** at 2048 bits (-16.0%) and 9434.652 -> **8564.325 ns** at 4096 bits (-9.2%). The B42 official benchmark verified 768 input pairs and 2688 operation results. All Unix, Windows and Coverage GitHub Actions passed (external Codecov statuses remained separate). The residual is not proven to lie in any single stage: decimal-power construction, multiprecision `divrem`, and leaf conversion were not isolated on the user's Apple M2.

B43 targets a separate, demonstrable source of allocation/copy overhead in the **shared** B22 recursive string assembly. Previously `to_decimal_recursive(value,powers)` returned a new string at every recursive node, constructing separate high and low strings, padding the low string to the exact divisor width and concatenating both. B43 keeps the same power table, selected splitting divisor, exact `divrem` quotient and remainder, 16-limb recursion cutoff and architecture-selected leaf conversion, but threads **one** `std::string& output` down the recursion. At each high branch it passes only the inherited minimum decimal field width not already owned by the lower branch; each low branch receives exactly its divisor's decimal width. Only leaf nodes format or emit leading zeroes, so root output has no leading zeros and all lower fields are exactly zero-padded. This handles recursively split *remainders* that themselves have many digits, including zero, small residuals and long internal decimal-zero runs. The root pre-reserves `bit_length/3+1` characters as a conservative, overflow-safe upper bound on decimal length (`log10(2)<1/3`), and shares its already computed bit length with power-table creation.

Unlike B42, this string-assembly optimization is **shared** across ARM64 and x86-64, using each architecture's existing radix-10^18 / radix-10^9 leaf formatter; it does not reenable the historically slower 128/64 decimal division on x86-64. Do not assume a cross-platform speedup until independently measured on each host. The change removes intermediate recursive string allocations, not the leaf string creation or the multiprecision division cost; the empirical 7000 ns (4096-bit M2) goal remains unverified.

All **880** randomized/adversarial 0..256-limb scenarios passed under Clang AddressSanitizer and UndefinedBehaviorSanitizer, including powers of two, zeros, saturated operands, random dense inputs and distant nonzero limbs. All earlier B22/B25/B31/B42 decimal scenarios and the official 84-operation benchmark are retained. B43 CI and same-M2 performance measurements are pending.

### B44 standalone production-decimal stage diagnostic (19 October 2026)


All B43 Unix, Windows and Coverage GitHub Actions passed; external Codecov patch/project checks failed separately.

B44 introduces only benchmarks/benchmark_integer_decimal_profile.cpp and its separate CMake benchmark target, without modifying production arithmetic, ISA dispatch, public APIs, or any existing benchmark file. It reuses integer_specifications(bits), the exact 32-case official signed Integer input pool, and profiles the corresponding magnitudes. The default sizes are 1024/2048/4096 bits, with 20000 iterations and 31 samples. Example: ./build/benchmarks/benchmark_integer_decimal_profile --bits 1024,2048,4096 --iterations 20000 --samples 31

The diagnostic measures complete native to_decimal, and independently replays B22 power-table construction via production square, recursion divisions via production divrem, leaf conversions via production to_decimal (its <=16-limb leaf dispatch), and append/padding with precomputed exact leaf strings. The 1024-bit size bypasses recursive splits. The recursions are reconstructed outside all timed regions. These are same-input non-additive stage PROXIES, not exact in-situ call-stack timings; do not sum their fractions or claim causal bottlenecks without repeatable same-host evidence. Apple M2 B44 measurements remain pending.

### B44.1 logging-compatible output and user diagnostic measurements (19 October 2026)


The B44.1 standalone benchmark fixes an output-only bug: although `logging_command_line()` already parses `-t/--theme`, `-v/--verbosity` and `-s/--scheduler`, B44 emitted its benchmark report through `std::cout`, outside the Ariadne logger. The B44.1 benchmark now uses `LOGGING_PRINTLN` for all normal report lines, including component measurements and the final verification banner. This restores normal logger theme/verbosity behavior, following all other benchmark executables. The standalone diagnostic's `--summary-only` and `--details` options are shared parser options and do not change its stage measurements. Fatal errors still go to `std::cerr`.

They suggest leaf conversion dominates the measured replays at 2048 bits and recursive division at 4096 bits, while append/padding is small at both. Measurements are user-provided, not independently reproduced; native performance improvement is not claimed. B44's Unix, Windows and Coverage GitHub Actions passed; the separate Codecov patch/project checks were failing.

### B45 remove exact whole-zero divisor limb prefixes (19 October 2026)


The user-supplied B44 4096-bit Apple M2 diagnostic measured production to_decimal 8236.562 ns; independently replayed recursive divrem 4440.488 ns, leaf conversion 3387.277 ns, decimal-power construction 361.692 ns, and append/padding 50.379 ns. Stage timings are separate non-additive proxies, not exact production CPU attribution. At 2048 bits leaf conversion was the largest measured replay (1649.181 ns).

B45 optimizes general shared BigUnsigned long division for a divisor with k>=1 complete zero low limbs. Let B=2^64, a=a_hi*B^k+a_lo with 0<=a_lo<B^k, and b=b_hi*B^k. Then the quotient equals floor(a_hi/b_hi), and the remainder equals (a_hi mod b_hi)*B^k+a_lo. After the original argument-size and scratch validations and output initialization, recurse into the same selected kernel on the canonical high input limb subspans, reusing the existing scratch area. Copy the original low dividend limbs into the remainder only when requested. The original joint-output, quotient-only and remainder-only results, public scratch contract, zero-filled surplus output limbs and no-alias requirements are all preserved. This exact optimization skips normalization, quotient-digit products and correction work for provably zero divisor limbs. Decimal powers 10^(9*2^j) contain many such low limbs because of their factor 2^(9*2^j).

No original tests or official benchmark workloads are changed. Both design journals describe the general arithmetic proof. Full CI and same-host Apple M2 benchmarking are pending. Do not claim an actual B45 speedup until those measurements.

### B46 ARM64 radix-10^18 exact reciprocal leaf division (19 October 2026)


B45's independent Apple M2 stage benchmark showed a reduction of the 4096-bit decimal recursion's divrem replay from 4440.488 to 3678.875 ns and full to_decimal from 8236.562 to 7520.004 ns. Correctness verification passed 768 paired inputs and 2688 operation results. All B45 Unix, Windows and Coverage Actions succeeded; external Codecov patch/project statuses remained failing. The remaining B45 stage replays at 4096 bits were divrem 3678.875 ns, leaf conversion 3464.458 ns, power construction 380.217 ns and append 51.048 ns; these separate timings are non-additive. At 2048 bits the leaf replay was 1674.200 ns and the full production conversion 2813.100 ns.

B46 changes ONLY the existing native ARM64 wide-leaf formatter to_decimal_blocks_wide in src/big_unsigned.cpp. It keeps the radix 10^18, >=9-limb ARM64 dispatch, work-block extraction, B42 pre-reserved blocks and exactly-sized two-digit output, B43 recursion, B45 division, x86-64/Windows base-10^9 fallback, and the unchanged 84-case official benchmark. ARM64 Clang currently lowers the wide 128/64 constant division of each leaf work limb to __udivti3. B46 replaces that division with an exact reciprocal calculation using 64x64->128 products and at most one quotient adjustment.

Let B=2^64 and d=10^18 (<B/2), and take R=floor(B^2/d), which equals floor((B^2-1)/d) because d does not divide B^2. For each leaf division step, N=h*B+l with 0<=h<d and 0<=l<B. Define q0=floor(N*R/B^2). Since 0<=N/d-N*R/B^2 < d/B < 1, q0 is either floor(N/d) or one less, never greater. Compute q0 using the high 128 bits of the product of two 128-bit values: R=r_hi*B+r_lo, N=h*B+l, q0=h*r_hi + floor((h*r_lo+l*r_hi+floor(l*r_lo/B))/B). The cross term fits in unsigned 128 bits with this d. The exact provisional remainder N-q0*d is in [0,2*d), strictly below B; thus the 64-bit modular subtraction l-q0*d is already exact. If provisional remainder>=d, increment q0 and subtract d once. No floating-point estimate, external arithmetic dependency, new public API, global cache, or ISA-specific assembly is introduced.

Before integration, a C++20 proof-of-implementation probe checked 2,643,208 independent hi/lo pairs against the compiler's exact unsigned 128-bit quotient and remainder under Clang ASan/UBSan; 26,303 cases required the single correction. An ARM64 Clang 17 assembly probe compiled the new recurrence into inline integer instructions without any __udivti3 call, whereas the previous implementation emitted such a call. The assembler probe does NOT establish Apple M2 runtime speedups.

All older test corpora and benchmarks are preserved. A regression would require revisiting or reverting the reciprocal path; no runtime speedup is claimed yet.

### B47 lower wide decimal parsing threshold and align input chunks (19 October 2026)


The direct cause is visible in BigUnsigned::from_decimal: the B24 unsigned-128 radix-10^18 parser was enabled only when input.size()>309, while the benchmark's 2^1023+offset values have exactly 308 decimal digits. They therefore used the B19 slower radix-10^9 two-half-products-per-limb path. This behavior is architectural, not a benchmark anomaly.

B47 makes three isolated changes to the existing C++20 parsing path on GCC/Clang compilers supporting unsigned __int128: (1) select the 18-digit wide radix when the input has >=32 digits; (2) pre-reserve ceil(input.size()/18) 64-bit limbs, a conservative bound because 10^18<2^64, with overflow-safe size arithmetic; (3) consume the first possibly short decimal group first, and then use the compile-time constant multiplier 10^18 on every remaining group. Inputs shorter than 32 digits still use the existing radix-10^9 portable path. MSVC and compilers without unsigned __int128 continue to use that portable path at all lengths, preserving cross-platform semantics. The ASCII validation is unchanged and precedes either implementation.

The isolated Clang x86-64 C++20 proxy (not an M2 or full Integer benchmark) measured the 308-digit, 1024-bit input parser at about 573.9 ns for the old 9-digit path and 263.9 ns for the aligned, reserved 18-digit path. At 128/256 bits the corresponding proxy was 48.8/100.8 ns (old) versus 33.8/55.9 ns (new). These are independent algorithm measurements and do not establish native production speedups. Preserve all historical decimal, parse, and signed Integer suites; no benchmark source is changed.

The primary target is construct_decimal at 1024 bits, not to_decimal; verify no regression at short sizes and on ARM64, x86-64, and portable toolchains. B47 integrated CI and native Apple M2 measurements are pending.

### B48 ARM64 small-divisor fast path in shared divide_wide (19 October 2026)


The original division-only benchmark pool has 32 fixed pairs: 16 with small, single-limb divisors, eight with half-width divisors and eight near-equal magnitudes. The benchmark and its operand pool must remain unchanged.

B48 optimizes the shared low-level divide_wide(high,low,divisor) when the existing precondition high<divisor holds and divisor<=2^32-1, on native ARM64 GCC/Clang only. Previously the unsigned-128 numerator quotient/remainder may call the compiler's __udivti3 helper; x86-64 already uses B26's hardware DIVQ. Write N=high*2^64+low. Divide U=high*2^32+(low>>32) by divisor to obtain the high 32 bits of the quotient and remainder rU. Divide L=rU*2^32+(low&0xffffffff) to obtain the low 32 quotient bits and final remainder. The invariants high<divisor<2^32 and rU<divisor imply U<2^64, L<2^64, and both quotient halves fit 32 bits, so both stages use exact native 64-bit division, never a 128-bit arithmetic runtime helper. The high==0 fast path, ARM64 larger-divisor existing implementation, x86-64 assembly DIVQ, Windows MSVC _udiv128 and portable fallback are unchanged. The improvement feeds all three joint, quotient-only and remainder-only BigUnsigned division modes, including the ordinary one-limb divisors used by Integer's official benchmark.

Standalone C++20 Clang ASan/UBSan checked 1,200,156 exact small quotient/remainder pairs against independent unsigned __int128 division, including saturated quotient, high=0 and divisor 1. A Clang 17 AArch64 assembly probe confirmed the two-stage narrow implementation emits native UDIV/MSUB and no __udivti3 call; the previous wide implementation emitted __udivti3. This is assembly verification, not an Apple M2 timing measurement. All original coverage and the official 84-operation workload are preserved. Benchmark on the same M2 before claiming a B48 speedup; check in particular quotient/remainder at 1024, 2048 and 4096 bits, as well as to_decimal to detect any unintended ARM64 regression.

### B49 ARM64 pre-inverted narrow quotient loop (20 October 2026)


B48 versus B47 reduced native remainder timings at 1024/2048/4096 bits from 121.002/257.183/546.194 ns to 103.625/221.458/467.669 ns (-14.36%/-13.89%/-14.38%). Quotient timings changed more modestly from 110.192/236.183/489.544 to 108.296/226.177/484.854 ns (-1.72%/-4.24%/-0.96%). All B48 Unix/Windows/Coverage Actions passed; external Codecov patch/project checks failed separately.

B49 targets quotient generation for single-limb divisors 2..2^32-1 on native ARM64 GCC/Clang. These account for 16 of the 32 official division-only workload cases (small divisors 3..35), so the goal is to accelerate the quotient-only path while RETAINING B48's remainder-only small-divisor fast path. The existing shared BigUnsigned::divrem_limb_into function used divide_wide per limb; B48's divide_wide requires two ARM64 UDIV operations per 64-bit limb for a 32-bit divisor. B49 selects a new exact pre-inverted path when the dividend contains >=8 limbs, the divisor is between 2 and 2^32-1 inclusive, and native ARM64 unsigned 128-bit multiplication is available. Outside that case the old loop, MSVC, x86-64 and portable paths are unchanged. The optimization is general for arbitrary 64-bit limb inputs, NOT specialized to the benchmark's operand values.

Let d be the divisor, B=2^32, and R=floor(2^64/d). The reciprocal is computed once using a 64-bit quotient/remainder operation: R=floor((2^64-1)/d) + [((2^64-1)%d)==d-1]. For each 64-bit dividend word, split its high/low 32-bit halves; each 64-bit partial numerator X is strictly below d*B<=2^64. Set q0=floor(X*R/2^64) from the high half of a 64x64 product. Since 0<=X/d-X*R/2^64<X/2^64<1, q0 is either the exact floor(X/d) or one smaller. The provisional residual X-q0*d is in [0,2*d), so one compare/subtraction (and increment of q0) yields an exact 32-bit quotient half and remainder. Each full quotient word combines both halves. The approach replaces repeated hardware divisions with UMULH/MSUB-style operations and requires no new public API, dependency or global table. The B48 remainder-only path and shared long-division digit estimator remain unchanged; the large-divisor paths are untouched. The returned full remainder, output zero-filling, public failure semantics and no-alias requirement remain unchanged.

An AArch64 Clang 17 assembly probe for the core quotient step showed native UMULH and MSUB, without division helper calls. Those checks establish arithmetic correctness for the isolated implementation, NOT performance on Apple M4.

Measure the existing official 84-operation benchmark on M4, 20,000 iterations and 31 samples, comparing B48/M4 with B49/M4. A performance improvement is NOT claimed until same-host measurements are supplied.

### B50 wider 19-digit input radix for construct_decimal (20 October 2026)


The B49 official 84-case benchmark is measured on the SAME Apple M4 ARM64 as B47 and B48. Although B49 reduced native quotient at 2048 bits to 198.350 ns, construct_decimal remains a large contributor to total time.

B50 is a narrowly scoped change to the wide 128-bit-supported BigUnsigned::from_decimal path that B47 selected for >=32-digit inputs. Each radix-10^18 block currently contains at most 18 decimal digits in one 64-bit limb. B50 uses the larger exact radix 10^19=10000000000000000000, still strictly below 2^64, to consume 19 digits per pass and reduce full multiprecision Horner scans. The first possibly shorter group is still consumed first, all subsequent groups use the constant radix, and the pre-reservation changes from ceil(decimal_digits/18) to ceil(decimal_digits/19). Each wide product has the form digit*radix+carry, where digit<2^64, radix<=10^19<2^64 and carry<radix. Hence the maximum product is (2^64-1)*radix+(radix-1)=2^64*radix-1<2^128, so the existing unsigned __int128 calculation stays exact. The historical base-10^9 portable path on Windows/MSVC and other compilers without unsigned 128-bit remains unchanged; the unsigned-wide path is still selected at >=32 digits. ASCII validation and signed Integer construction are unchanged.

A preliminary isolated Clang 17 x86-64 C++20 comparison (NOT an Apple M4 result) of otherwise identical Horner parsers showed radix-10^19 close to radix-10^18 at 1024 bits, roughly 3-5% faster at 2048 bits, and roughly 1-3% faster at 4096 bits. A prototype balanced parser and a 38-digit two-limb radix were slower in that environment and were NOT integrated. The numbers above are input-parser proxy measurements only, not independent Apple M4 evidence.

All standard CI matrices must pass. Keep B49 as the rollback baseline if the 19-digit wide path does not show a repeatable same-host benefit.

### B51 first x86-64 performance campaign: revisit B46 reciprocal leaves (20 October 2026)

The code review found 1) x86-64 production multiply dispatch depends on compiler: GCC enables the fused kernel and associated large balanced Karatsuba and symmetric-square specializations, while Clang generally chooses the B3 row-schoolbook path; 2) portable and x86-64 division kernels already include x86 DIVQ and several B28/B36 quotient fast paths; and 3) x86-64 decimal leaves remained on radix-10^9 in B50. Historical B25/B30 radix-10^18 experiments REGRESSED in the official x86-64 benchmark. B30 used a real 128/64 DIVQ loop, which is not the same algorithm as B46's exact pre-inverted radix-10^18, so these old measurements do not settle the current hypothesis.

B51 changes only the existing compile-time guards around to_decimal_blocks_wide and its 9+-limb selection in src/big_unsigned.cpp: native x86-64 GNU/Clang targets providing unsigned __int128 now reuse the *unchanged B46* quotient/remainder reciprocal and B42 18-digit pair formatter, just as ARM64 does. The B46 quotient uses 64x64->128 products and at most one remainder correction; it requires neither a runtime 128-bit divide nor an x86-specific instruction. The full B19 base-10^9 fallback remains selected for <=8-limb magnitudes, MSVC and other architectures; ARM64 behavior is unchanged. No arbitrary handwritten assembly is introduced before the available competing mulq/adcq, Comba, fused and Karatsuba kernels are benchmarked on the actual target host. This is a **controlled first x86-64 experiment**.

Clang 17 -O3 baseline-x86 medians (ns):
- 8 limbs: narrow 482.030, reciprocal wide 304.750;
- 9 limbs: narrow 548.132, reciprocal wide 336.512;
- 12 limbs: narrow 956.700, reciprocal wide 524.170;
- 16 limbs: narrow 1697.650, reciprocal wide 847.686;
- 32 limbs: narrow 7771.950, reciprocal wide 3069.380. GCC 14 -O3 baseline-x86 medians (ns):
- 8 limbs: narrow 551.198, reciprocal wide 346.259;
- 9 limbs: narrow 645.396, reciprocal wide 366.258;
- 12 limbs: narrow 1093.210, reciprocal wide 553.028;
- 16 limbs: narrow 1998.440, reciprocal wide 797.383;
- 32 limbs: narrow 8642.830, reciprocal wide 2632.700. These are **isolated local CPU timings**, not B51 official results and not a same-host comparison with the uploaded x86-64 baseline. Run CI across GCC, Clang, AppleClang, Windows/MSVC and ThreadSanitizer, and compare the official x86-64 benchmark on the user's machine; revert the x86 selector if the integrated experiment repeats B30's earlier regression.

Use diagnostic data to decide whether MULQ/ADCQ, BMI2/ADX gated kernels, square specialization or Karatsuba thresholds are justified for that target, rather than assuming M4 ARM64 instruction-level tuning transfers directly.

### B52 x86-64 GCC 32-limb assembly dispatch informed by native benchmarks (20 October 2026)


The user's official B51 x86-64 benchmark (20,000 iterations, 31 samples, unchanged 84 cases) demonstrated that enabling the B46 exact reciprocal base-10^18 leaf under x86-64 significantly accelerated Integer::to_decimal versus B50 on the SAME host: 1024 bits 1793.826->720.012 ns (-59.86%); 2048 bits 4526.997->3280.192 ns (-27.54%); 4096 bits 11715.378->9396.203 ns (-19.80%). This CONFIRMS that the B46 reciprocal implementation, unlike the historical pre-B46 DIVQ implementation rejected in B30.1, can beat the old x86 radix-10^9 output path in production.

The user also ran the existing benchmark_big_unsigned_multiply on that SAME x86-64 host, requesting all experimental variants, 20,000 iterations and 31 samples with Karatsuba cutoff 32. The selected production kernel was reported as fused (x86-64 GCC). For preallocated dense general products: bits fused MULQ/ADCQ relative assembly delta 128 16.409 16.709 +1.83% 256 22.367 23.956 +7.10% 1024 205.249 220.956 +7.66% 2048 954.267 810.044 -15.11% 4096 2964.676 3153.710 +6.38% At 2048 bits the historical portable and Comba alternatives were much worse; Karatsuba cutoff 32 was 961.873 ns, close to fused 954.267. At 4096 bits Karatsuba cutoff 32 was 2994.552 ns, again near fused; changing every general product to MULQ/ADCQ would REGRESS the observed 4096-bit result. The x86-64 hardware CPU model and target compiler flags are not supplied in the log: do not extrapolate these absolute values to a different CPU, and do not claim BMI2/ADX availability.

Other sizes, every truncated output below 64 limbs and Clang/AppleClang/MSVC/ARM64 all retain their previous selectors. The existing assembly has been benchmarked on the user's host and exercises the same full-product recurrence; no new assembly, CPU feature probing, public API or benchmark code is introduced. No value-dependent benchmark specialization exists: the 32-limb dispatch is valid for arbitrary limb values and applies equally to standalone owning products and Karatsuba base cases. The owning GCC Karatsuba 64-limb path uses 32x32 z0/z2 leaves, so it may also benefit; its 33x33 cross term remains fused. Only a new on-target integrated run can confirm either gain.

If this prediction fails, inspect build flags / path before adjusting dispatch.

Subsequent campaign stages should profile the existing DIVQ quotient estimator, reciprocal based small-divisor paths, and optional BMI2/ADX kernels only after checking actual target ISA and paired benchmark results; do not silently expand the narrow 32-limb selector without evidence.

### B53 port ARM64 B38/B40 sparse-row multiplication to x86-64 (20 October 2026)


On the user's SAME x86-64 host, the B52 native multiplication diagnostic confirmed the 32x32-limb MULQ/ADDQ/ADCQ dispatch: 2048-bit preallocated selected 953.096 ns (B51) to 853.847 ns (B52), 10.41% faster; at 4096 bits cutoff-32 Karatsuba 2994.552->2675.843 ns, 10.64% faster. In the official 84-operation benchmark, 2048-bit Integer mul fell 1030.988->861.190 ns (-16.47%), and 4096-bit Integer mul 3117.157->2845.357 ns (-8.72%). At 4096 bits pow3 remained approximately unchanged 8363.683->8362.002 ns, which exposes the cost of executing full schoolbook rows for sparse powers.

B53 transfers the *general mathematical* optimizations of ARM64 B38/B40 to native x86-64 multiplication, without special-casing the benchmark's pow3 operands. The x86 GCC selected fused unsigned __int128 row-schoolbook now has two templated implementations: sparse skips any exact-zero outer limb row after initially zero-filling the output; dense has no zero check inside the outer loop. The existing x86 GNU/Clang MULQ/ADDQ/ADCQ kernel gets the same two specializations, preserving instruction syntax and carry propagation. If a row multiplier ai==0, every partial product, carry and output update from that row is identically zero, including the final carry slot which earlier rows cannot touch. Thus skipping a row is exact for all full/truncated and noncanonical raw inputs.

As in ARM64 B40, when active outer rows = min(a.size(),output.size()) is <64 or b.size()<64, select the sparse variant without scanning; larger cases first scan up to a threshold max(2,floor(active_rows/8)) exact zero limbs and select sparse only upon reaching it, otherwise dense. The prescan exits early when sparse. This transfers ARM64's protection against costly dense-data branches at large sizes. The GCC 32x32 B52 assembly selector remains unchanged; it now benefits from zero-row elision too. The owning GCC Karatsuba size/balance selector is unchanged (a separate candidate for B54), and existing x86 symmetric square remains unchanged. Clang's default B3 schoolbook and MSVC are unchanged; only the explicit x86 carry benchmark on Clang uses the new template selection. ARM64 and portable kernels and all public exception, view ownership, no-alias, exactness and scratch contracts are unchanged. Neither BMI2 nor ADX is required.

The proxy is NOT the user's x86 target and cannot establish a performance claim. Production test corpus adds 16 operand-length pairs x24 adversarial/ random patterns = 384 source pairs. Coverage spans the B52 31/32/33-limb selector, B40 63/64/65 adaptive prescan, 128-limb unbalanced operands, saturated, zero, alternating, sparse and dense patterns.

Inspect dense 1024/2048/4096 preallocated multiply measurements for regressions and avoid treating Apple M4 ARM64 absolute timings as same-host comparisons. B53 has no on-target performance measurements yet. All standard CI checks should pass before acceptance.

### B53.1 fix doubled x86-64 inline-assembly escape sequences (20 October 2026)


The user's native x86-64 GCC -O3 build failed in src/big_unsigned_x86-64.cpp with assembler "invalid character '\\' in operand 1" near line 200. Root cause: B53 generated the five consecutive MULQ/ADDQ/ADCQ assembly string fragments with C++ source spelling `\\\\n\\\\t` (double backslashes), which produces literal backslash characters instead of the required newline/tab separators. B53.1 changes exactly these five source occurrences to `\\n\\t`. No instruction, register constraint, algorithm, dispatch, public API, test corpus, benchmark, or platform logic is modified. The fixed sequence was compiled locally on x86-64 with GCC 14.2 and Clang 17 (-O3 -std=gnu++20 -Wall -Wextra -Wpedantic -Werror) and each compiler passed 250000 exact products against unsigned-128 arithmetic. This isolated assembly smoke check does not replace full project CI or the user's native B53 benchmark. B53 official same-host measurements remain pending.

### B54 x86-64 GCC adapt square and Karatsuba to operand density; retune 32-limb dispatch (20 October 2026)


Versus B52 on the same x86-64 CPU, 1024/2048/4096-bit Integer::mul changed 252.938/861.190/2845.357 to 79.513/114.049/701.262 ns; pow3 changed 738.724/2193.092/8362.002 to 310.053/686.622/2176.326 ns. Dense preallocated selected multiplication was essentially stable at 1024/2048/4096 bits (204.242/853.847/2969.094 ns B52 versus 205.667/861.181/2934.027 ns B53). The B53.1 build escape fix restored correct MULQ/ADCQ asm newlines; all seven GitHub Actions eventually passed (external Codecov patch/project failures reported independently).

B54 addresses three remaining architecture-general opportunities, without operand-value-specific benchmark special cases:

1. **Retune the B52 32x32-limb dispatch**: on the SAME B53 x86 machine and compilation, explicit dense 2048-bit `prealloc_fused` measured 801.142 ns, whereas forced `prealloc_x86_64_carry` measured 858.091 ns and selected 861.181 ns. These new B53 fused measurements supersede B52's earlier preference for assembly (B52 fused 1008.084 ns versus assembly 854.808 ns). B54 removes the special-case 32x32 selected assembly routing, using B53's exact sparse-aware fused recurrence at all GCC input dimensions. The explicit ASM entrypoint remains available for comparison and optional testing. Sparse 32x32 performance is NOT YET established on the user's machine: only B54's official 84-case benchmark can validate that changing the 32x32 selection does not undo sparse gains.

2. **Transfer ARM64 B38/B40 adaptive symmetric square**: the existing B32 GCC x86-64 square forms the upper off-diagonal triangle and doubles it, then adds the diagonal products. B54 introduces sparse/dense compile-time instantiations. If an outer limb is zero, omit its entire off-diagonal row. Earlier triangular rows never write the row's final output slot `output[i+n]`, and the output is initialized to zero. All diagonal products and carry passes are preserved. For <=64 limbs, select the sparse variant without a prescan; for >64 limbs classify the first n-1 active outer rows with B53's exact-zero threshold max(2,floor((n-1)/8)) and select dense without per-row checks when the threshold is not reached. The public square() size >=8 dispatch, original arithmetic, allocation and output contracts remain. This should especially reduce cost for sparse 1024/2048/4096-bit squares, previously roughly 166.670/492.103/1821.579 ns on x86 B53.

3. **Transfer ARM64 B39 owning Karatsuba density gate**: the B35 GCC x86 owning mul() used cutoff-32 Karatsuba for all >=64-limb nearly balanced products, including sparse powers where B53's zero-row schoolbook is now cheap. B54 requires both operands to have <=floor(size/8) exact-zero limbs before choosing this same balanced Karatsuba path. Otherwise fall through to the current sparse-aware row-schoolbook selected kernel. Inputs of insufficient size or balance bypass density scanning as before; no new scratch allocation occurs for rejected Karatsuba. Dense products retain cutoff-32 Karatsuba. The same size/balance and zero-density criteria were independently used on ARM64 in B39.

All earlier test cases are preserved; the official benchmark and its 84 original workloads remain byte-for-byte unchanged.

Acceptance: run full CI including GCC, Clang, AppleClang, Windows/MSVC, Coverage and sanitizer. Repeat the unchanged 31-sample, 20,000-iteration official x86-64 benchmark and the explicit preallocated dense multiplication diagnostic with fused/carry/Karatsuba alternatives. If B54 harms sparse 2048 multiplication or dense 4096 Karatsuba, revert or retune only the responsible selection without changing any benchmark input or test expectations.

### B55 x86-64 independently verified MULQ/ADCQ multiply-subtract division step (20 October 2026)


B55 examined the first hypothesis, directly copying B49's preinverted narrow divisor strategy from ARM64 to x86-64. A standalone x86-64 AMD EPYC 9V74 GCC 14.2 and Clang 17 proxy compared exact scalar DIVQ 128/64 division with B49's two reciprocal multiply-high steps per 64-bit quotient digit, and full remainder-only equivalents, all using the unchanged divisor-independent input recurrence. On this proxy B49 reciprocal was *consistently slower*: GCC 32 limbs divisor=35 gave 119.167 ns for DIVQ and 181.581 ns for reciprocal quotient; the rem-only path 127.614 versus 206.651 ns. Clang 32 limbs divisor=35 gave 119.395 versus 221.048 ns for quotient. The exact CPU differs from the user's unknown x86 model, but the results disprove automatic portability and justify retaining B26's existing native DIVQ for narrow divisors on x86. The ARM64 B49 code and all related selection stay UNCHANGED.

The second hypothesis targets the genuine multi-limb divisor path, which normalizes the dividend and denominator, estimates each quotient digit with divide_wide, and subtracts digit*divisor from an m+1 limb window. B27 already fused the multiply and subtract in C++, but on x86-64 its 128-bit wide arithmetic still incurs separate borrow comparisons. B55 adds a baseline native MULQ -> ADDQ/ADCQ -> SUBQ/ADCQ *single carry chain* for the subtraction recurrence, selected only when m>=8 limbs on GNU/Clang x86-64, after all output and workspace contracts have been validated. Each step computes high:low = divisor[i]*digit+carry and then window[i]-=low while incorporating the borrow into the next carry. The exact algebra is unchanged. SparseZero==true retains B29's handling of normalized zero divisor limbs by subtracting only incoming carry. Final high-limb subtraction and underflow flag are identical.

A local AMD EPYC proxy of the exact native recurrence was compiled with -O3 -std=gnu++20 -Wpedantic -Werror using GCC 14.2 and Clang 17, plus Clang ASan/UBSan. In the GCC local timing proxy, the 16-limb kernel fell 48.10 -> 43.17 ns (~10.2%) and the 64-limb kernel 78.99 -> 71.40 ns (~9.6%); Clang 16-limb 19.31 -> 16.44 ns (~14.9%). Small 1/2 limb cases were similar/slower, motivating B55's >=8-limb selector. These are isolated local proxy timings, NOT evidence of the integrated user-host speedup.

No adjustments to the official 84-case workload or benchmark expectation. CI must pass on GCC, Clang, macOS AppleClang, Windows MSVC and all sanitizer/coverage configurations. Preserve B54 as the rollback baseline if division gets slower or no repeatable native improvement appears.

### B56 category-resolved x86-64 division profiling (20 October 2026)


The B55 benefits were repeatable at the granularity of these measurements; B55 remains the accepted production arithmetic baseline.

Before further instruction-level changes, distinguish which *algorithmic division category* dominates. The unchanged official workload mixes small single-limb divisors, half-width multiprecision denominators and near-equal operands. Its overall quotient/remainder medians cannot isolate the expensive subpath, nor can their medians be decomposed into additive stage costs. B56 intentionally makes NO PRODUCTION ARITHMETIC CHANGE.

This is **diagnostic instrumentation only**, with no numerical or speedup claims before the executable passes CI and is run on the user's same x86-64 host. A recommended initial diagnostic run is: cmake --build build --target benchmark_big_unsigned_division_profile -j8 ./build/benchmarks/benchmark_big_unsigned_division_profile -t dark \ --bits 1024,2048,4096 --iterations 4000 --samples 15 Repeat slowest categories at 20000 iterations/31 samples before choosing B57. Inspect whether multiprecision quotient/remainder paths still dominate, or if overhead for small/near-equal denominators is more consequential. Avoid mistaking preallocated and owning cost differences for arithmetic speedups. Do not change oracle coverage or adjust benchmark inputs to make a faster path win.

### B56.1 fix move-only independent exact-reference fixtures pool construction (20 October 2026)


Named-value return optimization is optional and the compiler must find a legal move constructor: GCC rejected it with "use of deleted function std::array<Fixture,32>::array". Both live objects remain valid for independent destruction; self move-assignment is guarded. There is no modification to any production library code, test source, official benchmark workload or diagnostic measurement logic. Run the complete benchmark and CI to establish integration and native timing; that has NOT yet been confirmed.

### B56.2 ensure exhaustive shape switch has required default (20 October 2026)


The B56.1 diagnostic CI recorded GCC 13 Release and AppleClang 21 Release build failures, while Clang 16 Release/TSan passed. CI enables WARNINGS_AS_ERRORS=ON. Review of the new diagnostic found a 12-case switch that intentionally enumerated every Shape but omitted the explicit default label required by the project's -Wswitch-default policy. B56.2 moves the existing logic_error into a switch default case without changing which inputs or cases can execute. An isolated C++20 reproducer confirmed -Wswitch-default -Werror rejects the old fully-enumerated switch on both GCC 14 and Clang 17. This is a warning-policy compliance fix, not a change to the arithmetic, oracles, measured operations, input fixtures or official benchmark. The exact CI error annotations were not accessible through the GitHub connector; confirmation requires new CI results. No user-host B56 timings have been recorded.

### B57 x86-64 long single-limb preinverse and universal power-of-two division (20 October 2026)


The preallocated comparison likewise showed high arithmetic costs. The unchanged official Integer benchmark's nontrivial quotient/remainder corpus has 16/32 one-limb divisors (50% of inputs), 8/32 ~half-width divisors and 8/32 near-equal operands, as verified in benchmarks/benchmark_exact_arithmetic_workload.hpp. B56 CLI and production arithmetic remain intact.

B57 evaluated a **different 64/64 normalized pi1 reciprocal** from ARM64 B49's two 32-bit stages. For B=2^64 and normalized d in [B/2,B), compute inv=floor((B^2-1)/d)-B. This fits 64 bits. Its exact construction is a SINGLE B26 DIVQ computing (hi=(~d),lo=UINT64_MAX)/d; the high limb is strictly <d, so DIVQ cannot overflow. Normalize the entire quotient stream virtually, without allocations: shift the divisor by countl_zero, start the remainder at the top discarded bits, and produce each normalized input limb from adjacent input limbs. The pi1 recurrence uses one 64x64->128 product plus modulo-2^64 subtraction/correction to give the exact digit and next normalized remainder, then shifts the final remainder down. The result is mathematically identical to division using DIVQ.

The revised single-DIVQ inverse setup was also verified with both compilers over the same cases. The median proxy times varied by operand size and compiler: in representative GCC runs at 64 limbs divisor=35, DIVQ quotient ~272 ns, normalized reciprocal quotient ~229 ns; at 32 limbs the reciprocal was not better, and for remainder-only it was often worse. Clang's 32/64-limb normalized reciprocal was generally not faster than DIVQ. These tests run on a proxy, not the user's x86 CPU, and are NOT integrated B57 performance claims. Therefore the reciprocal is selected **only** for >=64 raw quotient limbs on GNU x86-64 GCC, not Clang, not MSVC, not ARM64 and not remainder-only (which remains on B26 DIVQ); existing AMD and other x86 microarchitectures may differ. Small quotient loops retain B26 and ARM64 B49 paths as appropriate. No BMI2/ADX required.

B57 also recognizes exact *power-of-two divisors* and performs the quotient/remainder by shift and mask. BigUnsigned::divrem_limb_into uses the fast exact 2^s path, including 2^0=1; the multi-limb divrem_selected_into handles b=B^k*2^s and all three compile-time result modes. Its fast path is after all existing error, output-size and scratch-size checks, so it maintains their fail-before- mutation contracts; the results are fully zero-filled as before. For arbitrary dividends including a<b or noncanonical raw limb inputs on the single-limb API, the exact mathematical identity holds. The selected multi-limb algorithm does not allocate intermediate buffers or change public exception behavior; ARM64, Windows/MSVC and other architectures can benefit from this ISA-independent power-of-two optimization. The B45 trailing zero-limb and B36 near-equal optimizations continue to handle all non-power-of-two inputs unchanged.

The production suite adds 1440 raw-limb single-divisor pairs across 0/1/7/8/16/31/32/63/64/65/96/128 limbs, 15 different divisors and eight randomized/adversarial patterns, plus 2160 power-of-two multi-limb cases across nine dividend sizes, eight divisor sizes, five bit positions and six patterns. The official 84-case benchmark and the B56 diagnostic benchmark are UNCHANGED.

The main same-host acceptance test is to repeat the original 84-case benchmark at 20k iterations and 31 samples and the existing B56 profile at 4k iterations and 15 samples, checking separate limb_small, limb_wide, power_of_two and short_2 categories at 1024/2048/4096 bits for wins and regressions. Since GCC 64-limb reciprocal behavior is experimental on the user's unknown x86 CPU, retain B56/B55 as rollback baseline if results regress.

### B57.1 restore warnings-as-errors clean bit-count conversions (20 October 2026)


The user's native x86-64 GCC build of B57 succeeded but warned at three new conversion sites in src/big_unsigned.cpp, lines 502, 525 and 819: std::countr_zero/std::countl_zero return int whereas the destination variable is unsigned. All seven B57 GitHub Actions checks failed at their build steps with WARNINGS_AS_ERRORS=ON, including GCC, Clang, AppleClang and Windows. B57.1 makes these three conversions explicit via static_cast<unsigned>; the mathematical count range 0..64 guarantees they are nonnegative. GCC and Clang compiled isolated C++20 -Wsign-conversion -Werror reproducers. No algorithm, thresholds, test corpus or APIs change.

Re-run under controlled conditions before attributing unrelated decimal construction timing changes to B57.

B56/B57 division profile at 4096 bits: limb_small quot 1281.487->314.578 ns (-75.5%), limb_wide quot 1349.135->312.612 ns (-76.8%), power_of_two rem 749.174->87.441 ns (-88.3%). Generic single-limb remainder limb_wide rem 1322.548->1306.280 ns, nearly unchanged, and short_2 quot 1665.143->1706.745 ns (+2.5%). The 84-case official benchmark and existing tests remain unchanged. Full CI retest required.

### B57.2 fix AppleClang B57 single-limb test shadow warning (20 October 2026)


AppleClang arm64 Release CI failed to compile tests/test_big_unsigned.cpp with -Wshadow -Werror: the B57 test loop introduced a local `remainder` while test_division_limb() already held the structured- binding variable `remainder` from the fixed radix-boundary case. It is solely a warnings-as-errors source compatibility fix; full AppleClang and cross-platform CI validation is still pending.

### B58 diagnostic: decimal cutoff-16 versus cutoff-32 and recursive split costs (20 October 2026)


The independent B44 diagnostic measured, at 2048 bits, full 3115.608 ns, division replay 1212.801 ns and leaf replay 1411.688 ns; at 4096 bits, full 8595.959 ns, division replay 4789.633 ns and leaf replay 3056.365 ns. Replay fractions are not additive and are not exact in-situ attribution.

B58 is deliberately measurement-only, not a production threshold change. An internal benchmark-only symbol in src/big_unsigned.cpp duplicates the B22/B43 split recursion with a 32-limb leaf cutoff while invoking the SAME production to_decimal_blocks, square and divrem implementations. The B44 diagnostic now times production cutoff-16 against cutoff-32 on the same 32 official-input magnitudes at 1024/2048/4096 bits. The split results, like the original stage replays, are independent non-additive proxies, not literal in-situ call-stack measurements.

B58 source compilation, integrated execution, CI and user-host timings remain to be verified.

### B59 promote measured GCC x86-64 decimal cutoff-32 (20 October 2026)


B59 promotes the measured 32-limb leaf cutoff ONLY for __x86_64__ GCC (not Clang or MSVC) with unsigned __int128. ARM64, x86 Clang, Windows and other targets retain cutoff-16. The existing duplicated B58 cutoff-16/32 recursion and top-level formatter are consolidated into a single template parameterized by a compile-time cutoff: production selects the architecture constant while both diagnostic versions remain callable. The 84-operation official benchmark and original 32-input benchmark workload are unchanged. The standalone profile compares production against explicit diagnostic cutoff-16 AND cutoff-32, preserving repeatable same-host pre-promotion comparisons rather than letting the experimental comparator become identical to production.

The B58 diagnosis suggests optimizing multiprecision division for decimal-power divisors after this promotion; B59 itself does not alter division. Repeat the official 84-operation benchmark at 20000 iterations and 31 samples on the same host; validate all release, sanitizer, coverage and Windows CI before treating B59 as cross-platform verified. This documentation records the user's B58 results, not an independently executed B59 integration or on-target timing.

### B60 diagnose GCC x86-64 single-limb remainder-only pi1 (20 October 2026)


Thus the primary weighted-total parity objective is now met, though a single successful full run does not guarantee cross-run stability. Production 2048/4096 to_decimal was 2068.069/6834.012 ns, versus 3159.524/8717.960 ns in B57.2. Unrelated construct_decimal changes should not be attributed mechanically to the B59 selector. All B59 GitHub Actions workflows (Unix, Windows, Coverage) have completed successfully as checked on 20 October 2026.

The B59 production-versus-cutoff diagnostic independently reported at 1024 bits 684.641/684.808/685.323 ns (production/16/32); at 2048 bits 2036.237/3111.808/ 2031.133 ns; at 4096 bits 6685.362/8586.050/ 6686.932 ns. These non-additive diagnostics must NOT be treated as exact new production stack attribution.

B60 does NOT change the successful B59 production division/decimal selector or the immutable 84-operation benchmark. The B56 division profiler receives one additional owning-result `rem_pi1` row on GCC x86-64 for its existing exact 32-pair `limb_small` and `limb_wide` categories. The diagnostic-only function in src/big_unsigned.cpp reuses the B57 exact normalized preinverse digit recurrence, discards quotient digits, and returns an owning BigUnsigned remainder. Other hosts retain the existing production rem code; no public header/API or benchmark input pool has been changed.

No on-target B60 measurements or full CI results have yet been verified.

Reject the candidate if it loses, regardless of earlier isolated proxy assumptions. Keep B59 as rollback baseline for weighted-total performance.

### Post-B59 performance acceptance gates: per-size parity and worst-case cap (20 October 2026)


1. 2. 3. Preserve the production semantics, benchmark input corpus and result-lifetime matching. Recheck near-threshold results with repeated same-host runs; do not claim durable success from noise-level wins.

Avoid regressing ARM64 while tuning x86; if the two acceptance gates are applied cross-platform, ARM64's 4096-bit aggregate also needs improvement.

The B60 diagnostic commit 351ff6761711 leaves the B59 production arithmetic unchanged. B60 CI and same-host rem_pi1 timing are not yet confirmed. A change may improve one failure yet regress another: require an updated complete 84-case comparison before promotion, with a recorded rollback baseline. No algorithm or benchmark change is made by this acceptance-policy journal entry.

### B60.1 fix namespace of diagnostic declaration under GCC (20 October 2026)


The user's Ubuntu x86-64 GCC B60 build failed before running the division profile. The new forward declaration for `rem_limb_preinverse_diagnostic` in benchmarks/benchmark_big_unsigned_division_profile.cpp was mistakenly written inside the file's anonymous namespace, creating a *nested*, distinct `Ariadne` namespace. `BigUnsigned` belongs to the global `::Ariadne` namespace, so the compiler rejected the declaration and subsequent calls. B60.1 moves only that forward declaration into global `::Ariadne`, before opening the anonymous namespace. Native integration/CI must still confirm the fix; no performance result is claimed for B60.1.

### B61 promote measured single-limb remainder-only pi1 on GCC x86-64 (20 October 2026)


B61 accepts only the experimentally supported GCC x86-64 single-limb remainder path for >=8 canonical dividend limbs and a non-power-of-two divisor. The implementation extracts the validated B60 exact normalized pi1 recurrence into one shared internal helper, called by production rem-only selected division and by its diagnostic wrapper. No quotient buffer or scratch allocation is introduced. The B57 power-of-two fast path executes before the candidate. Quotient-only and joint division preserve B57's existing selection. Other compilers, short inputs, ARM64 and Windows retain their old paths.

The 32-pair pool, 12 shapes and official 84-case benchmark are untouched. Numerical comparisons remain independent; do not change test expectations to make the optimization pass.

The B60.1 diagnostic result alone does not establish that B61 meets either integrated gate. CI and full user-host benchmark remain to be checked before calling B61 an integrated success. Roll back to B59/B60.1 if the new production selector regresses independent correctness or aggregate performance.

### B62 diagnostic: quotient-only single-limb pi1 threshold at 16/32 limbs (20 October 2026)


The global worst is now quotient at 2048 bits, 3.386x, followed by quotient at 1024 bits, 3.282x. At 1024 and 2048 bits the remainder fell from 286.099/539.469 ns to 176.753/320.884 ns. The aggregate total worsened from B59 0.977x to 0.999x despite the remainder gain: construction and other operations moved as well, and the short sizes cannot invoke the new >=8-limb remainder kernel. Do not attribute unrelated shifts to B61; repeat full runs under comparable host conditions to estimate noise.

B57 already selects exact quotient-producing normalized pi1 on GCC x86-64 when the input has >=64 limbs. At 1024/2048 bits (16/32 limbs), it still uses per-digit native DIVQ. The pi1 comparator replays B57's normalized inverse estimation and digit correction for any input length, including 16 and 32 limbs. The DIVQ comparator uses the old exact recurrence independently of current production size selection. Both fill the same reusable quotient output span and return a 64-bit remainder. This kernel-only preallocated timing is not an owning `quot()` speedup claim and cannot be added to the original owning or divrem rows.

The new diagnostic is not yet compilation/CI-verified or independently measured on the user's host.

### B63 x86-64 GCC quotient-only normalized pi1 threshold 16 limbs (20 October 2026)


At 1024 bits: small 90.369 vs 300.657 ns (pi1 3.33x faster), wide 83.630 vs 318.015 ns (pi1 3.80x faster). At 2048 bits: small 147.365 vs 614.124 ns (4.17x), wide 144.513 vs 638.188 ns (4.42x). At 4096 bits: small 271.682 vs 1213.953 ns (4.47x), wide 266.891 vs 1293.786 ns (4.85x).

Crucially, the unchanged official 32-case integer division workload in benchmark_exact_arithmetic_workload.hpp uses single-limb divisors for i%4==0 and i%4==3, i.e. 16 of 32 inputs. For i%4==1 and i%4==2 it uses half-width and near-equal divisors. Consequently B63 can improve the real official `quotient` benchmark at 1024/2048 bits without modifying the workload, but cannot address the multi-limb half of those operation cases.

B63 promotes the measured `>=16`-limb pi1 threshold in `BigUnsigned::divrem_limb_into` ONLY for GCC x86-64 with __int128, replacing its earlier `>=64` threshold. It does not change the exact preinverse algorithm, the power-of-two fast path, the B61 rem-only >=8-limb selection, multi-limb Knuth division, ABI or public headers. Clang, ARM64, MSVC and short quotient inputs retain their existing selection. No test expectations or official benchmark inputs are changed.

Repeated on-target runs should confirm borderline wins. B61 is a baseline, not proof of B63 outcomes. All B62 Unix/Windows/Coverage GitHub Actions workflows completed successfully; B63 integrated CI and performance validation are pending.

### B69 test ARM64 carry-chain transfer against GCC x86 fused at short widths (20 October 2026)


Cross-machine absolute timing comparisons reflect different CPU and compiler performance and do not isolate an implementation cause.

The candidate signed multiplication was 0.949x/1.070x/0.984x and signed square 1.005x/1.006x/1.010x. This shows repeatable same-semantics wins for String adapter avoidance, while manual signed multiply/square construction offers no consistent benefit. Large x86 operands often show minimal wrapper differences, so blaming sign handling alone is unwarranted.

B69 makes the proposed ARM64-inspired kernel test on the SAME GCC x86-64 host, rather than comparing ARM64 and x86 absolute times. The existing B68 standalone diagnostic gains a third section ONLY when the production library reports both selected_multiply_uses_fused() and has_x86_64_carry_multiply(). The whole-result experiment constructs a canonical BigUnsigned from the eight-word stack buffer; that may allocate and copy differently from the current production owning mul(), so its timings are A/B between *forced variants*, not a direct claim that merely flipping the production kernel will yield exactly the same performance.

All B69 operations use exactly the existing 32 official Integer input pairs at each selected size. The candidate uses the existing BigUnsigned::mul_into_fused and BigUnsigned::mul_into_x86_64_carry entrypoints: no new architecture assembly, no production dispatch, ABI, public API, official benchmark workload/weights, oracle or other test behavior is changed.

Do not confuse the ARM64 success with proof x86 carry will win. B69 integrated build and host measurements remain unverified at commit authoring time.

### B71 test short owning stack construction inside production library TU (20 October 2026)


The near-unity and changing order do not support changing the GCC dispatch. These mixed figures further suggest compiler/code-layout or benchmark-run effects rather than a consistently better arithmetic kernel.

B70's noinline fused-owning variants incurred ~7.8%/7.5%/6.1% overhead for mul and ~9.9%/7.3%/7.1% for square at 53/128/256 bits. Noisy but repeatable same-run noinline production core wrappers only saved about 0.972x/0.972x/0.986x in nonaliased mul and 0.948x/0.950x/0.990x in self-square. They do NOT reconstruct all code generation or ABI effects across translation units. B70's simple production-vs-forced owning comparison has a larger gap, but cannot yet attribute it to instructions alone.

B71 explicitly moves a new diagnostic, exact short owning BigUnsigned product control INTO THE SAME src/big_unsigned.cpp translation unit as the unchanged production BigUnsigned::mul implementation. For input magnitude lengths at most four limbs, it calculates the selected-kernel product directly in a zeroed eight-limb local std::array, constructs a canonical BigUnsigned from the resulting limb view and returns the fully owned result. A second function does exactly the same but explicitly uses the fused kernel. Both fall back to normal mul() above the four-limb diagnostic limit or on other toolchains, so no production dispatch, public API or ARM64 algorithm is modified.

The existing B67/B68/B69/B70 benchmark gains side-by-side independent median timings of: (1) owning BigUnsigned::mul vs selected stack helper, (2) selected vs forced-fused stack helper, (3) signed owning Integer wrapper calling mul() vs the new stack helper, with equivalent non-inlining, and the same for actual self-square operands. All new helpers return the same exact owning type and are linked as production-library functions, rather than handwritten benchmark-local mul_into entrypoints. They may still exhibit different code placement or inlining across build options; this diagnostic does not create an exact controlled linker/inlining experiment.

No new public API or additional CMake build source is needed. Return of the user's B71 same-host results is required before proposing any short-product production change.

The short-size performance strategy remains: keep fused selected on GCC x86-64 pending robust contrary evidence; keep native ARM64 carry-chain; record direct signed decimal adapter improvements from B68 as a separate possible production change. No reweighted workload was introduced (per user). Run from build in one line: ./benchmarks/benchmark_integer_small_operations_profile -t dark --bits 53,128,256 --iterations 20000 --samples 31 B71 compilation/CI/host measurements are pending.

### B71.1 fix B71 benchmark wrapper declaration order (20 October 2026)


The user reported a native GCC x86-64 benchmark build failure in `benchmark_integer_small_operations_profile.cpp`: `verify_b69` referenced `b71_stack_selected_signed` and `b71_stack_fused_signed` before either was declared. Move the two existing diagnostic function definitions intact above `verify_b69` in the same anonymous namespace. No algorithm, oracle, test pool, production arithmetic, dispatch, build flag, official 84-case benchmark, or benchmark timing methodology changes. The fix does not claim performance or successful native build before actual CI/user verification.

### B72 production short owning stack product on GCC x86-64 (20 October 2026)


At 53 bits, unsigned multiply/square candidate/baseline were 0.950x/0.954x and signed owning 0.981x/0.978x. At 128 bits they were 0.978x/0.981x and 0.969x/1.006x. At 256 bits they were 0.928x/0.884x and 0.925x/0.902x. Selected/fused stack controls were nearly equivalent, sometimes reversing rank. Thus changing fused to carry-chain, which the prior B69 tests showed consistently slower, is not justified. The notable short owning gain at 256 bits appears linked to result construction rather than a faster multiplication kernel. The B70 inline/noinline controls have mixed results (including reversed rank at 128 bits in B71); do not treat these as a precise additive call cost.

B72 promotes the precisely measured B71 diagnostic short owning construction to **production BigUnsigned::mul** ONLY under the existing GCC x86-64 fused-build guard and only for nonzero canonical inputs of at most BigUnsigned::inline_capacity=4 limbs each. It uses an eight-limb zeroed std::array, writes the exact full product with the UNCHANGED selected mul_into fused kernel, then constructs a fully owned normalized BigUnsigned from a view of the local buffer. Results larger than four limbs still use the existing owned vector representation. If either input is >=5 limbs, the original owning construction/Karatsuba selection remains untouched. Zero returns at the original early check; size arithmetic overflow checks retain their original order. No global buffers, special benchmark input assumptions, numerical approximation, public header/ABI, FloatMP arithmetic, ARM64, Clang or MSVC selected algorithms are changed. The B71 stack selected/fused diagnostics remain available as controls, although the B71 production baseline at <=4 limbs will now use the promoted variant; compare with B71.1's pre-promotion commit for historical A/B reference.

The existing 6000-random multiplication oracle and all previous tests remain intact. The original official 84-operation benchmark and its 32-case input pools are unchanged. ARM64 baseline preservation requires its own build/performance verification. B72 CI, compilation and measured production performance are pending at publication time.

### B74 short owning arithmetic and generic buffered powers (20 October 2026)


The user's post-B73 x86-64 20000x31 official run verifies all 768 paired inputs and 2688 results. The remaining native-sum excesses are 4.352/122.464/178.269 ns at the three short sizes.

B74 targets owning intermediate work, retaining the existing kernels:

- One-limb unsigned sums return a scalar magnitude, with a separate two-limb carry result when needed. Subtraction whose significant range shrinks to one limb returns that exact scalar difference. Canonical size, not stale moved-from inline storage, decides zeros.
- Owning `quot`/`rem` with at most four dividend limbs use a private shared template. Shorter dividends return q=0/r=a; scalar inputs use exact unsigned division/remainder. Equal-width q=0/1 cases use B36's existing high-limb criterion before scratch preparation. The other cases call the UNCHANGED selected division template with local output of four limbs and scratch of nine limbs (n+m+1<=9). Zero-divisor checks still precede every dispatch. Joint divrem, low-level validation/scratch contracts, >4-limb owning division, normalization, preinverse thresholds and ISA kernels are unchanged. These storage/scalar paths are shared C++20 on all architectures.
- A new exact `pow(BigUnsigned const&,uint64_t)` uses generic MSB-first exponentiation. If base.bit_length()*exponent<=1024, two alternating local buffers replace owning intermediate results. The bound is checked by division, avoiding overflow; zero, exponent zero (including 0^0=1), and exponent one are handled separately. Every intermediate exponent is at most the requested exponent. Although each canonical intermediate fits sixteen limbs, ceil(x/64)+ceil(y/64) may exceed ceil((x+y)/64) by one. Each buffer therefore has SEVENTEEN limbs (272 bytes in total), allowing the full kernel output before normalization. The input is borrowed immutably, buffers/output are disjoint, and every selected kernel fully overwrites its output before it is read. No uninitialized limb is read and no stack view escapes: only the final BigUnsigned owns/copies storage. The >=8-limb symmetric-square selection and selected fused/carry/portable multiplication are preserved. Above the conservative budget, ordinary owning MSB-first square/multiply remains arbitrary-size. The budget is a bounded-storage policy, not a measured optimal crossover or an exponent-three special case.
- `Integer::pow` calls this unsigned engine only under the existing GCC x86-64 fused-build guard, keeping the same sign parity and canonical zero. Other signed front ends retain their old schedule. No public existing signature/layout or division convention changes; the unsigned pow overload is new.

Validation on Apple M2 Ultra / AppleClang 21 Release warnings-as-errors: 20/20 CTest pass; complete Integer and BigUnsigned Debug ASan/UBSan suites pass; complete GCC 16.2 native ARM64 BigUnsigned suite passes. Forced-portable (no native backend, __SIZEOF_INT128__ undefined) ASan/UBSan runs the actual power, short-division, arithmetic-boundary and common-prefix subtraction tests successfully. The forced build suppresses only the pre-existing unused fused-helper warning; normal native builds keep warnings-as-errors. No MSVC/x86 native build is claimed. The 65/341-bit dense powers explicitly exercise the seventeenth limb. All earlier independent and randomized tests remain intact.

Both repeated unchanged ARM64 official 20000x31 runs verify 768 paired inputs and 2688 results. Signed ARM64 power is intentionally unchanged; the separate same-signed profile verifies buffered powers before/after timing. Its MSB-control -> buffered ns are 27.156 -> 27.654 (53), 61.435 -> 56.331 (128), 90.904 -> 56.292 (256): neutral at 53, modest at 128, larger at 256. These ARM64 controls do not project native GCC x86-64 gains and are not additive components of the official total.

2048/4096 still fail strict ARM64 totals in these runs.

The unchanged FloatMP 1000x11 four-size benchmark verifies 256 paired inputs and 3456 results for both the saved B73 library and B74. These runs do not establish FloatMP performance acceptance; its pre-existing short-size limits remain unmet. No floating-point arithmetic/rounding algorithm changed.

Native GCC x86-64 and Windows CI/performance remain pending. The binding acceptance is still all SIX official per-size totals <1x and EVERY operation/size <2x, with stable same-host repeated runs. No benchmark corpus, weights or oracle was altered.

### FloatMP D143 — skip a known-zero leading quotient digit (23 October 2026)


In its ordinary division loop, a leading dividend limb smaller than the nonzero divisor has quotient digit zero and is already the exact initial remainder. Initialize that digit/remainder directly and divide only the remaining words. Empty inputs and trailing output zeroing retain their contracts; preinverted and power-of-two selections are unchanged. This portable shortcut reuses all existing architecture-specific divide_wide implementations and adds none. See floatmp_design.md D143 for the encompassing branch and its local performance/validation scope.

### B76 — freeze all 2,800 seeded BigUnsigned GCD/decimal oracles (23 October 2026)


Preserve its entire std::mt19937_64 seed 0xB16B00B5, conditional draws for both operand lengths, exact 0/20/32/96-limb width selections, fills and override order (zero limbs / UINT64_MAX), all 2800 cases, input immutability, exact unsigned GCD, canonical bit-length, decimal formatting and decimal roundtrip checks. No sampling.

The 14 TSV fixture files tests/data/big_unsigned_b76_gcd_decimal_ <chunk>.tsv contain all 2800 original source visits in ordered 200-row chunks. Every record freezes iteration, raw requested n/m lengths, both exact unsigned source magnitudes in hexadecimal, their exact Euclidean gcd in hexadecimal, and the independently derived source magnitude in base 10. The first 12 source visits were crosschecked against a separately compiled C++20 std::mt19937_64 generator, including the width override, zero/MAX fills and conditional random draws. Sources and outcomes are checked against TSV on every visit, with strict row, file boundary and extra-row validation. Also verify gcd(y,x) and source operands after arithmetic to retain and strengthen borrowed-source immutability/commutativity checking. The original to_decimal/ from_decimal canonicality and roundtrip assertions remain.

The subsequent 1800-visit random decimal parser test and all later seeded radix parse cases continue from the identical generator state. No change to production code, CMake or benchmarks. Cross-platform CI is required before B76 is considered validated.

### B77 — reproducible standalone C++20 oracle generators (23 October 2026)


Version tests/data/generator_big_unsigned_gcd.cpp, generator_big_unsigned_b76.cpp, oracle_exact_integer.hpp and README-generators.md. They independently regenerate all 15 GCD and 14 B76 GCD/decimal TSVs using independent base-2^64 limbs, binary GCD, exact decimal conversion and unchanged std::mt19937_64 streams. Both GCC and Clang C++20 warnings-as-errors builds reproduced all 29 committed TSV files byte for byte without changing them. By default --verify does not write files; --write explicitly regenerates them. tests/CMakeLists.txt compiles both programs as dependencies of the aggregate tests build target, without registering them for CTest or running them during compilation. Other TSVs, including certified FloatMP transcendental fixtures, still have no versioned generator. Full CI remains required.

### B77.1 — strict Clang signed-conversion fix for independent oracle helper (23 October 2026)


Fix `independent_oracle::Nat::trailing_zero_bits()` in `tests/data/oracle_exact_integer.hpp`: `std::countr_zero(uint64_t)` returns `int` by contract, while the bit offset is `std::size_t`. Cast the nonnegative count explicitly to `std::size_t` and use an unsigned 64-bit-word multiplier. The exact expression was compiled and executed with GCC and Clang under the reported warnings-as-errors flag set. Full repository CI should still be checked independently after this commit.

### B78 — version five independent exact data generators (23 October 2026)


Add standalone C++20 `tests/data/generator_exact_small.cpp` for the already frozen datasets `rational_r3.tsv` (6 ordered comparisons), `rational_r4.tsv` (54 exact floor/ceil samples), `rational_r8.tsv` (16 near-2^64 rational comparisons), `integer_power.tsv` (100 signed powers), and `big_unsigned_r14.tsv` (28 squares). R4 retains the canonical 1/1 case for base/base; R14 consumes all 12 source-RNG visits per width before selecting the original four structural classes. No extra sampling, seed changes or test modifications are introduced.

GCC 14 and Clang 17 with full strict warning flags (including -Wsign-conversion -Werror) compiled and verified the generator. All five regenerated datasets matched their original Git blob SHA-1 IDs exactly; `--verify` rejected an intentionally corrupted test TSV without writing it. Only the generator, its build-only registration, the generator README and this journal changed. The independent source corpus is now reproducible for 34 TSV files (29 B77 plus five B78). Full Unix/Windows/Sanitizers/Coverage CI remains to be checked.

### B79 — reproduce BigUnsigned D104 power and D103 short-division fixtures (23 October 2026)


Version standalone `tests/data/generator_big_unsigned_power_div.cpp` and register it as a build-only generator target. It independently reproduces `tests/data/big_unsigned_power.tsv` (373 existing source/result rows) and `tests/data/big_unsigned_short_div.tsv` (324 existing rows) using only standard C++20 and the independent unsigned-limb B77 oracle helper. Multiplication and powering use binary shift/add and exponentiation by squaring; quotient and remainder use bit-aligned repeated subtraction. For division, each reference verifies `a=b*q+r` and `0<=r<b`, with a separate explicit no-result representation for zero divisors.

For the D104 random corpus, replay each of the original `mt19937_64` (0xB7400A11) source visits across 10 widths and 12 modes, including unselected modes, before reproducing exactly the same five selected source modes and six selected exponents, followed by all 16 small sample and 57 1024-bit-boundary power rows. For D103, replay the entire seeded (0xB74D1A1DE) 6 x 6 x 24 source corpus and its original overwrite order; reproduce all nine selected modes per width pair, including the exceptional zero-divisor cases. No sampled row was discarded or added.

GCC 14 and Clang 17 compiled and passed `--verify` with the full strict CI-style warning list (`-Wsign-conversion -Wshadow -Werror` included). Both regenerated files exactly matched the existing Git blob hashes (c0dd421950bbdcfbc86adbda19d76c1773f84d61 and 5a0c2d9d5b2bd19f5ab6f0b44415ca10e6e27b8e). Corrupting a fixture causes read-only verification to fail. The generator README is updated; 36 committed TSVs now have versioned independent generators. Full repository CI (Unix, Windows, Sanitizers, Coverage) is still required.

### B82 — reproduce B72 and R17 source-keyed exact fixtures (23 October 2026)


Add standalone `tests/data/generator_exact_b72_r17.cpp` using the standard-library independent radix-2^64 Nat helper. Reproduce the existing `big_unsigned_b72.tsv` across all 36 (n,m) width pairs and all eight selected mode classes, 288 rows. Replay all 32 original `mt19937_64(0xB7204455)` source visits per pair to preserve the stream, including modes not frozen in the TSV. Compute all products and source squares with an independent binary shift-and-add limb multiplier. Reproduce the existing `rational_r17.tsv` across all 7 widths and all 6 x 6 signed rational orderings, 252 rows, by exact sign handling and cross-product comparison. Inputs are derived from the original test source; no TSV inputs are re-used as an oracle.

Compiled and generated both TSVs with GCC 14 and Clang 17 C++20 `-Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror`. All generated Git object IDs match the previous frozen blobs: `big_unsigned_b72.tsv` 87e5d405d7374f9a19963936b56761e089417f4a; `rational_r17.tsv` 3b53a0dc9a8c0908083ece7039778c2631c19696. Register generator as a build-only tests dependency and update README. 44 of 115 TSV files now have independent versioned generators. Full CI remains to be checked.

### B83 — dedicated generators and explicit generate-data CMake targets (23 October 2026)


Separate the independently compiled C++20 oracle programs from the aggregate `tests` build target in `tests/CMakeLists.txt`. Register a dedicated `generators` target with all six currently versioned oracle generators and set their executable targets EXCLUDE_FROM_ALL so that ordinary builds and `tests` do not compile them implicitly. Register an explicitly invoked `generate-data` target depending on `generators`, which executes each generator in a fixed sequence with `--write` and `--data-dir` set to the absolute source-tree `tests/data` directory. The commands are intentionally non-incremental and run on each target invocation, allowing the developer to measure overall execution time externally; no runtime measurement is embedded in the build logic. Single-command-sequence execution stops on the first generator failure, without parallel writes to source-tree TSVs.

The target covers only the 44 TSVs already reproduced by the six generators, leaving the remaining 71 versioned TSVs as they are. There is no automatic invocation during compilation or CTest, and no switch to build-tree fixture paths. Update `README-generators.md` with the explicit targets and warning that `--write` overwrites source-tree files.

Verified the target dependency and command-sequence behavior with a standalone CMake/Ninja probe: default `all` and `tests` build no generator; `generators` builds both sample executable targets; and `generate-data` builds dependencies and runs every `--write` command with the absolute source-data directory. This probe is not a full Numeric CI run.

### B85 — independent BigUnsigned B35/B37/B38 golden generators (23 October 2026)


Add `tests/data/generator_big_unsigned_mul.cpp` as a standalone C++20 source for `big_unsigned_b35.tsv` (14 x 5 = 70 cases), `big_unsigned_b37.tsv` (13 x 7 = 91 squares) and `big_unsigned_b38.tsv` (14 x 8 = 112 multiplication/truncation cases). The generator replays the complete original C++ mt19937_64 streams including all unselected visits: 32, 96 and 24 per original width pair, respectively. It reproduces all historic transformations and ordering, then derives products using independent binary shift-and-add unsigned-limb arithmetic and six low-limb results by exact truncation of the product.

Compiled with GCC and Clang under C++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror. Regenerated Git blob IDs exactly match the previously committed fixtures: B35 98e2be0ee5f55d17761a359e6501c9892a2259ce; B37 9ccb32059de5c5ae6664e6eaa7066ac2f87995a9; B38 6dbf2deb93bd07e41a57c44e05d599c4136113a4. The additional fifth B35 header/comment line is preserved exactly. Register the eighth executable through ORACLE_GENERATORS to make it available in the explicit generators and generate-data targets without coupling it to tests. The 53rd of 115 original frozen TSVs is now reproducible. Full CI remains pending.

### B86 — B39/B40/B54 source-reproducible multiplication fixtures (23 October 2026)


Add independent standard-C++20 `tests/data/generator_big_unsigned_dispatch.cpp` for B39 (19 x 8 = 152 cases), B40 (18 x 9 = 162), B54 products (17 x 8 = 136) and B54 squares (13 x 8 = 104). Original 20/16/32/64 mode cycles and source RNG draws are fully replayed; B54's product and square stages share the seeded mt19937_64(0xB54A6C4) stream, including all 17 x 32 product-source iterations before its first square.

The original four TSV object IDs were reproduced byte-for-byte: B39 82c2d5dfa3c171bc2d0588393aa1ee3a5a6dd671; B40 5f603317b54f2fd19178d8fcb061d0799a25eb38; B54 product 2bf1ee6c6259df30c627ad249ac963f16ea3641c; B54 square ed71c32efb60aea636775edc61f4131b7ac90c61. GCC and Clang compiled and verified the generator under stringent C++20 warnings including -Wsign-conversion -Wshadow -Werror. An intentionally corrupted fixture was rejected by read-only --verify. Register the ninth executable in the ORACLE_GENERATORS list for explicit `generators` and `generate-data`; retain decoupled `tests`. Coverage rises to 57 of 115 original versioned TSVs. Full multi-platform repository CI remains to be checked.

### B92 — complete seeded BigUnsigned randomized multiplication fixture (23 October 2026)


Add `tests/data/generator_big_unsigned_random_mul.cpp` as a standalone, standard-library-only C++20 generator. Reconstruct the initial shared mt19937_64(0x1282561024) stream state by replaying all 7 x 12 R14 source visits, then iterate all 6,000 original random multiplication source rounds. Preserve the two operand length draws, every word, all-ones and zero transforms keyed to original repetition, and the conditional extra RNG draw for low length when repetition%5 == 4. Freeze exactly the 600 original adaptive one-per-ten sample offsets, but compute all source operands, exact full products, squares and low products independently using radix-2^64 binary shift-and-add arithmetic. No frozen expected product or input data is embedded in the generator.

GCC 14 and Clang 17 both compiled with C++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror. Generated exactly 942,214 bytes with original Git SHA-1 `e06e8b7f8caf0eaaaa42bce30894561063226eac`. An intentionally corrupted copy is rejected by read-only --verify. Register the fifteenth generator in explicit `generators` and `generate-data`, leaving tests and default all decoupled. Total exact source-reproducible fixtures: 95/115. Full CI pending.

### B93 — independent FloatMP directed rational bounds generator (23 October 2026)


Add standalone standard-C++20 `tests/data/generator_floatmp_directed_bounds.cpp` which reproduces every original input and exact directed dyadic bound in `floatmp_directed_bounds.tsv`: six significand precisions (7, 32, 53, 64, 128, 192) times fifteen signed, zero or dyadic rational inputs. Determine floor(log2|q|) using exact unsigned limb bit lengths and shifts; independently floor/ceil the magnitude on its p-bit lattice using long quotient/remainder; reverse bounds for negative values; reduce resulting rational numerators/denominators using binary GCD and exact division. Source inputs are reconstructed from the original test's explicit 15-element input set, not copied from the frozen TSV.

GCC 14 and Clang 17 -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror compiled and ran successfully, producing identical 7,742-byte fixtures matching the original SHA-1 `ef8e337b8d4271371c919044d3eec16ca6fb731c` byte for byte. An optimization-specific GCC array-bounds warning from a temporary Nat copy was eliminated by constructing the scaled Nat from its original unsigned source rather than disabling warnings. Register the sixteenth generator in `ORACLE_GENERATORS` for the explicit `generators` and `generate-data`, not `tests`. Source-reproducible fixture coverage rises from 95 to 96 of 115. Full CI remains pending.

### B94 — independent certified integer-only FloatMP pi fixture (23 October 2026)


Add standalone C++20 `tests/data/generator_floatmp_pi.cpp` to reproduce `floatmp_pi.tsv` from Machin's formula pi=16*atan(1/5)-4*atan(1/239). Derive a rigorous integer-only fixed-point interval at scale 2^1200 using 400 Taylor terms per arctangent. Each series term is bracketed by its exact integer floor and ceiling, with a rigorous alternating next-term tail bound. Signed interval combination yields bounds for pi. Prove both endpoints lie in a common 384-bit directed-rounding cell; prove both give the same 320-bit nearest-even result. Emit the exact canonical reduced rational lower/upper and nearest-even expected values.

GCC 14 and Clang 17 compiled under -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror. Both generated the identical 1,250-byte original blob SHA-1 `33eccacd41fa17d8453c098ede5cab4c662b5fc5`, verified with Git hash-object and --verify. Register as generator #17 for explicit `generators` and `generate-data`; default builds and `tests` remain independent. The number of reproducible TSVs increases from 96 to 97 of 115. Full CI pending.

### B95 — independent two-stage FloatMP get_d binary64 oracle (23 October 2026)


Add `tests/data/generator_floatmp_get_d.cpp` as standalone C++20 source for all 272 original `floatmp_get_d.tsv` conversions. Recreate the original 17 input Rationals by their exact algebraic definitions: integers and thirds, binary halfway ties, extremely small dyadics, positive and negative exponent boundaries and the IEEE overflow tie. For precisions p=53,128,256,1024, first nearest-even round the original rational to the exact p-bit binary dyadic value. Then perform a separate correctly rounded conversion to IEEE-754 binary64 under N/D/U/Z, including biased normal exponent/mantissa, 2^-1074 subnormal grid, negative zero, positive/negative infinity, highest finite value and signed/toward-zero overflow saturation. All steps use independent radix-2^64 integer long division, exact remainders, bit lengths and explicit midpoint-even tie rules.

Independent Python BigInt and C++20 implementations both generate the same original 272 source-keyed records. GCC 14 and Clang 17 compiled under -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror; both emitted exactly 65,751 bytes matching original Git SHA-1 `7365a0f27accc0b364bf3f6de1db4dc06a0084c2`. The optimization-specific false-positive GCC array-bounds warning caused by returning a copied zero Rational was avoided by rejecting zero at the first-stage source as no zero source is present; all real zero/subnormal conversion paths remain tested at stage two. Register generator #18 in explicit `generators` and `generate-data`, without coupling `tests`. Source-reproducible coverage is now 98/115. Full platform CI remains to be checked.

### B96 — independent exact FloatMP integer-power generator (23 October 2026)


Add standalone `tests/data/generator_floatmp_power_exact.cpp` for all 2,400 original `floatmp_pow.tsv` keys: five p-bit precisions (7, 53, 128, 256, 1024); eight original input magnitudes with both signs, with duplicate signed zeros deduplicated; eight degrees (0,1,2,3,7, 8,13,31); four directed rounding modes. Inputs are reconstructed from the original source enumeration. Each signed source Rational first undergoes exact nearest-even input quantization to p binary digits. Compute the exact integer power using independent binary shift/add radix-2^64 multiplication, keeping the exact dyadic scale, followed by one final p-bit round with exact discarded-bit/tie comparison for N/D/U/Z modes. All expected numerators/denominators are canonicalized by stripping common factors of two.

GCC 14 and Clang 17 both compiled using C++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror. The generated TSV contains exactly 2,400 entries, 106,895 bytes, and matches the original Git blob SHA-1 `ae8049364ce05bfae1e381af2ba68127a7de1a21`. Clang read-only --verify accepted the independently generated GCC file; performance measured locally around 0.34-0.41 seconds exclusive of compilation. Register executable #19 for explicitly invoked `generators` and `generate-data`; leave `tests` decoupled. Coverage increases to 99 of 115 TSV fixtures. Full CI remains to be checked.

### B97 — exact certified FloatMP sqrt reference generator (23 October 2026)


Add standalone C++20 `tests/data/generator_floatmp_sqrt.cpp` to reproduce `floatmp_sqrt.tsv` from all original 14 positive Rational inputs, including source-specific p-bit steps around one/four, plus the four shifts (-2000,-127,127,2000) of the nearest-rounded source 5/4. Five significand precisions (7,53,128,256,1024) and all four output modes preserve all 360 original test visits. Exactly 300 source-keyed references remain after explicit duplicate rounded-input removal, matching the original fixture ordering and comments.

Independently nearest-even round source operands using exact radix-2^64 quotient/remainder. For the positive dyadic x at precision p, compute its square-root target unit from floor(log2(x)) with correct floor semantics for negative binary exponents. Certify both q^2*B<=A and (q+1)^2*B>A, compare exact 4*A against (2*q+1)^2*B for nearest-even, and perform outward rounding for upward while keeping downward and toward-zero exact. Normalize all output fractions by cancelling binary powers.

GCC 14 and Clang 17 strict C++20 warning builds and execution produced exactly 76,870 bytes and Git blob SHA-1 `b83956f9533a77c2dc8465f69ef0ebf47055a198`, matching the frozen fixture byte for byte. Clang read-only --verify accepted GCC's result. Register generator #20 in explicit `generators` and `generate-data`, not `tests`. Generator coverage reaches 100 of 115 TSV files.

### B98 — independent short integer-power boundary generator (23 October 2026)


Add standalone C++20 `tests/data/generator_floatmp_pow_short.cpp` for all 3,680 source-keyed exact signed power references in `floatmp_pow_short.tsv`. Reconstruct low and high binade edges as 2^(p-1)/2^p, 2^p/2^p, 2^(p+3), 2^(p+4) and each near-edge source offset derived by the original degree-dependent shortcut threshold. For degrees 2,3,7,8,13,16,17, choose offsets 0,1, t-1,t,t+1, excluding duplicates and invalid binade crossings. The lower threshold uses 2^floor((p-3)/2) for p>=3 while the upper threshold uses 2^floor((p-2)/2); at even p=64 these differ, which is required to preserve the exact original boundary corpus. For p=2 the lower power is clamped to 2^0 before dividing by the positive degree. Both source signs are represented. Sort by exact signed scaled operand magnitude, then exponent degree, under the original D/N/U/Z lexical mode sequence. All candidate input magnitudes, source fractions and expected results are recomputed from boundary formulas, not extracted from original TSV content.

Compute the original input's independent nearest-even p-bit quantization, exact dyadic integer powering via shift/add radix-2^64 arithmetic, then one p-bit final signed N/D/U/Z round using discarded integer remainder and midpoint ties. GCC 14 and Clang 17 compile with -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror. Generated 3,680 rows and 462,642 bytes, exactly matching the original Git blob SHA-1 `9b980ec173a529f89c78f8ddb83d6870da27d755`; Clang --verify matches GCC's generated bytes. Register standalone generator #21 in explicit `generators` and `generate-data`, not `tests`. Coverage is now 101/115. Full platform CI pending.

### B99 — independently reproducible RoundedFloat directed fixture (24 October 2026)


Add standalone standard C++20 `tests/data/generator_rounded_float_directed.cpp` for all 24 exact source rows / 48 strict comparisons in the existing `rounded_float_directed.tsv`. Reconstruct all original exp/log/atan rational sources algebraically. Use signed integer fixed-point intervals at B=2^1400, with outward rounding at each arithmetic operation; exp uses 241 Taylor terms on |x|/16 plus a twice-next-term tail and four squarings (negative arguments via reciprocal), log uses 271 positive atanh terms plus an explicit 4*y^543 tail, and atan uses 271 alternating terms and Machin-certified pi/4 reduction for inputs near +/-1. Certify that both interval endpoints agree when rounded independently downward and upward to 128 binary significant bits, including negative values and canonical reduced dyadic fractions.

Strict -std=c++20 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror builds with GCC and Clang regenerate 5,181 identical bytes; Git blob SHA-1 is `f20fd0de94317e1999cfe1c4d7ac6718c443a9db`, matching the versioned file. A corrupted TSV was rejected by --verify. Register generator #22 in the opt-in `generators` and `generate-data` targets. Add a separate opt-in `verify-data` target running all registered generators with --verify only, never overwriting source fixtures. Coverage reaches 102/115. Full platform CI generator verification remains to be added.

### B100 — exact FloatMP add/sub source reproduction (24 October 2026)


Add standalone C++20 `tests/data/generator_floatmp_add_sub.cpp`, reproducing all 1,780 canonical unique rows in the pre-existing `floatmp_add_sub.tsv` from the ORIGINAL complete five-precision, 10x10 ordered operand corpus and both add/sub operation visits. Source operands are +1/-1, +/-nearest-even(7/11), +/-2^-precision, 2^-200, 2^200, 2^-1400, and 2^1400. Form exact signed dyadic sums and differences with independent 64-bit limb arithmetic, preserving all large exponent gaps. Expected values and input literals are not read from the TSV. Preserve its original eight comment lines, original first-visit row ordering, and 1,780 unique rows. GCC/Clang strict warnings-as-errors builds both produce 1,104,075 bytes matching existing blob `b4552cfbdea659730486dd5bcdb3bb3a62be7984`; read-only --verify passes and an altered fixture is rejected without overwriting.

Register generator #23 in the existing explicit `generators`, `generate-data`, and `verify-data` targets. Coverage reaches 103/115. Full cross-platform CI still must be verified.

### B101 — certified exponential FloatMP fixture reconstruction (24 October 2026)


Add independent C++20 `tests/data/generator_floatmp_exp.cpp`, reconstructing the original 18 signed rational source keys in their precise original order and their 53, 128, 256, 1024-bit four-mode correctly rounded exp results (288 cases). Independent radix-2^32 multiplication within radix-2^64 Nat limbs, exact aligned-subtraction long division and outward fixed-point enclosures at B=2^1400 implement the Taylor series exp(|x|/16) through term 165. Twice term 166 bounds the positive remainder at max |x|/16=7/16; four outward squares and positive reciprocal enclosures reconstruct exp(+-|x|) without any floating-point transcendental oracle. Each endpoint is rounded using an exact independent binary grid and nearest-even parity/direct modes, then canonical dyadic rationals must agree before emission.

GCC/Clang strict C++20 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror compile, both read-only --verify and deliberate fixture corruption rejection passed on the original data. Exactly 68,122 bytes / 288 rows reproduce Git blob `92ce3c0eeda7f73f5e0385abf483f93367792891` byte-for-byte. Register generator #24 in opt-in `generators`, `generate-data`, `verify-data`; coverage 104/115. Previous B100 Unix, Windows, Coverage and Sanitizers workflow runs all completed success. CI for B101 is pending.

### B102 — certified logarithm fixture regeneration (24 October 2026)


Add standalone standard-library C++20 `generator_floatmp_log.cpp` covering all 444 source-keyed `floatmp_log.tsv` rows, including 284 centered-log unique references from six precision visits and 160 D96 additions. Preserve original source order, nearest-even input rounding for centered cases and unrounded original rational lookup keys in D96. Certify log via an integer fixed-point B=2^1800 enclosure of 2*atanh(|x-1|/(x+1)), with 1101 positive series terms plus rigorous 4*y^2205 remainder. At x=8 or 1/8 use 3*log(2) identity to obtain a sufficiently tight enclosure for 1024-bit output. Round both interval endpoints independently under N/D/U/Z and require exact agreement; negate with correct sign-aware directed rounding. Restore the D96 header in its original mid-file position. GCC and Clang strict warnings-as-errors builds and read-only --verify pass, corruption fails. Generated 122,768 bytes with SHA-1 `7cae466575432ff238dd73ecdeeb599b7e17ceed`, matching main's original Git blob exactly. Generator #25 joins opt-in CMake targets; coverage 105/115. No original fixtures or test/production code changed.

### B103 — certified independent sine and cosine fixture regeneration (24 October 2026)


Add one standard C++20 generator `tests/data/generator_floatmp_sin_cos.cpp` reproducing the complete original `floatmp_sin.tsv` (281,134 bytes / 1,080 rows) and `floatmp_cos.tsv` (278,016 bytes / 1,080 rows) with independent arbitrary-width integer arithmetic. Rebuild every source Rational in D91, D92 and D95 from original test rules and retain the exact precision schedule, sign order, original source-keyed map semantics, independent nearest-even input rounding, first-seen duplicate removal and all 4 output rounding modes. D91 adds 336 rows per file, D92 456, and D95 288. Preserve mid-file D92 and D95 explanatory comments.

Evaluate both sine and cosine through 240 signed alternating Taylor terms at fixed-point B=2^1800, computing all term enclosures using portable 32-bit digit multiplication, exact Nat division and directed outward integer rounding. Add a rigorous enclosure of the first omitted term; odd/even identities handle negative input. Independently round both signed interval endpoints N/D/U/Z at the requested output precision, rejecting any case for which they disagree. Caching on the exact nearest-even rounded input magnitude shares the certified sin and cos evaluation across all four modes, signs, precision visits and overlapping test groups.

Strict GCC/Clang C++20 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror compile and --verify pass for both files; corrupted fixture --verify fails. Exact original Git blob hashes: `9b75d8fd4992f62e288c151e8f930dedd66a778e` (sin) and `31475db715e25928b31b34d31e32735e2cbb97be` (cos). Register generator #26 in explicit opt-in CMake targets `generators`, `generate-data`, `verify-data`; coverage becomes 107/115 (93.0 percent). Full CI of B103 pending.

### B104 — correctly rounded centered FloatMP exponential fixture (24 October 2026)


Add independent standalone C++20 `generator_floatmp_centered_exp.cpp`. Reconstruct the original D85 source algorithm: six precision visits (53,7,128,256,1024,53), 15 magnitudes around 1/2 and at 1/4,3/4, 1/3,2/3,1,2, both signs, four rounding directions. Independently round every source Rational to nearest-even at p bits and deduplicate first-seen (p,mode,signed rounded dyadic) keys, giving 568 fixture rows while preserving 720 original checks. All values are rigorously certified by outward fixed-point integer enclosures (B=2^1800) of exp(|x|/16) through term 165, twice-next-term tail, four squarings, and directed reciprocals for negative inputs; interval endpoints must round to precisely the same p-bit N/D/U/Z rational. Strict GCC/Clang C++20 builds, --verify and corruption rejection pass. Regenerated 170,882 bytes have Git blob `a221dfa3cc3e8feff2b425aff8290180fe994d3a`, identical to the original TSV. Register generator #27 in explicit `generators`, `generate-data` and `verify-data`; coverage reaches 108/115.

### B105 — certified atan oracle regeneration (24 October 2026)


Add C++20 independent `generator_floatmp_atan.cpp`. Rebuild D87, D88 and D89 original source-keyed references (640/200/500 rows), including all precision visits, signed operands, exact original Rational keys, and independent nearest-even input quantization. Use signed outward fixed-point integer intervals at B=2^1800, 1,150 alternating atan Taylor terms at maximum, Machin's certified pi/4, direct Taylor at <=1/2, complementary pi/4 argument reduction for 1/2<x<1, and pi/2-atan(1/x) above 1. An early break when the next power interval is below 2^180 (at B=2^1800) preserves a rigorous symmetric omitted-tail enclosure, greatly improving runtime. Every p-bit N/D/U/Z result is accepted only if independently rounded lower/upper interval endpoints agree. Strict GCC and Clang C++20 compilation, --verify and deliberate corruption rejection succeed. Exact original 309,281-byte Git blob `9daa471b286fb3b3b5cd0d4a804a4977cfd16079` is reproduced. Register generator #28 in opt-in targets; coverage becomes 109/115. CI for B105 is pending.

### B106 — independently certified tangent fixture regeneration (24 October 2026)


Add standalone standard C++20 `tests/data/generator_floatmp_tan.cpp`, reproducing `floatmp_tan.tsv` exactly (303,143 bytes / 1,208 rows). Rebuild D92 (576 unique keys), D93 (240 new), D94 (200 new), and D95 (192 new) from their original test loops: local p-bit rational offsets from 1/2, centered p/2+16 offsets, wider binary rational sources, signed order, precision visits and first-seen deduplication. Independently round the original Rational source to nearest-even p-bit binary before evaluating tangent. Use 240 outward alternating integer-only Taylor terms for sin/cos at B=2^1800; enclose omitted term, require cosine interval exclude zero, divide signed numerator interval by signed denominator with four endpoint combinations, rigorously outward floor/ceil, and use oddness at negative sources. Support negative cosine for |x|>pi/2, not merely the near-1/2 window. Require identical rounded exact rational from both endpoint bounds in every N/D/U/Z mode.

Strict GCC and Clang -std=c++20 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror compile, and read-only --verify pass; a deliberately corrupted 303,143-byte fixture is rejected. The reproduced Git blob is exactly `0b73c09cd203961ace02c86c3af8c866db8c7ac4`, identical to the original version-controlled fixture. Add generator #29 to opt-in `generators`, `generate-data`, `verify-data` targets; coverage reaches 110/115, with five originals outstanding. Full CI for B106 is pending.

### B107 — independent certified asin/acos golden regeneration (24 October 2026)


Add `tests/data/generator_floatmp_asin_acos.cpp` using standard C++20-only arbitrary-width integer limbs. Reconstruct all D86 (520) and additional D90 (208) keys from the original signed Rational source and precision visitation rules, rounding each operand nearest even before evaluation. For each original magnitude independently certify sqrt(1-x^2) through floor integer square root of floor((d^2-n^2)*B^2/d^2), including exactness detection for the upper interval endpoint, where B=2^1800. Both interval endpoints must round to identical binary rationals for all N/D/U/Z modes.

GCC and Clang strict C++20 warnings-as-errors compiles and --verify of both fixtures passed; intentional corruption is detected. Exact original Git blobs: asin 189,909 bytes / 728 rows, `50d8070f7cc049ad76840f51f2d4042ce9f73cb4`; acos 189,361 bytes / 728 rows, `9bc0e3b101c737783565f4dc67c0033ec6fa7c91`. Register generator #30 in opt-in CMake targets. Generator coverage is 112/115, with three original TSVs outstanding. CI for B107 pending.

### B108 — certified independent tanh fixture (24 October 2026)


Add standard-C++20-only `generator_floatmp_tanh.cpp` covering the full 1,148 rows of `floatmp_tanh.tsv`, rebuilding 13 source precision families, all signed inputs sorted exactly as originally, and four rounding modes. Recover the original near-half cutoff selections, including the separate ceil((p+33)/2) rule for the wide precision regression cases, and original nearest-even input quantization. Compute tanh(x) independently with an outward integer interval for exp(2|x|), via exp(|x|/32) Taylor through 409 terms and six directed squarings at B=2^2304, followed by monotone interval fraction (E-1)/(E+1). For |x|=1000, use exp(2000)>2^1800 to prove 0<1-tanh(x)<2^-1799, preserving strict directed rounding at p<=1024, including correct 1 predecessor. Directed negative values use sign-aware mode swapping, and all other certified interval endpoints must round identically at requested p.

Strict GCC and Clang C++20 warnings-as-errors, --verify and corruption rejection passed. Exact original Git blob SHA-1 `7bc870f9e9f1006e3f5c0850c6e1eed1c5e68a2b` and 345,917 bytes. Generator #31 joins opt-in `generators`, `generate-data`, `verify-data` targets. Coverage reaches 113/115, leaving `big_unsigned_sqrt.tsv` and `floatmp_rounding.tsv` outstanding.

### B109 — independently reproduced exact FloatMP rational rounding (24 October 2026)


Add standard-C++20-only `tests/data/generator_floatmp_rounding.cpp`. The D79 section uses original exact source generation rules, including dividing operands after nearest-even input quantization and doing N/D/U/Z arithmetic with operands after the corresponding mode's quantization. Preserve original deduplication, five compound operations, midpoint ties-to-even and sign-aware directed rounding. No TSV data is read to construct source inputs or expected numbers. Implementation uses independently implemented exact radix-2^32 multiplication on Nat radix-2^64 limbs, signed exact rational arithmetic, binary GCD, quotient/remainder and p-bit rational canonicalization. Strict GCC and Clang builds and --verify pass and corrupted-file validation fails. Generated 4,794 rows / 1,225,939 bytes match Git blob exactly: `8d5c927485bc2c171853bb1502ed6b9d0d2aa51b`. Register #32 in opt-in generator and verify-data targets, raising coverage to 114/115. No originals, consumer tests or production algorithms were changed. Remaining big_unsigned_sqrt.tsv still requires recovery of the original random source algorithm.

### B110 — 115/115 fixtures with qualified historical entropy recovery (24 October 2026)


Add standalone C++20 `generator_big_unsigned_sqrt.cpp` plus the input-only `big_unsigned_sqrt_legacy_seeds.hpp` bank. Reconstruct 21 initial small edge integers, an extra five, 69 powers-of-two boundary inputs from the original widths, 37 historical random magnitudes, and 33 boundary cases from 11 original square-root magnitudes. The last 48 sources are frozen as seeds ONLY because the original 2026-10-23 random generator algorithm and seed were not retained in committed source/history. No expected roots/remainders are frozen in source. Recompute all exact floor root and remainder using independent binary-pair restoring arithmetic, and validate with a separate binary shift/add square and strict successor bound. JavaScript BigInt independently reconstructed the entire 161-row, 48,161-byte original fixture byte-for-byte, matching Git blob `686e648fbe32f70303b1b01c38629bd03e3672a3`. The standalone C++ generator builds in strict C++20 GCC/Clang warning-as-error mode; full C++ fixture verification on actual committed seeds and full CI still require confirmation. Do not misrepresent this as recovery of the unrecoverable pseudorandom stream.

Register generator #33 into opt-in `generators`, `generate-data` and `verify-data`; all 115 existing fixtures now have a registered generator, with the above 48-source provenance qualification.

### B111 — consolidated rounding generator, deterministic sqrt corpus (24 October 2026)


Reduce generator executables 33 -> 32 by merging D80 add/sub into `generator_floatmp_rounding.cpp`. The shared standalone C++20 source outputs both pre-existing exact TSVs without byte changes: rounding 4,794 rows SHA `8d5c927485bc2c171853bb1502ed6b9d0d2aa51b`, add/sub 1,780 rows SHA `b4552cfbdea659730486dd5bcdb3bb3a62be7984`. Delete the redundant generator executable/source but keep all original D80 visits, modes and first-visit deduplication.

Replace the 48 historic non-reproducible sqrt entropy inputs by a new SplitMix64 recurrence directly in `generator_big_unsigned_sqrt.cpp`. Document 64-bit wrap, seed `0x526f6f74534d3230`, limb generation, width masking and exact top bits, and sparse low-word-zero variant. Delete frozen seed header. Keep all 161 test rows and original size classes: 91 fixed edge/power-of-two cases, 37 generated numbers, 33 neighboring-square cases from 11 generated roots. Expected roots and remainders for 70 changed input magnitudes are regenerated from independent restoring arithmetic; Python's math.isqrt cross-checks all 161 triples and the PRNG is separately reproduced. TSV Git blob `5ed4a8d6705a2b39dd44d5ebcf9f5f666019d2bd`, 48,272 bytes. Both GCC/Clang strict C++20 generator builds, --verify, and corruption rejection pass; no consumer or production source changed.

### A1 — independently generated product-accumulation fixture (24 October 2026)


Preserve the original full mt19937_64(0x23ACC00) draw stream, original source transforms including empty operands, two directional accumulate calls, exact output-buffer lengths and exceptional unchanged-buffer assertions. Consumer reads, validates every original source key, result and escaped bit with an independent canonical BigUnsigned-from-hex parser, checks full buffer size, source immutability and no missing or extra records. The existing independent shift/add product routine of `generator_big_unsigned_mul.cpp` also regenerates this dataset to minimize number of executables. Independent JS BigInt + mt19937_64 reference matches all original 1,728 visits: 1,441,196 bytes and Git SHA `2e1d31c8b7e9339693ddbaa7f22d2741d213012e`. Full CI validation of consumer and generator is pending.

### A2 — B32 runtime external reference square numerical oracle migrated to TSV (25 October 2026)


The seed `mt19937_64(0xB3205A7E)`, 11 original widths, all 192 repetitions, all random draws (including overwritten draws), and the ordered 13/17/19/23 divisibility transforms are unchanged. Each row freezes original width, repetition, canonical input hex and exact square. The consumer checks every source key and input, both `square(x)` and `mul(x,x)` against the independent exact expected square, source immutability and canonical shape; it rejects missing or extra rows. The existing `generator_big_unsigned_mul.cpp` uses its independent exact binary shift/add multiplication to produce the 2,112-row B32 fixture alongside B35/B37/B38 and A1: 117 TSVs, still 32 generator executables. The independently reproduced B32 fixture is 3,328,858 bytes, Git blob SHA `a9afcd45514a36df0faf7167fa275bb7569d6d1f`. GCC 14 and Clang 17 strict C++20 standalone source builds regenerate identical B32 bytes. Full CI/consumer validation is reported separately from generator-level validation; do not count it as passed until the test jobs complete. No production code was modified.

### A3 — 128/64 edge and B26 independent exact-reference fixtures migrated to TSV (25 October 2026)


Migrate all 192 original 128/64 edge visits (8 divisors x 4 high-word choices x 6 low words) and all 2,048 B26 randomized DIVQ visits in `TestBigUnsigned::test_division_limb()`, without reducing coverage. The original `mt19937_64(0xB260D1B1)` source, conditional RNG draws, divisor/high/low selections and high<divisor precondition are preserved. The consumer verifies every ordered source key, canonical numerator and operands, exact full 64-bit remainder and quotient, unchanged source, original saturation assertions, and missing/extra TSV rows.

The existing standalone C++20 `tests/data/generator_big_unsigned_power_div.cpp` now also emits `big_unsigned_div128by64.tsv` using independent aligned-subtraction long division and binary shift/add multiplication for the identity numerator=divisor*quotient+remainder and remainder<divisor. Exactly 2,240 rows, 227,016 bytes, Git blob `f9084b5d88d01b2acf356e5f4b0cb9a0435ff7cd`. Strict GCC/Clang C++20 compilations agree byte-for-byte; a separately implemented Python mt19937_64 and built-in multiple-precision division reproduce all 2,240 source keys, quotients and remainders. 118 source TSV fixtures are registered across the same 32 generators. Consumer regression/CI results must be checked separately after publication; the opt-in full `verify-data` target is not part of standard CI.

### A4 — random 6,000-case independent exact-reference fixtures migrated (25 October 2026)


Preserve `mt19937_64(0xD1D1D1)`, every random draw, original 0/64/1..7-limb selection and ordered modulo-13/17/19 transforms, all divisor selection branches and both owning/non-owning division calls. Every row contains repetition, exact decimal divisor, canonical dividend, exact quotient and full 64-bit remainder. The consumer verifies input identity, ordered keys, both quotient paths, full-width remainder, scratch/output guard clearing, source immutability, and all-row exhaustion.

A second Python implementation of multiple-precision division and an independent mt19937_64 engine reproduces all 6,000 cases byte-for-byte; the first eight random draws also match a compiled std::mt19937_64 C++20 witness. Fixture 1,194,017 bytes, Git blob `5e9b7c6801c1fe9c6220a23235b50f592c8b1b84`. 119 TSVs across 32 existing generators. No production code or prior fixture changed. Standard CI does not run the opt-in all-fixture `verify-data` target.

### A5 — B48/B49 numerical oracle conversion (25 October 2026)


The consumer validates ordered width/divisor/mode, original source, canonical expected quotient/remainder, output guard clearing and TSV exhaustion.

The consolidated C++20 power/division generator independently computes all four new A5 fixture files (B48 560; B49 784; B57 1,800; B61 1,560 visits) using portable bit-at-a-time division with a 65th remainder carry and independent binary shift/add Euclidean identity checks. C++20 GCC and Clang mt19937_64 input witnesses and exact Python big-integer division independently reproduce all four files. SHA Git blobs: B48 `1a56559c3f5c529a5c222fe1136b69be1d3aed28`, B49 `4f82aac238df03238128ad5d41be21b9e2b24a22`, B57 `dc2a5205d8aa4a912f73f9e675522b1531f65794`, B61 `f9e8a720ece587c13d68802588cfe75ebb006e67`. Only B48/B49 consumer conversions are counted as migrated here. CI for A5 is pending. No production changes.

### A6 — complete runtime independent exact-reference fixtures of single-limb division (25 October 2026)


Use the existing independently certified `big_unsigned_b57.tsv` and `big_unsigned_b61.tsv` fixtures added in A5, retaining every original `mt19937_64` draw, all seven explicit transforms and the unchanged eighth random mode, original width and divisor sets and visited order.

The B57 consumer validates the full quotient and 64-bit remainder from `divrem_limb_into` against the exact TSV, retaining output guard zeroing, then checks `divrem_limb` quotient and remainder against the certified result. B61 validates `divrem` quotient and full remainder against the exact TSV, preserving the independent `rem` call and checking it matches that certified remainder. All expected values are checked through the canonical hex parser, including zero, and TSV rows must be complete and ordered. Original invalid-input, short-output and unchanged-buffer exception assertions are retained.

Total A1-A6 migrated test visits: 16,784. All existing 123 TSV fixtures retain their committed hashes; there are still 32 standalone generator executables. A5's independent Python and portable C++20 verification of B57/B61 fixture arithmetic is reused, not recharacterized as a new run. Full A6 consumer CI and opt-in `verify-data` remain separate checks; report CI as passed only when observed.

### A7 — exact single-limb division independent exact-reference fixtures (25 October 2026)


Preserve the full std::mt19937_64(0x22E1AC7) draw sequence, all original divisors, input-word transforms, both low-level `divexact_limb_into` and owning `quot_exact` result comparisons, canonical input/quotient shape, source immutability, output guards and all invalid/inexact/unchanged-buffer exception checks. The consumer rejects missing, reordered or extra records.

Independently replaying all 6,200 original std::mt19937_64 cases compiled with both GCC and Clang and Python multiple-precision multiplication reproduces every fixture byte-for-byte. Git blobs, respectively small (4,960 rows) `d4a68ff02f44a9273ca8bf0a30aa4874a57ab295`, length 65 (620 rows) `6a5c0cdee6d167bfedaff54917cc725bd73793ae`, length 129 first 78 divisors (312 rows) `0697660ce3945381db5840683f408ee73de55e1d` and remaining 77 divisors (308 rows) `c42a99815380e11bce4f3976773b8507483f47bb`. No production sources changed; the same 32 generator executables now cover 127 TSVs. A1–A7 cumulative migrated visits: 22,984. A7 consumer CI is to be checked after publication, separately from local fixture verification and opt-in `verify-data`.

### A8 — B57 power-of-two multiprecision division oracle (25 October 2026)


Convert all 2,160 original B57 visits in `TestBigUnsigned::test_division()`: nine original dividend widths x eight divisor widths x five bit positions x six input modes. Preserve the complete seeded std::mt19937_64(0xB570B172) input stream, including draws subsequently overwritten, all original ordered transforms, the power-of-two divisors, and all owning and caller-owned division paths. Each source-keyed TSV row records original width, divisor width, bit position and mode, canonical dividend and divisor, exact quotient and remainder. Consumer checks every key and both sources, `divrem`, `quot`, `rem`, `divrem_into`, zeroed quotient/remainder output guards, scratch guard integrity and absence of extra records.

The already registered independent standalone C++20 generator `generator_big_unsigned_power_div.cpp` produces three fixtures using Nat right shift and exact subtraction of the shifted quotient; it explicitly checks the Euclidean identity and remainder bound. GCC and Clang independently reproduce all 2,160 seeded input visits, while Python exact bit-shift and mask calculations produce the same source keys, dividends, divisors, quotients and remainders. Fixture blob SHAs and counts: short widths 1,2,8,16,32 (1,200 rows, 755,283 bytes) `35e1d41e2b862720d835a52fe477093aaeeb21a4`; medium 64,65 (480 rows, 826,491 bytes) `a4f0689a89079310a8bd85f066a9dfe69b687d45`; large 96,128 (480 rows, 1,301,002 bytes) `c5a0a0bf80b2de7a4cf1eecca40365822f0a1fe2`. The three committed Git blobs match Python-produced bytes exactly. The total is now 130 versioned TSVs across the unchanged 32 generator executables, 25,144 cumulative A1-A8 migrated visits. No production changes. Full consumer CI and opt-in generator-wide `verify-data` are separate validations.

### A9 — B55 multiprecision division independent exact-reference fixtures (25 October 2026)


Preserve the full `std::mt19937_64(0xB550DEC1)` draw stream, all ordered source transforms, the divisor high-bit normalization, and the 7/8/9 kernel dispatch boundary. Exact source-keyed fixtures check both operands, quotient and remainder on `divrem`, `quot`, `rem` and `divrem_into`; existing output-capacity, zero-surplus and scratch-guard assertions remain. All rows must appear in original order, with none missing or extra.

Separately regenerate source inputs with std::mt19937_64 compiled using GCC and Clang, and use Python multiple-precision divmod to certify every result. Exact bytewise Git blob identities: B55 short (216 visits, 178,716 bytes) `929963b14605a7a0eaa3f78aeacb8b398e615858`; medium (144 visits, 221,212 bytes) `3d8e5a3070ff80af9b01cd59a8b0a01a4da53882`; large (72 visits, 239,017 bytes) `7837342521a24a6cb672dbd364b37f0f4f44ab70`. Total: 133 versioned TSV fixtures and 32 existing generator executables; A1–A9 cumulative migrated test visits: 25,576. No production code changes. Full consumer CI and opt-in `verify-data` must be reported separately from oracle checks.

### A10 — B27 fused division oracle migration (25 October 2026)


Migrate both original fixed adversarial quotient-correction cases and all 2,048 original mt19937_64(0xB27071) randomized visits, preserving the original conditional RNG draws and ordered 13/19/23/29 transforms. The source-keyed consumer separately compares divrem, quot and rem against independent exact quotient and remainder, checks source immutability, scratch sizing and remainder bounds, and detects missing or extra records.

Reproduction through independently ported mt19937_64 and Python bigint division yields the committed blobs: first (1,026 rows, 847,088 bytes) `a4a2aa7f948a8bd9a673700f48bcd5efbb5ed101`; second (1,024 rows, 844,784 bytes) `18448ac0544c03bddb54c22488f507dd6778ffb2`. 135 fixtures, 32 standalone executables, 27,626 cumulatively migrated A1-A10 visits. Production code remains untouched. Consumer CI and optional verify-data are separate checks.


## Independent GMP-oracle retirement — archived A11–A34 record (25–26 October 2026)

The following inventory, milestone descriptions and validation outcomes were
preserved from the former standalone GMP oracle migration journal before its
removal. Its historical table is an audit of `tests/test_big_unsigned.cpp`,
not an assertion about other test suites; the original arithmetic semantics
and independence requirements are unchanged. Section heading depth was
adjusted to fit this design journal; the underlying descriptions and
measurements are preserved.


> **Current design perspective — October 2026.** This journal treats each native arithmetic component as an autonomous implementation. Progress is measured by exact, independently generated oracles, cross-platform regression tests, absolute operation latency, and speedups between revisions of the *same* algorithm. Historical cross-backend performance ratios and retired adapter details are omitted; mathematical proofs and native implementation decisions are retained.

Baseline: `main` at `60306926d7093bcf8ebee0867ec5ffe2600401c9` (A10), 25 October 2026.

### Scope and counting

Count **one test visit** per original input/mode combination whose expected
numeric output is supplied by GMP, regardless of the number of assertions.
Exclude GMP interoperability tests (`test_gmp_bridge()`), code/branch
coverage and migrations preceding A1. This inventory applies only to
`tests/test_big_unsigned.cpp`. The denominator remains **76,323** for
progress reporting, barring a documented correction to the audit.

### Current state — A34

| Method | Migrated | Remaining |
| --- | ---: | ---: |
| `test_division()` | 16,514 | 0 |
| `test_views_and_scratch()` | 11,600 | 0 |
| `test_multiply_kernel_paths()` | 11,432 | 0 |
| `test_randomized()` | 10,000 | 0 |
| `test_gcd_decimal()` | 3,793 | 0 |
| `test_division_limb()` | 12,944 | 0 |
| `test_exact_division_limb()` | 6,200 | 0 |
| `test_multiplication()` | 2,112 | 0 |
| `test_product_accumulation()` | 1,728 | 0 |
| **Total** | **76,323** | **0** |

Progress: **76,323 / 76,323 = 100.000%**.
A10 baseline was 27,626 migrated and 48,697 remaining (36.196%).

### Remaining inventory


### Commit milestones

| Phase | Newly migrated | Cumulative migrated |
| --- | ---: | ---: |
| A1–A10 | 27,626 | 27,626 |
| A11 — multiprecision B48 | 336 | 27,962 |
| A12 — randomized arithmetic, first segment | 1,250 | 29,212 |
| A13 — randomized arithmetic, second segment | 1,250 | 30,462 |
| A14 — randomized arithmetic, third segment | 1,250 | 31,712 |
| A15 — randomized arithmetic, fourth segment | 1,250 | 32,962 |
| A16 — randomized arithmetic, fifth segment | 1,250 | 34,212 |
| A17 — randomized arithmetic, sixth segment | 1,250 | 35,462 |
| A18 — randomized arithmetic, seventh segment | 1,250 | 36,712 |
| A19 — randomized arithmetic, final segment | 1,250 | 37,962 |
| A20 — randomized decimal parsing | 1,800 | 39,762 |
| A21 — structured decimal parsing | 484 | 40,246 |
| A22 — decimal formatting | 1,509 | 41,755 |
| A23 — scratch-buffer arithmetic | 5,000 | 46,755 |
| A24 — R38 shifted add/sub | 6,600 | 53,355 |
| A25 — B52 multiply-kernel | 288 | 53,643 |
| A26 — B53 sparse-row multiply | 2,304 | 55,947 |
| A27 — general multiply kernels | 6,000 | 61,947 |
| A28 — B8 + fixed Karatsuba | 2,840 | 64,787 |
| A29 — general division | 6,000 | 70,787 |
| A30 — B28 equal-high division | 2,048 | 72,835 |
| A31 — B36 high-ratio division | 1,280 | 74,115 |
| A32 — B29 sparse division | 1,024 | 75,139 |
| A33 — B45 zero-low/decimal division | 544 | 75,683 |
| A34 — divisor normalization | 640 | 76,323 |

A11 maintains the eight dividend widths, fourteen divisors, three original modes, all PRNG draws,
and separately validates exact `divrem`, `quot` and `rem` outputs and original scratch size
invariants against a new independent source-keyed fixture.
`tests/data/generator_big_unsigned_power_div.cpp` regenerates `big_unsigned_b48_multiprecision.tsv`
without GMP, Numeric or Boost. The fixture has 336 rows, 305,944 bytes, Git blob
`11042c7f930f66082d8a89d6670376c03786938d`. Other GMP division blocks remain unchanged.

### Post-commit reporting

After each commit report newly migrated test visits, cumulative
migrated/remaining numbers, percentage using the 76,323 baseline,
method-level remaining cases, and observed CI status separately from
fixture checks. Do not count multiple assertions on a single input as
multiple visits or infer `verify-data` success from ordinary CI.

### A12 implementation

A12 replaces GMP only for the first 1,250 seeded visits of
`test_randomized()`, using `big_unsigned_a12_randomized_first.tsv`.
The remaining 8,750 cases deliberately retain GMP. All 10,000 visits
continue to use the same single MT19937_64 engine and consume the same
shift draw each time; comparison, addition, conditional subtraction,
and both shifts have exact independent expectations. The C++20
standalone generator uses independent Nat arithmetic and the exact
fixture was cross-checked using a separately compiled GCC/Clang
random-input witness plus Python multiple-precision arithmetic.
The fixture is 1,202,560 bytes, Git blob
`487f06f3921d72d2de9caee06bf9c23957ee6012`.

### A13 — randomized arithmetic, second 1,250 visits (25 October 2026)

Migrate visits 1,250–2,499 of `test_randomized()` using the second
source-keyed fixture `tests/data/big_unsigned_a13_randomized_second.tsv`.
The entire original mt19937_64(0xBA17B00D) draw stream and original
n/m input selections, saturated/zero transformations, comparison,
addition, conditional subtraction, and left/right shifts are retained.
Both A12 and A13 are read sequentially from their own exact TSV; the
same C++20 random engine in the consumer serves all 10,000 visits.
The independent generator replays the first 1,250 full visits before
emitting the next 1,250, with GMP retained for visits 2,500–9,999.

GCC and Clang C++20 source witnesses generate identical 1,250 inputs.
Python built-in multiple-precision arithmetic separately generates
all expected values, confirmed byte-for-byte against the Git fixture
blob SHA `ae3a7c98ef3ac9dedea8e44d5a2d410daccb342e`
(1,217,436 bytes). Generator and consumer share no production arithmetic.
No production code or prior fixture was changed. Generator-wide
`verify-data` is not implied by normal CI.

### A14 — randomized visits 2,500–3,749 (25 October 2026)

Migrate 1,250 additional visits using an exact source-keyed TSV,
retaining the single original random engine, comparison/addition,
conditional subtraction and both shift operations. GMP still serves
visits 3,750–9,999. The standalone C++20 generator replays all 2,500
preceding visits and records results computed using independent Nat
arithmetic. The independently computed fixture has 1,250 records,
1,191,197 bytes and Git blob
`bdb0b5f9235005b2bc31d18a7e3b21c2af474abe`.

### A15 — randomized visits 3,750–4,999 (25 October 2026)

Migrate 1,250 additional visits with
`big_unsigned_a15_randomized_fourth.tsv`. The consumer still executes the
single original `mt19937_64(0xBA17B00D)` stream for all 10,000 visits and
retains comparison, addition, conditional subtraction and both shifts.
GMP remains the numerical oracle only for visits 5,000–9,999.

The standalone C++20 generator replays visits 0–3,749 before emitting A15.
A separate Python multiple-precision implementation reproduces the complete
A14 fixture exactly, including Git blob
`bdb0b5f9235005b2bc31d18a7e3b21c2af474abe`, then independently computes
A15. The A15 fixture has 1,250 records, 1,216,271 bytes and Git blob
`e23f97770356c28d0b76aef3c7d9d5c7dee8a0c5`. No production code or prior
fixture is changed. Generator-wide `verify-data` is not implied by normal CI.

### A16 — randomized visits 5,000–6,249 (25 October 2026)

Migrate 1,250 additional visits using
`big_unsigned_a16_randomized_fifth.tsv`. The consumer retains the original
single `mt19937_64(0xBA17B00D)` engine across all 10,000 visits, including
all operand draws, conditional full-array transformations and the shift draw.
Comparison, addition, conditional subtraction and both shifts continue to be
checked against exact external expectations. GMP remains the numerical oracle
for visits 6,250–9,999.

The standalone C++20 generator replays visits 0–4,999 before emitting A16.
A separate multiple-precision implementation reproduces A15 exactly, including
Git blob `e23f97770356c28d0b76aef3c7d9d5c7dee8a0c5`, before independently
computing A16. The A16 fixture contains 1,250 records, 1,204,292
bytes and Git blob `bbb0fc03a3dae88154205334fafdff5fb04f9619`. No production code or previous fixture
is changed. Generator-wide `verify-data` remains a separate check from CI.

### A17 — randomized visits 6,250–7,499 (25 October 2026)

Migrate 1,250 additional visits using
`big_unsigned_a17_randomized_sixth.tsv`. The consumer preserves the original
single `mt19937_64(0xBA17B00D)` engine across all 10,000 visits, including
the complete operand draw stream, conditional full-array transformations and
shift draw. Comparison, addition, conditional subtraction and both shifts are
checked against exact external expectations. GMP remains the numerical oracle
for visits 7,500–9,999.

The standalone C++20 generator replays visits 0–6,249 before emitting A17.
A separate multiple-precision implementation reproduces A16 exactly,
including Git blob `bbb0fc03a3dae88154205334fafdff5fb04f9619`, before
independently computing A17. The A17 fixture contains 1,250 records,
1,224,250 bytes and Git blob `2b595bab5b7a83fc2c9a27c188d2fba204b6c30e`. No production code or prior
fixture is changed. Generator-wide `verify-data` remains separate from CI.

### A18 — randomized visits 7,500–8,749 (25 October 2026)

Migrate 1,250 additional visits using
`big_unsigned_a18_randomized_seventh.tsv`. The consumer preserves the
original single `mt19937_64(0xBA17B00D)` engine across all 10,000 visits,
including the complete operand draw stream, conditional full-array
transformations and shift draw. Comparison, addition, conditional subtraction
and both shifts are checked against exact external expectations. GMP remains
the numerical oracle for visits 8,750–9,999.

The standalone C++20 generator replays visits 0–7,499 before emitting A18.
A separate multiple-precision implementation reproduces A17 exactly,
including Git blob `2b595bab5b7a83fc2c9a27c188d2fba204b6c30e`, before
independently computing A18. The A18 fixture contains 1,250 records,
1,212,953 bytes and Git blob `b1c7416dd60cc79511d545a706df56b18b70bf74`. No production code or prior
fixture is changed. Generator-wide `verify-data` remains separate from CI.

### A19 — randomized visits 8,750–9,999 (25 October 2026)

Migrate the final 1,250 `test_randomized()` visits using
`big_unsigned_a19_randomized_eighth.tsv`. The original single
`mt19937_64(0xBA17B00D)` engine, operand draws, conditional full-array
transformations and shift draws are unchanged. Comparison, addition,
conditional subtraction and both shifts now use exact external expectations
for all 10,000 visits; the runtime numerical GMP branch is removed from this
method.

The standalone C++20 generator replays visits 0–8,749 before emitting A19.
A separate multiple-precision implementation reproduces A18 exactly,
including Git blob `b1c7416dd60cc79511d545a706df56b18b70bf74`, before
independently computing A19. The A19 fixture contains 1,250 records,
1,202,168 bytes and Git blob `15f9d3656bf5fe77ac2b97e0ba6faee3b0883ea3`. No production code or prior
fixture is changed. Generator-wide `verify-data` remains separate from CI.

### A20 — randomized decimal parsing (25 October 2026)

Migrate all 1,800 randomized decimal parsing visits in
`test_gcd_decimal()` to the source-keyed
`big_unsigned_a20_decimal_random.tsv` fixture. The consumer retains the
single original `mt19937_64(0xB16B00B5)` engine after the complete B76
stream, so every input length and every decimal digit is unchanged.
Each row verifies repetition, length and the full source decimal before
checking the parsed value against an exact hexadecimal expectation and
its canonical decimal representation.

The standalone `generator_big_unsigned_decimal.cpp` replays the same
2,800 B76 visits only to advance the random engine, then computes A20
expectations with the independent `Nat` implementation in
`oracle_exact_integer.hpp`. It has no Numeric, GMP or Boost dependency.
The 1,800-row fixture is 1,011,827 bytes with Git blob
`5238d33947659f227a12efaaf8ffd63fa8770903`. An independently compiled
C++20 reproduction generated that exact blob. No production code or
previous fixture changes. The remaining numerical GMP oracle visits in
`test_gcd_decimal()` are the 484 structured parsing cases and 1,509
formatting cases.

### A21 — structured decimal parsing (25 October 2026)

Migrate the remaining 484 structured decimal parsing visits in
`test_gcd_decimal()` to `big_unsigned_a21_decimal_structured.tsv`.
The fixture preserves every original input string in execution order:
24 direct storage-boundary cases, 45 long B24 cases, 156 B47 cases,
234 B50 cases and 25 exact `2^(bits-1)+offset` cases. B47/B50 randomized
strings continue from the exact post-A20 state of the shared
`mt19937_64(0xB16B00B5)` stream.

The standalone decimal generator uses only the independent `Nat` helper
and emits 484 rows, 279,466 bytes, Git blob
`17a68c45126d418b00f9e984e95ea9bd486fef8e`. A separately compiled C++20
reproduction produced the same blob. Runtime GMP parsing and decimal
conversion are removed from these 484 cases. The only remaining numerical
GMP oracle visits in `test_gcd_decimal()` are 1,509 decimal-formatting
cases; production code and prior fixtures are unchanged.

### A22 — decimal formatting (25 October 2026)

Migrate the final 1,509 numerical GMP oracle visits in
`test_gcd_decimal()` to `big_unsigned_a22_decimal_format.tsv`.
Every row contains the canonical source magnitude in hexadecimal and an
independently generated decimal expectation. The consumer retains the
original random limb streams and dispatch-boundary constructions, checks
each native source against the frozen hex, then compares `to_decimal()`
with the TSV expectation and retains the native parse roundtrip assertion.

Coverage is B25 50, B25.1 18, B30 128, B42 208, B46 352, B51 288,
B59 216, B43 80, B31 119 and the final radix-power block 50, totaling
1,509 visits. The standalone decimal generator emits a 1,404,780-byte
fixture with Git blob `095b30979096b45222ebf019bbd5482588deb3dd`.
A separately compiled C++20 reproduction generated the same blob.
After A22, `test_gcd_decimal()` has no remaining numerical GMP oracle
visits; the fixed audit therefore moves this method to 3,793/3,793 migrated.

### A23 — scratch-buffer arithmetic (25 October 2026)

Migrate the first 5,000 visits of `test_views_and_scratch()` to five
1,000-row source-keyed TSVs. The original `mt19937_64(0xB22B22)` stream,
all generated widths, all-ones source transformations and the conditional
truncated-product length draw are unchanged. Each row independently fixes
the source magnitudes and exact add/sub/truncated-product expectations.

The standalone `generator_big_unsigned_views.cpp` uses independent `Nat`
arithmetic and shift-add multiplication with no Numeric, GMP or Boost
dependency. A separately compiled C++20 multiprecision witness reproduced
all five Git blobs exactly. No production code changes. The remaining
`test_views_and_scratch()` inventory is exactly the 6,600 R38 shifted
add/sub visits planned for A24.

### A24 — R38 shifted add/sub (25 October 2026)

Migrate the remaining 6,600 visits of `test_views_and_scratch()` to an
independent source TSV plus ten expectation TSV segments. The 264 original
base pairs preserve the exact `mt19937_64(0xB38038)` stream, generated
widths and all-ones transformations. Every one of the 25 shift pairs per
source is keyed by width, case index and both shifts, with independently
computed sum, comparison direction and absolute difference.

The standalone views generator uses only `Nat::shl`, add/sub and source
serialization; a separate C++20 multiprecision witness reproduced all 11
A24 blobs exactly. No production code changes. After A24,
`test_views_and_scratch()` is fully migrated at 11,600/11,600 visits, so
the fixed audit has only `test_division()` and
`test_multiply_kernel_paths()` remaining.

### A25 — B52 multiply-kernel oracle (25 October 2026)

Migrate all 288 B52 source visits in `test_multiply_kernel_paths()` to
`big_unsigned_a25_b52.tsv`. The exact seeded source stream, structural
transforms, canonicalized operands, full-product checks, fused-kernel
agreement, oversized destination guards and `full-1` low-product checks
are preserved. Only the runtime numerical GMP oracle is replaced.

The existing standalone multiplication generator uses independent shift/add
`Nat` arithmetic and emits source hex, exact full product and exact truncated
product. The fixture has 288 rows, 774,359 bytes and Git blob
`db88450f727a185a1b5543e1b2323081d6893088`. No production code or prior
fixture changes. Remaining `test_multiply_kernel_paths()` inventory is
2,304 B53 visits, 6,000 general-kernel visits, 2,800 B8 cutoff visits and
40 fixed Karatsuba cutoff visits, for 11,144 total.

### A26 — B53 sparse-row multiply oracle (25 October 2026)

Migrate all 2,304 B53 truncation visits to a 384-row source-keyed fixture,
with six expected low products per original source pair. Preserve all 16
width pairs, 24 modes, exact seeded draws, zero-row density transforms,
selected/fused/x86 carry agreement and original truncation lengths.

The standalone multiplication generator uses independent shift/add arithmetic;
`big_unsigned_a26_b53.tsv` is 2,112,618 bytes with Git blob
`a22d13dc7e098789870a45cd38625a0f00078ebd`. No production changes.
Remaining `test_multiply_kernel_paths()` inventory is 6,000 general-kernel
visits, 2,800 B8 cutoff visits and 40 fixed Karatsuba visits = 8,840.

### A27 — general multiply kernels (25 October 2026)

Migrate all 6,000 original general multiply-kernel visits to six sequential
1,000-row source-keyed TSVs. Preserve the single `mt19937_64(0xB3B3B3)`
stream, operand widths, every source transformation and the conditional
low-product-length draw. The test retains selected, portable, Comba, fused,
ARM64 carry, x86 carry, schoolbook and Karatsuba agreement checks; only the
full and truncated numerical oracle moves from GMP to the fixtures.

The standalone multiplication generator computes exact products using
independent shift/add `Nat` arithmetic. No production code changes.
Remaining `test_multiply_kernel_paths()` inventory is 2,800 B8 cutoff
visits and 40 fixed Karatsuba cutoff visits = 2,840.

### A28 — B8 and fixed Karatsuba oracle (25 October 2026)

Migrate the final 2,840 visits of `test_multiply_kernel_paths()` using two
source-keyed exact fixtures. B8 contributes 700 original sources x four
cutoffs; the fixed stress corpus contributes eight sources x five cutoffs.
The latter continues the exact B8 random engine state, preserving every
historical draw. All scratch/guard/output invariants remain in the consumer.

The standalone multiplication generator computes products with independent
shift/add arithmetic. The B8 fixture blob is
`16848b2efecd3cdc5535e03e0e4694df1890e09f`; the fixed fixture blob is
`51233e8a709ad097bfd1435c0a9171fdef7dcdb6`. No production changes.
After A28, `test_multiply_kernel_paths()` is 11,432/11,432 migrated and the
only remaining audited numerical GMP-oracle method is `test_division()`
with 11,536 visits.

### A29 — general division oracle (25 October 2026)

Migrate all 6,000 general randomized division visits to six sequential
source-keyed TSVs while preserving the exact `mt19937_64(0xB15D1A1D)`
stream, operand-width branches and all deterministic transforms. The consumer
retains divrem/quot/rem agreement, `divrem_into` output sizing, scratch guards,
remainder bounds and the native Euclidean identity assertion.

The new standalone remaining-division generator uses aligned-subtraction
long division and shift/add multiplication, with no Numeric/GMP/Boost
runtime dependency. No production code changes. Remaining audited numerical
GMP-oracle visits: B28 2,048, B36 1,280, B29 1,024, B45 544 and normalization
640 = 5,536 total.

### A30 — B28 equal-high division (25 October 2026)

Migrate all 2,048 B28 equal-high-limb visits to two exact source-keyed TSVs.
The original seeded stream, equality/borrow source transforms, joint and
single-output paths, zero-tail checks and the special guarantee that this
fast path leaves its scratch buffer untouched are all preserved. Only the
runtime GMP numerical oracle is replaced.

The standalone remaining-division generator emits two 1,024-row fixtures
using independent long division. No production code changes. Remaining
numerical GMP-oracle inventory is B36 1,280, B29 1,024, B45 544 and
normalization 640 = 3,488 visits.

### A31 — B36 high-ratio division (25 October 2026)

Migrate all 1,280 B36 visits to exact source-keyed TSVs while preserving all
widths, seeded draws and quotient-ratio classes. The consumer still executes
joint, quotient-only and remainder-only paths and verifies remainder bounds;
only the numerical expected values move off GMP.

No production changes. Remaining audited numerical GMP-oracle visits are
B29 1,024, B45 544 and normalization 640 = 2,208 total.

### A32 — B29 sparse division (25 October 2026)

Migrate all 1,024 B29 visits to exact source-keyed TSVs while preserving the
conditional random draw schedule and sparse/dense divisor shapes. The test
continues to exercise joint and single-output APIs and remainder bounds;
only the numerical GMP oracle is removed.

No production changes. Remaining audited numerical GMP-oracle visits are
B45 shifted 504, B45 decimal 40, and normalization 640 = 1,184 total.

### A33 — B45 zero-low/decimal division (25 October 2026)

Migrate all 544 B45 numerical GMP-oracle visits: 504 zero-low-limb shifted
divisor cases plus 40 decimal-power divisor cases. The original seeded
`mt19937_64(0xB450D1D1)` stream, source transforms, divrem/quot/rem checks,
`divrem_into` sizing, zero-tail assertions and scratch guards are preserved.
The source-keyed fixture records numerator, denominator, quotient and remainder.

The standalone remaining-division generator emits
`big_unsigned_a33_b45.tsv` without Numeric/GMP/Boost. No production code
changes. Remaining audited numerical GMP-oracle visits: the 640 divisor
normalization cases only.

### A34 — divisor normalization (25 October 2026)

Migrate the final 640 audited numerical GMP-oracle visits in
`test_division()`: all 64 divisor-normalization shifts crossed with the ten
original dividend modes. Inputs remain deterministically constructed in the
consumer; the source-keyed TSV verifies both operands before checking exact
quotient and remainder, `quot`/`rem` agreement, remainder bounds and the
native Euclidean identity.

The standalone remaining-division generator emits
`big_unsigned_a34_normalization.tsv` without Numeric/GMP/Boost. No production
code changes. This completes the fixed BigUnsigned audit at **76,323 / 76,323
= 100.000%** migrated numerical GMP-oracle visits. GMP interoperability tests
remain outside this denominator by definition.


### MSVC checked-span range fix (26 October 2026)

Windows Debug with `/RTC1` exposed a deterministic stall in the first
unequal-width subtraction (`2` limbs minus `1` limb). The arithmetic loop
itself was finite; the bug was the tail clearing range passed to
`std::fill`: its begin iterator came from `output.subspan(count)` while its
end iterator came from the original `output` span. That happens to behave like
a pointer range on some standard libraries, but it is not a valid checked
iterator pair on MSVC Debug.

All three occurrences in `BigUnsigned` (add, subtract and Karatsuba tail
clearing) now bind the subspan to a local `tail` and pass
`tail.begin(), tail.end()` to the algorithm. Regression coverage includes
explicit 2-limb/1-limb addition and subtraction. The temporary Windows tracing
used to isolate the problem was removed after the fix.


A follow-up Windows Debug run exposed five additional instances of the same
invalid iterator pairing in the single-limb and exact-division tail clearing:
`quotient.subspan(a.size()).begin()` was paired with `quotient.end()`.
Those sites are now fixed identically by naming the tail subspan and using
both iterators from that same view. Existing division tests already exercise
oversized quotient buffers, so no synthetic coverage reduction or replacement
was needed.
