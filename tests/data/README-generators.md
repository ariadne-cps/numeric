# Reproducible exact golden-data generators

The files `generator_big_unsigned_gcd.cpp`, `generator_big_unsigned_b76.cpp`,
`generator_big_unsigned_decimal.cpp`, `generator_big_unsigned_views.cpp`, `generator_exact_small.cpp`, `generator_big_unsigned_power_div.cpp`, `generator_big_unsigned_division_remaining.cpp`, `generator_rational_exact.cpp`,
`generator_exact_b72_r17.cpp`, `generator_integer_structural.cpp`, `generator_big_unsigned_mul.cpp`,
`generator_big_unsigned_dispatch.cpp`, `generator_integer_random.cpp`,
`generator_rational_binary64.cpp`, `generator_rational_shared_stream.cpp`,
`generator_rational_r31.cpp`, `generator_rational_random.cpp`,
`generator_big_unsigned_random_mul.cpp`, `generator_floatmp_directed_bounds.cpp`,
`generator_floatmp_pi.cpp`, `generator_floatmp_get_d.cpp`, `generator_floatmp_power_exact.cpp`,
`generator_floatmp_sqrt.cpp`, `generator_floatmp_pow_short.cpp`,
`generator_rounded_float_directed.cpp`, `generator_floatmp_exp.cpp`, `generator_floatmp_log.cpp`,
`generator_floatmp_sin_cos.cpp`, `generator_floatmp_centered_exp.cpp`, `generator_floatmp_atan.cpp`,
`generator_floatmp_tan.cpp`, `generator_floatmp_asin_acos.cpp`, `generator_floatmp_tanh.cpp`,
`generator_floatmp_rounding.cpp`, `generator_big_unsigned_sqrt.cpp`, and `oracle_exact_integer.hpp`
are **standalone C++20** sources. They use only the standard library, with an independent
unsigned-limb implementation for exact arithmetic. In particular, **do not link them to Numeric, GMP
or Boost**: expected values must not come from the code being tested.

### CMake targets

After configuring the project, invoke these targets explicitly:

```sh
cmake --build <build-dir> --target tests
cmake --build <build-dir> --target generators
cmake --build <build-dir> --target generate-data
cmake --build <build-dir> --target verify-data
```

`tests` builds the test suite without any generator dependency.
`generators` builds all standalone C++20 generator executables; the
generator executables are also excluded from the default `all` target.
`generate-data` first builds the generators if necessary, then runs
**every registered generator in sequence** with `--write` and an absolute
`--data-dir` pointing at the existing **source-tree** `tests/data`.
This target executes on each invocation, even if no input changed; it
overwrites only the currently covered TSV files (204 of 204 after the
27 October 2026 coverage iteration).
`verify-data` instead builds and runs every registered generator with
`--verify` and never rewrites fixtures. It fails on any missing or different
byte. The other TSVs remain committed and untouched. The generators' own
messages show which files were written. A failed generator stops the
sequence. Neither `tests` nor CTest runs the generators automatically.

The TSVs remain version-controlled: review `git diff -- tests/data`
after regeneration. Removing versioned fixtures or switching consumers
to build-tree paths is deliberately deferred until all generators are
implemented and the remaining GMP oracle migration is complete.

For a manual, read-only comparison rather than regeneration, use the
following standalone `--verify` commands.

From the repository root, using a C++20 compiler:

```sh
c++ -std=c++20 -O2 -Wall -Wextra -Werror tests/data/generator_big_unsigned_gcd.cpp -o /tmp/numeric-oracle-gcd
c++ -std=c++20 -O2 -Wall -Wextra -Werror tests/data/generator_big_unsigned_b76.cpp -o /tmp/numeric-oracle-b76
c++ -std=c++20 -O2 -Wall -Wextra -Werror tests/data/generator_big_unsigned_views.cpp -o /tmp/numeric-oracle-views
c++ -std=c++20 -O2 -Wall -Wextra -Werror tests/data/generator_big_unsigned_decimal.cpp -o /tmp/numeric-oracle-decimal
c++ -std=c++20 -O2 -Wall -Wextra -Wsign-conversion -Werror tests/data/generator_exact_small.cpp -o /tmp/numeric-oracle-exact-small
c++ -std=c++20 -O2 -Wall -Wextra -Wsign-conversion -Werror tests/data/generator_big_unsigned_power_div.cpp -o /tmp/numeric-oracle-power-div
c++ -std=c++20 -O2 -Wall -Wextra -Wsign-conversion -Werror tests/data/generator_big_unsigned_division_remaining.cpp -o /tmp/numeric-oracle-division-remaining
c++ -std=c++20 -O2 -Wall -Wextra -Wsign-conversion -Werror tests/data/generator_rational_exact.cpp -o /tmp/numeric-oracle-rational-exact
c++ -std=c++20 -O2 -Wall -Wextra -Wsign-conversion -Wshadow -Werror tests/data/generator_exact_b72_r17.cpp -o /tmp/numeric-oracle-b72-r17
c++ -std=c++20 -O2 -Wall -Wextra -Wsign-conversion -Wshadow -Werror tests/data/generator_integer_structural.cpp -o /tmp/numeric-oracle-structural
c++ -std=c++20 -O2 -Wall -Wextra -Wsign-conversion -Wshadow -Werror tests/data/generator_big_unsigned_mul.cpp -o /tmp/numeric-oracle-b35-b37-b38
c++ -std=c++20 -O2 -Wall -Wextra -Wsign-conversion -Wshadow -Werror tests/data/generator_big_unsigned_dispatch.cpp -o /tmp/numeric-oracle-dispatch
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_integer_random.cpp -o /tmp/numeric-oracle-integer-random
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_rational_binary64.cpp -o /tmp/numeric-oracle-binary64
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_rational_shared_stream.cpp -o /tmp/numeric-oracle-rational-stream
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_rational_r31.cpp -o /tmp/numeric-oracle-r31
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_rational_random.cpp -o /tmp/numeric-oracle-r0
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_big_unsigned_random_mul.cpp -o /tmp/numeric-oracle-random-mul
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatmp_directed_bounds.cpp -o /tmp/numeric-oracle-floatmp-bounds
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatmp_pi.cpp -o /tmp/numeric-oracle-pi
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatmp_get_d.cpp -o /tmp/numeric-oracle-get-d
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatmp_power_exact.cpp -o /tmp/numeric-oracle-power
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatmp_sqrt.cpp -o /tmp/numeric-oracle-floatmp-sqrt
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatmp_pow_short.cpp -o /tmp/numeric-oracle-short-pow
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_rounded_float_directed.cpp -o /tmp/numeric-oracle-rounded-directed
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatmp_exp.cpp -o /tmp/numeric-oracle-floatmp-exp
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatmp_log.cpp -o /tmp/numeric-oracle-floatmp-log
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatmp_sin_cos.cpp -o /tmp/numeric-oracle-floatmp-sin-cos
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatmp_centered_exp.cpp -o /tmp/numeric-oracle-floatmp-centered-exp
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatmp_atan.cpp -o /tmp/numeric-oracle-floatmp-atan
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatmp_tan.cpp -o /tmp/numeric-oracle-floatmp-tan
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatmp_asin_acos.cpp -o /tmp/numeric-oracle-floatmp-asin-acos
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatmp_tanh.cpp -o /tmp/numeric-oracle-floatmp-tanh
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatmp_rounding.cpp -o /tmp/numeric-oracle-floatmp-rounding
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_big_unsigned_sqrt.cpp -o /tmp/numeric-oracle-bigunsigned-sqrt
/tmp/numeric-oracle-gcd --verify --data-dir tests/data
/tmp/numeric-oracle-b76 --verify --data-dir tests/data
/tmp/numeric-oracle-exact-small --verify --data-dir tests/data
/tmp/numeric-oracle-power-div --verify --data-dir tests/data
/tmp/numeric-oracle-rational-exact --verify --data-dir tests/data
/tmp/numeric-oracle-b72-r17 --verify --data-dir tests/data
/tmp/numeric-oracle-structural --verify --data-dir tests/data
/tmp/numeric-oracle-b35-b37-b38 --verify --data-dir tests/data
/tmp/numeric-oracle-dispatch --verify --data-dir tests/data
/tmp/numeric-oracle-integer-random --verify --data-dir tests/data
/tmp/numeric-oracle-binary64 --verify --data-dir tests/data
/tmp/numeric-oracle-rational-stream --verify --data-dir tests/data
/tmp/numeric-oracle-r31 --verify --data-dir tests/data
/tmp/numeric-oracle-r0 --verify --data-dir tests/data
/tmp/numeric-oracle-random-mul --verify --data-dir tests/data
/tmp/numeric-oracle-floatmp-bounds --verify --data-dir tests/data
/tmp/numeric-oracle-pi --verify --data-dir tests/data
/tmp/numeric-oracle-get-d --verify --data-dir tests/data
/tmp/numeric-oracle-power --verify --data-dir tests/data
/tmp/numeric-oracle-floatmp-sqrt --verify --data-dir tests/data
/tmp/numeric-oracle-short-pow --verify --data-dir tests/data
/tmp/numeric-oracle-rounded-directed --verify --data-dir tests/data
/tmp/numeric-oracle-floatmp-exp --verify --data-dir tests/data
/tmp/numeric-oracle-floatmp-log --verify --data-dir tests/data
/tmp/numeric-oracle-floatmp-sin-cos --verify --data-dir tests/data
/tmp/numeric-oracle-floatmp-centered-exp --verify --data-dir tests/data
/tmp/numeric-oracle-floatmp-atan --verify --data-dir tests/data
/tmp/numeric-oracle-floatmp-tan --verify --data-dir tests/data
/tmp/numeric-oracle-floatmp-asin-acos --verify --data-dir tests/data
/tmp/numeric-oracle-floatmp-tanh --verify --data-dir tests/data
/tmp/numeric-oracle-floatmp-rounding --verify --data-dir tests/data
/tmp/numeric-oracle-bigunsigned-sqrt --verify --data-dir tests/data
```

`--verify` is the default and never writes to the data directory. It
regenerates every covered fixture from the test's original source-generation
rules and checks **byte-for-byte identity** against the committed TSV,
including comments, ordering, exact source operands, and expected results.
It returns nonzero on missing or different content. The optional `--write`
regenerates the complete covered fixture set and **overwrites** these TSVs;
review every diff before committing. No sampling is used.

### Coverage of versioned C++ generators

| Generator | Committed TSVs regenerated | Original test cases |
|---|---|---:|
| `generator_big_unsigned_gcd.cpp` | `big_unsigned_r6.tsv`, `r8`, `r10`, `r11`, `r21`, all nine `r26_*`, `r27` | 5,561 source-and-GCD checks, including 3,240 R26 variants |
| `generator_big_unsigned_b76.cpp` | all fourteen `big_unsigned_b76_gcd_decimal_*.tsv` | 2,800 source-keyed GCD/decimal cases |
| `generator_big_unsigned_decimal.cpp` | `big_unsigned_a20_decimal_random.tsv`, `big_unsigned_a21_decimal_structured.tsv`, `big_unsigned_a22_decimal_format.tsv` | 1,800 randomized parsing + 484 structured parsing + 1,509 formatting cases |
| `generator_big_unsigned_views.cpp` | five `big_unsigned_a23_views_*.tsv`, `big_unsigned_a24_shifted_sources.tsv`, ten `big_unsigned_a24_shifted_*.tsv` | 5,000 scratch-buffer arithmetic + 6,600 source-keyed shifted add/sub cases |
| `generator_exact_small.cpp` | `rational_r3.tsv`, `rational_r4.tsv`, `rational_r8.tsv`, `integer_power.tsv`, `big_unsigned_r14.tsv` | 6 + 54 + 16 + 100 + 28 = 204 exact cases |
| `generator_big_unsigned_power_div.cpp` | `big_unsigned_power.tsv`, `big_unsigned_short_div.tsv` | 373 exact power cases, 324 exact division cases; full input RNG replay |
| `generator_big_unsigned_division_remaining.cpp` | six `big_unsigned_a29_division_*.tsv`, two `big_unsigned_a30_b28_*.tsv`, two `big_unsigned_a31_b36_*.tsv`, two `big_unsigned_a32_b29_*.tsv` files | 6,000 general, 2,048 B28, 1,280 B36, and 1,024 B29 visits |
| `generator_rational_exact.cpp` | `rational_r28.tsv`, `rational_r30.tsv`, `rational_r22.tsv`, `rational_r16.tsv`, `rational_r19.tsv`, `rational_r20.tsv` | 36 signed add/sub, 1,689 normalization, 160 scalar mul/div, 84 R16 alias/add/sub, 102 R19 product/division, 160 R20 cross-add/sub |
| `generator_exact_b72_r17.cpp` | `big_unsigned_b72.tsv`, `rational_r17.tsv` | 288 B72 products/squares and 252 R17 comparisons |
| `generator_integer_structural.cpp` | `integer_b66.tsv`, `integer_sub.tsv`, `big_unsigned_sub.tsv`, `integer_b33.tsv`, `integer_b37.tsv`, `integer_b38.tsv` | 360 power, 192 signed subtract, 100 unsigned subtract, 72 signed squares, 72 signed powers, 30 signed cubes |
| `generator_big_unsigned_mul.cpp` | `big_unsigned_b35.tsv`, `big_unsigned_b37.tsv`, `big_unsigned_b38.tsv`, `big_unsigned_a25_b52.tsv`, `big_unsigned_a26_b53.tsv`, six `big_unsigned_a27_general_*.tsv` files, `big_unsigned_a28_b8.tsv`, `big_unsigned_a28_karatsuba_fixed.tsv` | 70 Karatsuba products, 91 squares, 112 products with six truncated variants, 288 B52, 2,304 B53, 6,000 general, 2,800 B8 and 40 fixed Karatsuba visits |
| `generator_big_unsigned_dispatch.cpp` | `big_unsigned_b39.tsv`, `big_unsigned_b40.tsv`, `big_unsigned_b54_product.tsv`, `big_unsigned_b54_square.tsv` | 152 B39, 162 B40, 136 B54 products, 104 B54 squares |
| `generator_integer_random.cpp` | `integer_b17.tsv`, `integer_b18.tsv` | 250 signed arithmetic/shift cases and 170 signed arithmetic/division cases |
| `generator_rational_binary64.cpp` | `rational_r7.tsv`, `rational_r12.tsv` | 384 R7 rational probes; 2,048 R12 one-limb seeded conversions (both signs) |
| `generator_rational_shared_stream.cpp` | `rational_r29.tsv`, all sixteen `rational_r13_mode*_limbs*.tsv`, `rational_r5.tsv` | 216 R29, 1,536 R13, 48 R5 rows; exact shared original RNG state |
| `generator_rational_r31.cpp` | 14 `rational_r31_<bits>_<chunk>.tsv` fixtures | All 960 R31/R32 cases, 10 widths, 24 source indices, four signs; preserve 2048/4096-bit chunks |
| `generator_rational_random.cpp` | `rational_random.tsv` | Full 1,200-visit RNG replay and all 120 original stratified exact R0 cases |
| `generator_big_unsigned_random_mul.cpp` | `big_unsigned_random_mul.tsv` | All 6,000 original seeded multiplication visits, 600 selected source-keyed full, square and low products |
| `generator_floatmp_directed_bounds.cpp` | `floatmp_directed_bounds.tsv` | 90 certified exact rational lower/upper dyadic bounds, 6 precisions × 15 signed inputs |
| `generator_floatmp_pi.cpp` | `floatmp_pi.tsv` | Machin's formula: 400-term rigorous 1,200-bit integer interval, directed 384-bit and nearest-even 320-bit π |
| `generator_floatmp_get_d.cpp` | `floatmp_get_d.tsv` | 272 exact p-bit FloatMP-to-binary64 conversions, four precisions and four IEEE rounding modes |
| `generator_floatmp_power_exact.cpp` | `floatmp_pow.tsv` | 2,400 original signed integer-power references at five precisions, eight degrees and four rounding modes |
| `generator_floatmp_sqrt.cpp` | `floatmp_sqrt.tsv` | 300 unique input/mode exact rounded square roots; 360 source checks retained by consumer |
| `generator_floatmp_pow_short.cpp` | `floatmp_pow_short.tsv` | 3,680 exact boundary-focused integer powers, 5 precisions, 7 degrees, 4 directed modes |
| `generator_rounded_float_directed.cpp` | `rounded_float_directed.tsv` | 24 rigorous 128-bit exp/log/atan directed pairs, 48 original strict bound comparisons |
| `generator_floatmp_exp.cpp` | `floatmp_exp.tsv` | 320 certified, correctly rounded exp results across 20 signed inputs, 4 precisions and 4 rounding modes |
| `generator_floatmp_log.cpp` | `floatmp_log.tsv`, `floatmp_log_native.tsv` | 452 certified logs plus 32 native exponent limit results at 53/128/256/1024 bits in all four modes |
| `generator_floatmp_sin_cos.cpp` | `floatmp_sin.tsv`, `floatmp_cos.tsv` | 1,080 exact source-keyed independently certified outputs each, across D91, D92, D95 original test families |
| `generator_floatmp_centered_exp.cpp` | `floatmp_centered_exp.tsv` | 568 exact certified reference rows from all 720 original centered-exp visits |
| `generator_floatmp_atan.cpp` | `floatmp_atan.tsv` | 1,340 certified atan references reproducing all original D87, D88, D89 lookup keys |
| `generator_floatmp_tan.cpp` | `floatmp_tan.tsv` | 1,208 certified tangent values across D92, D93, D94 and D95 original sources |
| `generator_floatmp_asin_acos.cpp` | `floatmp_asin.tsv`, `floatmp_acos.tsv` | 728 certified inverse-trigonometric references per TSV, preserving D86/D90 original sources |
| `generator_floatmp_tanh.cpp` | `floatmp_tanh.tsv` | 1,164 certified tanh references from 13 precisions, 4 modes, signed/cutoff/saturation probes and binary64 discarded-word boundaries |
| `generator_floatmp_rounding.cpp` | `floatmp_rounding.tsv`, `floatmp_add_sub.tsv` | 4,794 exact-rational rounding rows and 1,780 independent D80 add/sub rows; both existing TSVs unchanged |
| `generator_big_unsigned_sqrt.cpp` | `big_unsigned_sqrt.tsv` | 166 exact sqrtrem triples; fully specified seeded SplitMix64 input generation, 37 random width profiles, 33 square-neighbor cases and 5 buffered-division/correction boundaries |

**Current reproducibility scope (27 October 2026):** all 203 versioned TSVs
are covered by the 42 registered standalone generators, including the 16
new fixtures described below. The complete set is checked byte-for-byte by
`verify-data`; ordinary CTest still consumes committed fixtures without
running generators. Historical counts and implementation notes below record
the migration sequence rather than the current inventory. Original frozen
header text, including historical FloatAP names and provenance descriptions,
is preserved as part of the byte-level reproducibility contract.
Transcendental oracles continue to require certified outward bounds that
round to the same exact result; approximate host floating-point outputs
are not accepted.

The exact-limb helper is a small isolated oracle implementation, **not** a
new Numeric production arithmetic backend. It uses a binary GCD, radix-2^64
limbs, and long division by 10^9 to format decimals; its algorithms differ
from Numeric's production GCD and formatting implementations. Generator
code is compiled independently and does not become a Numeric library target.

Initial adoption was checked with GCC and Clang C++20 warnings-as-errors.
The generated 29 files matched the Git blob SHA-1 object IDs of all 29
previously committed TSVs, **without modifying the reference data**.
This is reproducibility evidence, not a claim of cross-platform compiler
validation or formal proof of all independent-oracle algorithms.

B78 independently regenerates five additional existing exact fixture families,
for 34 versioned TSV files in total. R4 preserves canonical `1/1` rather
than the unreduced `2^bits / 2^bits` representation; R14 consumes every
original seeded source draw even though only four structurally distinct
cases per limb width are retained. GCC 14 and Clang 17 builds with strict
warnings (including `-Wsign-conversion -Werror`) reproduced the five
unchanged committed fixture blob IDs byte for byte. An intentionally
altered TSV was correctly rejected by `--verify` without rewriting it.
This is not yet a full-platform CI result.

B79 reproduces the D104 BigUnsigned power and D103 short-owned division
fixtures without changing them (36 TSV files now have versioned generators).
The power fixture replays all source draws from the original 10 limb widths
and 12 source modes, retaining the same five selected modes and six
exponents; the deterministic sample and boundary power cases remain intact.
The division fixture replays all 864 source draws and retains the original
nine distinct tested mode indices per width pair, including zero divisors.
Both files exactly matched their previously committed Git blob IDs with
GCC 14 and Clang 17, including the strict Clang/CI warning set with
`-Wsign-conversion -Wshadow -Werror`. `--verify` rejects a modified fixture
without altering it. The independent long-division generator additionally
checks `a == b*q+r` and `r < b` for every nonzero divisor. CI remains pending.

B80 adds independent signed-rational and fraction-normalization sources for R28, R30, and R22 (39
reproducible TSVs overall). Signed source and exact result fractions are derived without linking or
executing Numeric, GMP or Boost. The R28 source columns intentionally preserve the original
unreduced numerator/denominator in cases with a common factor, while mathematical result columns are
reduced canonically. R30 includes all nine dense fallback rows at original 128..8192-bit boundaries.
The three files were reproduced byte-for-byte against their pre-existing Git blob SHA-1 hashes on
GCC and Clang with strict CI-style warnings, including `-Wsign-conversion -Wshadow -Werror`; a
deliberately modified fixture was rejected in `--verify` mode. No checked-in TSV or consumer
changed. Keep TSVs versioned until every generator and GMP oracle migration is complete; no
automatic `generate-data` target is introduced.

B81 extends the same independent rational fixture generator to all R16,
R19 and R20 original cases. Each source pair is explicitly reconstructed
from the historical pre-TSV test source, not read back from a TSV; exact
signed rational outputs are canonicalized through independent C++20 limb
arithmetic. All 346 additional rows were also independently reproduced
using exact JavaScript BigInt arithmetic and their complete Git blob hashes
matched the previous TSVs, including headers and row ordering. For R20,
the fourth pair has identical operands, essential to preserve cancellation
and aliasing; it must not be replaced with the first unequal pair.
No TSVs are removed, and `generate-data` is still deferred.

B82 independently reproduces B72 and R17. The B72 source consumes all
1,152 original mt19937_64 source visits (including 24 nonselected modes
per width pair) before emitting exactly 36 x 8 = 288 source-keyed
product/square rows. R17 regenerates all 7 x 6 x 6 = 252 exact signed
cross-comparisons using separate shift/add limb multiplication. GCC 14
and Clang 17 strict warning builds generated both TSVs with their
original Git blob IDs, byte-for-byte. No fixtures are changed, no GMP
oracles are replaced, and automatic data generation is still deferred.

B84 adds six C++20 exact structural fixtures while keeping every original
TSV versioned. The historical source streams are replayed across all original
unselected RNG visits, not just retained samples. The signed and unsigned
expected outputs are calculated by an independent radix-2^64 implementation,
using binary shift/add powers and exact canonical sign handling. Both GCC
and Clang strict C++20 builds regenerated all six TSVs byte-for-byte; their
Git blob SHA-1 IDs match the original committed fixtures. `generators` now
compiles seven executables, while `generate-data` runs all seven in sequence
and overwrites only the 50 reproducible TSVs. `tests` remains decoupled.

B85 adds exact independent B35/B37/B38 multiplication, squaring and truncated
multiplication fixtures, while preserving the complete seeded RNG advance
through discarded modes. The three generated TSVs match the original Git
blob SHA-1 values, verified after both GCC and Clang strict-warning builds.
The eighth program is wired through `generators` and `generate-data`; the
source-tree TSVs remain versioned and `tests` does not depend on generators.

B86 adds the four independent B39/B40/B54 multiplication and square fixtures.
All original source-stream draws are replayed, including discarded modes.
B54's square fixture intentionally shares the exact RNG state after all
17 x 32 product iterations. With strict GCC and Clang C++20 warnings,
all four regenerated TSVs match their original Git blob IDs byte-for-byte.
The ninth generator participates in explicit `generators` and
`generate-data` targets, not `tests` or the default build.

B87 restores the full original std::mt19937_64 visits for B17 (2,500)
and B18 (1,700), freezes their original 1-in-10 stratified sample
positions, and derives exact signed add/subtract/multiply, truncating
quotient/remainder and B17 truncating right shift using independent
base-2^64 integer arithmetic. GCC and Clang C++20 strict-warning
compilation and execution matched both committed Git blob SHA-1 IDs.
The tenth program is available from explicit `generators` and
`generate-data`; existing TSVs remain unchanged and tracked.

B88 reproduces R7 and R12 conversion references with exact integer mantissa/exponent extraction and
round-toward-zero behavior, independent of the host FP rounding mode, runtime GMP. The original R12
mt19937_64 source sequence is fully replayed for each of four host modes. Both generated TSVs match
the original checked-in blobs byte for byte, with GCC and Clang strict C++20 warnings; the read-only
--verify also rejects deliberate TSV corruption. The eleventh generator is registered in the
explicit generators and generate-data targets.

B89 regenerates all 18 R29/R13/R5 files from the exact shared Rational mt19937_64 stream. Replaying
the 1,200 original R0 source visits is mandatory, followed by interleaving R29, four R13 limb-width
groups and R5 within EACH host-mode iteration; generating R29 for all modes first would break the
reproducibility of downstream inputs. Exact integer arithmetic yields GMP-style binary64 truncation,
without GMP. GCC and Clang strict compilation and --verify reproduce all 18 original Git blobs
byte-for-byte. Registered as the twelfth explicit generator; all original files remain versioned and
unchanged.

B90 adds the complete deterministic R31/R32 corpus: 960 exact rational
source combinations in 14 original TSV chunks. Chunk indices and exact
header ranges are preserved at widths 2048 and 4096 bits. Both GCC and
Clang compile and run the independent C++20 generator under strict warnings,
and every output matches the original Git blob SHA-1, with no copied data.
The generator is registered in explicit `generators` and `generate-data`;
`tests` and the default build remain decoupled. Data stays versioned.

B91 adds the complete independent R0 random Rational fixture generator.
It replays all 1,200 original seeded mt19937_64 iterations and freezes
exactly the same 120 adaptive one-per-ten sample offsets, without reading
TSV input. Independent signed radix-2^64 limbs reproduce raw integer
sources, canonical fractions, arithmetic, floor/ceil and GMP-compatible
round-toward-zero binary64 bit patterns. GCC and Clang strict builds match
the existing `rational_random.tsv` Git SHA exactly; read-only --verify
rejects corruption. It is the fourteenth explicitly registered generator.

B92 reproduces the full original-source BigUnsigned multiplication corpus:
replay the initial 7 x 12 R14 seeded visits, then all 6,000 random
multiplication source iterations with their four structural transforms and
conditional fifth-class low-product RNG draws. Freeze exactly the 600
original adaptive one-per-ten selected indices; they are selection metadata,
not copied expected data. Independent shift/add limb multiplication and
truncation yield the original 942,214-byte TSV and its exact Git SHA-1.
GCC and Clang strict C++20 builds succeed. Explicit `generators` and
`generate-data` include the fifteenth program; test builds remain separate.

B93 reproduces `floatmp_directed_bounds.tsv` from original rational inputs, computing exact
floor/ceil on a power-of-two lattice and canonical dyadic fractions at precisions 7, 32, 53, 64, 128
and 192. Both GCC 14 and Clang 17 compile and regenerate all 90 rows with strict C++20 warnings,
producing the original 7,742-byte Git blob. The independent generator uses no host floating-point
calculations or runtime GMP/FloatMP. Its addition to the explicit `generators` and `generate-data`
targets does not change tests, existing fixtures or the default build.

B94 independently recreates `floatmp_pi.tsv` using exact integer fixed-point Taylor intervals for
Machin's formula, 400 terms per atan at scale 2^1200. Explicit integer interval widths and next-term
remainders certify directed 384-bit endpoints and the same nearest-even 320-bit value from BOTH ends
of the interval. GCC and Clang C++20 strict builds reproduce the original 1,250-byte Git blob
without GMP, Numeric or host floating point. The seventeenth generator is available from
`generators` and `generate-data`; all data remains versioned and `tests` does not depend on
generators.

B95 reproduces `floatmp_get_d.tsv`: independently round each of the 17
original exact rational sources to an intermediate p-bit FloatMP number
(precision 53, 128, 256 or 1024) with nearest-even ties, then convert
the exact dyadic result to IEEE-754 binary64 in each of four directed
modes. Signed zero, normal/subnormal, carry, midpoint-even and saturating
overflow use pure integer quotient/remainder. A separate exact Python
integer implementation reproduced all 272 outputs before the C++20
GCC/Clang strict build and 65,751-byte original Git blob comparison.
This eighteenth generator belongs to `generators`/`generate-data`, not
`tests`, and neither GMP nor original fixture values are consulted.

B96 reproduces all 2,400 exact FloatMP integer-power fixture keys from 15 distinct original signed
input rationals, five binary precisions, eight exponent choices and all four output rounding modes.
The input is rounded to p-bit nearest/even using independent integer quotient and remainder, its
exact dyadic significand raised to the original exponent, and the final result rounded ONCE to
precision p using nearest/even or signed directed rules. All outputs are canonical dyadic rationals.
GCC 14 and Clang 17 strict C++20 builds and independent execution match the original 106,895-byte
Git blob; neither Numeric nor GMP/Boost are used. The nineteenth generator participates in explicit
`generators` and `generate-data`, not `tests`; all TSVs remain versioned.

B97 recreates FloatMP square root expectations from the ORIGINAL 14 positive
rational sources and four shifted nearest-rounded inputs at each of five
precisions, preserving source order and deduplicating identical rounded keys.
Exact binary integer division, pairwise restoring isqrt, signed scale
bookkeeping, and strict midpoint-even comparisons generate 300 reference
rows matching the original 76,870-byte `floatmp_sqrt.tsv` Git blob.
GCC and Clang C++20 strict builds and cross-compiler read-only --verify
pass. The twentieth generator is included only in explicit `generators`
and `generate-data`, while all original fixtures remain versioned.

B98 regenerates 3,680 source-keyed short integer powers from exact
low/high binade-edge dyadics, preserving five precision sets (2,7,53,
64,65) and seven degree-dependent binade distances. Near each lower
and upper binade boundary the offset is derived by exact integer
quotients of powers of two, including the distinct lower-vs-upper
threshold at even precision p=64. All signed distinct source inputs
are ordered by numeric rational value and degree, nearest-even rounded
at p bits before exponentiation, and finally independently rounded under
N/D/U/Z. GCC 14 and Clang 17 strict builds reproduce all 462,642 bytes
of the original `floatmp_pow_short.tsv`; no frozen operand values or
expected outputs are copied into source. Generator #21 participates only
in the explicit `generators`/`generate-data` targets.

B99 independently regenerates the exact 24-source 128-bit directed
`rounded_float_directed.tsv` certificate with standard C++20 integers.
Its exp, log and atan interval enclosures use an outward-rounded
fixed-point scale of 2^1400; Taylor remainders, Machin's pi/4 reduction,
negative-argument signs and reciprocal bounds are explicitly enclosed.
Both enclosures must lead to the same 128-bit downward and upward binary
values. The source contains no expected TSV values and reads no TSV
inputs. GCC and Clang C++20 warnings-as-errors builds regenerated all
5,181 bytes with identical original Git blob SHA-1
`f20fd0de94317e1999cfe1c4d7ac6718c443a9db`; read-only --verify
also rejected a deliberately modified file. The 22nd generator is
registered in `generators`, `generate-data` and the new non-writing
`verify-data` CMake target. Generator coverage is 102 of 115 TSVs.
No existing fixture, consumer, runtime oracle or production source changed.

B100 independently regenerates `floatmp_add_sub.tsv`: the complete
10 x 10 ordered source pair corpus, both exact additions/subtractions,
all four directed/nearest output modes, and precision sets 7, 53, 128,
256 and 1024 bits. The original 7/11 source is first nearest-even
quantized to the selected precision, and the other source dyadics are
constructed from their original powers of two. Every exact signed sum
and difference is formed in arbitrary-width integer limbs before a
single p-bit output rounding, with canonical rational serialization.
Keys are deduplicated only after reproducing their source order, leaving
all 1,780 distinct rows while existing test consumers retain all
10 x 10 x 2 x 4 x 5 original checks. GCC and Clang strict C++20
compilation, read-only --verify, and deliberate fixture corruption
rejection pass. The reproduced 1,104,075-byte source-tree TSV has Git
blob SHA-1 `b4552cfbdea659730486dd5bcdb3bb3a62be7984`, exactly the
original committed blob. This adds generator #23 to the explicit
`generators`, `generate-data` and `verify-data` targets, without
modifying a TSV, test consumer, production source or runtime GMP oracle.
Coverage is now 103 of 115 TSVs; all-platform CI is pending.

B101 independently regenerates `floatmp_exp.tsv` with correctly rounded, rigorously certified
exponential values. Recreate nine exact rational magnitudes (1/4, 3/4, 1/2, 1, 3, 7, 2^-16, 2^-64
and 1023/1024) and both signs, in the original order, then emit four bit widths (53/128/256/1024)
and N/D/U/Z modes for all 288 original checks. At B=2^1400, compute outward integer interval Taylor
terms k=0..165 for exp(|x|/16) and bound the positive remainder by twice term 166. Four directed
outward squares enclose exp(|x|), and exp(-|x|) is the outward reciprocal of this positive interval.
Round BOTH interval endpoints independently in each mode and require matching exact dyadic rational
outputs before accepting any row. C++20 standard-library-only with independently implemented 32-bit
digit multiplication and exact limb division; no GMP, FloatMP, Numeric, Boost or host floating exp.
Strict GCC and Clang warnings-as-errors builds passed; `--verify` passed on both and rejected
deliberately corrupted source data. Regenerated exactly 288 rows / 68,122 bytes, identical to Git
blob SHA-1 `92ce3c0eeda7f73f5e0385abf483f93367792891` on main. The 24th independent generator
participates only in opt-in CMake targets `generators`, `generate-data`, `verify-data`. No original
TSV, test consumer, production source or runtime GMP oracle was changed. Generator coverage reaches
104/115. Full CI pending for this commit.

B102 independently reproduces all 444 rows (122,768 bytes) in
`floatmp_log.tsv`. Recreate the original six centered-precision visits,
all 15 centered rational source operands, nearest-even input rounding,
unique-key de-duplication, and the original D96 input sequence at four
precisions. Use exact integer limb arithmetic at B=2^1800 and outward
atanh series through k=1100 with a bounded positive remainder; for
1/8 and 8 use the exact identity log(8)=3*log(2) to avoid an
insufficient direct Taylor remainder at precision 1024. Both enclosure
endpoints must yield identical N/D/U/Z p-bit rounded dyadic rationals.
D96's comment is emitted at its original in-file position, not in
the first header block. GCC/Clang strict C++20 compiles and both
--verify runs pass; deliberately corrupted fixtures fail verification.
The regenerated blob SHA-1 `7cae466575432ff238dd73ecdeeb599b7e17ceed`
is identical to main's original fixture. Generator #25 is opt-in through
`generators`, `generate-data`, `verify-data`. Coverage now 105/115.
No fixture, consumer, runtime GMP oracle or production source changed.

B103 independently regenerates both `floatmp_sin.tsv` and `floatmp_cos.tsv` without changing their
original contents. The standalone C++20 generator reconstructs precisely the three input families
(D91 centered offsets p/2+16, D92 local offsets p, D95 wider exact binary dyadics), preserving the
full precision-visit order, signed source lookup keys, original nearest-even FloatMP input
quantization and first-seen key deduplication. There are exactly 336 D91, 456 D92 and 288 D95 unique
source/mode rows for each function (1,080 each). Evaluate 240 alternating integer Taylor terms for
sin and cos with outward integer interval products/quotients at B=2^1800; separately enclose the
first omitted tail term. Negative sine uses exact oddness and cosine uses exact evenness. All four
N/D/U/Z roundings require both enclosures to agree in the target precision; a shared cache reuses
exact rounded input magnitudes. The generator does not read TSV files to obtain either input
operands or expected answers, and uses no Numeric, GMP, Boost or host floating-point sine/cosine.
GCC and Clang strict C++20 builds and --verify pass on both original files; deliberate corruption
fails. Original Git blob SHAs are unchanged: `floatmp_sin.tsv`:
`9b75d8fd4992f62e288c151e8f930dedd66a778e`, 281,134 bytes; `floatmp_cos.tsv`:
`31475db715e25928b31b34d31e32735e2cbb97be`, 278,016 bytes. Generator #26 participates only in the
opt-in targets `generators`, `generate-data`, `verify-data` and increases coverage to 107/115,
leaving eight original TSV files without generators. No tests, production source or runtime GMP
oracle were changed.

B104 adds the standard-library-only `generator_floatmp_centered_exp.cpp`
covering all 568 distinct rows of `floatmp_centered_exp.tsv`. Recreate
all six precision visits, all 15 rational magnitudes and both signs,
quantizing the source independently to nearest-even p-bit binary before
forming each original lookup key. Deduplicate only repeated keys,
retaining all 720 original consumer checks. Certify exp(|x|) by
165 Taylor terms at B=2^1800, twice-next-term remainder, four
outward interval squarings and for negative arguments an outward
reciprocal. Require identical N/D/U/Z p-bit rounding of both interval
endpoints. GCC and Clang C++20 warnings-as-errors builds and
read-only --verify pass, deliberate fixture corruption is detected.
Generated 170,882 bytes identical to the original blob SHA-1
`a221dfa3cc3e8feff2b425aff8290180fe994d3a`. This is generator
#27 and coverage becomes 108/115; CMake targets are opt-in.
No original TSV, test consumer, runtime GMP oracle or production code
has been modified. Cross-platform CI for B104 is pending.

B105 independently regenerates all 1,340 original `floatmp_atan.tsv`
source-keyed references (309,281 bytes). Preserve D87 centered-source
precision visits, D88 Dyadic regression probes and D89 cache/centered
argument-reduction cases in exact visit order; independently round
original inputs nearest-even at each precision. Use portable unsigned
integer limbs at B=2^1800, an outward alternating Taylor recurrence,
and certified Machin pi/4 identity with interval argument reduction.
Small-divisor long division and an early termination at a rigorously
bounded tiny omitted tail permit full generation efficiently without
weakening the 1,150-term maximum. Independently round both interval
endpoints under N/D/U/Z and reject ambiguous values. GCC/Clang strict
C++20 builds and --verify pass; corrupted fixtures are rejected.
The unchanged original TSV Git blob SHA is
`9daa471b286fb3b3b5cd0d4a804a4977cfd16079`.
Generator #28 joins only opt-in CMake targets, reaching 109/115.
No input fixture, test consumer, runtime GMP oracle or production code
has changed.

B106 independently regenerates all 1,208 original `floatmp_tan.tsv`
references: 576 D92 local sin/cos/tan, 240 D93 centered tangent,
200 D94 one-square-root cosine tangent and 192 D95 exact dyadic
trigonometric recurrence/regression references. Reconstruct every
original source value algebraically in original visitation order,
including all signed cases, precision repeats, and first-seen-key
deduplication. Convert each source independently to nearest-even
p-bit binary before evaluating its exact tangent.

The generator uses standard C++20 unsigned integer limbs and outward B=2^1800 fixed-point
alternating sine and cosine Taylor series (240 terms and omitted-term enclosure). It encloses the
signed quotient sin(x)/cos(x) with rigorous outward rational interval division and correctly handles
negative cosine for arguments past pi/2. Oddness handles negative arguments; all four rounding modes
N/D/U/Z require both interval endpoints to agree as exact binary rational outputs. No Numeric, GMP,
Boost or host floating-point transcendental is used to calculate expectations. GCC and Clang strict
C++20 warnings-as-errors builds and --verify pass; deliberate fixture corruption fails. Original TSV
byte count 303,143 and blob `0b73c09cd203961ace02c86c3af8c866db8c7ac4` are reproduced exactly.
Generator #29 joins opt-in CMake targets, raising coverage to 110 of 115 TSVs. No fixture, test
consumer, production code, or runtime GMP oracle has been modified.

B107 adds an independent C++20 generator for both asin and acos
fixtures. It reconstructs the original D86 centered source visits
and D90 signed endpoints, with exact source keys and independently
nearest-even rounded binary inputs. Rigorous fixed-point B=2^1800
integer enclosures calculate sqrt(1-x^2) through digit-pair integer
square root of a rational-scaled integer, then asin(x)=2*atan(x/(1+
sqrt(1-x^2))), and acos(x)=pi/2-asin(x) using a certified Machin pi/4.
All interval endpoints must round to identical p-bit N/D/U/Z values.
728 rows per file: 520 D86 and 208 new D90 references. GCC and
Clang strict C++20 compiles, --verify for both TSVs, and deliberate
fixture corruption rejection pass. The independently regenerated
189,909-byte asin blob is `50d8070f7cc049ad76840f51f2d4042ce9f73cb4`;
189,361-byte acos blob is `9bc0e3b101c737783565f4dc67c0033ec6fa7c91`.
Both hashes match the original version-controlled fixtures exactly.
Generator #30 joins the opt-in targets, reaching 112/115 TSVs.
No existing TSVs, test consumers, runtime GMP oracles, or production
sources are changed.

B108 adds independent C++20 `generator_floatmp_tanh.cpp`, recreating
all 1,148 source-sorted `floatmp_tanh.tsv` rows. All 13 original
precision sets and mode groups are generated algebraically (17, 29,
or 18 source rationals per mode); exact source keys are independently
rounded nearest-even before applying a certified tanh. Use B=2^2304,
an outward Taylor enclosure of exp(|x|/32) through term 409 with
next-term remainder and six outward squarings to enclose exp(2|x|).
Directed division bounds (exp(2x)-1)/(exp(2x)+1). Large |x|=1000
is certified by the strict nonattainment of 1 and the bound
1-tanh(1000)<2^-1799, so directed rounding returns the correct
predecessor at p<=1024 rather than spuriously returning 1.
Preserve all signed nearest-even and N/D/U/Z modes and the distinct
p-even/p-odd cutoff rules. Original comment lines stay unchanged.
Strict GCC and Clang C++20 warning-as-error builds, read-only
--verify, and deliberate fixture-corruption rejection all passed.
The generated 345,917 bytes have exactly Git blob SHA-1
`7bc870f9e9f1006e3f5c0850c6e1eed1c5e68a2b`, identical to
the original fixture. Generator #31 joins opt-in targets, reaching
113/115 TSV files. No original TSV, test consumer, runtime GMP
oracle or production code changes.

B109 independently regenerates all 4,794 references in `floatmp_rounding.tsv` without consulting
source TSV operands or expected values. Rebuild 232 base rows, 832 additional multiprecision
quotient rows and 3,730 D79 references in precisely their original first-visit order, preserving the
original mid-file section comments. The D79 inputs include 320 source ratios, 320 dyadics, 2,520
small-limb numerator/denominator and quotient visits, and 570 new integer arithmetic cases, with
original duplicate removal. Form compound quotient inputs only after independent nearest-even source
rounding, and compound arithmetic operands only after their own N/D/U/Z conversion. Round each
resulting exact signed ratio with independent 32-bit-product limb arithmetic, exact division and
nearest-even midpoint parity. No expected values, source operands, Numeric, GMP or Boost are taken
from TSV files. GCC and Clang strict C++20 warnings-as-errors builds and read-only --verify both
pass. An intentionally corrupted fixture is rejected. The regenerated 1,225,939 bytes match the
original Git blob `8d5c927485bc2c171853bb1502ed6b9d0d2aa51b` byte-for-byte. Generator #32 joins
opt-in CMake targets, reaching 114/115 original TSV files. Only `big_unsigned_sqrt.tsv` remains.

B110 completes the generator registration for all 115 original TSVs.
`generator_big_unsigned_sqrt.cpp` independently computes all 161
integer floor-square-root/remainder triples with radix-4 restoring
arithmetic and separately checks n=root^2+remainder and the strict
upper bound (root+1)^2>n using a binary shift-and-add square.
The first 22 small edge inputs and 69 powers-of-two +/-1 cases are
recreated algebraically. **Important provenance limitation:** the
original 37 random input magnitudes and the 11 original random square
roots used to build 33 adjacent-square cases do not have a recoverable
PRNG/seed in repository history. To preserve every original case, 48
*input-only entropy seeds* are frozen in
`big_unsigned_sqrt_legacy_seeds.hpp`. These are historical input
identities, NOT precalculated oracle outputs. All 161 root/remainder
results are independently recomputed; the original fixture is not
read during generation. Therefore this last fixture is reproducible
from committed seeds, but not reconstructible from the original random
stream. This qualification must remain explicit when reporting 115/115.
An independent JavaScript BigInt radix-4 restoring implementation
reproduced all 161 triples, exactly 48,161 bytes and original Git
blob `686e648fbe32f70303b1b01c38629bd03e3672a3`.
The C++ source compiled in strict-warning C++20 with GCC and Clang;
full compiled execution against the committed input bank and CI
of B110 must still be confirmed. No original TSV, test consumer,
production source, or runtime GMP oracle changed. Generator #33
joins opt-in targets `generators`, `generate-data`, `verify-data`.

B111 consolidation: `generator_floatmp_rounding.cpp` now emits both
`floatmp_rounding.tsv` (4,794 rows, SHA `8d5c927485bc2c171853bb1502ed6b9d0d2aa51b`)
and `floatmp_add_sub.tsv` (1,780 rows, SHA
`b4552cfbdea659730486dd5bcdb3bb3a62be7984`). These two files
remain byte-for-byte unchanged; remove `generator_floatmp_add_sub.cpp`
and its CMake executable. Generator count becomes 32 for all 115 TSVs.

The previously unrecoverable 48 legacy input-only seeds in
`big_unsigned_sqrt_legacy_seeds.hpp` are no longer needed. The
existing `generator_big_unsigned_sqrt.cpp` now specifies the entire
input stream via SplitMix64 seeded with `0x526f6f74534d3230`,
64-bit wraparound, low-limb-first word order, exact-width top-bit masks,
and sparse zero-low-limb cases. Keep the same 161 cases, including
37 random magnitudes in their original size classes and 33 boundary
cases around squares of 11 newly generated roots. The first 91
small and power-of-two input cases remain identical; 70 random-derived
cases intentionally change, and all sqrt/rem expectations are
recomputed independently. The new TSV is 48,272 bytes and Git SHA
`5ed4a8d6705a2b39dd44d5ebcf9f5f666019d2bd`. Strict GCC/Clang
C++20 builds, independent Python math.isqrt + PRNG cross-checks, and
both generator --verify/corruption checks passed. No test consumer or
production implementation was modified. Standard CI tests do not yet
build or run the opt-in `verify-data` target.

### A1 — first runtime GMP oracle removal (24 October 2026)

Replace the complete runtime GMP numerical reference in
`TestBigUnsigned::test_product_accumulation()` by a new independent
source-keyed `big_unsigned_accumulation.tsv`. The original nested loops
are retained with all mt19937_64(0x23ACC00) draws, exact width choices,
zero/all-ones/sparse source transforms, all 9 x 32 x 3 x 2 = 1,728
operation visits and overflow/borrow escape checks. The TSV records
length, trial, extra capacity, sign of accumulation, both operands,
the original buffer, exact modular result and escaped flag. The
consumer verifies source keys and expected values without GMP,
checks all operand/buffer immutability, canonical limbs and buffer
length, rejects missing/out-of-order/extra rows, and retains original
exception/unchanged destination assertions. The generator now lives
in existing `generator_big_unsigned_mul.cpp` alongside B35/B37/B38,
not in an additional executable: 116 TSVs and still 32 generators.
Independently regenerated 1,728 cases with original C++ mt19937_64
and a second JavaScript BigInt mt19937_64 implementation; both agree
byte-for-byte on 1,441,196 bytes, Git blob
`2e1d31c8b7e9339693ddbaa7f22d2741d213012e`. Standalone GCC and
Clang strict C++20 generator compilation plus --verify passed on the
fixture. This is one converted GMP runtime test method; other GMP
numerical oracle methods have deliberately not been touched.

### A2 — complete B32 runtime GMP square oracle migration (25 October 2026)

Migrate all 11 x 192 = 2,112 original B32 squared-input visits in
`TestBigUnsigned::test_multiplication()` from GMP to a source-keyed
independent `big_unsigned_b32.tsv` oracle. The seed
`mt19937_64(0xB3205A7E)`, 11 original widths, all 192 repetitions,
all random draws (including overwritten draws), and the ordered
13/17/19/23 divisibility transforms are unchanged. Each row freezes
original width, repetition, canonical input hex and exact square.
The consumer checks every source key and input, both `square(x)` and
`mul(x,x)` against the independent exact expected square, source
immutability and canonical shape; it rejects missing or extra rows.
Only this B32 numerical GMP oracle is removed; other GMP checks
remain. The existing `generator_big_unsigned_mul.cpp` uses its
independent exact binary shift/add multiplication to produce the
2,112-row B32 fixture alongside B35/B37/B38 and A1: 117 TSVs,
still 32 generator executables. The independently reproduced B32
fixture is 3,328,858 bytes, Git blob SHA
`a9afcd45514a36df0faf7167fa275bb7569d6d1f`.
GCC 14 and Clang 17 strict C++20 standalone source builds regenerate
identical B32 bytes. Full CI/consumer validation is reported
separately from generator-level validation; do not count it as passed
until the test jobs complete. No production code was modified.

### A3 — 128/64 edge and B26 runtime GMP migration (25 October 2026)

Migrate all 192 original 128/64 edge visits (8 divisors x 4 high-word
choices x 6 low words) and all 2,048 B26 randomized DIVQ visits in
`TestBigUnsigned::test_division_limb()`, without reducing coverage.
The original `mt19937_64(0xB260D1B1)` source, conditional RNG draws,
divisor/high/low selections and high<divisor precondition are preserved.
The consumer verifies every ordered source key, canonical numerator
and operands, exact full 64-bit remainder and quotient, unchanged
source, original saturation assertions, and missing/extra TSV rows.

The existing standalone C++20
`tests/data/generator_big_unsigned_power_div.cpp` now also emits
`big_unsigned_div128by64.tsv` using independent aligned-subtraction
long division and binary shift/add multiplication for the identity
numerator=divisor*quotient+remainder and remainder<divisor. Neither
GMP nor production division is used in the generator. Exactly
2,240 rows, 227,016 bytes, Git blob
`f9084b5d88d01b2acf356e5f4b0cb9a0435ff7cd`.
Strict GCC/Clang C++20 compilations agree byte-for-byte; a separately
implemented Python mt19937_64 and built-in multiple-precision
division reproduce all 2,240 source keys, quotients and remainders.
118 source TSV fixtures are registered across the same 32 generators.
Other division/GMP blocks, all production sources and all previously
committed fixtures remain unchanged. Consumer regression/CI results
must be checked separately after publication; the opt-in full
`verify-data` target is not part of standard CI.

### A4 — seeded random single-limb division numerical oracle (25 October 2026)

Convert the original first 6,000 visits of `TestBigUnsigned::test_division_limb()`
from runtime GMP to source-keyed `big_unsigned_div_limb_random.tsv`.
Preserve `mt19937_64(0xD1D1D1)`, every random draw, original
0/64/1..7-limb selection and ordered modulo-13/17/19 transforms,
all divisor selection branches and both owning/non-owning division calls.
Every row contains repetition, exact decimal divisor, canonical dividend,
exact quotient and full 64-bit remainder. The consumer verifies input
identity, ordered keys, both quotient paths, full-width remainder,
scratch/output guard clearing, source immutability, and all-row exhaustion.

Extend the already registered C++20
`generator_big_unsigned_power_div.cpp` with independent aligned-subtraction
long division and binary shift/add multiplication; certify Euclidean
identity n=q*d+r and remainder<divisor independently, not from Numeric
or GMP. A second Python implementation of multiple-precision division
and an independent mt19937_64 engine reproduces all 6,000 cases
byte-for-byte; the first eight random draws also match a compiled
std::mt19937_64 C++20 witness. Fixture 1,194,017 bytes, Git blob
`5e9b7c6801c1fe9c6220a23235b50f592c8b1b84`.
119 TSVs across 32 existing generators. B48, B49, B57/B63 and B61
GMP blocks remain for later steps, as do GMP tests in other methods.
No production code or prior fixture changed. Standard CI does not
run the opt-in all-fixture `verify-data` target.

### A5 — B48/B49 numerical oracle conversion (25 October 2026)

Convert every original B48 and B49 numerical GMP division comparison
(560 + 784 = 1,344 visits) to source-keyed TSV results while preserving
original mt19937_64 random draws, all dividend transforms, both
owning/non-owning quotient paths and exact full-width remainder checks.
The consumer validates ordered width/divisor/mode, original source,
canonical expected quotient/remainder, output guard clearing and
TSV exhaustion. Other GMP blocks B57/B63 and B61 are not yet migrated.

The consolidated C++20 power/division generator independently computes
all four new A5 fixture files (B48 560; B49 784; B57 1,800; B61
1,560 visits) using portable bit-at-a-time division with a 65th
remainder carry and independent binary shift/add Euclidean identity
checks. C++20 GCC and Clang mt19937_64 input witnesses and exact
Python big-integer division independently reproduce all four files.
SHA Git blobs: B48 `1a56559c3f5c529a5c222fe1136b69be1d3aed28`,
B49 `4f82aac238df03238128ad5d41be21b9e2b24a22`,
B57 `dc2a5205d8aa4a912f73f9e675522b1531f65794`,
B61 `f9e8a720ece587c13d68802588cfe75ebb006e67`.
Only B48/B49 consumer conversions are counted as migrated here.
CI for A5 is pending. No production changes.

### A6 — B57/B63 and B61 numerical GMP oracle removal (25 October 2026)

Complete the remaining source-keyed replacement in
`TestBigUnsigned::test_division_limb()`: B57/B63 contributes
15 widths x 15 divisors x 8 modes = 1,800 original seeded
visits, and B61 contributes 13 widths x 15 divisors x 8 modes =
1,560 visits, for 3,360 new GMP-free oracle visits.
Use the existing independently certified
`big_unsigned_b57.tsv` and `big_unsigned_b61.tsv` fixtures
added in A5, retaining every original `mt19937_64` draw,
all seven explicit transforms and the unchanged eighth random
mode, original width and divisor sets and visited order.

The B57 consumer validates the full quotient and 64-bit remainder
from `divrem_limb_into` against the exact TSV, retaining output
guard zeroing, then checks `divrem_limb` quotient and remainder
against the certified result. B61 validates `divrem` quotient
and full remainder against the exact TSV, preserving the
independent `rem` call and checking it matches that certified
remainder. All expected values are checked through the canonical
hex parser, including zero, and TSV rows must be complete and
ordered. Original invalid-input, short-output and unchanged-buffer
exception assertions are retained.

No numerical GMP oracle remains inside `test_division_limb()`;
other GMP tests and production implementation are unchanged.
Total A1-A6 migrated test visits: 16,784.
All existing 123 TSV fixtures retain their committed hashes;
there are still 32 standalone generator executables. A5's
independent Python and portable C++20 verification of B57/B61
fixture arithmetic is reused, not recharacterized as a new run.
Full A6 consumer CI and opt-in `verify-data` remain separate
checks; report CI as passed only when observed.

### A7 — exact single-limb division GMP oracle migration (25 October 2026)

Migrate all 6,200 visits of `test_exact_division_limb()` (10 original
limb lengths x 155 original divisors x 4 modes) from runtime GMP
generation and numerical comparisons to four source-keyed exact TSV
fixtures. Preserve the full std::mt19937_64(0x22E1AC7) draw sequence,
all original divisors, input-word transforms, both low-level
`divexact_limb_into` and owning `quot_exact` result comparisons,
canonical input/quotient shape, source immutability, output guards and
all invalid/inexact/unchanged-buffer exception checks. The consumer
rejects missing, reordered or extra records.

Existing independent C++20 `generator_big_unsigned_power_div.cpp` regenerates the datasets via
Nat::mul_small and independently checks the exact inverse using the bitwise single-limb divide
implementation; no Numeric, GMP or Boost dependency. Independently replaying all 6,200 original
std::mt19937_64 cases compiled with both GCC and Clang and Python multiple-precision multiplication
reproduces every fixture byte-for-byte. Git blobs, respectively small (4,960 rows)
`d4a68ff02f44a9273ca8bf0a30aa4874a57ab295`, length 65 (620 rows)
`6a5c0cdee6d167bfedaff54917cc725bd73793ae`, length 129 first 78 divisors (312 rows)
`0697660ce3945381db5840683f408ee73de55e1d` and remaining 77 divisors (308 rows)
`c42a99815380e11bce4f3976773b8507483f47bb`. No production sources changed; the same 32 generator
executables now cover 127 TSVs. A1–A7 cumulative migrated visits: 22,984. A7 consumer CI is to be
checked after publication, separately from local fixture verification and opt-in `verify-data`.

### A8 — B57 power-of-two multiprecision division oracle (25 October 2026)

Convert all 2,160 original B57 visits in
`TestBigUnsigned::test_division()`: nine original dividend widths x
eight divisor widths x five bit positions x six input modes. Preserve
the complete seeded std::mt19937_64(0xB570B172) input stream, including
draws subsequently overwritten, all original ordered transforms,
the power-of-two divisors, and all owning and caller-owned division paths.
Each source-keyed TSV row records original width, divisor width,
bit position and mode, canonical dividend and divisor, exact
quotient and remainder. Consumer checks every key and both sources,
`divrem`, `quot`, `rem`, `divrem_into`, zeroed quotient/remainder
output guards, scratch guard integrity and absence of extra records.
Other GMP oracles in `test_division()` are intentionally unchanged.

The already registered independent standalone C++20 generator
`generator_big_unsigned_power_div.cpp` produces three fixtures
using Nat right shift and exact subtraction of the shifted quotient;
it explicitly checks the Euclidean identity and remainder bound.
GCC and Clang independently reproduce all 2,160 seeded input
visits, while Python exact bit-shift and mask calculations produce
the same source keys, dividends, divisors, quotients and remainders.
Fixture blob SHAs and counts: short widths 1,2,8,16,32 (1,200 rows,
755,283 bytes) `35e1d41e2b862720d835a52fe477093aaeeb21a4`;
medium 64,65 (480 rows, 826,491 bytes)
`a4f0689a89079310a8bd85f066a9dfe69b687d45`;
large 96,128 (480 rows, 1,301,002 bytes)
`c5a0a0bf80b2de7a4cf1eecca40365822f0a1fe2`.
The three committed Git blobs match Python-produced bytes exactly.
The total is now 130 versioned TSVs across the unchanged 32 generator
executables, 25,144 cumulative A1-A8 migrated visits.
No production changes. Full consumer CI and opt-in generator-wide
`verify-data` are separate validations.

### A9 — B55 multiprecision division GMP oracle migration (25 October 2026)

Replace all 6 divisor widths (7,8,9,16,32,64) x 3 length
offsets (2,9,32) x 24 original input modes = 432 original B55
GMP oracle visits in `test_division()`. Preserve the full
`std::mt19937_64(0xB550DEC1)` draw stream, all ordered source
transforms, the divisor high-bit normalization, and the 7/8/9
kernel dispatch boundary. Exact source-keyed fixtures check both
operands, quotient and remainder on `divrem`, `quot`, `rem` and
`divrem_into`; existing output-capacity, zero-surplus and
scratch-guard assertions remain. All rows must appear in original
order, with none missing or extra. GMP blocks outside B55 remain.

Extend the existing independent C++20 power/division generator with B55 input replay, independent
Nat aligned-subtraction long division and shift/add Euclidean identity checks; no Numeric, GMP or
Boost is used. Separately regenerate source inputs with std::mt19937_64 compiled using GCC and
Clang, and use Python multiple-precision divmod to certify every result. Exact bytewise Git blob
identities: B55 short (216 visits, 178,716 bytes) `929963b14605a7a0eaa3f78aeacb8b398e615858`; medium
(144 visits, 221,212 bytes) `3d8e5a3070ff80af9b01cd59a8b0a01a4da53882`; large (72 visits, 239,017
bytes) `7837342521a24a6cb672dbd364b37f0f4f44ab70`. Total: 133 versioned TSV fixtures and 32 existing
generator executables; A1–A9 cumulative migrated test visits: 25,576. No production code changes.
Full consumer CI and opt-in `verify-data` must be reported separately from oracle checks.

### A10 — B27 fused division oracle migration (25 October 2026)

Migrate both original fixed adversarial quotient-correction cases and
all 2,048 original mt19937_64(0xB27071) randomized visits, preserving
the original conditional RNG draws and ordered 13/19/23/29 transforms.
The source-keyed consumer separately compares divrem, quot and rem
against independent exact quotient and remainder, checks source
immutability, scratch sizing and remainder bounds, and detects
missing or extra records. The shared GMP oracle remains for other
B28/B36/B29/B45 tests, which are unchanged.

The existing independent C++20 generator regenerates the two B27 TSVs via aligned-subtraction Nat
division, independently checks the Euclidean identity using binary shift/add multiplication, without
GMP/Boost/Numeric. Reproduction through independently ported mt19937_64 and Python bigint division
yields the committed blobs: first (1,026 rows, 847,088 bytes)
`a4a2aa7f948a8bd9a673700f48bcd5efbb5ed101`; second (1,024 rows, 844,784 bytes)
`18448ac0544c03bddb54c22488f507dd6778ffb2`. 135 fixtures, 32 standalone executables, 27,626
cumulatively migrated A1-A10 visits. Production code remains untouched. Consumer CI and optional
verify-data are separate checks.

### A11 — B48 multiprecision division source-keyed oracle (25 October 2026)

Migrate all 8 widths x 14 divisors x 3 modes = 336 original B48
numerical oracle visits in `test_division()`. Preserve seeded
`mt19937_64(0xB48A11D1)` draws, transformation precedence, and the
three independent result APIs. The existing independent C++20 generator
now emits `big_unsigned_b48_multiprecision.tsv` with 336 source-keyed
exact results, checking the Euclidean identity. Git blob
`11042c7f930f66082d8a89d6670376c03786938d` (305,944 bytes).
The GMP oracle for B28 and other division groups remains unchanged.
The audit is maintained separately at `doc/designs/gmp_migration.md`.

### A12 — first 1,250 randomized BigUnsigned GMP-oracle visits

Use an independent source-keyed TSV for visits 0..1249 in
`test_randomized()`, retaining GMP for visits 1250..9999. The
original seeded random stream is shared across both branches.
Exact comparison, addition, conditional subtraction and both shifts
are verified; fixture blob `487f06f3921d72d2de9caee06bf9c23957ee6012`.
Progress tracked in `doc/designs/gmp_migration.md`.

### A13 — second 1,250 randomized arithmetic oracle visits

Add `big_unsigned_a13_randomized_second.tsv` for original seeded visits
1,250–2,499; preserve all original random draws, including the first
1,250 discarded only by the independent generator. The consumer compares
all original arithmetic and shift outputs against independently
certified exact results and continues using GMP for visits 2,500–9,999.
Exact fixture 1,217,436 bytes, Git blob
`ae3a7c98ef3ac9dedea8e44d5a2d410daccb342e`.
The audit is updated in `doc/designs/gmp_migration.md`.

### A14 — third 1,250 randomized arithmetic oracle visits

Migrate original visits 2,500–3,749 to independent
`big_unsigned_a14_randomized_third.tsv`, with source-keyed
input, shift and expected arithmetic. The generator replays the full
original source stream. Fixture Git SHA
`bdb0b5f9235005b2bc31d18a7e3b21c2af474abe`, 1,250 rows.
Migration inventory: `doc/designs/gmp_migration.md`.

### A15 — fourth 1,250 randomized arithmetic oracle visits

Migrate original visits 3,750–4,999 to independent
`big_unsigned_a15_randomized_fourth.tsv`, retaining the single original
random engine and all source/input/shift transformations. GMP remains only
for visits 5,000–9,999 in `test_randomized()`. The fixture has 1,250 rows,
1,216,271 bytes and Git blob
`e23f97770356c28d0b76aef3c7d9d5c7dee8a0c5`.
The independent Python multiple-precision witness reproduces the preceding
A14 blob byte-for-byte before generating A15. Migration inventory:
`doc/designs/gmp_migration.md`.

### A16 — fifth 1,250 randomized arithmetic oracle visits

Migrate original visits 5,000–6,249 to independent
`big_unsigned_a16_randomized_fifth.tsv`, preserving the single original
`mt19937_64(0xBA17B00D)` stream, all input transforms and the shift draw.
GMP remains only for visits 6,250–9,999 in `test_randomized()`.
The fixture has 1,250 rows, 1,204,292 bytes and Git blob
`bbb0fc03a3dae88154205334fafdff5fb04f9619`. A separate multiple-precision witness first reproduces
the complete A15 fixture byte-for-byte before computing A16. Migration
inventory: `doc/designs/gmp_migration.md`.

### A17 — sixth 1,250 randomized arithmetic oracle visits

Migrate original visits 6,250–7,499 to independent
`big_unsigned_a17_randomized_sixth.tsv`, preserving the original single
`mt19937_64(0xBA17B00D)` stream, conditional input transformations and shift
draw. GMP remains only for visits 7,500–9,999 in `test_randomized()`.
The fixture has 1,250 rows, 1,224,250 bytes and Git blob
`2b595bab5b7a83fc2c9a27c188d2fba204b6c30e`. A separate multiple-precision witness first reproduces
the A16 fixture exactly before independently computing A17. Migration
inventory: `doc/designs/gmp_migration.md`.

### A18 — seventh 1,250 randomized arithmetic oracle visits

Migrate original visits 7,500–8,749 to independent
`big_unsigned_a18_randomized_seventh.tsv`, preserving the original single
`mt19937_64(0xBA17B00D)` stream, conditional input transformations and shift
draw. GMP remains only for visits 8,750–9,999 in `test_randomized()`.
The fixture has 1,250 rows, 1,212,953 bytes and Git blob
`b1c7416dd60cc79511d545a706df56b18b70bf74`. A separate multiple-precision witness first reproduces
the A17 fixture exactly before independently computing A18. Migration
inventory: `doc/designs/gmp_migration.md`.

### A19 — final 1,250 randomized arithmetic oracle visits

Migrate original visits 8,750–9,999 to independent
`big_unsigned_a19_randomized_eighth.tsv`, preserving the original single
`mt19937_64(0xBA17B00D)` stream, conditional input transformations and shift
draw. This removes the remaining numerical GMP branch from
`test_randomized()`; all 10,000 original visits now use source-keyed exact
external expectations. The fixture has 1,250 rows, 1,202,168 bytes and Git
blob `15f9d3656bf5fe77ac2b97e0ba6faee3b0883ea3`. A separate multiple-precision witness first
reproduces the A18 fixture exactly before independently computing A19.
Migration inventory: `doc/designs/gmp_migration.md`.

### A20 — randomized decimal parsing oracle migration

Migrate the original 1,800 randomized decimal parsing visits in
`test_gcd_decimal()` to `big_unsigned_a20_decimal_random.tsv`.
The consumer preserves the existing `mt19937_64(0xB16B00B5)` state by
executing the complete preceding 2,800-visit B76 stream before generating
the same decimal strings. Each TSV row source-keys repetition, length and
full decimal input, then supplies exact hexadecimal and canonical decimal
expectations.

`generator_big_unsigned_decimal.cpp` is standalone C++20 and uses the
independent `Nat` arithmetic in `oracle_exact_integer.hpp`; it does not
link Numeric, GMP or Boost. The fixture contains 1,800 rows, 1,011,827
bytes and Git blob `5238d33947659f227a12efaaf8ffd63fa8770903`.
An independently compiled C++20 reproduction produced the same Git blob,
providing a byte-for-byte check of the random stream and expected values.
Migration inventory: `doc/designs/gmp_migration.md`.

### A21 — structured decimal parsing oracle migration

Migrate the remaining 484 structured parsing visits in
`test_gcd_decimal()` to `big_unsigned_a21_decimal_structured.tsv`.
The fixture preserves all original boundary strings, B47/B50 randomized
strings and the 25 exact `2^(bits-1)+offset` sources. The consumer checks
sequential source identity for every explicitly constructed input; for the
25 power cases it additionally checks the independently frozen parsed value
against the original native shift/add construction.

`generator_big_unsigned_decimal.cpp` replays B76 and A20 to preserve the
shared `mt19937_64(0xB16B00B5)` state before generating B47/B50 random
strings. Exact expected values use only the standalone `Nat` implementation.
The fixture contains 484 rows, 279,466 bytes and Git blob
`17a68c45126d418b00f9e984e95ea9bd486fef8e`. A separately compiled C++20
reproduction generated the identical blob. After A21 no numerical GMP
oracle remains in decimal parsing; only decimal formatting cases remain in
`test_gcd_decimal()`.

### A22 — decimal formatting oracle migration

Migrate all 1,509 remaining decimal-formatting visits in
`test_gcd_decimal()` to `big_unsigned_a22_decimal_format.tsv`.
Each row freezes the exact source magnitude in canonical hexadecimal and
the independently generated decimal rendering. Source construction in the
consumer preserves the original limb patterns and seeded streams; decimal
power cases are rebuilt natively from the same mathematical inputs before
being checked against the frozen source hex.

The fixture covers B25 (50), B25.1 (18), B30 (128), B42 (208), B46 (352),
B51 (288), B59 (216), B43 (80), B31 (119), and the final radix-power block
(50), for 1,509 visits total. `generator_big_unsigned_decimal.cpp` uses
only standalone `Nat` arithmetic and exact decimal conversion. The fixture
is 1,404,780 bytes with Git blob
`095b30979096b45222ebf019bbd5482588deb3dd`; a separately compiled C++20
reproduction generated the identical blob. After A22, `test_gcd_decimal()`
contains no numerical GMP oracle calls.

### A23 — scratch-buffer arithmetic oracle migration

Migrate the 5,000 original buffer arithmetic visits in
`test_views_and_scratch()` to five source-keyed TSV segments of 1,000
visits each. The consumer preserves the single original
`mt19937_64(0xB22B22)` stream, all-ones transformations and the conditional
`low` draw performed only for repetitions congruent to 4 modulo 5. Each
row freezes generated widths, `low`, both source magnitudes, exact sum,
conditional nonnegative difference and truncated product.

`generator_big_unsigned_views.cpp` is standalone C++20 and computes the
expectations with independent `Nat` add/sub plus shift-add multiplication.
The five fixture blob SHAs are
`596e3b10f9ba28a70b36ed6da6fce94080839380`,
`6a51e17a1a70c6247a54c7652b881b2ee357bb5d`,
`cb3d6a863d56721b6dc6f959141f539adbd4a829`,
`9f16abdf7e4cff8845e0b1ef13ae91cc4d0e5738`, and
`d2cc565c0eadd9b8fe34fda785c072891bb40d81`. A separate C++20
multiprecision witness reproduced all five files byte-for-byte. R38 shifted
operations remain on their existing oracle for A24.

### A24 — shifted add/sub oracle migration

Migrate all 6,600 R38 shifted add/sub visits in
`test_views_and_scratch()`. The 264 base source pairs are frozen once in
`big_unsigned_a24_shifted_sources.tsv`; ten expectation files cover the
25 shift-pair combinations per base pair. The first seven expectation files
contain 825 rows each and the three 129-limb-heavy tail files contain 275
rows each, avoiding duplication of the large base operands while preserving
source keys `(width, case_index, sa, sb)` for every visit.

The generator replays the exact `mt19937_64(0xB38038)` stream and all
all-ones transformations, then computes shifted values, sums and absolute
differences with standalone `Nat` arithmetic. A separate C++20
multiprecision witness reproduced the source TSV and all ten expectation
TSVs byte-for-byte. Source blob: `5fd306317e0dc7f71bae8af3b5d894d991101589`.
Expectation blobs: `ed88323a46cd24ce210a77b2299fda63c7aaae08`,
`e35de3e0fc8e5c1d499653a87c777be53a2cd398`,
`0c760d3f73e91f0b8a56c26d2ad317cf4b38a947`,
`b3ea7798c97b96eaf8c473b01c915b85756e2397`,
`71f88fc348dab5a9869b0ffa7fe3ce4a56c27f74`,
`b05c52c00db7e6fffbd16a630873f541292f9ac5`,
`ad4f6c8420a9ce90ea72cec26a972b6e3904e561`,
`ada36cb9605cf4f70b0026ad8ebbd58f89d93155`,
`0d684e99ccad45c483eda803b8dab2e9be86ab99`, and
`2fc3a3aa79ac055a4daf38f786664481d0525ed2`.
After A24, `test_views_and_scratch()` contains no numerical GMP oracle.

### A25 — B52 multiply-kernel oracle migration

Migrate all 288 B52 n/m/mode visits in `test_multiply_kernel_paths()`
to `big_unsigned_a25_b52.tsv`. The consumer preserves the exact
`mt19937_64(0xB5203232)` stream and all five structural transforms,
verifies both source magnitudes, then checks the full product and the
`full-1` truncated product from the fixture while retaining all kernel
agreement, oversized-output and zero-tail assertions.

The existing standalone `generator_big_unsigned_mul.cpp` now emits A25
using its independent shift/add `Nat` multiplication and low-product
helper, with no Numeric, GMP or Boost dependency. The fixture contains
288 rows, 774,359 bytes and Git blob
`db88450f727a185a1b5543e1b2323081d6893088`. Normal CI and
`verify-data` remain distinct validation gates.

### A26 — B53 sparse-row multiply oracle migration

Migrate all 2,304 B53 low-product visits in `test_multiply_kernel_paths()`.
The fixture stores 384 original source pairs and six exact low products per
pair, preserving the full `mt19937_64(0xB530B040)` stream, all density
transforms and every original truncation length. The consumer still compares
selected, fused and x86-64 carry kernels word-for-word before checking the
independent low-product expectation.

The existing standalone multiplication generator emits
`big_unsigned_a26_b53.tsv` using independent shift/add `Nat` arithmetic.
The fixture contains 384 source rows / 2,304 visits, 2,112,618 bytes, Git
blob `a22d13dc7e098789870a45cd38625a0f00078ebd`. Normal CI does not imply
`verify-data` success.

### A27 — general multiply-kernel oracle migration

Migrate all 6,000 general-kernel visits in `test_multiply_kernel_paths()`
to six source-keyed TSV segments of 1,000 visits each. The consumer keeps
one continuous `mt19937_64(0xB3B3B3)` stream, all original width draws,
all-ones/zero transforms and the conditional low-product length draw. Each
row freezes both source magnitudes, the exact full product and the exact
truncated product for the original `low` length.

`generator_big_unsigned_mul.cpp` emits the six A27 files using independent
shift/add `Nat` multiplication. Segmenting is only for repository/tooling
size; execution order and RNG state are unchanged. No production code is
modified. Normal CI and `verify-data` remain separate gates.

### A28 — B8 and fixed Karatsuba oracle migration

Migrate the final 2,840 visits in `test_multiply_kernel_paths()`: 700 B8
sources exercised at four cutoffs (2,800 visits) and eight fixed sources at
five cutoffs (40 visits). The fixed corpus continues the exact
`mt19937_64(0xB8B8B8)` stream after all 700 B8 sources. The consumer keeps
all scratch sizing, high-zero output checks, guard-limb assertions and
cutoff-specific execution; only the exact product oracle moves to TSV.

`generator_big_unsigned_mul.cpp` emits `big_unsigned_a28_b8.tsv` and
`big_unsigned_a28_karatsuba_fixed.tsv` with independent shift/add `Nat`
multiplication. The B8 fixture is 2,445,194 bytes (blob
`16848b2efecd3cdc5535e03e0e4694df1890e09f`); the fixed fixture is 43,579
bytes (blob `51233e8a709ad097bfd1435c0a9171fdef7dcdb6`).

### A29 — general division oracle migration

Migrate all 6,000 general randomized visits in `test_division()` to six
1,000-row source-keyed TSV segments. The consumer preserves the single
`mt19937_64(0xB15D1A1D)` stream, all width choices and every source
transformation. Each row freezes dividend, divisor, quotient and remainder;
all existing joint/single-output, scratch, guard, zero-tail and Euclidean
identity checks remain active.

`generator_big_unsigned_division_remaining.cpp` is a standalone C++20
generator using independent aligned-subtraction long division and shift/add
multiplication. The six A29 fixture blobs are
`7b6b55f4727df0f458c922df46365d8d2f92dff7`,
`9c29bfd502d7abc0ebb27853464170c96d0d311d`,
`388adcbd50795a61c795e32b8403b5858681f76f`,
`7789485fe33bc3a7c31de069537ae3ddcf4353f4`,
`432a61307b1af89ff8548289857d447d6629841b` and
`5ee40a6f52a8c94c765da8dcbad76c07bf1e6ac7`.

### A30 — B28 equal-high division oracle migration

Migrate all 2,048 B28 visits to two source-keyed TSV segments while
preserving the exact `mt19937_64(0xB280C27)` stream, equal-high transforms,
all three public division modes, surplus-output clearing and the fast-path
requirement that scratch/guards remain untouched. Each row freezes numerator,
denominator, quotient and remainder.

The standalone remaining-division generator emits two 1,024-row files using
independent exact long division. Their blobs are
`1f44f222aab02b942935c47b43bf6307577ba47c` and
`d41a4d633262ace82442664c05daf0fb9227f4cd`.

### A31 — B36 high-ratio division oracle migration

Migrate all 1,280 B36 equal-size high-ratio visits to two exact TSV
segments. Preserve all ten operand widths, 128 repetitions per width, the
single `mt19937_64(0xB360A11)` stream, saturated/zero transforms and every
q=0/q=1/q>=2 boundary class. Joint and single-output division paths remain
checked against the frozen exact quotient/remainder.

The remaining-division generator emits two independently computed fixtures;
blobs: `d7f1e6ddeccf60d4085baf5eaedf46b7bb920c21` and
`9bf497e9817766610c9ff96b7511b58c92fd04b4`.

### A32 — B29 sparse division oracle migration

Migrate all 1,024 B29 sparse/dense divisor visits to two exact TSV segments.
The consumer preserves the exact conditional `mt19937_64(0xB290D1F5)`
draw sequence, sparse midpoint injection, partially dense divisor draws,
all-ones dividends and high-limb overrides. Joint, quotient-only and
remainder-only paths remain checked against frozen exact results.

The remaining-division generator emits two 512-row fixtures; blobs:
`f56baf8176a40e6b0f51967d653754a0d5451735` and
`0b74ec4c834107fe055776195f80b8e9f2e5c721`.

### A33 — B45 zero-low and decimal-power division

`generator_big_unsigned_division_remaining.cpp` now also emits
`big_unsigned_a33_b45.tsv`, covering 504 shifted zero-low-limb divisor cases
and 40 decimal-power divisor cases. The original seeded stream and all source
transforms are preserved; each row is source-keyed and contains exact quotient
and remainder computed by the standalone long-division oracle. The fixture has
544 rows, 743,901 bytes and Git blob
`a29d884326073e4972fae8a47c4610721de10003`.

### A34 — divisor-normalization division

`generator_big_unsigned_division_remaining.cpp` also emits
`big_unsigned_a34_normalization.tsv`, covering all 64 normalization shifts x
10 modes = 640 exact division cases. Each row freezes both source operands and
exact quotient/remainder; the fixture has 640 rows, 86,823 bytes and Git blob
`da602174c15f52d33838659fe2596e63cdcbdc24`. This completes the fixed
76,323-visit BigUnsigned numerical GMP-oracle migration audit.

### Decimal exact coverage oracle

`generator_decimal_exact.cpp` emits `decimal_exact.tsv`: 64 source-keyed
pairs covering equal values with different scales, unequal values, signs,
zero, sums and products. Expectations use bounded standard `int64_t`
arithmetic (input coefficients <=1250 and scales <=3), independently of
Numeric, GMP and Boost. It supports the usual `--write` / `--verify` modes.

### Integer conversion boundaries

`generator_integer_conversion.cpp` emits `integer_conversion.tsv`: 27 signed
sources around 2^31, 2^32, 2^63 and 2^64 with independently computed signed
and unsigned 32/64-bit representability. The consumer checks native int plus all four
long/long-long destination types already present in the coverage report,
selecting the frozen expectation by the platform's destination width.
The generator uses only the standalone exact Nat helper.

### Scalar comparisons and binary64 wrappers

`generator_scalar_comparison.cpp` emits 27 Int/Rational comparison rows,
using bounded exact cross-products and explicit infinity/NaN ordering.
The consumer tests the native Int comparison and generic Integer comparisons,
including the explicit public generic inequality overload on unordered values.

`generator_floatdp_wrappers.cpp` emits eight IEEE binary64 bit patterns
(both zero signs, 1/2, 3/2, 2^60 and their negatives) and their exact unary
zero, identity, negation and half results. Expectations are integer bit-field
transformations; no host floating arithmetic or Numeric oracle is involved.
Every result is checked bit-for-bit in all four builtin rounding modes.
Both standalone generators support write and byte-for-byte verify.

### Dyadic signed addition/subtraction dispatch

`generator_dyadic_add_sub.cpp` emits eight `dyadic_add_sub_<limbs>.tsv`
fixtures, 256 rows each. Widths 1/5/16/17/31/32/64/65 bracket both platform
crossovers between shifted, stack-aligned and owned algorithms. Exponent
gaps 0/1/63/64 bracket virtual-shift dispatch; both exponent orders, signs,
zero, equal and adjacent mantissas cover cancellation, carry and borrow.
The second exponent base also exercises comparison of fractions below one.
Expectations use independent signed Nat arithmetic with a common power-of-two
denominator and canonicalization. No production code is used or linked.

### Dyadic binary64 conversion and rounding

`generator_dyadic_conversion.cpp` emits 168 source-keyed dyadics with exact
truncating binary64 bit patterns, representability flags, power-of-two
denominators, floor and ceil. Expected IEEE results are assembled with
integer bit fields, without host floating arithmetic. Sources cover both
signs, zero, 53/64/117/128-bit extraction boundaries, overflow and subnormal
or underflow boundaries. Exactly representable rows also verify FloatDP
round trips; every source verifies an exact 2048-bit FloatMP round trip and
Approximation conversion. The fixture has an independent Python witness
checking the decoded result's adjacent-representable interval and exact
floor/ceil. Standalone --write/--verify are supported.

### Mixed exact interval operations and powers

`generator_dyadic_mixed_bounds.cpp` emits dyadic_mixed_bounds.tsv (285
arithmetic rows), dyadic_mixed_comparisons.tsv (300 three-valued comparison
rows), and dyadic_bounds_power.tsv (25 power rows). Small signed rational
cross-products fit int64_t throughout. Bounds are computed from endpoint
vertices; even powers include zero when the input crosses it. Comparisons
certify every result as true, false or indeterminate. Both operand orders,
singletons, disjoint/overlapping intervals and positive/negative domains
exercise the public mixed Dyadic/Rational wrappers. Quotients use divisors
excluding zero. A separate Python Fraction witness verifies arithmetic and
powers. The standalone generator supports byte-for-byte --verify.

### FloatMP rounded min/max

`generator_floatmp_minmax.cpp` emits `floatmp_minmax.tsv` with 648 exact
source-keyed cases covering NaN, infinities, signed zeros, finite order,
precision reduction and all four rounding modes. Small bounded signed
integer fractions provide all numerical expectations, including the first
nearest/even quantization to each input precision. A separate enumerated
Python Fraction witness reproduces every row. The generator supports
--write and byte-for-byte --verify independently of Numeric/GMP/MPFR.

### FloatMP exact construction

`generator_floatmp_construction.cpp` emits 42 exact constructor values in
`floatmp_construction.tsv`, reducing small signed fractions with the standard
integer GCD. Inputs fit eight significant bits, so all results are exact at
8/53/512 bits. The consumer checks default and directed constructor wrappers,
precision preservation, empty construction and self-move with inline and heap
storage. No Numeric, GMP or floating-point oracle is used.

### FloatMP native exponent limits

`generator_floatmp_exponent_limits.cpp` emits 76 rows in
`floatmp_exponent_limits.tsv`. Signed arbitrary-size exponent arithmetic
and bounded exact significand rounding compute finite native fields or
prove an overflow. Cases use native min/max values at 2/7 bits, four rounding
modes, and neighboring shifts. The consumer checks significand, exponent,
precision and kind directly, avoiding enormous power-of-two allocations.
A Python arbitrary-integer/Fraction witness independently verifies all rows.

The coverage extension adds exp(+256) and exp(-256) at the same four precisions and
rounding modes, for 320 total rows. The independent fixed scale is now B=2^1800;
Taylor degree 400 covers |x|/16<=16 while degree 165 remains sufficient for the
original inputs. The tail ratio after the first omitted term is below 1/2.
Both outward endpoints must certify the same result. All original 288 numerical
rows remain unchanged. The added negative input exercises a zero intermediate
lower endpoint and subsequent precision refinement during reconstruction.

The native logarithm fixture uses `min` and `max` source metadata instead of
materializing the input as a rational. For precision p these are exactly
`2^(INT64_MIN+p-1)` and `(2^p-1)*2^INT64_MAX`. The generator encloses the
logarithm of the normalized significand and adds the exact integer multiple
of its independently bounded log(2). The native power is computed with the
standalone unsigned limbs, including `INT64_MAX+p-1`, which can exceed a
signed 64-bit integer. Both interval endpoints must round to the same signed
p-bit value in each of the four modes. All original 444 log rows are preserved.

The tanh fixture now also includes both signs of two exact p=128 inputs,
`1/2+162899406*2^-128` and `1/2+562322472*2^-128`, in all four modes.
Their affine center corrections lie at a directed rounding boundary and
near a nearest-even midpoint. These 16 rows exercise the genuine ambiguity
fallback. Expected results use the same independent outward exp certificate
at 2304 bits; the original 1164 rows are retained (1180 rows in total).

Two additional exact logarithm sources `1-2^-98` and `1+2^-98`, at 128 bits
in all four modes, exercise zero quantized lower coefficients and an exact
upper coefficient square. The original 444 rows remain unchanged; the
standard log fixture now has 452 rows, in addition to the 32 native limit rows.
