# Reproducible exact golden-data generators

The files `generator_big_unsigned_gcd.cpp`, `generator_big_unsigned_b76.cpp`,
`generator_exact_small.cpp`, `generator_big_unsigned_power_div.cpp`, `generator_rational_exact.cpp`,
`generator_exact_b72_r17.cpp`, `generator_integer_structural.cpp`, `generator_big_unsigned_mul.cpp`,
`generator_big_unsigned_dispatch.cpp`, `generator_integer_random.cpp`,
`generator_rational_binary64.cpp`, `generator_rational_shared_stream.cpp`, `generator_rational_r31.cpp`, `generator_rational_random.cpp`, `generator_big_unsigned_random_mul.cpp`,
`generator_floatap_directed_bounds.cpp`, `generator_floatap_pi.cpp`, `generator_floatap_get_d.cpp`, `generator_floatap_power_exact.cpp`, `generator_floatap_sqrt.cpp`, `generator_floatap_pow_short.cpp`, `generator_rounded_float_directed.cpp`,
`generator_floatap_add_sub.cpp`, `generator_floatap_exp.cpp`,
`generator_floatap_log.cpp`, `generator_floatap_sin_cos.cpp`,
`generator_floatap_centered_exp.cpp`, `generator_floatap_atan.cpp`,
and `oracle_exact_integer.hpp` are **standalone C++20** sources. They use only
the standard library, with an independent unsigned-limb implementation for
exact arithmetic. In particular, **do not link them to Numeric, GMP, MPFR,
FloatMP or Boost**: expected values must not come from the code being tested.

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
overwrites only the currently covered TSV files (109 of 115 at B105).
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
c++ -std=c++20 -O2 -Wall -Wextra -Wsign-conversion -Werror tests/data/generator_exact_small.cpp -o /tmp/numeric-oracle-exact-small
c++ -std=c++20 -O2 -Wall -Wextra -Wsign-conversion -Werror tests/data/generator_big_unsigned_power_div.cpp -o /tmp/numeric-oracle-power-div
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
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatap_directed_bounds.cpp -o /tmp/numeric-oracle-floatap-bounds
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatap_pi.cpp -o /tmp/numeric-oracle-pi
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatap_get_d.cpp -o /tmp/numeric-oracle-get-d
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatap_power_exact.cpp -o /tmp/numeric-oracle-power
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatap_sqrt.cpp -o /tmp/numeric-oracle-floatap-sqrt
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatap_pow_short.cpp -o /tmp/numeric-oracle-short-pow
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_rounded_float_directed.cpp -o /tmp/numeric-oracle-rounded-directed
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatap_add_sub.cpp -o /tmp/numeric-oracle-floatap-add-sub
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatap_exp.cpp -o /tmp/numeric-oracle-floatap-exp
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatap_log.cpp -o /tmp/numeric-oracle-floatap-log
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatap_sin_cos.cpp -o /tmp/numeric-oracle-floatap-sin-cos
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatap_centered_exp.cpp -o /tmp/numeric-oracle-floatap-centered-exp
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wsign-conversion -Wshadow -Werror tests/data/generator_floatap_atan.cpp -o /tmp/numeric-oracle-floatap-atan
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
/tmp/numeric-oracle-floatap-bounds --verify --data-dir tests/data
/tmp/numeric-oracle-pi --verify --data-dir tests/data
/tmp/numeric-oracle-get-d --verify --data-dir tests/data
/tmp/numeric-oracle-power --verify --data-dir tests/data
/tmp/numeric-oracle-floatap-sqrt --verify --data-dir tests/data
/tmp/numeric-oracle-short-pow --verify --data-dir tests/data
/tmp/numeric-oracle-rounded-directed --verify --data-dir tests/data
/tmp/numeric-oracle-floatap-add-sub --verify --data-dir tests/data
/tmp/numeric-oracle-floatap-exp --verify --data-dir tests/data
/tmp/numeric-oracle-floatap-log --verify --data-dir tests/data
/tmp/numeric-oracle-floatap-sin-cos --verify --data-dir tests/data
/tmp/numeric-oracle-floatap-centered-exp --verify --data-dir tests/data
/tmp/numeric-oracle-floatap-atan --verify --data-dir tests/data
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
| `generator_exact_small.cpp` | `rational_r3.tsv`, `rational_r4.tsv`, `rational_r8.tsv`, `integer_power.tsv`, `big_unsigned_r14.tsv` | 6 + 54 + 16 + 100 + 28 = 204 exact cases |
| `generator_big_unsigned_power_div.cpp` | `big_unsigned_power.tsv`, `big_unsigned_short_div.tsv` | 373 exact power cases, 324 exact division cases; full input RNG replay |
| `generator_rational_exact.cpp` | `rational_r28.tsv`, `rational_r30.tsv`, `rational_r22.tsv`, `rational_r16.tsv`, `rational_r19.tsv`, `rational_r20.tsv` | 36 signed add/sub, 1,689 normalization, 160 scalar mul/div, 84 R16 alias/add/sub, 102 R19 product/division, 160 R20 cross-add/sub |
| `generator_exact_b72_r17.cpp` | `big_unsigned_b72.tsv`, `rational_r17.tsv` | 288 B72 products/squares and 252 R17 comparisons |
| `generator_integer_structural.cpp` | `integer_b66.tsv`, `integer_sub.tsv`, `big_unsigned_sub.tsv`, `integer_b33.tsv`, `integer_b37.tsv`, `integer_b38.tsv` | 360 power, 192 signed subtract, 100 unsigned subtract, 72 signed squares, 72 signed powers, 30 signed cubes |
| `generator_big_unsigned_mul.cpp` | `big_unsigned_b35.tsv`, `big_unsigned_b37.tsv`, `big_unsigned_b38.tsv` | 70 Karatsuba products, 91 squares, 112 products with six truncated variants |
| `generator_big_unsigned_dispatch.cpp` | `big_unsigned_b39.tsv`, `big_unsigned_b40.tsv`, `big_unsigned_b54_product.tsv`, `big_unsigned_b54_square.tsv` | 152 B39, 162 B40, 136 B54 products, 104 B54 squares |
| `generator_integer_random.cpp` | `integer_b17.tsv`, `integer_b18.tsv` | 250 signed arithmetic/shift cases and 170 signed arithmetic/division cases |
| `generator_rational_binary64.cpp` | `rational_r7.tsv`, `rational_r12.tsv` | 384 R7 rational probes; 2,048 R12 one-limb seeded conversions (both signs) |
| `generator_rational_shared_stream.cpp` | `rational_r29.tsv`, all sixteen `rational_r13_mode*_limbs*.tsv`, `rational_r5.tsv` | 216 R29, 1,536 R13, 48 R5 rows; exact shared original RNG state |
| `generator_rational_r31.cpp` | 14 `rational_r31_<bits>_<chunk>.tsv` fixtures | All 960 R31/R32 cases, 10 widths, 24 source indices, four signs; preserve 2048/4096-bit chunks |
| `generator_rational_random.cpp` | `rational_random.tsv` | Full 1,200-visit RNG replay and all 120 original stratified exact R0 cases |
| `generator_big_unsigned_random_mul.cpp` | `big_unsigned_random_mul.tsv` | All 6,000 original seeded multiplication visits, 600 selected source-keyed full, square and low products |
| `generator_floatap_directed_bounds.cpp` | `floatap_directed_bounds.tsv` | 90 certified exact rational lower/upper dyadic bounds, 6 precisions × 15 signed inputs |
| `generator_floatap_pi.cpp` | `floatap_pi.tsv` | Machin's formula: 400-term rigorous 1,200-bit integer interval, directed 384-bit and nearest-even 320-bit π |
| `generator_floatap_get_d.cpp` | `floatap_get_d.tsv` | 272 exact p-bit FloatAP-to-binary64 conversions, four precisions and four IEEE rounding modes |
| `generator_floatap_power_exact.cpp` | `floatap_pow.tsv` | 2,400 original signed integer-power references at five precisions, eight degrees and four rounding modes |
| `generator_floatap_sqrt.cpp` | `floatap_sqrt.tsv` | 300 unique input/mode exact rounded square roots; 360 source checks retained by consumer |
| `generator_floatap_pow_short.cpp` | `floatap_pow_short.tsv` | 3,680 exact boundary-focused integer powers, 5 precisions, 7 degrees, 4 directed modes |
| `generator_rounded_float_directed.cpp` | `rounded_float_directed.tsv` | 24 rigorous 128-bit exp/log/atan directed pairs, 48 original strict bound comparisons |
| `generator_floatap_add_sub.cpp` | `floatap_add_sub.tsv` | 1,780 unique source-keyed exact add/sub roundings from all 10×10 operand pairs, both operations, 4 modes and 5 precisions |
| `generator_floatap_exp.cpp` | `floatap_exp.tsv` | 288 certified, correctly rounded exp results across 18 signed inputs, 4 precisions and 4 rounding modes |
| `generator_floatap_log.cpp` | `floatap_log.tsv` | 444 certified correctly rounded logs, including 284 original centered inputs and 160 D96 cases |
| `generator_floatap_sin_cos.cpp` | `floatap_sin.tsv`, `floatap_cos.tsv` | 1,080 exact source-keyed independently certified outputs each, across D91, D92, D95 original test families |
| `generator_floatap_centered_exp.cpp` | `floatap_centered_exp.tsv` | 568 exact certified reference rows from all 720 original centered-exp visits |
| `generator_floatap_atan.cpp` | `floatap_atan.tsv` | 1,340 certified atan references reproducing all original D87, D88, D89 lookup keys |

**Scope limitation:** the repository contains many other TSVs, including
`floatap_*.tsv`, other `rational_*.tsv`, and other `integer_*.tsv`, whose
independent generator sources are **not yet versioned**. A descriptive TSV
header is not a substitute for an executable generator. These sources do
not claim to reproduce or validate those other fixtures. In particular,
floating-point transcendental oracles require rigorously outward-rounded
intervals with a proof that both endpoints round to the identical binary
result; approximate `double`/`long double` outputs are not accepted.

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

B80 adds independent signed-rational and fraction-normalization sources for
R28, R30, and R22 (39 reproducible TSVs overall). Signed source and exact
result fractions are derived without linking or executing Numeric, GMP,
MPFR, FloatMP or Boost. The R28 source columns intentionally preserve the
original unreduced numerator/denominator in cases with a common factor,
while mathematical result columns are reduced canonically. R30 includes
all nine dense fallback rows at original 128..8192-bit boundaries.
The three files were reproduced byte-for-byte against their pre-existing
Git blob SHA-1 hashes on GCC and Clang with strict CI-style warnings,
including `-Wsign-conversion -Wshadow -Werror`; a deliberately modified
fixture was rejected in `--verify` mode. No checked-in TSV or consumer
changed. Keep TSVs versioned until every generator and GMP oracle
migration is complete; no automatic `generate-data` target is introduced.

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

B88 reproduces R7 and R12 conversion references with exact integer
mantissa/exponent extraction and round-toward-zero behavior, independent
of the host FP rounding mode, runtime GMP or MPFR. The original R12
mt19937_64 source sequence is fully replayed for each of four host modes.
Both generated TSVs match the original checked-in blobs byte for byte,
with GCC and Clang strict C++20 warnings; the read-only --verify also
rejects deliberate TSV corruption. The eleventh generator is registered
in the explicit generators and generate-data targets.

B89 regenerates all 18 R29/R13/R5 files from the exact shared Rational
mt19937_64 stream. Replaying the 1,200 original R0 source visits is
mandatory, followed by interleaving R29, four R13 limb-width groups and
R5 within EACH host-mode iteration; generating R29 for all modes first
would break the reproducibility of downstream inputs. Exact integer
arithmetic yields GMP-style binary64 truncation, without GMP or MPFR.
GCC and Clang strict compilation and --verify reproduce all 18 original
Git blobs byte-for-byte. Registered as the twelfth explicit generator;
all original files remain versioned and unchanged.

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

B93 reproduces `floatap_directed_bounds.tsv` from original rational inputs,
computing exact floor/ceil on a power-of-two lattice and canonical dyadic
fractions at precisions 7, 32, 53, 64, 128 and 192. Both GCC 14 and
Clang 17 compile and regenerate all 90 rows with strict C++20 warnings,
producing the original 7,742-byte Git blob. The independent generator
uses no host floating-point calculations or runtime GMP/MPFR/FloatAP.
Its addition to the explicit `generators` and `generate-data` targets
does not change tests, existing fixtures or the default build.

B94 independently recreates `floatap_pi.tsv` using exact integer fixed-point
Taylor intervals for Machin's formula, 400 terms per atan at scale 2^1200.
Explicit integer interval widths and next-term remainders certify directed
384-bit endpoints and the same nearest-even 320-bit value from BOTH ends
of the interval. GCC and Clang C++20 strict builds reproduce the original
1,250-byte Git blob without GMP, MPFR, Numeric or host floating point.
The seventeenth generator is available from `generators` and `generate-data`;
all data remains versioned and `tests` does not depend on generators.

B95 reproduces `floatap_get_d.tsv`: independently round each of the 17
original exact rational sources to an intermediate p-bit FloatAP number
(precision 53, 128, 256 or 1024) with nearest-even ties, then convert
the exact dyadic result to IEEE-754 binary64 in each of four directed
modes. Signed zero, normal/subnormal, carry, midpoint-even and saturating
overflow use pure integer quotient/remainder. A separate exact Python
integer implementation reproduced all 272 outputs before the C++20
GCC/Clang strict build and 65,751-byte original Git blob comparison.
This eighteenth generator belongs to `generators`/`generate-data`, not
`tests`, and neither GMP nor original fixture values are consulted.

B96 reproduces all 2,400 exact FloatAP integer-power fixture keys from
15 distinct original signed input rationals, five binary precisions,
eight exponent choices and all four output rounding modes. The input
is rounded to p-bit nearest/even using independent integer quotient and
remainder, its exact dyadic significand raised to the original exponent,
and the final result rounded ONCE to precision p using nearest/even or
signed directed rules. All outputs are canonical dyadic rationals. GCC 14
and Clang 17 strict C++20 builds and independent execution match the
original 106,895-byte Git blob; neither Numeric nor GMP/MPFR/Boost are
used. The nineteenth generator participates in explicit `generators`
and `generate-data`, not `tests`; all TSVs remain versioned.

B97 recreates FloatAP square root expectations from the ORIGINAL 14 positive
rational sources and four shifted nearest-rounded inputs at each of five
precisions, preserving source order and deduplicating identical rounded keys.
Exact binary integer division, pairwise restoring isqrt, signed scale
bookkeeping, and strict midpoint-even comparisons generate 300 reference
rows matching the original 76,870-byte `floatap_sqrt.tsv` Git blob.
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
of the original `floatap_pow_short.tsv`; no frozen operand values or
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

B100 independently regenerates `floatap_add_sub.tsv`: the complete
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

B101 independently regenerates `floatap_exp.tsv` with correctly rounded,
rigorously certified exponential values. Recreate nine exact rational
magnitudes (1/4, 3/4, 1/2, 1, 3, 7, 2^-16, 2^-64 and 1023/1024)
and both signs, in the original order, then emit four bit widths
(53/128/256/1024) and N/D/U/Z modes for all 288 original checks.
At B=2^1400, compute outward integer interval Taylor terms k=0..165
for exp(|x|/16) and bound the positive remainder by twice term 166.
Four directed outward squares enclose exp(|x|), and exp(-|x|) is the
outward reciprocal of this positive interval. Round BOTH interval
endpoints independently in each mode and require matching exact dyadic
rational outputs before accepting any row. C++20 standard-library-only
with independently implemented 32-bit digit multiplication and exact
limb division; no GMP, MPFR, FloatAP, FloatMP, Numeric, Boost or host
floating exp. Strict GCC and Clang warnings-as-errors builds passed;
`--verify` passed on both and rejected deliberately corrupted source
data. Regenerated exactly 288 rows / 68,122 bytes, identical to Git
blob SHA-1 `92ce3c0eeda7f73f5e0385abf483f93367792891` on main.
The 24th independent generator participates only in opt-in CMake
targets `generators`, `generate-data`, `verify-data`. No original TSV,
test consumer, production source or runtime GMP oracle was changed.
Generator coverage reaches 104/115. Full CI pending for this commit.

B102 independently reproduces all 444 rows (122,768 bytes) in
`floatap_log.tsv`. Recreate the original six centered-precision visits,
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

B103 independently regenerates both `floatap_sin.tsv` and
`floatap_cos.tsv` without changing their original contents. The
standalone C++20 generator reconstructs precisely the three input
families (D91 centered offsets p/2+16, D92 local offsets p, D95 wider
exact binary dyadics), preserving the full precision-visit order,
signed source lookup keys, original nearest-even FloatAP input
quantization and first-seen key deduplication. There are exactly
336 D91, 456 D92 and 288 D95 unique source/mode rows for each
function (1,080 each). Evaluate 240 alternating integer Taylor terms
for sin and cos with outward integer interval products/quotients at
B=2^1800; separately enclose the first omitted tail term. Negative
sine uses exact oddness and cosine uses exact evenness. All four
N/D/U/Z roundings require both enclosures to agree in the target
precision; a shared cache reuses exact rounded input magnitudes.
The generator does not read TSV files to obtain either input operands
or expected answers, and uses no Numeric, GMP, MPFR, FloatMP, Boost
or host floating-point sine/cosine. GCC and Clang strict C++20 builds
and --verify pass on both original files; deliberate corruption fails.
Original Git blob SHAs are unchanged:
`floatap_sin.tsv`: `9b75d8fd4992f62e288c151e8f930dedd66a778e`,
281,134 bytes; `floatap_cos.tsv`:
`31475db715e25928b31b34d31e32735e2cbb97be`, 278,016 bytes.
Generator #26 participates only in the opt-in targets
`generators`, `generate-data`, `verify-data` and increases coverage
to 107/115, leaving eight original TSV files without generators.
No tests, production source or runtime GMP oracle were changed.

B104 adds the standard-library-only `generator_floatap_centered_exp.cpp`
covering all 568 distinct rows of `floatap_centered_exp.tsv`. Recreate
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

B105 independently regenerates all 1,340 original `floatap_atan.tsv`
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
