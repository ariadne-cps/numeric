# Reproducible exact golden-data generators

The files `generator_big_unsigned_gcd.cpp`, `generator_big_unsigned_b76.cpp`,
`generator_exact_small.cpp`, `generator_big_unsigned_power_div.cpp`, `generator_rational_exact.cpp`,
and `oracle_exact_integer.hpp` are **standalone C++20** sources. They use only
the standard library, with an independent unsigned-limb implementation for
exact arithmetic. In particular, **do not link them to Numeric, GMP, MPFR,
FloatMP or Boost**: expected values must not come from the code being tested.

From the repository root, using a C++20 compiler:

```sh
c++ -std=c++20 -O2 -Wall -Wextra -Werror tests/data/generator_big_unsigned_gcd.cpp -o /tmp/numeric-oracle-gcd
c++ -std=c++20 -O2 -Wall -Wextra -Werror tests/data/generator_big_unsigned_b76.cpp -o /tmp/numeric-oracle-b76
c++ -std=c++20 -O2 -Wall -Wextra -Wsign-conversion -Werror tests/data/generator_exact_small.cpp -o /tmp/numeric-oracle-exact-small
c++ -std=c++20 -O2 -Wall -Wextra -Wsign-conversion -Werror tests/data/generator_big_unsigned_power_div.cpp -o /tmp/numeric-oracle-power-div
c++ -std=c++20 -O2 -Wall -Wextra -Wsign-conversion -Werror tests/data/generator_rational_exact.cpp -o /tmp/numeric-oracle-rational-exact
/tmp/numeric-oracle-gcd --verify --data-dir tests/data
/tmp/numeric-oracle-b76 --verify --data-dir tests/data
/tmp/numeric-oracle-exact-small --verify --data-dir tests/data
/tmp/numeric-oracle-power-div --verify --data-dir tests/data
/tmp/numeric-oracle-rational-exact --verify --data-dir tests/data
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
| `generator_rational_exact.cpp` | `rational_r28.tsv`, `rational_r30.tsv`, `rational_r22.tsv` | 36 signed add/sub, 1,689 normalization, 160 scalar mul/div cases |

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
