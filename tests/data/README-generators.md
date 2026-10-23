# Reproducible exact golden-data generators

The files `generator_big_unsigned_gcd.cpp`, `generator_big_unsigned_b76.cpp`,
`generator_exact_small.cpp`, and `oracle_exact_integer.hpp` are **standalone C++20** sources. They use only
the standard library, with an independent unsigned-limb implementation for
exact arithmetic. In particular, **do not link them to Numeric, GMP, MPFR,
FloatMP or Boost**: expected values must not come from the code being tested.

From the repository root, using a C++20 compiler:

```sh
c++ -std=c++20 -O2 -Wall -Wextra -Werror tests/data/generator_big_unsigned_gcd.cpp -o /tmp/numeric-oracle-gcd
c++ -std=c++20 -O2 -Wall -Wextra -Werror tests/data/generator_big_unsigned_b76.cpp -o /tmp/numeric-oracle-b76
c++ -std=c++20 -O2 -Wall -Wextra -Wsign-conversion -Werror tests/data/generator_exact_small.cpp -o /tmp/numeric-oracle-exact-small
/tmp/numeric-oracle-gcd --verify --data-dir tests/data
/tmp/numeric-oracle-b76 --verify --data-dir tests/data
/tmp/numeric-oracle-exact-small --verify --data-dir tests/data
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
