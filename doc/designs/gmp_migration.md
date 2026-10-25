# GMP numerical-oracle migration — BigUnsigned

Baseline: `main` at `60306926d7093bcf8ebee0867ec5ffe2600401c9` (A10), 25 October 2026.

## Scope and counting

Count **one test visit** per original input/mode combination whose expected
numeric output is supplied by GMP, regardless of the number of assertions.
Exclude GMP interoperability tests (`test_gmp_bridge()`), code/branch
coverage and migrations preceding A1. This inventory applies only to
`tests/test_big_unsigned.cpp`. The denominator remains **76,323** for
progress reporting, barring a documented correction to the audit.

## Current state — A19

| Method | Migrated | Remaining |
| --- | ---: | ---: |
| `test_division()` | 4,978 | 11,536 |
| `test_views_and_scratch()` | 0 | 11,600 |
| `test_multiply_kernel_paths()` | 0 | 11,432 |
| `test_randomized()` | 10,000 | 0 |
| `test_gcd_decimal()` | 0 | 3,793 |
| `test_division_limb()` | 12,944 | 0 |
| `test_exact_division_limb()` | 6,200 | 0 |
| `test_multiplication()` | 2,112 | 0 |
| `test_product_accumulation()` | 1,728 | 0 |
| **Total** | **37,962** | **38,361** |

Progress: **37,962 / 76,323 = 49.739%**.
A10 baseline was 27,626 migrated and 48,697 remaining (36.196%).

## Remaining inventory

- `test_division()`: generalized random 6,000; B28, B36, B29, B45,
  normalization and remaining special division cases 5,536.
- `test_views_and_scratch()`: 5,000 buffer arithmetic plus 6,600
  R38 shifted operation visits.
- `test_multiply_kernel_paths()`: 11,432, comprising B52/B53,
  general kernels, B8 and Karatsuba.
- `test_gcd_decimal()`: 2,284 decimal parsing visits and 1,509
  decimal formatting visits.

## Commit milestones

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

A11 maintains the eight dividend widths, fourteen divisors, three
original modes, all PRNG draws, and separately validates exact
`divrem`, `quot` and `rem` outputs and original scratch size
invariants against a new independent source-keyed fixture.
`tests/data/generator_big_unsigned_power_div.cpp` regenerates
`big_unsigned_b48_multiprecision.tsv` without GMP, Numeric, MPFR
or Boost. The fixture has 336 rows, 305,944 bytes, Git blob
`11042c7f930f66082d8a89d6670376c03786938d`.
Other GMP division blocks remain unchanged.

## Post-commit reporting

After each commit report newly migrated test visits, cumulative
migrated/remaining numbers, percentage using the 76,323 baseline,
method-level remaining cases, and observed CI status separately from
fixture checks. Do not count multiple assertions on a single input as
multiple visits or infer `verify-data` success from ordinary CI.

## A12 implementation

A12 replaces GMP only for the first 1,250 seeded visits of
`test_randomized()`, using `big_unsigned_a12_randomized_first.tsv`.
The remaining 8,750 cases deliberately retain GMP. All 10,000 visits
continue to use the same single MT19937_64 engine and consume the same
shift draw each time; comparison, addition, conditional subtraction,
and both shifts have exact independent expectations. The C++20
standalone generator uses independent Nat arithmetic and the exact
fixture was cross-checked using a separately compiled GCC/Clang
random-input witness plus Python arbitrary-precision arithmetic.
The fixture is 1,202,560 bytes, Git blob
`487f06f3921d72d2de9caee06bf9c23957ee6012`.

## A13 — randomized arithmetic, second 1,250 visits (25 October 2026)

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
Python built-in arbitrary-precision arithmetic separately generates
all expected values, confirmed byte-for-byte against the Git fixture
blob SHA `ae3a7c98ef3ac9dedea8e44d5a2d410daccb342e`
(1,217,436 bytes). Generator and consumer share no production arithmetic.
No production code or prior fixture was changed. Generator-wide
`verify-data` is not implied by normal CI.

## A14 — randomized visits 2,500–3,749 (25 October 2026)

Migrate 1,250 additional visits using an exact source-keyed TSV,
retaining the single original random engine, comparison/addition,
conditional subtraction and both shift operations. GMP still serves
visits 3,750–9,999. The standalone C++20 generator replays all 2,500
preceding visits and records results computed using independent Nat
arithmetic. The independently computed fixture has 1,250 records,
1,191,197 bytes and Git blob
`bdb0b5f9235005b2bc31d18a7e3b21c2af474abe`.

## A15 — randomized visits 3,750–4,999 (25 October 2026)

Migrate 1,250 additional visits with
`big_unsigned_a15_randomized_fourth.tsv`. The consumer still executes the
single original `mt19937_64(0xBA17B00D)` stream for all 10,000 visits and
retains comparison, addition, conditional subtraction and both shifts.
GMP remains the numerical oracle only for visits 5,000–9,999.

The standalone C++20 generator replays visits 0–3,749 before emitting A15.
A separate Python arbitrary-precision implementation reproduces the complete
A14 fixture exactly, including Git blob
`bdb0b5f9235005b2bc31d18a7e3b21c2af474abe`, then independently computes
A15. The A15 fixture has 1,250 records, 1,216,271 bytes and Git blob
`e23f97770356c28d0b76aef3c7d9d5c7dee8a0c5`. No production code or prior
fixture is changed. Generator-wide `verify-data` is not implied by normal CI.

## A16 — randomized visits 5,000–6,249 (25 October 2026)

Migrate 1,250 additional visits using
`big_unsigned_a16_randomized_fifth.tsv`. The consumer retains the original
single `mt19937_64(0xBA17B00D)` engine across all 10,000 visits, including
all operand draws, conditional full-array transformations and the shift draw.
Comparison, addition, conditional subtraction and both shifts continue to be
checked against exact external expectations. GMP remains the numerical oracle
for visits 6,250–9,999.

The standalone C++20 generator replays visits 0–4,999 before emitting A16.
A separate arbitrary-precision implementation reproduces A15 exactly, including
Git blob `e23f97770356c28d0b76aef3c7d9d5c7dee8a0c5`, before independently
computing A16. The A16 fixture contains 1,250 records, 1,204,292
bytes and Git blob `bbb0fc03a3dae88154205334fafdff5fb04f9619`. No production code or previous fixture
is changed. Generator-wide `verify-data` remains a separate check from CI.

## A17 — randomized visits 6,250–7,499 (25 October 2026)

Migrate 1,250 additional visits using
`big_unsigned_a17_randomized_sixth.tsv`. The consumer preserves the original
single `mt19937_64(0xBA17B00D)` engine across all 10,000 visits, including
the complete operand draw stream, conditional full-array transformations and
shift draw. Comparison, addition, conditional subtraction and both shifts are
checked against exact external expectations. GMP remains the numerical oracle
for visits 7,500–9,999.

The standalone C++20 generator replays visits 0–6,249 before emitting A17.
A separate arbitrary-precision implementation reproduces A16 exactly,
including Git blob `bbb0fc03a3dae88154205334fafdff5fb04f9619`, before
independently computing A17. The A17 fixture contains 1,250 records,
1,224,250 bytes and Git blob `2b595bab5b7a83fc2c9a27c188d2fba204b6c30e`. No production code or prior
fixture is changed. Generator-wide `verify-data` remains separate from CI.

## A18 — randomized visits 7,500–8,749 (25 October 2026)

Migrate 1,250 additional visits using
`big_unsigned_a18_randomized_seventh.tsv`. The consumer preserves the
original single `mt19937_64(0xBA17B00D)` engine across all 10,000 visits,
including the complete operand draw stream, conditional full-array
transformations and shift draw. Comparison, addition, conditional subtraction
and both shifts are checked against exact external expectations. GMP remains
the numerical oracle for visits 8,750–9,999.

The standalone C++20 generator replays visits 0–7,499 before emitting A18.
A separate arbitrary-precision implementation reproduces A17 exactly,
including Git blob `2b595bab5b7a83fc2c9a27c188d2fba204b6c30e`, before
independently computing A18. The A18 fixture contains 1,250 records,
1,212,953 bytes and Git blob `b1c7416dd60cc79511d545a706df56b18b70bf74`. No production code or prior
fixture is changed. Generator-wide `verify-data` remains separate from CI.

## A19 — randomized visits 8,750–9,999 (25 October 2026)

Migrate the final 1,250 `test_randomized()` visits using
`big_unsigned_a19_randomized_eighth.tsv`. The original single
`mt19937_64(0xBA17B00D)` engine, operand draws, conditional full-array
transformations and shift draws are unchanged. Comparison, addition,
conditional subtraction and both shifts now use exact external expectations
for all 10,000 visits; the runtime numerical GMP branch is removed from this
method.

The standalone C++20 generator replays visits 0–8,749 before emitting A19.
A separate arbitrary-precision implementation reproduces A18 exactly,
including Git blob `b1c7416dd60cc79511d545a706df56b18b70bf74`, before
independently computing A19. The A19 fixture contains 1,250 records,
1,202,168 bytes and Git blob `15f9d3656bf5fe77ac2b97e0ba6faee3b0883ea3`. No production code or prior
fixture is changed. Generator-wide `verify-data` remains separate from CI.
