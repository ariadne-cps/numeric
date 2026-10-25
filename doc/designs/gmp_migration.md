# GMP numerical-oracle migration — BigUnsigned

Baseline: `main` at `60306926d7093bcf8ebee0867ec5ffe2600401c9` (A10), 25 October 2026.

## Scope and counting

Count **one test visit** per original input/mode combination whose expected
numeric output is supplied by GMP, regardless of the number of assertions.
Exclude GMP interoperability tests (`test_gmp_bridge()`), code/branch
coverage and migrations preceding A1. This inventory applies only to
`tests/test_big_unsigned.cpp`. The denominator remains **76,323** for
progress reporting, barring a documented correction to the audit.

## Current state — A11

| Method | Migrated | Remaining |
| --- | ---: | ---: |
| `test_division()` | 4,978 | 11,536 |
| `test_views_and_scratch()` | 0 | 11,600 |
| `test_multiply_kernel_paths()` | 0 | 11,432 |
| `test_randomized()` | 0 | 10,000 |
| `test_gcd_decimal()` | 0 | 3,793 |
| `test_division_limb()` | 12,944 | 0 |
| `test_exact_division_limb()` | 6,200 | 0 |
| `test_multiplication()` | 2,112 | 0 |
| `test_product_accumulation()` | 1,728 | 0 |
| **Total** | **27,962** | **48,361** |

Progress: **27,962 / 76,323 = 36.635%**.
A10 baseline was 27,626 migrated and 48,697 remaining (36.196%).

## Remaining inventory

- `test_division()`: generalized random 6,000; B28, B36, B29, B45,
  normalization and remaining special division cases 5,536.
- `test_views_and_scratch()`: 5,000 buffer arithmetic plus 6,600
  R38 shifted operation visits.
- `test_multiply_kernel_paths()`: 11,432, comprising B52/B53,
  general kernels, B8 and Karatsuba.
- `test_randomized()`: 10,000 seeded arithmetic/shift visits.
- `test_gcd_decimal()`: 2,284 decimal parsing visits and 1,509
  decimal formatting visits.

## Commit milestones

| Phase | Newly migrated | Cumulative migrated |
| --- | ---: | ---: |
| A1–A10 | 27,626 | 27,626 |
| A11 — multiprecision B48 | 336 | 27,962 |

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
