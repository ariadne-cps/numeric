# Reduced frozen oracle sampling (remove-gmp)

The original files were certified independently. Sampling is performed only on those existing frozen rows: expected results are never recomputed from the implementation under test. The original mt19937_64 stream is consumed for every original iteration; the C++ checks skip excluded cases only after constructing and mutating the random words. Both sequential TSV readers reject missing, out-of-order, and extra rows.

* BigUnsigned power: all ten width counts, modes 0/1/2/3/7 (zero, saturation, top-bit, dense, alternate saturation) and exponents 0/1/2/8/16/17 (identities, squaring, dispatch/high-bit boundaries). Every S (irregular base) and B (1024-bit limit) row is retained. 1513 -> 373 rows (4.06x).
* BigUnsigned short division: all 36 limb-width pairs (including empty dividend/divisor, equal widths, and the 4/5-limb edge); modes 0/1/2/3/4/5/11/17/23 preserve each distinct mode%6 category plus multiple dense high modes. 864 -> 324 rows (2.67x).

A tenfold reduction is not justified for these two datasets without dropping valuable operand widths or power/boundary classes. Future independent GMP/MPFR replacements should prefer source-keyed TSV fixtures and preserve algorithm thresholds, rounding modes, operand signs, equal/unequal widths, aliasing and exceptional inputs; record per-dataset counts. This step does not claim to remove all remaining runtime GMP/MPFR uses.
