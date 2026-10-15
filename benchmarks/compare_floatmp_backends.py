#!/usr/bin/env python3
"""Compare FloatVP to the MPFR reference, including aggregate slowdown ratios.

Every row is a median nanoseconds/operation from benchmark_floatmp_backends.
The geometric mean gives each (operation, precision) equal influence;
the total-time ratio weights measurements by the observed operation cost.
"""
import argparse
import csv
import math
import sys
from collections import defaultdict


def summarize(pairs):
    """Return geomean, total-time ratio and worst ratio for (operation,bits,mpfr,vp)."""
    ratios = [(vp / mpfr, operation, bits) for operation, bits, mpfr, vp in pairs]
    geomean = math.exp(math.fsum(math.log(ratio) for ratio, _, _ in ratios) / len(ratios))
    total_ratio = math.fsum(vp for _, _, _, vp in pairs) / math.fsum(
        mpfr for _, _, mpfr, _ in pairs
    )
    worst_ratio, worst_operation, worst_bits = max(ratios)
    return geomean, total_ratio, worst_ratio, worst_operation, worst_bits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_file", help="dual-backend CSV measured on one machine")
    parser.add_argument("--summary-only", action="store_true",
                        help="show just per-precision and overall slowdown, not every operation")
    parser.add_argument("--max-slowdown", type=float,
                        help="fail if ANY individual FloatVP/MPFR ratio exceeds this threshold")
    args = parser.parse_args()
    if args.max_slowdown is not None and (
        not math.isfinite(args.max_slowdown) or args.max_slowdown <= 0
    ):
        parser.error("--max-slowdown must be positive and finite")

    rows = {}
    with open(args.csv_file, newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            try:
                key = (row["backend"], row["operation"], int(row["bits"]))
                ns = float(row["median_ns"])
            except (KeyError, TypeError, ValueError) as error:
                parser.error("invalid benchmark CSV row: " + str(error))
            if key in rows:
                parser.error("duplicate benchmark result: " + repr(key))
            if not math.isfinite(ns) or ns <= 0:
                parser.error("benchmark timing must be positive and finite: " + repr(key))
            rows[key] = (ns, row["iterations"], row["samples"])

    if not any(key[0] == "floatvp" for key in rows):
        parser.error("missing FloatVP rows: enable BUILD_FLOATVP_COMPARISON "
                     "once FloatVP is implemented")
    if not args.summary_only:
        print("operation,bits,mpfr_ns,floatvp_ns,slowdown_floatvp_over_mpfr")

    pairs_by_precision = defaultdict(list)
    violations = 0
    for operation, bits in sorted({(operation, bits) for _, operation, bits in rows}):
        reference = rows.get(("mpfr", operation, bits))
        candidate = rows.get(("floatvp", operation, bits))
        if reference is None or candidate is None:
            parser.error("unpaired result: " + operation + " / " + str(bits))
        if reference[1:] != candidate[1:]:
            parser.error("mismatched iterations/samples: " + operation + " / " + str(bits))

        mpfr_ns, vp_ns = reference[0], candidate[0]
        ratio = vp_ns / mpfr_ns
        if not args.summary_only:
            print(f"{operation},{bits},{mpfr_ns:.6f},{vp_ns:.6f},{ratio:.3f}")
        pairs_by_precision[bits].append((operation, bits, mpfr_ns, vp_ns))
        if args.max_slowdown is not None and ratio > args.max_slowdown:
            violations += 1

    print()
    print("# FloatVP / FloatMP aggregate slowdown (x; >1 means slower)")
    print("precision,operations,geomean_x,total_time_x,worst_x,worst_operation,worst_bits")
    overall = []
    for bits, pairs in sorted(pairs_by_precision.items()):
        gm, total, worst, operation, worst_bits = summarize(pairs)
        print(f"{bits},{len(pairs)},{gm:.3f},{total:.3f},{worst:.3f},{operation},{worst_bits}")
        overall.extend(pairs)
    gm, total, worst, operation, worst_bits = summarize(overall)
    print(f"ALL,{len(overall)},{gm:.3f},{total:.3f},{worst:.3f},{operation},{worst_bits}")

    if violations:
        print(f"{violations} per-operation regressions exceed threshold", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
