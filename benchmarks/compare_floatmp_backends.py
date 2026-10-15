#!/usr/bin/env python3
"""Compare MPFR and FloatVP benchmark CSV results from the same executable."""
import argparse
import csv
import sys

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_file", help="dual-backend CSV measured on one machine")
    parser.add_argument("--max-slowdown", type=float,
                        help="fail when any FloatVP/MPFR ratio exceeds this threshold")
    args = parser.parse_args()
    rows = {}
    with open(args.csv_file, newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            key = (row["backend"], row["operation"], int(row["bits"]))
            if key in rows:
                parser.error("duplicate result: "+repr(key))
            rows[key] = (float(row["median_ns"]), row["iterations"], row["samples"])
    if not any(key[0]=="floatvp" for key in rows):
        parser.error("missing FloatVP rows: enable BUILD_FLOATVP_COMPARISON once FloatVP is implemented")
    print("operation,bits,mpfr_ns,floatvp_ns,slowdown_floatvp_over_mpfr")
    violations = 0
    for operation,bits in sorted({(op,bits) for _,op,bits in rows}):
        a = rows.get(("mpfr",operation,bits))
        b = rows.get(("floatvp",operation,bits))
        if a is None or b is None:
            parser.error("unpaired result: "+operation+" / "+str(bits))
        if a[1:] != b[1:]:
            parser.error("mismatched iteration/sample counts: "+operation)
        if a[0]<=0:
            parser.error("invalid reference timing for "+operation)
        ratio = b[0]/a[0]
        print(f"{operation},{bits},{a[0]:.3f},{b[0]:.3f},{ratio:.3f}")
        if args.max_slowdown is not None and ratio > args.max_slowdown:
            violations += 1
    if violations:
        print(f"{violations} performance regressions exceed threshold",file=sys.stderr)
        return 1
    return 0

if __name__=="__main__":
    sys.exit(main())
