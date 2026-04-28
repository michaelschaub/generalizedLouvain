#!/usr/bin/env python3
"""Read a TSV produced by run_bench.sh and emit a markdown table."""
import argparse
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv", help="path to bench_results.tsv")
    args = ap.parse_args()

    rows = []
    with open(args.tsv) as f:
        header = next(f).strip().split("\t")
        for line in f:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))

    if not rows:
        print("No rows.", file=sys.stderr)
        sys.exit(1)

    fams = {}
    for r in rows:
        fams.setdefault(r["label"], []).append(r)

    print("| family | n | m | t (master) | t (PR) | speed-up | Q match |")
    print("|---|---|---|---|---|---|---|")
    for fam in fams:
        rs = sorted(fams[fam], key=lambda r: int(r["n"]))
        for r in rs:
            n = int(r["n"]); m = int(r["m"])
            old_t = float(r["old_s"]); new_t = float(r["new_s"])
            sp = float(r["speedup"]) if r["speedup"] not in ("", "?") else None
            sp_s = f"{sp:.1f}×" if sp is not None else "—"
            qc = r["Q_check"]
            if qc == "exact":
                qc_s = "exact"
            elif qc.startswith("rel="):
                qc_s = qc.replace("rel=", "rel Δ ")
            else:
                qc_s = qc
            print(f"| {fam} | {n:,} | {m:,} | {old_t:.3f} s | {new_t:.3f} s | **{sp_s}** | {qc_s} |")


if __name__ == "__main__":
    main()
