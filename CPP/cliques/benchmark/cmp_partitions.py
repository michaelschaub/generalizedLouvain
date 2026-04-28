#!/usr/bin/env python3
"""Compare two optimal_partitions.dat files for partition equivalence.

Each file has one row per Louvain hierarchy level, with space-separated
community IDs (column = node id). Two partitions are equivalent if they
induce the same set of node groupings, regardless of label numbering.
We canonicalise each level by sorting nodes within each community and
sorting the resulting tuples lexicographically, then compare.

Usage:
    cmp_partitions.py file_a file_b
"""
import sys


def read_levels(path):
    levels = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            labels = [int(x) for x in line.split()]
            levels.append(labels)
    return levels


def canonical(labels):
    groups = {}
    for node, comm in enumerate(labels):
        groups.setdefault(comm, []).append(node)
    return sorted(tuple(sorted(g)) for g in groups.values())


def main():
    if len(sys.argv) != 3:
        print("usage: cmp_partitions.py file_a file_b", file=sys.stderr)
        sys.exit(2)
    a = read_levels(sys.argv[1])
    b = read_levels(sys.argv[2])
    if len(a) != len(b):
        print(f"DIFF: hierarchy depths differ: {len(a)} vs {len(b)}")

    last_a = canonical(a[-1])
    last_b = canonical(b[-1])
    if last_a == last_b:
        print(f"FINAL_LEVEL: identical ({len(last_a)} communities, "
              f"{sum(len(g) for g in last_a)} nodes)")
    else:
        print(f"FINAL_LEVEL: DIFFERENT ({len(last_a)} vs {len(last_b)} communities)")

    matches = 0
    common = min(len(a), len(b))
    for i in range(common):
        ca = canonical(a[i])
        cb = canonical(b[i])
        ok = ca == cb
        matches += int(ok)
        print(f"  level {i}: {'EQUAL' if ok else 'DIFF'}  comms_a={len(ca)} comms_b={len(cb)}")
    print(f"\nLevels matching: {matches}/{common}")
    return 0 if (matches == common and len(a) == len(b)) else 1


if __name__ == "__main__":
    sys.exit(main())
