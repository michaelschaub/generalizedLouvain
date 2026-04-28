#!/usr/bin/env python3
"""Generate edge-list (.edj) + null-model (.edj) inputs for the Louvain
benchmark. Two families:

  - sbm: stochastic block model with k equally-sized blocks, intra-block
         edge prob `p_in`, inter-block prob `p_out`.
  - er:  Erdos-Renyi random graph with `m` distinct undirected edges.

The null-model file pairs the degree-based stationary distribution with
itself (suitable for normalised-stability runs at Markov time 1.0).
"""

import argparse
import random
import sys


def gen_sbm(n, k, p_in, p_out, seed):
    rng = random.Random(seed)
    block_size = n // k
    edges = []
    seen = set()
    for c in range(k):
        base = c * block_size
        nodes = list(range(base, base + block_size))
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                if rng.random() < p_in:
                    u, v = nodes[i], nodes[j]
                    edges.append((u, v, 1.0))
                    seen.add((u, v))
    inter_target = int(p_out * n * (n - 1) / 2)
    inter_added = 0
    tries = 0
    max_tries = max(inter_target * 10, 10000)
    while inter_added < inter_target and tries < max_tries:
        u = rng.randrange(n)
        v = rng.randrange(n)
        tries += 1
        if u == v or u // block_size == v // block_size:
            continue
        a, b = (u, v) if u < v else (v, u)
        if (a, b) in seen:
            continue
        edges.append((a, b, 1.0))
        seen.add((a, b))
        inter_added += 1
    return edges


def gen_er(n, m, seed):
    rng = random.Random(seed)
    edges = []
    seen = set()
    while len(edges) < m:
        u = rng.randrange(n)
        v = rng.randrange(n)
        if u == v:
            continue
        a, b = (u, v) if u < v else (v, u)
        if (a, b) in seen:
            continue
        seen.add((a, b))
        edges.append((a, b, 1.0))
    return edges


def write_graph(edges, n, edj_path, null_path):
    with open(edj_path, "w") as f:
        for u, v, w in edges:
            f.write(f"{u} {v} {w:.6e}\n")
    deg = [0.0] * n
    for u, v, w in edges:
        deg[u] += w
        deg[v] += w
    two_m = sum(deg)
    with open(null_path, "w") as f:
        for i in range(n):
            v = deg[i] / two_m if two_m > 0 else 0.0
            f.write(f"{v} {v}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("kind", choices=["sbm", "er"])
    ap.add_argument("n", type=int, help="number of nodes")
    ap.add_argument("--k", type=int, default=20, help="SBM blocks (default 20)")
    ap.add_argument("--p_in", type=float, default=0.3, help="SBM intra-block prob")
    ap.add_argument("--p_out", type=float, default=0.001, help="SBM inter-block prob")
    ap.add_argument("--m", type=int, help="ER edge count (default 10*n)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--edj", required=True, help="output edge-list path")
    ap.add_argument("--null", required=True, help="output null-model path")
    args = ap.parse_args()

    if args.kind == "sbm":
        edges = gen_sbm(args.n, args.k, args.p_in, args.p_out, args.seed)
    else:
        m = args.m if args.m is not None else 10 * args.n
        edges = gen_er(args.n, m, args.seed)
    write_graph(edges, args.n, args.edj, args.null)
    print(f"{args.kind} n={args.n} edges={len(edges)} -> {args.edj}", file=sys.stderr)


if __name__ == "__main__":
    main()
