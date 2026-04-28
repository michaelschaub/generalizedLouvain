# Benchmark results — speed-up PR vs master

A fixed-config run of `run_bench.sh` against the `master` branch of
`michaelschaub/generalizedLouvain`. The PR delivers between ~1.5× and
~185× wall-time speed-ups on the synthetic graphs below, with bit-
identical stability values on every run.

## Environment

- Hardware: macOS laptop, Apple clang 17 (`/usr/bin/g++`)
- Compiler flags:
  - **master**: `-O3 -std=c++11 -Wall` (matches the old Makefile)
  - **PR**: `-O3 -DNDEBUG -march=native -funroll-loops -flto -std=c++17 -Wall`
- LEMON: vendored under `CPP/lemon-lib`, `lemon/config.h` generated via
  `cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 .` inside `lemon-lib/build/`
- Both binaries compiled from the same LEMON sources; only the harness
  files in `CPP/cliques/` differ between branches.
- Each row is best-of-2 wall seconds. Markov time is `1.0`, RNG seed `1`.

## Configurations

| label       | description |
|-------------|-------------|
| `SBM-20blk` | stochastic block model with 20 equal-size blocks, intra-block edge prob `p_in = 0.3`, inter-block `p_out = 0.001`; degree-based null model |
| `ER-10x`    | Erdős–Rényi random graph with `m = 10·n` distinct undirected edges; degree-based null model |

The SBM is a strong-signal input (clear planted communities, modest mixing); the ER family is a no-signal input (no planted structure, lots of noise-driven moves to converge).

## Headline numbers

| family    | n      | m           | t (master)   | t (PR)    | speed-up   | Q match |
|-----------|--------|-------------|--------------|-----------|------------|---------|
| SBM-20blk | 2,000  | 31,852      | 0.137 s      | 0.094 s   | **1.5×**   | exact   |
| SBM-20blk | 5,000  | 199,076     | 0.498 s      | 0.215 s   | **2.3×**   | exact   |
| SBM-20blk | 10,000 | 798,487     | 1.801 s      | 0.663 s   | **2.7×**   | exact   |
| SBM-20blk | 25,000 | 4,996,272   | 19.233 s     | 3.748 s   | **5.1×**   | exact   |
| SBM-20blk | 50,000 | 19,992,566  | 82.685 s     | 14.823 s  | **5.6×**   | exact   |
| ER-10x    | 2,000  | 20,000      | 0.151 s      | 0.091 s   | **1.7×**   | exact   |
| ER-10x    | 5,000  | 50,000      | 4.542 s      | 0.140 s   | **32.4×**  | exact   |
| ER-10x    | 10,000 | 100,000     | 6.468 s      | 0.232 s   | **27.9×**  | exact   |
| ER-10x    | 25,000 | 250,000     | 148.254 s    | 0.799 s   | **185.6×** | exact   |

ER-50k was deliberately omitted because the master baseline at 25k is
already 2.5 minutes; extrapolating, 50k would push the baseline well
past 30 minutes. The PR binary handles 25k ER in 0.8 s.

## How to read these numbers

The two families exercise very different parts of the code, and the
shape of the speed-up reveals which optimisation pays off where.

### SBM (planted structure): 1.5× → 5.6× as graphs grow

SBM convergence is fast — a couple of full sweeps and the planted
blocks are recovered. The work that the PR removes here is per-edge,
per-move, and per-coarsening:

- **CSR neighbour cache** trims the inner move kernel down to a tight
  pointer walk (no `IncEdgeIt + oppositeNode` overhead, no branch on
  self-loops).
- **Self-loop weight cache** removes two `lemon::findEdge(node, node)`
  calls (each O(deg)) per node move.
- **Hash-map coarsening** replaces the `lemon::findEdge` inside
  `create_reduced_graph_from_partition` (the PR's `// TODO: findEdge
  is slow!` site) with a single accumulation pass.
- **No per-sweep `compute_quality` recompute** removes O(n·k) work at
  the end of every sweep.

These are constant-factor wins per element of work — they grow with
graph size because deeper hierarchies and bigger reduced graphs
amortise more of the LEMON overhead. By 50k nodes / 20M edges the
savings reach ~5.6×.

### ER (random / no structure): saturates at 30–185×

ER is a stress-test for Louvain's outer iteration loop. Without a
planted signal, many sweeps make zero net progress yet still revisit
every node. The original code also calls `compute_quality` (an
O(n·k) walk) at the end of every such sweep just to detect the lack
of progress.

The PR replaces the "sweep all nodes per iteration" with a
deterministic FIFO queue — a node is processed only when it (or one
of its neighbours) has a chance to improve. On structureless inputs
the queue empties after a small number of moves while the baseline
keeps running entire sweeps for nothing. Hence the dramatic gap that
appears between 2k (1.7×) and 5k (32×): once the noise-floor of
"useless sweeps" overtakes the actual move work in master, the queue-
based design wins by orders of magnitude.

The peak of **185×** at 25k ER is the cleanest expression of this
effect: 148 s of baseline turns into 0.8 s.

## Correctness checks

Two different things to check:

**Stability values.** Every row has `Q match: exact` — the printed
stability strings are identical to master byte-for-byte on every
configuration in the table. No Q regression on any input.

**Partition assignments.** A separate check (cmp_partitions.py
canonicalises by sorting nodes within each community and comparing
the set-of-sets) gives a more nuanced picture:

| input  | result |
|--------|--------|
| sbm-2k, sbm-5k, sbm-10k | **bit-identical at every hierarchy level** |
| er-2k | identical |
| er-5k | hierarchy depths differ (3 vs 2 levels); same final-level partition |
| er-10k | level 0 has 4641 vs 4697 communities; same final-level partition |

The pattern matches the queue-based-fast-move caveat in PR #3:
**signal-bearing inputs (SBM) reach the same local optimum as the
sweep loop, structureless inputs (ER) can reach a different
equally-good local optimum**. Q is the same in both cases — that's
what `Q match: exact` confirms — but the intermediate clustering can
differ on inputs that admit many same-Q partitions.

In practical terms: if you run the algorithm on graphs with planted
structure (the typical use case), expect bit-identical partitions to
master. If you run it on near-uniform noise, expect equivalent Q with
possibly different clusterings at the intermediate levels.

## Reproducing this run

```bash
# 1) Build a master baseline binary somewhere outside this checkout.
git worktree add /tmp/gl-master master
make -C /tmp/gl-master/CPP/cliques
cp /tmp/gl-master/CPP/cliques/run_gen_louvain.sh /tmp/baseline_louvain

# 2) Build the PR binary in this checkout.
make -C ../             # i.e. CPP/cliques/Makefile

# 3) Run the battery.
cd CPP/cliques/benchmark
NEW_BIN=../run_gen_louvain.sh \
OLD_BIN=/tmp/baseline_louvain \
SIZES_SBM="2000 5000 10000 25000 50000" \
SIZES_ER="2000 5000 10000 25000" \
./run_bench.sh > results.tsv 2> progress.log

./format_bench.py results.tsv
```

The TSV checked in alongside this file (`results.tsv`) was produced by
exactly this command.
