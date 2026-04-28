# Benchmark harness

Times the new Louvain binary against a master-built baseline on synthetic
SBM and Erdős–Rényi graphs. Confirms the per-run stability matches and
prints a markdown table.

## Files

- `gen_bench_graph.py` — emits an edge-list (`.edj`) and a degree-based
  null-model file (`.edj`) for either an SBM or an ER graph.
- `run_bench.sh` — drives the battery; writes a TSV on stdout and
  progress markers on stderr.
- `format_bench.py` — turns the TSV into a markdown table.
- `cmp_partitions.py` — compares two `optimal_partitions.dat` files by
  canonicalising each hierarchy level (sorted node sets per community)
  and reporting per-level equality. Useful when Q values match but you
  want to know whether the partitions themselves agree.
- `RESULTS.md` — detailed write-up of one run.

## Running

You need two binaries: the new one (built from this branch via the
top-level `Makefile`) and a master-built baseline. Build the baseline
into a fixed path before running, e.g.

```bash
git worktree add /tmp/gl-master master
make -C /tmp/gl-master/CPP/cliques        # produces run_gen_louvain.sh there
cp /tmp/gl-master/CPP/cliques/run_gen_louvain.sh /tmp/baseline_louvain
```

Then from this folder:

```bash
NEW_BIN=../run_gen_louvain.sh \
OLD_BIN=/tmp/baseline_louvain \
SIZES_SBM="2000 5000 10000 25000" \
SIZES_ER="2000 5000 10000" \
./run_bench.sh > results.tsv 2> progress.log

./format_bench.py results.tsv
```

Defaults are `SIZES_SBM="2000 5000 10000"` and `SIZES_ER="2000 5000 10000"`.
Generated graphs are dropped under `OUTDIR` (default `/tmp`).

## Configuration knobs

- SBM defaults: 20 equal-size blocks, `p_in=0.3`, `p_out=0.001`, seed 0.
- ER defaults: `m = 10·n` undirected edges, no self-loops, seed 0.
- Louvain runs at Markov time `1.0`, RNG seed `1`, best-of-2 wall seconds.

To override SBM parameters, edit the `--k / --p_in / --p_out` arguments
at the bottom of `run_bench.sh`, or call `gen_bench_graph.py` directly.

## What the table reports

| column     | meaning                                             |
|------------|-----------------------------------------------------|
| `t (master)` | best wall time across two runs of the master binary |
| `t (PR)`   | same, for the new binary                            |
| `speed-up` | `t (master) / t (PR)`                               |
| `Q match`  | `exact` if reported stability strings are identical, otherwise relative diff |
