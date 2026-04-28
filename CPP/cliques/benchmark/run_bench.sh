#!/usr/bin/env bash
# Time the master Louvain binary against the new one across a small
# battery of synthetic graphs (SBM + ER). Reports best-of-2 wall seconds
# and confirms the stability values match.
#
# Required environment:
#   NEW_BIN  path to the new run_gen_louvain.sh (default: ../run_gen_louvain.sh)
#   OLD_BIN  path to a master-built run_gen_louvain.sh (no default;
#            build master separately, e.g. into /tmp/baseline_louvain)
#
# Optional:
#   SIZES_SBM  space-separated SBM node counts (default "2000 5000 10000")
#   SIZES_ER   space-separated ER  node counts (default "2000 5000 10000")
#   OUTDIR     where to drop generated graphs (default /tmp)
#
# Output: TSV on stdout, progress markers on stderr.

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEW_BIN="${NEW_BIN:-$HERE/../run_gen_louvain.sh}"
OLD_BIN="${OLD_BIN:-}"
GEN="$HERE/gen_bench_graph.py"
SIZES_SBM="${SIZES_SBM:-2000 5000 10000}"
SIZES_ER="${SIZES_ER:-2000 5000 10000}"
OUTDIR="${OUTDIR:-/tmp}"

if [[ -z "$OLD_BIN" ]]; then
    echo "set OLD_BIN to a master-built run_gen_louvain.sh" >&2
    exit 1
fi
for f in "$NEW_BIN" "$OLD_BIN" "$GEN"; do
    if [[ ! -e "$f" ]]; then
        echo "missing: $f" >&2
        exit 1
    fi
done

best_of_two() {
    local bin="$1"; shift
    local edj="$1"; shift
    local null="$1"; shift
    local mt="$1"; shift
    local seed="$1"; shift
    local best=999999
    local stab=""
    for _ in 1 2; do
        local t0 t1 elapsed cmp_result out
        t0=$(python3 -c 'import time;print(time.perf_counter())')
        out=$("$bin" "$edj" "$null" "$mt" "$seed" 2>/dev/null | tail -5)
        t1=$(python3 -c 'import time;print(time.perf_counter())')
        elapsed=$(python3 -c "print($t1-$t0)")
        cmp_result=$(python3 -c "print(1 if $elapsed < $best else 0)")
        if [[ "$cmp_result" == "1" ]]; then
            best="$elapsed"
        fi
        stab=$(echo "$out" | grep -oE 'Stability:  [-0-9.eE+]+' | tail -1 | awk '{print $2}')
    done
    printf "%.3f\t%s" "$best" "$stab"
}

bench_row() {
    local label="$1"; shift
    local edj="$1"; shift
    local null="$1"; shift
    local mt="$1"; shift
    local seed="$1"; shift
    local n="$1"; shift
    local m="$1"; shift

    echo ">>> $label  (n=$n m=$m mt=$mt seed=$seed)" >&2
    local old_row new_row old_t old_q new_t new_q speedup q_match
    old_row=$(best_of_two "$OLD_BIN" "$edj" "$null" "$mt" "$seed")
    new_row=$(best_of_two "$NEW_BIN" "$edj" "$null" "$mt" "$seed")
    old_t=$(echo "$old_row" | cut -f1); old_q=$(echo "$old_row" | cut -f2)
    new_t=$(echo "$new_row" | cut -f1); new_q=$(echo "$new_row" | cut -f2)
    speedup=$(python3 -c "print(f'{$old_t/$new_t:.2f}')" 2>/dev/null || echo "?")
    if [[ "$old_q" == "$new_q" ]]; then
        q_match="exact"
    else
        q_match=$(python3 -c "
old=float('$old_q'); new=float('$new_q')
rel=abs(old-new)/max(abs(old),abs(new),1e-300)
print(f'rel={rel:.2e}')
" 2>/dev/null || echo "diff")
    fi
    printf "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$label" "$n" "$m" "$mt" "$old_t" "$new_t" "$speedup" "$q_match" "$old_q"
}

printf "label\tn\tm\tmt\told_s\tnew_s\tspeedup\tQ_check\tQ_old\n"

for n in $SIZES_SBM; do
    edj="$OUTDIR/bench_sbm_${n}.edj"
    nullf="$OUTDIR/bench_sbm_${n}_null.edj"
    python3 "$GEN" sbm "$n" --k 20 --p_in 0.3 --p_out 0.001 --seed 0 --edj "$edj" --null "$nullf" 2>/dev/null
    m=$(wc -l < "$edj" | tr -d ' ')
    bench_row "SBM-20blk" "$edj" "$nullf" 1.0 1 "$n" "$m"
done

for n in $SIZES_ER; do
    edj="$OUTDIR/bench_er_${n}.edj"
    nullf="$OUTDIR/bench_er_${n}_null.edj"
    m=$((10 * n))
    python3 "$GEN" er "$n" --m "$m" --seed 0 --edj "$edj" --null "$nullf" 2>/dev/null
    bench_row "ER-10x" "$edj" "$nullf" 1.0 1 "$n" "$m"
done
