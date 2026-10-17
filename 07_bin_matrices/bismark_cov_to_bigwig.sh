#!/bin/bash

# ---------------------------- SLURM job parameters ---------------------------
# Edit --account before submitting. Ignored when run with 'bash'.
#SBATCH --account=<your_account>
#SBATCH --job-name=bismark_cov_to_bigwig
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%j.log

# =============================================================================
# bismark_cov_to_bigwig.sh
#
# Bismark coverage file(s) (*.cov / *.cov.txt; columns: chr, start, end,
# % methylation, count methylated, count unmethylated) -> one BigWig of CpG
# methylation (fraction 0-1), for build_methylation_matrix.sh.
# Several files (replicates) are merged by summing the methylated and
# unmethylated counts per CpG. Chromosome names without "chr" get the prefix;
# only chromosomes in the chrom.sizes file are kept.
#
# Usage:
#   sbatch 07_bin_matrices/bismark_cov_to_bigwig.sh \
#       -g 07_bin_matrices/mm10_noMY.chrom.sizes -o P10_GO_5mC.bw \
#       GSM2690945_d10_GO_input_BS_1.cov.txt GSM2690946_d10_GO_input_BS_2.cov.txt
#
#   -g  chrom.sizes                                               [required]
#   -o  output BigWig                                             [required]
#
# Software: bedGraphToBigWig (UCSC Kent tools) if available, otherwise Python
# pyBigWig from the deepTools module (Saga: deepTools/3.5.0-foss-2021a).
# =============================================================================

set -o errexit
set -o nounset
set -o pipefail

if command -v module >/dev/null 2>&1; then
    module --quiet purge
    module load Kent_tools 2>/dev/null || true
    command -v bedGraphToBigWig >/dev/null 2>&1 || module load deepTools/3.5.0-foss-2021a 2>/dev/null || true
fi

SIZES=""; OUT=""
usage() { sed -n '/^# =====/,/^# =====/p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
while getopts ":g:o:h" opt; do
    case $opt in
        g) SIZES="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        h) usage ;;
        \?) echo "ERROR: invalid option -$OPTARG" >&2; usage ;;
        :)  echo "ERROR: option -$OPTARG requires an argument" >&2; usage ;;
    esac
done
shift $((OPTIND - 1))
[[ -z "$SIZES" || -z "$OUT" || $# -eq 0 ]] && usage
for f in "$SIZES" "$@"; do [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }; done
if command -v bedGraphToBigWig >/dev/null 2>&1; then
    WRITER="bedGraphToBigWig"
elif python3 -c "import pyBigWig" 2>/dev/null; then
    WRITER="pyBigWig"
else
    echo "ERROR: neither bedGraphToBigWig nor Python pyBigWig found (load the deepTools module)" >&2; exit 1
fi
echo "[$(date '+%F %T')] BigWig writer: $WRITER"

echo "[$(date '+%F %T')] Merging $# file(s): $*"
tmp="$(mktemp "${OUT}.XXXXXX.bedGraph")"
trap 'rm -f "$tmp"' EXIT

# Sum counts per CpG over all files; fraction = M / (M + U); bedGraph is 0-based
cat "$@" | awk -v sizes="$SIZES" 'BEGIN { OFS = "\t"
        while ((getline line < sizes) > 0) { split(line, a, "\t"); keep[a[1]] = 1 } }
    { c = ($1 ~ /^chr/) ? $1 : "chr" $1
      if (!(c in keep)) next
      k = c "\t" $2; m[k] += $5; u[k] += $6 }
    END { for (k in m) if (m[k] + u[k] > 0) {
              split(k, p, "\t"); print p[1], p[2] - 1, p[2], m[k] / (m[k] + u[k]) } }' \
    | LC_ALL=C sort -k1,1 -k2,2n > "$tmp"

echo "[$(date '+%F %T')] $(wc -l < "$tmp") CpGs"
if [[ "$WRITER" == "bedGraphToBigWig" ]]; then
    bedGraphToBigWig "$tmp" "$SIZES" "$OUT"
else
    python3 - "$tmp" "$SIZES" "$OUT" <<'PY'
import sys, pyBigWig
bg, sizes, out = sys.argv[1:4]
chroms = sorted((l.split()[0], int(l.split()[1])) for l in open(sizes) if l.strip())
bw = pyBigWig.open(out, "w")
bw.addHeader(chroms)               # same order as the bedGraph (sort -k1,1, C locale)
buf = {"c": [], "s": [], "e": [], "v": []}
def flush():
    if buf["c"]:
        bw.addEntries(buf["c"], buf["s"], ends=buf["e"], values=buf["v"])
        for k in buf: buf[k] = []
for line in open(bg):
    c, s, e, v = line.split()
    buf["c"].append(c); buf["s"].append(int(s)); buf["e"].append(int(e)); buf["v"].append(float(v))
    if len(buf["c"]) >= 100000: flush()
flush()
bw.close()
PY
fi
echo "[$(date '+%F %T')] Done: $OUT"
