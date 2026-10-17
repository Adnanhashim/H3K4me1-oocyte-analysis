#!/bin/bash

# ---------------------------- SLURM job parameters ---------------------------
# Edit --account before submitting. Ignored when run with 'bash'.
#SBATCH --account=<your_account>
#SBATCH --job-name=build_methylation_matrix
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/%x_%j.log

# =============================================================================
# build_methylation_matrix.sh
#
# DNA methylation matrix (bins x samples): mean methylation per bin from
# BigWig files (deepTools multiBigwigSummary BED-file). Same bins as the
# read-count matrix, so figure scripts can combine both by bin name.
#
# Input : -s sample sheet, tab-separated, header line, columns:
#              label   bigwig
#         BigWigs may be in % (0-100) or fractions (0-1): each column is checked
#         and divided by 100 if its maximum is above 1, so the matrix is 0-1.
#         -b bins BED from make_bins.sh
# Output: -o matrix, tab-separated: Bins (chr_start_end) + one column per label.
#         Bins without CpG coverage are NA (use -z to write 0 instead).
#
# Usage:
#   sbatch 07_bin_matrices/build_methylation_matrix.sh \
#       -s 07_bin_matrices/samples/mouse_methylation_samples.tsv \
#       -b data/bins/mm10_5kb_bins.bed \
#       -o data/matrices/mouse_5kb_MII_5mC.txt
#
#   -s  sample sheet                                              [required]
#   -b  bins BED                                                  [required]
#   -o  output matrix                                             [required]
#   -d  directory prepended to relative BigWig paths (default: none)
#   -z  write 0 instead of NA for bins without data
#
# Software: deepTools 3.5.0 (Saga HPC module: deepTools/3.5.0-foss-2021a)
# =============================================================================

set -o errexit
set -o nounset
set -o pipefail

if command -v module >/dev/null 2>&1; then
    module --quiet purge
    module load deepTools/3.5.0-foss-2021a
fi

SHEET=""; BINS=""; OUT=""; DIR=""; MISSING="NA"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
usage() { sed -n '/^# =====/,/^# =====/p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
while getopts ":s:b:o:d:zh" opt; do
    case $opt in
        s) SHEET="$OPTARG" ;;
        b) BINS="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        d) DIR="$OPTARG" ;;
        z) MISSING="0" ;;
        h) usage ;;
        \?) echo "ERROR: invalid option -$OPTARG" >&2; usage ;;
        :)  echo "ERROR: option -$OPTARG requires an argument" >&2; usage ;;
    esac
done
[[ -z "$SHEET" || -z "$BINS" || -z "$OUT" ]] && usage
for f in "$SHEET" "$BINS"; do [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }; done
command -v multiBigwigSummary >/dev/null 2>&1 || { echo "ERROR: multiBigwigSummary (deepTools) not found" >&2; exit 1; }

labels=(); bws=()
while IFS=$'\t' read -r label bw _; do
    [[ -z "$label" || "$label" == \#* || "$label" == "label" ]] && continue
    [[ -n "$DIR" && "$bw" != /* ]] && bw="${DIR%/}/$bw"
    [[ -f "$bw" ]] || { echo "ERROR: BigWig not found for $label: $bw" >&2; exit 1; }
    labels+=("$label"); bws+=("$bw")
done < "$SHEET"
[[ ${#labels[@]} -gt 0 ]] || { echo "ERROR: no samples in $SHEET" >&2; exit 1; }

echo "[$(date '+%F %T')] ${#labels[@]} samples, $(wc -l < "$BINS") bins"
echo "[$(date '+%F %T')] deepTools $(multiBigwigSummary --version 2>&1)"

mkdir -p "$(dirname "$OUT")"
tmpdir="$(mktemp -d "${OUT}.XXXXXX")"
trap 'rm -rf "$tmpdir"' EXIT

# bins as a real file (worker processes cannot read a <(...) pipe)
cut -f1-3 "$BINS" > "$tmpdir/bins.bed"
multiBigwigSummary BED-file -p "$THREADS" \
    -b "${bws[@]}" -l "${labels[@]}" \
    --BED "$tmpdir/bins.bed" \
    -o "$tmpdir/summary.npz" --outRawCounts "$tmpdir/summary.tab"

# Header from the sample sheet; rows: bin name + values (nan -> NA or 0)
{
    printf "Bins"; printf "\t%s" "${labels[@]}"; printf "\n"
    grep -v '^#' "$tmpdir/summary.tab" > "$tmpdir/values.tab"
    # per column: range check; % (max > 1) -> fraction
    awk -v labels="${labels[*]}" 'BEGIN { n = split(labels, lab, " ") }
        { for (i = 4; i <= NF; i++) if (tolower($i) != "nan") {
              v = $i + 0; j = i - 3
              if (!(j in mx) || v > mx[j]) mx[j] = v
              if (!(j in mn) || v < mn[j]) mn[j] = v } }
        END { for (j = 1; j <= n; j++) {
                  d = (mx[j] > 1) ? 100 : 1
                  printf "%s\t%s\n", j, d
                  printf "  %-28s min %s  max %s  -> %s\n", lab[j], mn[j], mx[j],
                         (d == 100 ? "% (divided by 100)" : "fraction (unchanged)") > "/dev/stderr" } }' \
        "$tmpdir/values.tab" > "$tmpdir/divisors.tab"
    awk -v miss="$MISSING" 'BEGIN { OFS = "\t" }
        NR == FNR { d[$1] = $2; next }
        { line = $1"_"$2"_"$3
          for (i = 4; i <= NF; i++) line = line OFS (tolower($i) == "nan" ? miss : $i / d[i - 3])
          print line }' "$tmpdir/divisors.tab" "$tmpdir/values.tab"
} > "$OUT"

echo "[$(date '+%F %T')] Done: $OUT ($(($(wc -l < "$OUT") - 1)) bins x ${#labels[@]} samples)"
