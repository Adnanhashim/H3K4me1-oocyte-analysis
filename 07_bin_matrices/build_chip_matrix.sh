#!/bin/bash

# ---------------------------- SLURM job parameters ---------------------------
# Edit --account before submitting. Ignored when run with 'bash'.
#SBATCH --account=<your_account>
#SBATCH --job-name=build_chip_matrix
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/%x_%j.log

# =============================================================================
# build_chip_matrix.sh
#
# Read-count matrix (bins x samples) for ChIP-seq / picoChIP-seq / RNA-seq /
# input BAM files with bedtools multicov. One matrix per genome and bin size
# holds every sample used in the manuscript; figure scripts select the columns
# they need.
#
# Input : -s sample sheet, tab-separated, header line, columns:
#              label   bam
#           label = column name in the matrix (e.g. MII_H3K4me3_rep1)
#           bam   = BAM file (indexed; an index is created if missing)
#         -b bins BED from make_bins.sh
# Output: -o matrix, tab-separated: Bins (chr_start_end) + one column per label
#
# Usage:
#   sbatch 07_bin_matrices/build_chip_matrix.sh \
#       -s 07_bin_matrices/samples/mouse_chip_samples.tsv \
#       -b data/bins/mm10_5kb_bins.bed \
#       -o data/matrices/mouse_5kb_chip_counts.txt
#
#   -s  sample sheet                                              [required]
#   -b  bins BED                                                  [required]
#   -o  output matrix                                             [required]
#   -d  directory prepended to relative BAM paths (default: none)
#
# Counting: bedtools multicov defaults (reads overlapping a bin; reads flagged
# as duplicates are not counted). BAMs are counted in parallel (one bedtools
# multicov per BAM, up to --cpus-per-task at a time) and pasted in sheet order;
# the result is identical to one multicov call with all BAMs.
#
# Software: BEDTools 2.30.0, SAMtools 1.13
# (Saga HPC modules: BEDTools/2.30.0-GCC-10.3.0, SAMtools/1.13-GCC-10.3.0)
# =============================================================================

set -o errexit
set -o nounset
set -o pipefail

if command -v module >/dev/null 2>&1; then
    module --quiet purge
    module load BEDTools/2.30.0-GCC-10.3.0 SAMtools/1.13-GCC-10.3.0
fi

SHEET=""; BINS=""; OUT=""; DIR=""
THREADS="${SLURM_CPUS_PER_TASK:-4}"
usage() { sed -n '/^# =====/,/^# =====/p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
while getopts ":s:b:o:d:h" opt; do
    case $opt in
        s) SHEET="$OPTARG" ;;
        b) BINS="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        d) DIR="$OPTARG" ;;
        h) usage ;;
        \?) echo "ERROR: invalid option -$OPTARG" >&2; usage ;;
        :)  echo "ERROR: option -$OPTARG requires an argument" >&2; usage ;;
    esac
done
[[ -z "$SHEET" || -z "$BINS" || -z "$OUT" ]] && usage
for f in "$SHEET" "$BINS"; do [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }; done
command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found in PATH" >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH" >&2; exit 1; }

# ---------------------------- read sample sheet ------------------------------
labels=(); bams=(); missing=()
while IFS=$'\t' read -r label bam _; do
    [[ -z "$label" || "$label" == \#* || "$label" == "label" ]] && continue
    [[ -n "$DIR" && "$bam" != /* ]] && bam="${DIR%/}/$bam"
    if [[ ! -f "$bam" ]]; then missing+=("$label: $bam"); continue; fi
    labels+=("$label"); bams+=("$bam")
done < "$SHEET"
if [[ ${#missing[@]} -gt 0 ]]; then
    echo "ERROR: ${#missing[@]} BAM file(s) not found:" >&2
    printf '  %s\n' "${missing[@]}" >&2
    exit 1
fi
[[ ${#labels[@]} -gt 0 ]] || { echo "ERROR: no samples in $SHEET" >&2; exit 1; }
for bam in "${bams[@]}"; do
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
        echo "[$(date '+%F %T')] Indexing $bam"
        samtools index -@ "$THREADS" "$bam"
    fi
done
printf '[%s] Samples:\n' "$(date '+%F %T')"
for i in "${!labels[@]}"; do printf '  %-32s %s\n' "${labels[$i]}" "${bams[$i]}"; done

# ---------------------------- count ------------------------------------------
mkdir -p "$(dirname "$OUT")"
tmpdir="$(mktemp -d "${OUT}.XXXXXX")"
trap 'rm -rf "$tmpdir"' EXIT
cut -f1-3 "$BINS" > "$tmpdir/bins.bed"

count_one() {   # $1 = index, $2 = BAM -> one column of counts
    bedtools multicov -bams "$2" -bed "$tmpdir/bins.bed" | cut -f4 > "$tmpdir/$1.count"
    echo "[$(date '+%F %T')] counted $2"
}
export -f count_one; export tmpdir
for i in "${!bams[@]}"; do printf '%s\t%s\n' "$i" "${bams[$i]}"; done \
    | xargs -P "$THREADS" -L 1 bash -c 'count_one "$0" "$1"'

for i in "${!bams[@]}"; do
    [[ $(wc -l < "$tmpdir/$i.count") -eq $(wc -l < "$tmpdir/bins.bed") ]] \
        || { echo "ERROR: counting failed for ${labels[$i]}" >&2; exit 1; }
done

{
    printf "Bins"; printf "\t%s" "${labels[@]}"; printf "\n"
    paste <(awk 'BEGIN{OFS="\t"} {print $1"_"$2"_"$3}' "$tmpdir/bins.bed") \
          $(for i in "${!bams[@]}"; do printf '%s ' "$tmpdir/$i.count"; done)
} > "$OUT"

echo "[$(date '+%F %T')] Done: $OUT ($(($(wc -l < "$OUT") - 1)) bins x ${#labels[@]} samples)"
