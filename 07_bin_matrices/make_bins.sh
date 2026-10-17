#!/bin/bash

# ---------------------------- SLURM job parameters ---------------------------
# Edit --account before submitting. Ignored when run with 'bash'.
#SBATCH --account=<your_account>
#SBATCH --job-name=make_bins
#SBATCH --time=0:15:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%j.log

# =============================================================================
# make_bins.sh
#
# Genome-wide fixed-size bins (bedtools makewindows) used for all bin-based
# count matrices.
#
#   mm10, chr1-19 + chrX (no chrM, chrY):  5 kb -> 526,765 bins (Fig. 3)
#                                          10 kb -> 263,389 bins (Fig. 2f,
#                                          Supp. Fig. 3b)
#
# Input : chromosome sizes file (07_bin_matrices/mm10_noMY.chrom.sizes)
# Output: BED (chr, start, end), one line per bin, in chrom.sizes order
#
# Usage:
#   bash 07_bin_matrices/make_bins.sh -g 07_bin_matrices/mm10_noMY.chrom.sizes \
#        -w 5000 -o data/bins/mm10_5kb_bins.bed
#
#   -g  chromosome sizes file (chrom <tab> length)                [required]
#   -w  bin size in bp                                            [required]
#   -o  output BED file                                           [required]
#
# Software: BEDTools 2.30.0 (Saga HPC module: BEDTools/2.30.0-GCC-10.3.0)
# =============================================================================

set -o errexit
set -o nounset
set -o pipefail

if command -v module >/dev/null 2>&1; then
    module --quiet purge
    module load BEDTools/2.30.0-GCC-10.3.0
fi

GENOME=""; WIDTH=""; OUT=""
usage() { sed -n '/^# =====/,/^# =====/p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
while getopts ":g:w:o:h" opt; do
    case $opt in
        g) GENOME="$OPTARG" ;;
        w) WIDTH="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        h) usage ;;
        \?) echo "ERROR: invalid option -$OPTARG" >&2; usage ;;
        :)  echo "ERROR: option -$OPTARG requires an argument" >&2; usage ;;
    esac
done
[[ -z "$GENOME" || -z "$WIDTH" || -z "$OUT" ]] && usage
[[ -f "$GENOME" ]] || { echo "ERROR: genome file not found: $GENOME" >&2; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found in PATH" >&2; exit 1; }

mkdir -p "$(dirname "$OUT")"
bedtools makewindows -g "$GENOME" -w "$WIDTH" > "$OUT"
echo "[$(date '+%F %T')] $(wc -l < "$OUT") bins of $WIDTH bp written to $OUT"
