#!/bin/bash

# ---------------------------- SLURM job parameters ---------------------------
# Edit --account before submitting. Ignored when run with 'bash'.
#SBATCH --account=<your_account>
#SBATCH --job-name=genome_coverage
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%j.log

# =============================================================================
# genome_coverage.sh
#
# Percentage of the genome covered by broad peaks (bedtools coverage).
# Used for Fig. 1c (H3K4me1 genome coverage, P10 GOs vs ENCODE tissues).
# Plot the result with plot_genome_coverage.R.
#
# For every peak file, bedtools coverage is run against a BED file spanning
# each chromosome (chrom, 0, size). Genome coverage (%) is the mean of the
# per-chromosome covered fraction (bedtools coverage column 7) x 100.
#
# Input : -g  genome file: chrom.sizes (2 columns) or BED (chrom, 0, size)
#         peak files (.broadPeak / .bed) as positional arguments
# Output (in <outdir>):
#   genome_coverage.txt                       one row per peak file (sample, coverage %)
#   per_chrom/<sample>_coverage_per_chrom.tsv raw bedtools coverage output
#
# Usage (SLURM; submit from the repository root, the logs/ folder must exist):
#   mkdir -p logs
#   sbatch genome_coverage.sh -g <mm10.chrom.sizes> -o <outdir> peaks/*.broadPeak
# Usage (interactive / local):
#   bash genome_coverage.sh -g <mm10.chrom.sizes> -o <outdir> peaks/*.broadPeak
#
#   -g  Genome sizes file (chrom.sizes or BED)                     [required]
#   -o  Output directory (default: ./genome_coverage)
#
# Example:
#   bash genome_coverage.sh -g mm10.chrom.sizes -o results/genome_coverage/ peaks/*.broadPeak
#
# Software: BEDTools 2.30.0
# (Saga HPC module: BEDTools/2.30.0-GCC-10.3.0)
# =============================================================================

# Enable error checking
set -o errexit
set -o nounset
set -o pipefail

# ---------------------------- load modules -----------------------------------
# Module names are those of the Saga HPC cluster; adjust to your system.
# Skipped automatically when the 'module' command is not available.
if command -v module >/dev/null 2>&1; then
    module --quiet purge
    module load BEDTools/2.30.0-GCC-10.3.0
fi

# ---------------------------- defaults ---------------------------------------
OUTDIR="./genome_coverage"
GENOME=""

usage() {
    sed -n '/^# =====/,/^# =====/p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
}

# ---------------------------- arguments --------------------------------------
while getopts ":g:o:h" opt; do
    case $opt in
        g) GENOME="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        h) usage ;;
        \?) echo "ERROR: invalid option -$OPTARG" >&2; usage ;;
        :)  echo "ERROR: option -$OPTARG requires an argument" >&2; usage ;;
    esac
done
shift $((OPTIND - 1))

[[ -z "$GENOME" ]] && { echo "ERROR: -g <genome sizes file> is required" >&2; usage; }
[[ -f "$GENOME" ]] || { echo "ERROR: genome file not found: $GENOME" >&2; exit 1; }
[[ $# -lt 1 ]]     && { echo "ERROR: no peak files given" >&2; usage; }

command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found in PATH" >&2; exit 1; }

mkdir -p "${OUTDIR%/}/per_chrom"

# ---------------------------- genome BED -------------------------------------
# Accept chrom.sizes (chrom, size) or BED (chrom, start, end); keep 3 columns
# so that bedtools coverage output column 7 is the covered fraction.
genome_bed="$(mktemp "${OUTDIR%/}/genome.XXXXXX.bed")"
trap 'rm -f "$genome_bed"' EXIT

awk 'BEGIN{OFS="\t"} /^#/ || NF==0 {next}
     NF==2 {print $1, 0, $2; next}
     {print $1, $2, $3}' "$GENOME" > "$genome_bed"

echo "[$(date '+%F %T')] Genome file: $GENOME ($(wc -l < "$genome_bed") chromosomes/contigs)"
echo "[$(date '+%F %T')] bedtools version: $(bedtools --version)"

# ---------------------------- coverage per peak file -------------------------
summary="${OUTDIR%/}/genome_coverage.txt"
printf "sample\tpeak_file\tgenome_coverage_pct\n" > "$summary"

for peaks in "$@"; do
    if [[ ! -f "$peaks" ]]; then
        echo "WARNING: $peaks not found, skipping" >&2
        continue
    fi

    sample="$(basename "$peaks")"
    sample="${sample%.gz}"
    sample="${sample%.broadPeak}"
    sample="${sample%.narrowPeak}"
    sample="${sample%.bed}"
    per_chrom="${OUTDIR%/}/per_chrom/${sample}_coverage_per_chrom.tsv"

    echo "[$(date '+%F %T')] $sample"
    bedtools coverage -a "$genome_bed" -b "$peaks" > "$per_chrom"

    # Mean covered fraction over chromosomes, as %
    awk -v s="$sample" -v f="$peaks" '
        { sum += $7; n++ }
        END {
            if (n == 0) { print "ERROR: empty coverage for " f > "/dev/stderr"; exit 1 }
            printf "%s\t%s\t%.2f\n", s, f, (sum / n) * 100
        }' "$per_chrom" >> "$summary"
done

echo "[$(date '+%F %T')] Done. Summary: $summary"
column -t -s $'\t' "$summary" 2>/dev/null || cat "$summary"
