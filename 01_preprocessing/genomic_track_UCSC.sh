#!/bin/bash

# ---------------------------- SLURM job parameters ---------------------------
# Edit --account before submitting. Ignored when run with 'bash'.
#SBATCH --account=<your_account>
#SBATCH --job-name=bigwig_UCSC
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --cpus-per-task=20
#SBATCH --output=logs/%x_%j.log

# =============================================================================
# genomic_track_UCSC.sh
#
# Generate RPKM-normalized BigWig coverage tracks for UCSC Genome Browser
# visualization (deepTools bamCoverage).
#
# This script was used to generate all genome browser tracks shown in the
# manuscript (Fig. 1a, 3d, 4a, 4c; Supp. Fig. 1a, 1f, 2a, 2f, 4a, 4b).
#
# Input : deduplicated, coordinate-sorted, indexed BAM file(s) (*_dedup.bam)
# Output: <prefix>_dedup_RPKM.bw  (10 bp bins, RPKM normalized)
#
# Normalization settings:
#   --binSize 10                      10 bp resolution
#   --normalizeUsing RPKM             reads per kilobase per million mapped reads
#   --ignoreForNormalization chrX chrM  sex chromosome / mitochondrial reads
#                                     excluded when computing the scaling factor
#   --ignoreDuplicates                duplicate reads are not counted
#   --extendReads                     reads extended to fragment length
#                                     (paired-end: fragment size from mate pairs)
#
# Usage (SLURM; submit from the repository root, the logs/ folder must exist):
#   mkdir -p logs
#   sbatch genomic_track_UCSC.sh -g <genome> -o <outdir> sample1_dedup.bam [sample2_dedup.bam ...]
# Usage (interactive / local):
#   bash genomic_track_UCSC.sh -g <genome> -t <threads> -o <outdir> sample1_dedup.bam [...]
#
#   -g  Genome: mm10 | hg38 | bosTau8, or an explicit effective genome size (integer)
#   -t  Number of threads (default: SLURM --cpus-per-task, otherwise 20)
#   -o  Output directory (default: ./bigwig)
#
# Example:
#   sbatch genomic_track_UCSC.sh -g mm10 -o bigwig/ bam/P10_GO_H3K4me1_rep1_dedup.bam
#
# Software: deepTools 3.5.0, SAMtools 1.13
# (Saga HPC modules: deepTools/3.5.0-foss-2021a, SAMtools/1.13-GCC-10.3.0)
# =============================================================================

# Enable error checking
set -o errexit
set -o nounset
set -o pipefail

# ---------------------------- load modules -----------------------------------
# Module names are those of the Saga HPC cluster; adjust to your system.
# Skipped automatically when the 'module' command is not available (e.g. when
# using the conda environment from environment.yml).
if command -v module >/dev/null 2>&1; then
    module --quiet purge
    module load deepTools/3.5.0-foss-2021a SAMtools/1.13-GCC-10.3.0
fi

# ---------------------------- defaults ---------------------------------------
THREADS="${SLURM_CPUS_PER_TASK:-20}"
OUTDIR="./bigwig"
GENOME=""

# Effective genome sizes (non-N bases) used in the analysis.
declare -A EFFECTIVE_GENOME_SIZE=(
    [mm10]=2652783500     # GRCm38 (deepTools documentation)
    [hg38]=2913022398     # GRCh38 (deepTools documentation)
    [bosTau8]=2640168518  # UMD3.1 (not in deepTools table)
)

usage() {
    sed -n '/^# =====/,/^# =====/p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
}

# ---------------------------- arguments --------------------------------------
while getopts ":g:t:o:h" opt; do
    case $opt in
        g) GENOME="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        h) usage ;;
        \?) echo "ERROR: invalid option -$OPTARG" >&2; usage ;;
        :)  echo "ERROR: option -$OPTARG requires an argument" >&2; usage ;;
    esac
done
shift $((OPTIND - 1))

[[ -z "$GENOME" ]] && { echo "ERROR: -g <genome> is required" >&2; usage; }
[[ $# -lt 1 ]]     && { echo "ERROR: no BAM files given" >&2; usage; }

# Resolve effective genome size (preset name or explicit integer)
if [[ "$GENOME" =~ ^[0-9]+$ ]]; then
    effective_genome_size="$GENOME"
elif [[ -n "${EFFECTIVE_GENOME_SIZE[$GENOME]+x}" ]]; then
    effective_genome_size="${EFFECTIVE_GENOME_SIZE[$GENOME]}"
else
    echo "ERROR: unknown genome '$GENOME' (use mm10, hg38, bosTau8 or an integer)" >&2
    exit 1
fi
if [[ "$effective_genome_size" -le 0 ]]; then
    echo "ERROR: effective genome size for '$GENOME' is not set. Pass it explicitly with -g <integer>." >&2
    exit 1
fi

command -v bamCoverage >/dev/null 2>&1 || { echo "ERROR: bamCoverage (deepTools) not found in PATH" >&2; exit 1; }
command -v samtools    >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH" >&2; exit 1; }

mkdir -p "$OUTDIR"

echo "[$(date '+%F %T')] Genome: $GENOME | effective genome size: $effective_genome_size | threads: $THREADS"
echo "[$(date '+%F %T')] deepTools version: $(bamCoverage --version 2>&1)"

# ---------------------------- main loop --------------------------------------
for bam in "$@"; do
    if [[ ! -f "$bam" ]]; then
        echo "WARNING: $bam not found, skipping" >&2
        continue
    fi

    # Index BAM if needed
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
        echo "[$(date '+%F %T')] Indexing $bam"
        samtools index -@ "$THREADS" "$bam"
    fi

    # Output prefix: strip directory and the _dedup.bam / .bam suffix
    out="$(basename "$bam")"
    out="${out%.bam}"
    out="${out%_dedup}"
    out="${OUTDIR%/}/${out}"

    echo "[$(date '+%F %T')] Processing $bam -> ${out}_dedup_RPKM.bw"

    bamCoverage -p "$THREADS" \
        --bam "$bam" \
        --outFileName "${out}_dedup_RPKM.bw" \
        --binSize 10 \
        --normalizeUsing RPKM \
        --effectiveGenomeSize "$effective_genome_size" \
        --ignoreForNormalization chrX chrM \
        --ignoreDuplicates \
        --extendReads

    echo "[$(date '+%F %T')] Done: ${out}_dedup_RPKM.bw"
done

echo "[$(date '+%F %T')] All BigWig tracks generated in $OUTDIR"
