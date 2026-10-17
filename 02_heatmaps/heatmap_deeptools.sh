#!/bin/bash

# ---------------------------- SLURM job parameters ---------------------------
# Edit --account before submitting. Ignored when run with 'bash'.
#SBATCH --account=<your_account>
#SBATCH --job-name=heatmap_deeptools
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --cpus-per-task=20
#SBATCH --output=logs/%x_%j.log

# =============================================================================
# heatmap_deeptools.sh
#
# Heatmap of ChIP-seq signal centred on broad peaks/domains, sorted by domain
# size (deepTools computeMatrix reference-point + plotHeatmap).
#
# Used for the size-sorted domain heatmaps in the manuscript (e.g. Fig. 1b, 2b).
#
# Input : -r  regions (broadPeak / BED file, e.g. broad H3K4me1 domains)
#         -s  RPKM BigWig track (from 01_preprocessing/genomic_track_UCSC.sh)
# Output (in <outdir>):
#   <name>_matrix_all_broadPeaks_RP.gz         deepTools matrix (input for plotHeatmap)
#   <name>_matrix_all_broadPeaks_RP_values.tab count matrix, plain text (regions x bins)
#   <name>_IndividualValues_domain_RP.tab      regions kept in the matrix (BED-like)
#   <name>_WB_HM_AP_sorted_regions.bed         regions in heatmap order (top to bottom)
#   <name>_WB_HM_AP.svg                        heatmap (white-black)
#
# computeMatrix settings:
#   reference-point, --referencePoint center   centred on the domain midpoint
#   -b 15000 -a 15000                          +/- 15 kb window (default, see -w)
#   --binSize 10                               10 bp bins
#   --skipZeros, --missingDataAsZero           drop all-zero regions; NaN -> 0
#   --blackListFileName                        ENCODE blacklist regions excluded
# plotHeatmap settings:
#   white-black colour scale, zMin 0, zMax 22 (default, see -z)
#   regions sorted by region length (domain size, descending)
#
# Usage (SLURM; submit from the repository root, the logs/ folder must exist):
#   mkdir -p logs
#   sbatch heatmap_deeptools.sh -r <peaks.broadPeak> -s <sample_RPKM.bw> -n <name> [options]
# Usage (interactive / local):
#   bash heatmap_deeptools.sh -r <peaks.broadPeak> -s <sample_RPKM.bw> -n <name> [options]
#
#   -r  Regions file (.broadPeak or .bed)                         [required]
#   -s  Signal BigWig file                                        [required]
#   -n  Sample name (output prefix and heatmap label)             [required]
#   -l  Blacklist BED file (default: none)
#   -z  Colour scale maximum, zMax (default: 22)
#   -w  Window up-/downstream of the domain centre in bp (default: 15000)
#   -L  Regions label on the heatmap (default: broadPeak)
#   -o  Output directory (default: ./heatmaps)
#   -t  Number of threads (default: SLURM --cpus-per-task, otherwise 20)
#
# Example:
#   sbatch heatmap_deeptools.sh -r peaks/P10_GO_H3K4me1.broadPeak \
#       -s bigwig/P10_GO_H3K4me1_rep1_dedup_RPKM.bw -n P10_GO_H3K4me1_rep1 \
#       -l mm10-blacklist.v2.bed -z 22 -o heatmaps/
#
# Software: deepTools 3.5.0
# (Saga HPC module: deepTools/3.5.0-foss-2021a)
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
    module load deepTools/3.5.0-foss-2021a
fi

# ---------------------------- defaults ---------------------------------------
THREADS="${SLURM_CPUS_PER_TASK:-20}"
OUTDIR="./heatmaps"
PEAKS=""
BIGWIG=""
NAME=""
BLACKLIST=""
ZMAX=22
WINDOW=15000
REGIONS_LABEL="broadPeak"

usage() {
    sed -n '/^# =====/,/^# =====/p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
}

# ---------------------------- arguments --------------------------------------
while getopts ":r:s:n:l:z:w:L:o:t:h" opt; do
    case $opt in
        r) PEAKS="$OPTARG" ;;
        s) BIGWIG="$OPTARG" ;;
        n) NAME="$OPTARG" ;;
        l) BLACKLIST="$OPTARG" ;;
        z) ZMAX="$OPTARG" ;;
        w) WINDOW="$OPTARG" ;;
        L) REGIONS_LABEL="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        \?) echo "ERROR: invalid option -$OPTARG" >&2; usage ;;
        :)  echo "ERROR: option -$OPTARG requires an argument" >&2; usage ;;
    esac
done

[[ -z "$PEAKS"  ]] && { echo "ERROR: -r <regions file> is required" >&2; usage; }
[[ -z "$BIGWIG" ]] && { echo "ERROR: -s <BigWig file> is required" >&2; usage; }
[[ -z "$NAME"   ]] && { echo "ERROR: -n <name> is required" >&2; usage; }
[[ -f "$PEAKS"  ]] || { echo "ERROR: regions file not found: $PEAKS" >&2; exit 1; }
[[ -f "$BIGWIG" ]] || { echo "ERROR: BigWig file not found: $BIGWIG" >&2; exit 1; }
if [[ -n "$BLACKLIST" && ! -f "$BLACKLIST" ]]; then
    echo "ERROR: blacklist file not found: $BLACKLIST" >&2; exit 1
fi

command -v computeMatrix >/dev/null 2>&1 || { echo "ERROR: computeMatrix (deepTools) not found in PATH" >&2; exit 1; }
command -v plotHeatmap   >/dev/null 2>&1 || { echo "ERROR: plotHeatmap (deepTools) not found in PATH" >&2; exit 1; }

mkdir -p "$OUTDIR"
prefix="${OUTDIR%/}/${NAME}"

# Blacklist is optional
blacklist_arg=()
[[ -n "$BLACKLIST" ]] && blacklist_arg=(--blackListFileName "$BLACKLIST")

echo "[$(date '+%F %T')] Sample: $NAME | regions: $PEAKS ($(grep -vc '^#' "$PEAKS") lines) | threads: $THREADS"
echo "[$(date '+%F %T')] Window: +/-${WINDOW} bp | zMax: $ZMAX | blacklist: ${BLACKLIST:-none}"
echo "[$(date '+%F %T')] deepTools version: $(computeMatrix --version 2>&1)"

# ---------------------------- 1. count matrix --------------------------------
echo "[$(date '+%F %T')] Computing matrix"
computeMatrix reference-point \
    -p "$THREADS" \
    --referencePoint center \
    -b "$WINDOW" -a "$WINDOW" \
    -R "$PEAKS" \
    -S "$BIGWIG" \
    --skipZeros \
    -o "${prefix}_matrix_all_broadPeaks_RP.gz" \
    --outFileNameMatrix "${prefix}_matrix_all_broadPeaks_RP_values.tab" \
    --outFileSortedRegions "${prefix}_IndividualValues_domain_RP.tab" \
    --missingDataAsZero \
    ${blacklist_arg[@]+"${blacklist_arg[@]}"} \
    --binSize 10

# ---------------------------- 2. heatmap -------------------------------------
echo "[$(date '+%F %T')] Plotting heatmap"
plotHeatmap \
    -m "${prefix}_matrix_all_broadPeaks_RP.gz" \
    -out "${prefix}_WB_HM_AP.svg" \
    --outFileSortedRegions "${prefix}_WB_HM_AP_sorted_regions.bed" \
    --colorList 'white, black' \
    --dpi 300 \
    --heatmapHeight 21 \
    --heatmapWidth 7 \
    --plotFileFormat "svg" \
    --samplesLabel "$NAME" \
    --regionsLabel "$REGIONS_LABEL" \
    --zMin 0 \
    --zMax "$ZMAX" \
    --sortUsing region_length \
    --missingDataColor "#FFFFFF"

echo "[$(date '+%F %T')] Done. Outputs in $OUTDIR:"
ls -1 "${prefix}"_* | sed 's/^/    /'
