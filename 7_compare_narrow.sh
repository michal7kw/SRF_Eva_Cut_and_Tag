#!/bin/bash

#===============================================================================
# SCRIPT: 7_compare.sh
# PURPOSE: Peak comparison and visualization analysis for Cut&Tag data
#
# DESCRIPTION:
# This script performs comprehensive peak comparison analysis between different
# conditions (TES, TESmut, TEAD1) to identify overlapping, unique, and 
# condition-specific peaks. It generates statistical summaries and creates
# visualization heatmaps and profiles for comparative analysis.
#
# KEY OPERATIONS:
# 1. Peak overlap analysis using bedtools intersect
# 2. Identification of condition-specific peaks
# 3. Generation of peak statistics and counts
# 4. Creation of signal heatmaps using deepTools
# 5. Generation of signal profiles around peak centers
#
# METHODOLOGY:
# - Uses bedtools for efficient genomic interval operations
# - Employs deepTools for matrix computation and visualization
# - Computes signal matrices centered on peak regions (±3kb)
# - Generates both heatmaps and line profiles for visual comparison
#
# IMPORTANT PARAMETERS:
# - Memory: 32GB (for large matrix computations)
# - Time: 4 hours (sufficient for visualization generation)
# - Threads: 16 (parallel processing for matrix computation)
# - Window: ±3kb around peak centers
# - Reference point: Peak center
#
# INPUTS:
# - Narrow peak files from MACS2 (05_peaks directory)
# - Normalized BigWig files (06_bigwig directory)
# - All sample replicates for comprehensive comparison
#
# OUTPUTS:
# - Overlap BED files (TES_TEAD1_overlap.bed, condition-specific peaks)
# - Peak statistics summary (peak_stats.txt)
# - Signal matrix (TES_peaks_matrix.gz)
# - Heatmap visualization (TES_peaks_heatmap.pdf)
# - Signal profile plot (TES_peaks_profile.pdf)
#
# DEPENDENCIES:
# - bedtools (genomic interval operations)
# - deepTools (computeMatrix, plotHeatmap, plotProfile)
# - Conda environment: bigwig
#
# USAGE:
# sbatch 7_compare.sh
#
# NOTES:
# - Focuses on TES peaks as reference for comparison
# - Includes all three conditions and their replicates
# - Generates publication-ready visualizations
# - Statistics help quantify overlap relationships
#===============================================================================

#SBATCH --job-name=7_analysis
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/7_analysis.err"
#SBATCH --output="./logs/7_analysis.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate bigwig

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG"
PEAKS_DIR="${BASE_DIR}/results/05_peaks_narrow"
CONSENSUS_DIR="${BASE_DIR}/results/11_combined_replicates_narrow"
BIGWIG_DIR="${BASE_DIR}/results/06_bigwig"
OUTPUT_DIR="${BASE_DIR}/results/07_analysis_narrow"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Determine which peak sets to use
# Priority: 1) Consensus peaks (if available), 2) Combined MACS2 peaks
echo "Determining peak sets to use for analysis..."

# Function to select best available peak file
select_peak_file() {
    local condition=$1
    local consensus_scored="${CONSENSUS_DIR}/peaks/${condition}_consensus_peaks_scored.bed"
    local consensus_full="${CONSENSUS_DIR}/peaks/${condition}_consensus_peaks.bed"
    local idr_union="${CONSENSUS_DIR}/idr/${condition}_idr_union.bed"
    local macs2_combined="${PEAKS_DIR}/${condition}_peaks.narrowPeak"

    if [ -f "$consensus_scored" ]; then
        echo "$consensus_scored"
        echo "  Using quality-scored consensus peaks for $condition" >&2
    elif [ -f "$consensus_full" ]; then
        echo "$consensus_full"
        echo "  Using full consensus peaks for $condition" >&2
    elif [ -f "$idr_union" ]; then
        echo "$idr_union"
        echo "  Using IDR union peaks for $condition" >&2
    elif [ -f "$macs2_combined" ]; then
        echo "$macs2_combined"
        echo "  Using MACS2 combined peaks for $condition (consensus not available)" >&2
    else
        echo ""
        echo "  ERROR: No peak file found for $condition" >&2
    fi
}

# Select peak files for each condition
TES_PEAKS=$(select_peak_file "TES")
TESMUT_PEAKS=$(select_peak_file "TESmut")
TEAD1_PEAKS=$(select_peak_file "TEAD1")

# Verify all peak files are available
if [ -z "$TES_PEAKS" ] || [ -z "$TESMUT_PEAKS" ] || [ -z "$TEAD1_PEAKS" ]; then
    echo "ERROR: Missing required peak files. Exiting."
    exit 1
fi

# 1. Comprehensive peak overlap analysis
# Now using consensus/high-quality peaks instead of raw MACS2 peaks
echo "Performing comprehensive peak overlap analysis using high-quality peak sets..."

# TES vs TEAD1 overlaps
bedtools intersect \
    -a $TES_PEAKS \
    -b $TEAD1_PEAKS \
    -u > ${OUTPUT_DIR}/TES_TEAD1_overlap.bed

bedtools intersect \
    -a $TES_PEAKS \
    -b $TEAD1_PEAKS \
    -v > ${OUTPUT_DIR}/TES_specific.bed

bedtools intersect \
    -a $TEAD1_PEAKS \
    -b $TES_PEAKS \
    -v > ${OUTPUT_DIR}/TEAD1_specific.bed

# TES vs TESmut overlaps
bedtools intersect \
    -a $TES_PEAKS \
    -b $TESMUT_PEAKS \
    -u > ${OUTPUT_DIR}/TES_TESmut_overlap.bed

bedtools intersect \
    -a $TES_PEAKS \
    -b $TESMUT_PEAKS \
    -v > ${OUTPUT_DIR}/TES_not_in_TESmut.bed

bedtools intersect \
    -a $TESMUT_PEAKS \
    -b $TES_PEAKS \
    -v > ${OUTPUT_DIR}/TESmut_specific.bed

# TESmut vs TEAD1 overlaps
bedtools intersect \
    -a $TESMUT_PEAKS \
    -b $TEAD1_PEAKS \
    -u > ${OUTPUT_DIR}/TESmut_TEAD1_overlap.bed

bedtools intersect \
    -a $TESMUT_PEAKS \
    -b $TEAD1_PEAKS \
    -v > ${OUTPUT_DIR}/TESmut_not_in_TEAD1.bed

bedtools intersect \
    -a $TEAD1_PEAKS \
    -b $TESMUT_PEAKS \
    -v > ${OUTPUT_DIR}/TEAD1_not_in_TESmut.bed

# Three-way overlaps
echo "Finding three-way overlaps..."
bedtools intersect \
    -a $TES_PEAKS \
    -b $TESMUT_PEAKS \
    -u | bedtools intersect \
    -a stdin \
    -b $TEAD1_PEAKS \
    -u > ${OUTPUT_DIR}/TES_TESmut_TEAD1_overlap.bed

# Peaks unique to each condition (not in any other)
bedtools intersect \
    -a $TES_PEAKS \
    -b $TESMUT_PEAKS \
    -v | bedtools intersect \
    -a stdin \
    -b $TEAD1_PEAKS \
    -v > ${OUTPUT_DIR}/TES_unique.bed

bedtools intersect \
    -a $TESMUT_PEAKS \
    -b $TES_PEAKS \
    -v | bedtools intersect \
    -a stdin \
    -b $TEAD1_PEAKS \
    -v > ${OUTPUT_DIR}/TESmut_unique.bed

bedtools intersect \
    -a $TEAD1_PEAKS \
    -b $TES_PEAKS \
    -v | bedtools intersect \
    -a stdin \
    -b $TESMUT_PEAKS \
    -v > ${OUTPUT_DIR}/TEAD1_unique.bed

# 2. Generate comprehensive peak overlap statistics
echo "Generating comprehensive statistics..."
{
    echo "=== ANALYSIS USING HIGH-QUALITY PEAK SETS ==="
    echo "Peak files used:"
    echo "  TES: $TES_PEAKS"
    echo "  TESmut: $TESMUT_PEAKS"
    echo "  TEAD1: $TEAD1_PEAKS"
    echo ""

    echo "=== PEAK COUNTS ==="
    wc -l $TES_PEAKS | awk '{print "TES total peaks: "$1}'
    wc -l $TESMUT_PEAKS | awk '{print "TESmut total peaks: "$1}'
    wc -l $TEAD1_PEAKS | awk '{print "TEAD1 total peaks: "$1}'

    echo ""
    echo "=== PAIRWISE OVERLAPS ==="
    wc -l ${OUTPUT_DIR}/TES_TEAD1_overlap.bed | awk '{print "TES-TEAD1 overlapping peaks: "$1}'
    wc -l ${OUTPUT_DIR}/TES_TESmut_overlap.bed | awk '{print "TES-TESmut overlapping peaks: "$1}'
    wc -l ${OUTPUT_DIR}/TESmut_TEAD1_overlap.bed | awk '{print "TESmut-TEAD1 overlapping peaks: "$1}'

    echo ""
    echo "=== CONDITION-SPECIFIC PEAKS ==="
    wc -l ${OUTPUT_DIR}/TES_specific.bed | awk '{print "TES-specific (vs TEAD1): "$1}'
    wc -l ${OUTPUT_DIR}/TEAD1_specific.bed | awk '{print "TEAD1-specific (vs TES): "$1}'
    wc -l ${OUTPUT_DIR}/TESmut_specific.bed | awk '{print "TESmut-specific (vs TES): "$1}'
    wc -l ${OUTPUT_DIR}/TES_not_in_TESmut.bed | awk '{print "TES not in TESmut: "$1}'
    wc -l ${OUTPUT_DIR}/TESmut_not_in_TEAD1.bed | awk '{print "TESmut not in TEAD1: "$1}'
    wc -l ${OUTPUT_DIR}/TEAD1_not_in_TESmut.bed | awk '{print "TEAD1 not in TESmut: "$1}'

    echo ""
    echo "=== THREE-WAY COMPARISONS ==="
    wc -l ${OUTPUT_DIR}/TES_TESmut_TEAD1_overlap.bed | awk '{print "Common to all three conditions: "$1}'
    wc -l ${OUTPUT_DIR}/TES_unique.bed | awk '{print "TES unique (not in any other): "$1}'
    wc -l ${OUTPUT_DIR}/TESmut_unique.bed | awk '{print "TESmut unique (not in any other): "$1}'
    wc -l ${OUTPUT_DIR}/TEAD1_unique.bed | awk '{print "TEAD1 unique (not in any other): "$1}'

    echo ""
    echo "=== OVERLAP PERCENTAGES ==="
    TES_TOTAL=$(wc -l < $TES_PEAKS)
    TESMUT_TOTAL=$(wc -l < $TESMUT_PEAKS)
    TEAD1_TOTAL=$(wc -l < $TEAD1_PEAKS)
    TES_TEAD1_OVERLAP=$(wc -l < ${OUTPUT_DIR}/TES_TEAD1_overlap.bed)
    TES_TESMUT_OVERLAP=$(wc -l < ${OUTPUT_DIR}/TES_TESmut_overlap.bed)
    TESMUT_TEAD1_OVERLAP=$(wc -l < ${OUTPUT_DIR}/TESmut_TEAD1_overlap.bed)

    echo "TES-TEAD1 overlap %: $(echo "scale=2; $TES_TEAD1_OVERLAP * 100 / $TES_TOTAL" | bc)% of TES peaks"
    echo "TES-TESmut overlap %: $(echo "scale=2; $TES_TESMUT_OVERLAP * 100 / $TES_TOTAL" | bc)% of TES peaks"
    echo "TESmut-TEAD1 overlap %: $(echo "scale=2; $TESMUT_TEAD1_OVERLAP * 100 / $TESMUT_TOTAL" | bc)% of TESmut peaks"

    echo ""
    echo "NOTE: These statistics use high-quality consensus peaks (when available)"
    echo "      which provide better reproducibility than raw MACS2 combined peaks."
} > ${OUTPUT_DIR}/peak_stats.txt

# 3. Create comprehensive heatmaps and profiles
echo "Creating comprehensive visualizations..."

# Function to create matrix and plots for a given peak set
create_visualization() {
    local PEAK_FILE=$1
    local OUTPUT_PREFIX=$2
    local TITLE=$3

    echo "Creating visualization for $TITLE..."

    # Compute matrix
    computeMatrix reference-point \
        --referencePoint center \
        -b 3000 -a 3000 \
        -R $PEAK_FILE \
        -S ${BIGWIG_DIR}/TES-1_CPM.bw \
           ${BIGWIG_DIR}/TES-2_CPM.bw \
           ${BIGWIG_DIR}/TES-3_CPM.bw \
           ${BIGWIG_DIR}/TESmut-1_CPM.bw \
           ${BIGWIG_DIR}/TESmut-2_CPM.bw \
           ${BIGWIG_DIR}/TESmut-3_CPM.bw \
           ${BIGWIG_DIR}/TEAD1-1_CPM.bw \
           ${BIGWIG_DIR}/TEAD1-2_CPM.bw \
           ${BIGWIG_DIR}/TEAD1-3_CPM.bw \
        --skipZeros \
        -o ${OUTPUT_DIR}/${OUTPUT_PREFIX}_matrix.gz \
        -p 16

    # Plot heatmap
    plotHeatmap \
        -m ${OUTPUT_DIR}/${OUTPUT_PREFIX}_matrix.gz \
        -out ${OUTPUT_DIR}/${OUTPUT_PREFIX}_heatmap.pdf \
        --colorMap RdBu_r \
        --whatToShow 'heatmap and colorbar' \
        --samplesLabel TES-1 TES-2 TES-3 TESmut-1 TESmut-2 TESmut-3 TEAD1-1 TEAD1-2 TEAD1-3 \
        --plotTitle "$TITLE Heatmap"

    # Plot profile
    plotProfile \
        -m ${OUTPUT_DIR}/${OUTPUT_PREFIX}_matrix.gz \
        -out ${OUTPUT_DIR}/${OUTPUT_PREFIX}_profile.pdf \
        --samplesLabel TES-1 TES-2 TES-3 TESmut-1 TESmut-2 TESmut-3 TEAD1-1 TEAD1-2 TEAD1-3 \
        --plotTitle "Signal at $TITLE"

    # Convert to PNG
    convert -density 300 ${OUTPUT_DIR}/${OUTPUT_PREFIX}_heatmap.pdf ${OUTPUT_DIR}/${OUTPUT_PREFIX}_heatmap.png
    convert -density 300 ${OUTPUT_DIR}/${OUTPUT_PREFIX}_profile.pdf ${OUTPUT_DIR}/${OUTPUT_PREFIX}_profile.png
}

# Create visualizations for all peak sets (using selected high-quality peaks)
create_visualization "$TES_PEAKS" "TES_peaks" "TES High-Quality Peaks"
create_visualization "$TESMUT_PEAKS" "TESmut_peaks" "TESmut High-Quality Peaks"
create_visualization "$TEAD1_PEAKS" "TEAD1_peaks" "TEAD1 High-Quality Peaks"

# Create visualizations for overlap sets (if they have sufficient peaks)
if [ -s ${OUTPUT_DIR}/TES_TEAD1_overlap.bed ]; then
    create_visualization "${OUTPUT_DIR}/TES_TEAD1_overlap.bed" "TES_TEAD1_overlap" "TES-TEAD1 Overlapping Peaks"
fi

if [ -s ${OUTPUT_DIR}/TES_TESmut_overlap.bed ]; then
    create_visualization "${OUTPUT_DIR}/TES_TESmut_overlap.bed" "TES_TESmut_overlap" "TES-TESmut Overlapping Peaks"
fi

if [ -s ${OUTPUT_DIR}/TES_TESmut_TEAD1_overlap.bed ]; then
    create_visualization "${OUTPUT_DIR}/TES_TESmut_TEAD1_overlap.bed" "TES_TESmut_TEAD1_overlap" "Three-way Overlapping Peaks"
fi

echo "Analysis complete!"