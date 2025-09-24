#!/bin/bash

#===============================================================================
# SCRIPT: 12_pairwise_overlap.sh
# PURPOSE: Comprehensive pairwise binding site overlap analysis
#
# DESCRIPTION:
# This script performs detailed pairwise overlap analysis between different
# experimental conditions (TES, TESmut, TEAD1) to identify shared and
# condition-specific binding sites. It generates Venn diagrams, overlap
# statistics, and comprehensive visualizations to understand binding
# specificity and relationships between different experimental conditions.
#
# KEY OPERATIONS:
# 1. Automatic peak file detection (consensus, broad, or narrow peaks)
# 2. Pairwise overlap analysis using bedtools intersect
# 3. Identification of unique and shared binding sites
# 4. Venn diagram generation for visual representation
# 5. Statistical analysis of overlap percentages
# 6. Comprehensive summary report generation
#
# METHODOLOGY:
# - Prioritizes consensus peaks, then broad peaks, then narrow peaks.
# - Uses bedtools for genomic interval overlap calculations
# - Employs R/VennDiagram for publication-quality visualizations
# - Calculates multiple overlap metrics (intersection, unique regions)
# - Generates both numerical and visual summaries
# - Provides biological interpretation of results
#
# IMPORTANT PARAMETERS:
# - Memory: 16GB (for large peak file operations)
# - Time: 2 hours (adequate for comprehensive analysis)
# - Threads: 8 (parallel processing for efficiency)
# - Overlap threshold: Minimum 1bp overlap
# - Peak format: Supports consensus BED, broadPeak, and narrowPeak
#
# INPUTS:
# - Peak files from step 11 (consensus), 5b (broad), or 5 (narrow)
# - Automatically detects and uses optimal peak type
# - Three experimental conditions: TES, TESmut, TEAD1
#
# OUTPUTS:
# - Pairwise overlap BED files for each comparison
# - Condition-specific unique peak sets
# - Venn diagrams for visual overlap representation
# - Comprehensive overlap statistics
# - Detailed summary report with biological interpretation
# - R scripts for custom visualization
#
# DEPENDENCIES:
# - bedtools for genomic interval operations
# - R with VennDiagram, ggplot2, dplyr packages
# - Conda environment: overlap
# - Peak calling results from previous steps
#
# USAGE:
# sbatch 12_pairwise_overlap.sh
#
# NOTES:
# - Automatically selects optimal peak type (consensus preferred)
# - Provides guidance for peak type selection
# - Essential for understanding binding specificity
# - Generates publication-ready visualizations
# - Includes biological interpretation and next steps
#===============================================================================

#SBATCH --job-name=pairwise_overlap
#SBATCH --output=./logs/12_pairwise_overlap.out
#SBATCH --error=./logs/12_pairwise_overlap.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00

echo "=== Pairwise Binding Site Overlap Analysis ==="
echo "Start time: $(date)"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate overlap

# Set up directories
OUTDIR="results/12_pairwise_overlap"
mkdir -p "$OUTDIR"/{overlaps,unique,venn,matrices,plots,stats}

# Check required directories and determine which peak files to use
CONSENSUS_PEAKS_AVAILABLE=false
BROAD_PEAKS_AVAILABLE=false
NARROW_PEAKS_AVAILABLE=false

# Check for consensus peaks (highest priority)
if [ -d "results/11_combined_replicates_broad/peaks" ]; then
    if find "results/11_combined_replicates_broad/peaks" -mindepth 1 -print -quit 2>/dev/null | grep -q .; then
        CONSENSUS_PEAKS_AVAILABLE=true
        echo "Consensus broad peaks directory found - using for overlap analysis."
    fi
elif [ -d "results/11_combined_replicates_narrow/peaks" ]; then
    if find "results/11_combined_replicates_narrow/peaks" -mindepth 1 -print -quit 2>/dev/null | grep -q .; then
        CONSENSUS_PEAKS_AVAILABLE=true
        echo "Consensus narrow peaks directory found - using for overlap analysis."
    fi
fi

# Check for broad peaks (second priority)
if [ -d "results/05_peaks_broad" ]; then
    BROAD_PEAKS_AVAILABLE=true
    echo "Broad peaks directory found."
fi

# Check for narrow peaks (fallback)
if [ -d "results/05_peaks_narrow" ]; then
    NARROW_PEAKS_AVAILABLE=true
    echo "Narrow peaks directory found."
fi

# Determine which peaks to use based on priority
if [ "$CONSENSUS_PEAKS_AVAILABLE" = true ]; then
    echo "Using CONSENSUS PEAKS for overlap analysis (recommended)"
    if [ -d "results/11_combined_replicates_broad/peaks" ]; then
        PEAK_DIR="results/11_combined_replicates_broad/peaks"
        ANALYSIS_TYPE="consensus_broad"
    else
        PEAK_DIR="results/11_combined_replicates_narrow/peaks"
        ANALYSIS_TYPE="consensus_narrow"
    fi
    PEAK_SUFFIX="_consensus_peaks.bed"
elif [ "$BROAD_PEAKS_AVAILABLE" = true ]; then
    echo "Using BROAD PEAKS for overlap analysis. Note: For best results, run 11b_combine_replicates_broad.sh to generate consensus peaks."
    PEAK_DIR="results/05_peaks_broad"
    PEAK_SUFFIX="_peaks.broadPeak"
    ANALYSIS_TYPE="broad"
elif [ "$NARROW_PEAKS_AVAILABLE" = true ]; then
    echo "Using NARROW PEAKS for overlap analysis (suboptimal). Note: For best results, generate consensus or broad peaks."
    PEAK_DIR="results/05_peaks_narrow"
    PEAK_SUFFIX="_peaks.narrowPeak"
    ANALYSIS_TYPE="narrow"
else
    echo "Error: No peak calling results found in results/11_combined_replicates_*/peaks, results/05_peaks_broad, or results/05_peaks_narrow."
    exit 1
fi

# Define conditions and peak files
CONDITIONS=("TES" "TESmut" "TEAD1")
declare -A PEAK_FILES=(
    ["TES"]="${PEAK_DIR}/TES${PEAK_SUFFIX}"
    ["TESmut"]="${PEAK_DIR}/TESmut${PEAK_SUFFIX}"
    ["TEAD1"]="${PEAK_DIR}/TEAD1${PEAK_SUFFIX}"
)

# Check if peak files exist
echo "=== Checking peak files ==="
for condition in "${CONDITIONS[@]}"; do
    if [ -f "${PEAK_FILES[$condition]}" ]; then
        peak_count=$(wc -l < "${PEAK_FILES[$condition]}")
        echo "$condition: $peak_count peaks (${PEAK_FILES[$condition]})
"
    else
        echo "Error: Peak file not found for $condition: ${PEAK_FILES[$condition]}"
        exit 1
    fi
done

echo ""
echo "=== Performing pairwise overlap analysis ==="

# Function to perform pairwise comparison
perform_pairwise_comparison() {
    local cond1=$1
    local cond2=$2
    local file1="${PEAK_FILES[$cond1]}"
    local file2="${PEAK_FILES[$cond2]}"
    
    echo "Analyzing $cond1 vs $cond2 using $ANALYSIS_TYPE peaks..." >&2

    # Create comparison directory
    mkdir -p "$OUTDIR/overlaps/${cond1}_vs_${cond2}"

    # Convert to BED3 format for consistent comparison, removing potential headers
    grep -vE '^#|^track' "$file1" | cut -f1-3 > "$OUTDIR/${cond1}_temp.bed"
    grep -vE '^#|^track' "$file2" | cut -f1-3 > "$OUTDIR/${cond2}_temp.bed"

    # Find overlapping peaks (minimum 1bp overlap)
    echo "  Finding overlapping regions..." >&2
    bedtools intersect -a "$OUTDIR/${cond1}_temp.bed" -b "$OUTDIR/${cond2}_temp.bed" -wa > "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond1}_overlapping_${cond2}.bed"
    bedtools intersect -a "$OUTDIR/${cond2}_temp.bed" -b "$OUTDIR/${cond1}_temp.bed" -wa > "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond2}_overlapping_${cond1}.bed"

    # Find actual overlapping regions (intersection)
    echo "  Finding intersection regions..." >&2
    bedtools intersect -a "$OUTDIR/${cond1}_temp.bed" -b "$OUTDIR/${cond2}_temp.bed" > "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond1}_${cond2}_intersection.bed"

    # Find unique peaks (no overlap)
    echo "  Finding unique regions..." >&2
    bedtools intersect -a "$OUTDIR/${cond1}_temp.bed" -b "$OUTDIR/${cond2}_temp.bed" -v > "$OUTDIR/unique/${cond1}_unique_vs_${cond2}.bed"
    bedtools intersect -a "$OUTDIR/${cond2}_temp.bed" -b "$OUTDIR/${cond1}_temp.bed" -v > "$OUTDIR/unique/${cond2}_unique_vs_${cond1}.bed"
    
    # Count results
    overlap_cond1=$(wc -l < "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond1}_overlapping_${cond2}.bed")
    overlap_cond2=$(wc -l < "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond2}_overlapping_${cond1}.bed")
    intersection_count=$(wc -l < "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond1}_${cond2}_intersection.bed")
    unique_cond1=$(wc -l < "$OUTDIR/unique/${cond1}_unique_vs_${cond2}.bed")
    unique_cond2=$(wc -l < "$OUTDIR/unique/${cond2}_unique_vs_${cond1}.bed")
    total_cond1=$(wc -l < "$OUTDIR/${cond1}_temp.bed")
    total_cond2=$(wc -l < "$OUTDIR/${cond2}_temp.bed")
    
    # Calculate statistics
    overlap_percent_cond1=$(awk "BEGIN {printf \"%.1f\", ($overlap_cond1/$total_cond1)*100}")
    overlap_percent_cond2=$(awk "BEGIN {printf \"%.1f\", ($overlap_cond2/$total_cond2)*100}")
    intersection_percent_cond1=$(awk "BEGIN {printf \"%.1f\", ($intersection_count/$total_cond1)*100}")
    intersection_percent_cond2=$(awk "BEGIN {printf \"%.1f\", ($intersection_count/$total_cond2)*100}")
    
    echo "  Results:" >&2
    echo "    $cond1 peaks overlapping with $cond2: $overlap_cond1 / $total_cond1 (${overlap_percent_cond1}%)" >&2
    echo "    $cond2 peaks overlapping with $cond1: $overlap_cond2 / $total_cond2 (${overlap_percent_cond2}%)" >&2
    echo "    Intersection (common regions): $intersection_count" >&2
    echo "    $cond1 unique peaks: $unique_cond1" >&2
    echo "    $cond2 unique peaks: $unique_cond2" >&2
    
    # Clean up temp files
    rm "$OUTDIR/${cond1}_temp.bed" "$OUTDIR/${cond2}_temp.bed"
    
    # Return statistics for summary
    echo "${cond1}_vs_${cond2}:${overlap_cond1}:${overlap_cond2}:${unique_cond1}:${unique_cond2}:${total_cond1}:${total_cond2}:${intersection_count}:${intersection_percent_cond1}:${intersection_percent_cond2}"
}

# Perform all pairwise comparisons
echo ""
declare -A COMPARISON_STATS
declare -ATOTAL_PEAKS
declare -AINTERSECTION_COUNTS
declare -AUNIQUE_PEAKS_1
declare -AUNIQUE_PEAKS_2

CONDITIONS=("TES" "TESmut" "TEAD1")
for condition in "${CONDITIONS[@]}"; do
    TOTAL_PEAKS[$condition]=$(wc -l < "${PEAK_FILES[$condition]}")
done

comparisons=("TES:TESmut" "TES:TEAD1" "TESmut:TEAD1")

for comparison in "${comparisons[@]}"; do
    IFS=':' read -r cond1 cond2 <<< "$comparison"
    echo "Processing comparison: ${cond1} vs ${cond2}"
    stats_output=$(perform_pairwise_comparison "$cond1" "$cond2")
    echo "Stats output: $stats_output"

    # Parse the stats output
    IFS=':' read -r comp_name overlap1 overlap2 unique1 unique2 total1 total2 intersection_count intersection_percent1 intersection_percent2 <<< "$stats_output"

    COMPARISON_STATS["${cond1}_vs_${cond2}"]="$stats_output"
    INTERSECTION_COUNTS["${cond1}_${cond2}"]="$intersection_count"
    UNIQUE_PEAKS_1["${cond1}_${cond2}"]="$unique1"
    UNIQUE_PEAKS_2["${cond1}_${cond2}"]="$unique2"

    echo "Parsed intersection_count: $intersection_count"
    echo ""
done

echo "=== Creating Venn diagram data ==="

# Create R script for Venn diagrams and advanced statistics
cat > "$OUTDIR/create_venn_diagrams.R" << EOF
library(VennDiagram)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Pass counts from shell script, ensuring they are not empty
tes_total=${TOTAL_PEAKS[TES]:-0}
tesmut_total=${TOTAL_PEAKS[TESmut]:-0}
tead1_total=${TOTAL_PEAKS[TEAD1]:-0}

tes_tesmut_intersection=${INTERSECTION_COUNTS[TES_TESmut]:-0}
tes_tead1_intersection=${INTERSECTION_COUNTS[TES_TEAD1]:-0}
tesmut_tead1_intersection=${INTERSECTION_COUNTS[TESmut_TEAD1]:-0}

tes_unique_tesmut=${UNIQUE_PEAKS_1[TES_TESmut]:-0}
tesmut_unique_tes=${UNIQUE_PEAKS_2[TES_TESmut]:-0}
tes_unique_tead1=${UNIQUE_PEAKS_1[TES_TEAD1]:-0}
tead1_unique_tes=${UNIQUE_PEAKS_2[TES_TEAD1]:-0}
tesmut_unique_tead1=${UNIQUE_PEAKS_1[TESmut_TEAD1]:-0}
tead1_unique_tesmut=${UNIQUE_PEAKS_2[TESmut_TEAD1]:-0}


cat(paste("Peak counts:\n"))
cat(paste("TES:", tes_total, "\n"))
cat(paste("TESmut:", tesmut_total, "\n"))
cat(paste("TEAD1:", tead1_total, "\n\n"))

# Create Venn diagrams for all pairwise comparisons
output_dir <- "results/12_pairwise_overlap/venn"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. TES vs TESmut
cat("Creating TES vs TESmut Venn diagram...\n")
png(file.path(output_dir, "TES_vs_TESmut_${ANALYSIS_TYPE}_venn.png"),
    width = 3000, height = 3000, res = 300) # Set resolution and dimensions
venn.plot1 <- draw.pairwise.venn(
    area1 = tes_total,
    area2 = tesmut_total,
    cross.area = tes_tesmut_intersection,
    category = c("TES", "TESmut"),
    col = "transparent",
    fill = c("lightblue", "lightcoral"),
    alpha = 0.6,
    cex = 1.5,
    fontfamily = "serif",
    main = "TES vs TESmut $ANALYSIS_TYPE Peak Overlap",
    main.cex = 1.8
)
grid.draw(venn.plot1)
dev.off()

# 2. TES vs TEAD1
cat("Creating TES vs TEAD1 Venn diagram...\n")
png(file.path(output_dir, "TES_vs_TEAD1_${ANALYSIS_TYPE}_venn.png"),
    width = 3000, height = 3000, res = 300) # Set resolution and dimensions
venn.plot2 <- draw.pairwise.venn(
    area1 = tes_total,
    area2 = tead1_total,
    cross.area = tes_tead1_intersection,
    category = c("TES", "TEAD1"),
    col = "transparent",
    fill = c("lightblue", "lightgreen"),
    alpha = 0.6,
    cex = 1.5,
    fontfamily = "serif",
    main = "TES vs TEAD1 $ANALYSIS_TYPE Peak Overlap",
    main.cex = 1.8
)
grid.draw(venn.plot2)
dev.off()

# 3. TESmut vs TEAD1
cat("Creating TESmut vs TEAD1 Venn diagram...\n")
png(file.path(output_dir, "TESmut_vs_TEAD1_${ANALYSIS_TYPE}_venn.png"),
    width = 3000, height = 3000, res = 300) # Set resolution and dimensions
venn.plot3 <- draw.pairwise.venn(
    area1 = tesmut_total,
    area2 = tead1_total,
    cross.area = tesmut_tead1_intersection,
    category = c("TESmut", "TEAD1"),
    col = "transparent",
    fill = c("lightcoral", "lightgreen"),
    alpha = 0.6,
    cex = 1.5,
    fontfamily = "serif",
    main = "TESmut vs TEAD1 $ANALYSIS_TYPE Peak Overlap",
    main.cex = 1.8
)
grid.draw(venn.plot3)
dev.off()

# Calculate detailed overlap statistics
cat("Calculating overlap statistics...\n")

# Create summary statistics table
stats_table <- data.frame(
    Comparison = c("TES vs TESmut", "TES vs TEAD1", "TESmut vs TEAD1"),
    Overlap = c(tes_tesmut_intersection, tes_tead1_intersection, tesmut_tead1_intersection),
    First_Only = c(tes_unique_tesmut, tes_unique_tead1, tesmut_unique_tead1),
    Second_Only = c(tesmut_unique_tes, tead1_unique_tes, tead1_unique_tesmut),
    First_Total = c(tes_total, tes_total, tesmut_total),
    Second_Total = c(tesmut_total, tead1_total, tead1_total),
    Overlap_Percent_First = round(c(tes_tesmut_intersection/tes_total*100,
                                   tes_tead1_intersection/tes_total*100,
                                   tesmut_tead1_intersection/tesmut_total*100), 1),
    Overlap_Percent_Second = round(c(tes_tesmut_intersection/tesmut_total*100,
                                    tes_tead1_intersection/tead1_total*100,
                                    tesmut_tead1_intersection/tead1_total*100), 1)
)

# Write statistics to file
write.table(stats_table, file.path(output_dir, "${ANALYSIS_TYPE}_peak_overlap_statistics.txt"), 
            sep="\t", row.names=FALSE, quote=FALSE)

cat("Venn diagram analysis completed.\n")
cat("Files created in:", file.path(output_dir), "\n")

# Print summary to console
print(stats_table)
EOF

# Run R script for Venn diagrams
echo "=== Creating Venn diagrams with R ==="
Rscript "$OUTDIR/create_venn_diagrams.R"

echo ""
echo "=== Creating signal comparison matrices ==="

# Check if BigWig files exist from previous analysis or create basic comparison
BIGWIG_DIR="results/06_bigwig"
if [ -d "$BIGWIG_DIR" ]; then
    
    # Find available BigWig files
    bigwig_files=""
    for condition in "${CONDITIONS[@]}"; do
        # Try combined BigWig first, then individual replicates
        if [ -f "$BIGWIG_DIR/${condition}.bw" ]; then
            bigwig_files="$bigwig_files $BIGWIG_DIR/${condition}.bw"
        elif [ -f "$BIGWIG_DIR/${condition}_combined.bw" ]; then
            bigwig_files="$bigwig_files $BIGWIG_DIR/${condition}_combined.bw"
        fi
    done
    
    if [ -n "$bigwig_files" ]; then
        echo "Creating signal comparison matrices for overlap regions..."
        
        # For each pairwise comparison, create matrices showing signals at overlap regions
        for comparison in "${comparisons[@]}"; do
            IFS=':' read -r cond1 cond2 <<< "$comparison"
            
            if [ -f "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond1}_overlapping_${cond2}.bed" ]; then
                echo "  Processing ${cond1} vs ${cond2} overlap regions..."
                
                computeMatrix reference-point \
                    --referencePoint center \
                    -b 2000 -a 2000 \
                    -R "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond1}_overlapping_${cond2}.bed" \
                    -S $bigwig_files \
                    --skipZeros \
                    -o "$OUTDIR/matrices/${cond1}_${cond2}_overlap_matrix.gz" \
                    -p 4 2>/dev/null || echo "    Warning: Matrix creation failed for ${cond1}_${cond2}"
            fi
        done
    fi
fi

echo ""
echo "=== Generating comprehensive summary ==="

# Create comprehensive summary file
SUMMARY_FILE="$OUTDIR/PAIRWISE_OVERLAP_SUMMARY.md"

cat > "$SUMMARY_FILE" << EOF
# Pairwise Binding Site Overlap Analysis Summary

**Generated on:** $(date)  
**Peak Type Used:** ${ANALYSIS_TYPE} peaks ($(echo ${ANALYSIS_TYPE} | tr '[:lower:]' '[:upper:]'))  
**Peak Files Source:** ${PEAK_DIR}

## Overview

This analysis compares binding sites between TES, TESmut, and TEAD1 conditions using pairwise overlap analysis with ${ANALYSIS_TYPE} peaks. $(if [ "$ANALYSIS_TYPE" = "consensus" ]; then echo "Consensus peaks represent the most reliable binding sites derived from replicate data."; elif [ "$ANALYSIS_TYPE" = "broad" ]; then echo "Broad peaks provide good overlap detection for Cut&Tag data by capturing extended binding domains."; else echo "Note: Consider using broad or consensus peaks for better overlap sensitivity."; fi)

## Peak Counts by Condition

EOF

for condition in "${CONDITIONS[@]}"; do
    if [ -f "${PEAK_FILES[$condition]}" ]; then
        peak_count=$(wc -l < "${PEAK_FILES[$condition]}")
        echo "- **$condition**: $peak_count peaks" >> "$SUMMARY_FILE"
    fi
done

cat >> "$SUMMARY_FILE" << EOF

## Pairwise Overlap Results

### TES vs TESmut
EOF

if [ -f "$OUTDIR/overlaps/TES_vs_TESmut/TES_overlapping_TESmut.bed" ]; then
    tes_overlap_tesmut=$(wc -l < "$OUTDIR/overlaps/TES_vs_TESmut/TES_overlapping_TESmut.bed")
    tesmut_overlap_tes=$(wc -l < "$OUTDIR/overlaps/TES_vs_TESmut/TESmut_overlapping_TES.bed")
    tes_unique=$(wc -l < "$OUTDIR/unique/TES_unique_vs_TESmut.bed")
    tesmut_unique=$(wc -l < "$OUTDIR/unique/TESmut_unique_vs_TES.bed")
    
    cat >> "$SUMMARY_FILE" << EOF
- TES peaks overlapping with TESmut: **$tes_overlap_tesmut**
- TESmut peaks overlapping with TES: **$tesmut_overlap_tes**
- TES-specific peaks (not in TESmut): **$tes_unique**
- TESmut-specific peaks (not in TES): **$tesmut_unique**

EOF
fi

cat >> "$SUMMARY_FILE" << EOF
### TES vs TEAD1
EOF

if [ -f "$OUTDIR/overlaps/TES_vs_TEAD1/TES_overlapping_TEAD1.bed" ]; then
    tes_overlap_tead1=$(wc -l < "$OUTDIR/overlaps/TES_vs_TEAD1/TES_overlapping_TEAD1.bed")
    tead1_overlap_tes=$(wc -l < "$OUTDIR/overlaps/TES_vs_TEAD1/TEAD1_overlapping_TES.bed")
    tes_unique_tead1=$(wc -l < "$OUTDIR/unique/TES_unique_vs_TEAD1.bed")
    tead1_unique_tes=$(wc -l < "$OUTDIR/unique/TEAD1_unique_vs_TES.bed")
    
    cat >> "$SUMMARY_FILE" << EOF
- TES peaks overlapping with TEAD1: **$tes_overlap_tead1**
- TEAD1 peaks overlapping with TES: **$tead1_overlap_tes**
- TES-specific peaks (not in TEAD1): **$tes_unique_tead1**
- TEAD1-specific peaks (not in TES): **$tead1_unique_tes**

EOF
fi

cat >> "$SUMMARY_FILE" << EOF
### TESmut vs TEAD1
EOF

if [ -f "$OUTDIR/overlaps/TESmut_vs_TEAD1/TESmut_overlapping_TEAD1.bed" ]; then
    tesmut_overlap_tead1=$(wc -l < "$OUTDIR/overlaps/TESmut_vs_TEAD1/TESmut_overlapping_TEAD1.bed")
    tead1_overlap_tesmut=$(wc -l < "$OUTDIR/overlaps/TESmut_vs_TEAD1/TEAD1_overlapping_TESmut.bed")
    tesmut_unique_tead1=$(wc -l < "$OUTDIR/unique/TESmut_unique_vs_TEAD1.bed")
    tead1_unique_tesmut=$(wc -l < "$OUTDIR/unique/TEAD1_unique_vs_TESmut.bed")
    
    cat >> "$SUMMARY_FILE" << EOF
- TESmut peaks overlapping with TEAD1: **$tesmut_overlap_tead1**
- TEAD1 peaks overlapping with TESmut: **$tead1_overlap_tesmut**
- TESmut-specific peaks (not in TEAD1): **$tesmut_unique_tead1**
- TEAD1-specific peaks (not in TESmut): **$tead1_unique_tesmut**

EOF
fi

cat >> "$SUMMARY_FILE" << EOF
## Files Generated

### Overlap Files
- 
`overlaps/`: Peak regions that overlap between conditions
- 
`unique/`: Condition-specific peaks with no overlap

### Visualization Files
- 
`venn/`: Venn diagrams showing overlap relationships
- 
`plots/`: Additional visualization plots

### Statistics Files
- 
`venn/${ANALYSIS_TYPE}_peak_overlap_statistics.txt
`: Detailed numerical overlap statistics
- 
`PAIRWISE_OVERLAP_SUMMARY.md
`: This comprehensive summary

### Analysis Files
- 
`matrices/`: Signal matrices for overlap regions (if BigWig files available)

## Interpretation Notes

1. **TES vs TESmut**: Comparison shows binding specificity of TES construct
2. **TES vs TEAD1**: Reveals similarity between TES and endogenous TEAD1 binding
3. **TESmut vs TEAD1**: Control comparison to validate TESmut background levels

## Next Steps

1. Review Venn diagrams for visual overlap assessment
2. Analyze signal intensity at overlap regions using provided matrices
3. Perform functional annotation of condition-specific peak sets
4. Consider motif analysis in overlap vs unique regions

EOF

echo "Files created:"
echo "- Overlap regions: $OUTDIR/overlaps/"
echo "- Unique regions: $OUTDIR/unique/"
echo "- Venn diagrams: $OUTDIR/venn/"
echo "- Summary report: $SUMMARY_FILE"
echo ""
# Optional: Add comparison with the other peak type if available
if [ "$BROAD_PEAKS_AVAILABLE" = true ] && [ "$NARROW_PEAKS_AVAILABLE" = true ]; then
    echo ""
    echo "=== Peak Type Comparison Available ==="
    
    if [ "$ANALYSIS_TYPE" = "broad" ]; then
        echo "Analysis used broad peaks. Narrow peaks also available for comparison."
        echo "Broad peaks typically show higher overlap rates due to extended binding domains."
    else
        echo "Analysis used narrow peaks. Broad peaks also available and recommended for Cut&Tag."
        echo "To run with broad peaks: sbatch 5b_broad_peak_calling.sh && sbatch 12_pairwise_overlap.sh"
    fi
    
    # Add a brief comparison to the summary file
    cat >> "$SUMMARY_FILE" << EOF

## Peak Type Notes

- **Current Analysis**: Used ${ANALYSIS_TYPE} peaks from ${PEAK_DIR}
$(if [ "$ANALYSIS_TYPE" = "broad" ]; then echo "- **Advantage**: Broad peaks capture extended binding domains, leading to more sensitive overlap detection"; else echo "- **Recommendation**: Consider running broad peak analysis (5b_broad_peak_calling.sh) for improved overlap sensitivity"; fi)
- **Alternative Available**: $(if [ "$ANALYSIS_TYPE" = "broad" ]; then echo "Narrow peaks available in results/05_peaks/"; else echo "Run 5b_broad_peak_calling.sh to generate broad peaks"; fi)

EOF
fi

echo ""
echo "=== Pairwise overlap analysis completed ==="
echo "Analysis type: $ANALYSIS_TYPE peaks"
echo "Results saved to: $OUTDIR"
echo "End time: $(date)"
