#!/bin/bash

#===============================================================================
# SCRIPT: 11b_combine_replicates_broad.sh
# PURPOSE: Combine biological replicates for average signal analysis using BROAD PEAKS
#
# DESCRIPTION:
# This script combines biological replicates to create average signal tracks
# and consensus peak sets for improved signal-to-noise ratio and statistical
# power using broad peaks. It merges BAM files, generates average BigWig tracks, creates
# high-confidence consensus peaks, and performs cross-condition signal analysis
# to compare binding patterns across different experimental conditions.
#
# KEY OPERATIONS:
# 1. BAM file merging per condition (TES, TESmut, TEAD1)
# 2. Average BigWig track generation from merged BAMs
# 3. Consensus peak set creation (peaks in ≥2 replicates)
# 4. Signal matrix computation around consensus peaks
# 5. Cross-condition signal analysis and visualization
# 6. Summary statistics generation
#
# METHODOLOGY:
# - Uses samtools for BAM merging and indexing
# - Employs bamCoverage for normalized average signal tracks
# - Applies bedtools for peak merging and consensus calling
# - Utilizes deepTools for signal matrix computation
# - Creates heatmaps and profile plots for visualization
# - Implements stringent consensus criteria (≥2/3 replicates)
#
# IMPORTANT PARAMETERS:
# - Memory: 16GB (for large BAM file operations)
# - Time: 2 hours (adequate for merging and analysis)
# - Threads: 8 (parallel processing for efficiency)
# - Bin size: 10bp (high-resolution signal tracks)
# - Window: ±3kb around peak centers
# - Consensus threshold: ≥2 replicates
#
# INPUTS:
# - Filtered BAM files from step 4 (per replicate)
# - Individual broad peak files from step 5b (broadPeak format)
# - BigWig files from step 6 (per replicate)
# - Three biological replicates per condition
#
# OUTPUTS:
# - Merged BAM files per condition
# - Average BigWig tracks (CPM normalized)
# - Consensus peak sets (high-confidence)
# - Signal matrices around consensus peaks
# - Cross-condition heatmaps and profiles
# - Comprehensive summary statistics
#
# DEPENDENCIES:
# - samtools for BAM manipulation
# - bamCoverage (deepTools) for signal tracks
# - bedtools for genomic interval operations
# - computeMatrix, plotHeatmap (deepTools) for visualization
# - Conda environment: combine_results
#
# USAGE:
# sbatch 11b_combine_replicates_broad.sh
#
# NOTES:
# - Improves signal quality through replicate averaging
# - Increases statistical confidence in peak calls
# - Enables robust cross-condition comparisons
# - Essential for publication-quality analysis
# - Provides foundation for downstream comparative studies
# - Uses BROAD PEAKS for enhanced consensus calling
#===============================================================================

#SBATCH --job-name=combine_replicates_broad
#SBATCH --output=./logs/11b_combine_replicates_broad.out
#SBATCH --error=./logs/11b_combine_replicates_broad.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00

echo "=== Combining Replicates for Average Signal Analysis (BROAD PEAKS) ==="
echo "Start time: $(date)"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate combine_results

# Set up directories
OUTDIR="results/11_combined_replicates_broad"
mkdir -p "$OUTDIR"/{bigwig,peaks,matrices}

# Check required files
REQUIRED_DIRS=("results/04_filtered" "results/05_peaks_broad" "results/06_bigwig")
for dir in "${REQUIRED_DIRS[@]}"; do
    if [ ! -d "$dir" ]; then
        echo "Error: Required directory $dir not found"
        exit 1
    fi
done

echo "=== Combining BAM files for average signal computation ==="

# Define conditions and replicates
declare -A CONDITIONS=(
    ["TES"]="TES-1 TES-2 TES-3"
    ["TESmut"]="TESmut-1 TESmut-2 TESmut-3"
    ["TEAD1"]="TEAD1-1 TEAD1-2 TEAD1-3"
)

# Combine BAM files per condition for average signal
for condition in "${!CONDITIONS[@]}"; do
    echo "Processing $condition condition..."

    # Get replicate names
    replicates=(${CONDITIONS[$condition]})

    # Check if BAM files exist
    bam_files=""
    for replicate in "${replicates[@]}"; do
        bam_file="results/04_filtered/${replicate}_filtered.bam"
        if [ -f "$bam_file" ]; then
            bam_files="$bam_files $bam_file"
            echo "  Found: $bam_file"
        else
            echo "  Warning: Missing $bam_file"
        fi
    done

    if [ -n "$bam_files" ]; then
        echo "  Merging BAM files for $condition..."
        samtools merge -f -@ 8 "$OUTDIR/${condition}_merged.bam" $bam_files

        # Index merged BAM
        echo "  Indexing merged BAM for $condition..."
        samtools index "$OUTDIR/${condition}_merged.bam"

        # Generate average BigWig from merged BAM
        echo "  Creating average BigWig for $condition..."
        bamCoverage \
            -b "$OUTDIR/${condition}_merged.bam" \
            -o "$OUTDIR/bigwig/${condition}_average.bw" \
            --binSize 10 \
            --normalizeUsing CPM \
            --effectiveGenomeSize 2913022398 \
            --numberOfProcessors 8 \
            --extendReads

        echo "  Completed $condition average signal track"

        # Calculate read counts for normalization info
        total_reads=$(samtools view -c "$OUTDIR/${condition}_merged.bam")
        echo "  $condition total reads: $total_reads"

    else
        echo "  Error: No BAM files found for $condition"
    fi
done

echo "=== Creating consensus peak sets from broad peaks ==="

# Use broad peaks for consensus calling
for condition in "${!CONDITIONS[@]}"; do
    echo "Processing consensus broad peaks for $condition..."

    # Get individual replicate peak files
    replicates=(${CONDITIONS[$condition]})
    peak_files=""
    existing_peaks=0

    for replicate in "${replicates[@]}"; do
        peak_file="results/05_peaks_broad/${replicate}_peaks.broadPeak"
        if [ -f "$peak_file" ]; then
            peak_files="$peak_files $peak_file"
            existing_peaks=$((existing_peaks + 1))
            echo "  Found: $peak_file"
        fi
    done

    if [ $existing_peaks -ge 2 ]; then
        echo "  Creating high-confidence consensus broad peaks for $condition..."

        # Merge all individual peaks and find overlaps
        cat $peak_files | sort -k1,1 -k2,2n > "$OUTDIR/peaks/${condition}_all_peaks.bed"

        # Use bedtools to find peaks present in at least 2 out of 3 replicates
        bedtools merge -i "$OUTDIR/peaks/${condition}_all_peaks.bed" -c 4 -o count > "$OUTDIR/peaks/${condition}_merged_raw.bed"

        # Filter for high-confidence peaks (present in >= 2 replicates)
        awk '$4 >= 2' "$OUTDIR/peaks/${condition}_merged_raw.bed" > "$OUTDIR/peaks/${condition}_consensus_peaks.bed"

        # Copy existing combined peak file as well
        if [ -f "results/05_peaks_broad/${condition}_peaks.broadPeak" ]; then
            cp "results/05_peaks_broad/${condition}_peaks.broadPeak" "$OUTDIR/peaks/${condition}_combined_peaks.broadPeak"
        fi

        consensus_count=$(wc -l < "$OUTDIR/peaks/${condition}_consensus_peaks.bed")
        echo "  $condition consensus broad peaks (≥2 replicates): $consensus_count"

    else
        echo "  Warning: Insufficient peak files for $condition consensus"
    fi
done

echo "=== Creating signal matrices around consensus broad peaks ==="

# Create signal matrices for each condition using their own consensus peaks
for condition in "${!CONDITIONS[@]}"; do
    if [ -f "$OUTDIR/peaks/${condition}_consensus_peaks.bed" ] && [ -f "$OUTDIR/bigwig/${condition}_average.bw" ]; then
        echo "Computing signal matrix for $condition around its consensus broad peaks..."

        computeMatrix reference-point \
            --referencePoint center \
            -b 3000 -a 3000 \
            -R "$OUTDIR/peaks/${condition}_consensus_peaks.bed" \
            -S "$OUTDIR/bigwig/${condition}_average.bw" \
            --skipZeros \
            -o "$OUTDIR/matrices/${condition}_self_matrix.gz" \
            --outFileSortedRegions "$OUTDIR/matrices/${condition}_self_regions.bed" \
            -p 8
    fi
done

echo "=== Cross-condition signal analysis ==="

# Create matrices showing each condition's signal on all other conditions' peaks
conditions=("TES" "TESmut" "TEAD1")
for peak_condition in "${conditions[@]}"; do
    if [ -f "$OUTDIR/peaks/${peak_condition}_consensus_peaks.bed" ]; then

        # Collect all BigWig files
        bigwig_files=""
        bigwig_labels=""

        for signal_condition in "${conditions[@]}"; do
            if [ -f "$OUTDIR/bigwig/${signal_condition}_average.bw" ]; then
                bigwig_files="$bigwig_files $OUTDIR/bigwig/${signal_condition}_average.bw"
                bigwig_labels="$bigwig_labels ${signal_condition}_avg"
            fi
        done

        if [ -n "$bigwig_files" ]; then
            echo "Computing cross-condition matrix for ${peak_condition} broad peaks..."

            computeMatrix reference-point \
                --referencePoint center \
                -b 3000 -a 3000 \
                -R "$OUTDIR/peaks/${peak_condition}_consensus_peaks.bed" \
                -S $bigwig_files \
                --samplesLabel $bigwig_labels \
                --skipZeros \
                -o "$OUTDIR/matrices/${peak_condition}_peaks_all_signals.gz" \
                --outFileSortedRegions "$OUTDIR/matrices/${peak_condition}_peaks_all_regions.bed" \
                -p 8

            # Create heatmap visualization
            plotHeatmap \
                -m "$OUTDIR/matrices/${peak_condition}_peaks_all_signals.gz" \
                -out "$OUTDIR/${peak_condition}_peaks_signal_heatmap.pdf" \
                --plotTitle "${peak_condition} broad peaks: Signal from all conditions" \
                --whatToShow 'heatmap and colorbar' \
                --colorMap RdBu_r \
                --zMin -2 --zMax 2

            # Convert heatmap PDF to PNG
            convert -density 300 "$OUTDIR/${peak_condition}_peaks_signal_heatmap.pdf" "$OUTDIR/${peak_condition}_peaks_signal_heatmap.png"

            # Create profile plot
            plotProfile \
                -m "$OUTDIR/matrices/${peak_condition}_peaks_all_signals.gz" \
                -out "$OUTDIR/${peak_condition}_peaks_signal_profile.pdf" \
                --plotTitle "${peak_condition} broad peaks: Average signal profiles" \
                --plotType lines \
                --perGroup

            # Convert profile PDF to PNG
            convert -density 300 "$OUTDIR/${peak_condition}_peaks_signal_profile.pdf" "$OUTDIR/${peak_condition}_peaks_signal_profile.png"
        fi
    fi
done

echo "=== Generating summary statistics ==="

# Create summary file
SUMMARY_FILE="$OUTDIR/REPLICATE_COMBINATION_SUMMARY_BROAD.txt"
echo "Cut&Tag Replicate Combination Analysis Summary (BROAD PEAKS)" > "$SUMMARY_FILE"
echo "Generated on: $(date)" >> "$SUMMARY_FILE"
echo "=============================================" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

for condition in "${!CONDITIONS[@]}"; do
    echo "=== $condition Condition ===" >> "$SUMMARY_FILE"

    # Read counts from merged BAM
    if [ -f "$OUTDIR/${condition}_merged.bam" ]; then
        reads=$(samtools view -c "$OUTDIR/${condition}_merged.bam")
        echo "Total reads in merged BAM: $reads" >> "$SUMMARY_FILE"
    fi

    # Consensus peak counts
    if [ -f "$OUTDIR/peaks/${condition}_consensus_peaks.bed" ]; then
        consensus_peaks=$(wc -l < "$OUTDIR/peaks/${condition}_consensus_peaks.bed")
        echo "Consensus broad peaks (≥2 replicates): $consensus_peaks" >> "$SUMMARY_FILE"
    fi

    # Combined peaks (original)
    if [ -f "$OUTDIR/peaks/${condition}_combined_peaks.broadPeak" ]; then
        combined_peaks=$(wc -l < "$OUTDIR/peaks/${condition}_combined_peaks.broadPeak")
        echo "Combined broad peaks (original): $combined_peaks" >> "$SUMMARY_FILE"
    fi

    echo "" >> "$SUMMARY_FILE"
done

echo "Files created:"
echo "- Average BigWig tracks: $OUTDIR/bigwig/"
echo "- Consensus broad peak sets: $OUTDIR/peaks/"
echo "- Signal matrices: $OUTDIR/matrices/"
echo "- Cross-condition heatmaps: $OUTDIR/*_heatmap.pdf"
echo "- Summary report: $SUMMARY_FILE"

echo "=== Broad peak replicate combination completed ==="
echo "End time: $(date)"