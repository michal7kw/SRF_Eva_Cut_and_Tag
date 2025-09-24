#!/bin/bash

#===============================================================================
# SCRIPT: 5b_broad_peak_calling.sh
# PURPOSE: Broad peak calling for Cut&Tag data using MACS2 for overlap analysis
#
# DESCRIPTION:
# This script performs broad peak calling on filtered BAM files from Cut&Tag
# experiments using MACS2 with broad peak mode. Broad peaks are particularly
# useful for overlap analysis as they capture broader binding domains and
# provide better sensitivity for detecting shared binding regions between
# conditions. This complements the narrow peak analysis from step 5.
#
# KEY OPERATIONS:
# 1. BAM file validation and indexing
# 2. Automatic detection of paired-end vs single-end sequencing
# 3. MACS2 broad peak calling with relaxed parameters
# 4. Processing of combined replicates for main overlap analysis
# 5. Peak count summary and comparison with narrow peaks
# 6. Optional comparative diagnostic plot generation
#
# METHODOLOGY:
# - Uses MACS2 callpeak with --broad flag for broad peak detection
# - Applies relaxed q-value threshold of 0.05 (vs 0.01 for narrow peaks)
# - Uses broad-cutoff of 0.1 for defining broad peak boundaries
# - Uses appropriate control samples (IggMs for TES/TESmut, IggRb for TEAD1)
# - Keeps all duplicate reads (--keep-dup all) as appropriate for Cut&Tag
# - Automatically selects BAMPE (paired-end) or BAM (single-end) format
#
# IMPORTANT PARAMETERS:
# - Memory: 16GB (sufficient for broad peak calling)
# - Time: 4 hours (adequate for multiple samples)
# - Threads: 8 (MACS2 is mostly single-threaded but allows parallel processing)
# - Q-value threshold: 0.05 (more permissive than narrow peaks)
# - Broad cutoff: 0.1 (threshold for broad peak boundaries)
# - Genome: hs (human genome size for MACS2)
#
# EXPECTED INPUTS:
# - Filtered BAM files from step 4 (04_filtered directory)
# - BAM files should be properly indexed
# - Control samples: IggMs_filtered.bam, IggRb_filtered.bam
# - Treatment samples: TES-*, TESmut-*, TEAD1-* replicates
#
# EXPECTED OUTPUTS:
# - *_peaks.broadPeak: Broad peak coordinates in ENCODE broadPeak format
# - *_peaks.gappedPeak: Gapped peak regions (if generated)
# - *_summits.bed: Peak summit positions
# - *_model.r: R script for MACS2 model visualization
# - *_peaks.xls: Detailed peak information in Excel format
# - *_macs2_broad.log: MACS2 execution logs for troubleshooting
# - broad_vs_narrow_peaks_comparison.pdf: Optional comparative plot
#
# DEPENDENCIES:
# - MACS2 (peak calling algorithm)
# - samtools (BAM file operations)
# - R and ggplot2 (optional, for comparative plots)
# - Conda environment: peak_calling_new
#
# USAGE:
# sbatch 5b_broad_peak_calling.sh
#
# NOTES:
# - Broad peaks are better suited for overlap analysis than narrow peaks
# - Uses more permissive parameters to capture broader binding domains
# - Focuses on combined replicates for main analysis
# - Provides comparison with narrow peaks from step 5
# - Error handling includes detailed logging and validation
#===============================================================================

#SBATCH --job-name=5b_broad_peaks
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/5b_broad_peaks.err"
#SBATCH --output="./logs/5b_broad_peaks.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva"
BAM_DIR="${BASE_DIR}/results/04_filtered"
OUTPUT_DIR="${BASE_DIR}/results/05_peaks_broad"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

echo "=========================================="
echo "Starting broad peak calling analysis"
echo "Purpose: Better suited for overlap analysis"
echo "Date: $(date)"
echo "=========================================="

# Function to check BAM file
check_bam_file() {
    local bam_file=$1
    local sample_name=$(basename $bam_file .bam)
    
    echo "Checking $sample_name..."
    
    # Check if file exists
    if [ ! -f "$bam_file" ]; then
        echo "ERROR: File $bam_file does not exist!"
        return 1
    fi
    
    # Check file size
    local file_size=$(ls -lh $bam_file | awk '{print $5}')
    echo "  File size: $file_size"
    
    # Count total reads
    local read_count=$(samtools view -c $bam_file 2>/dev/null)
    echo "  Total reads: $read_count"
    
    if [ "$read_count" -eq 0 ]; then
        echo "  WARNING: BAM file is empty!"
        return 1
    fi
    
    # Check for paired reads
    local paired_count=$(samtools view -f 1 -c $bam_file 2>/dev/null)
    local paired_pct=0
    if [ "$read_count" -gt 0 ]; then
        paired_pct=$((paired_count * 100 / read_count))
    fi
    echo "  Paired reads: $paired_count ($paired_pct%)"
    
    # Return paired status (0 = paired, 1 = single)
    if [ "$paired_pct" -gt 50 ]; then
        echo "  Status: PAIRED-END"
        return 0
    else
        echo "  Status: SINGLE-END"
        return 2
    fi
}

# Function to index BAM file if needed
index_bam_if_needed() {
    local bam_file=$1
    if [ ! -f "${bam_file}.bai" ]; then
        echo "  Creating index for $(basename $bam_file)..."
        samtools index $bam_file
    fi
}

echo ""
echo "Step 1: Checking and indexing BAM files"
echo "------------------------------------------"

# Check all BAM files and determine read type
PAIRED_END=true
for bam in ${BAM_DIR}/*.bam; do
    check_bam_file $bam
    status=$?
    if [ $status -eq 1 ]; then
        echo "ERROR: Problem with $bam"
        exit 1
    elif [ $status -eq 2 ]; then
        PAIRED_END=false
    fi
    index_bam_if_needed $bam
    echo ""
done

# Set format based on paired-end status
if [ "$PAIRED_END" = true ]; then
    FORMAT="BAMPE"
    echo "Using PAIRED-END mode (BAMPE) for MACS2"
else
    FORMAT="BAM"
    echo "Using SINGLE-END mode (BAM) for MACS2"
fi

echo ""
echo "Step 2: Running MACS2 broad peak calling"
echo "------------------------------------------"

# Function to run MACS2 with broad peak calling and error handling
run_macs2_broad() {
    local treatment_files="$1"
    local control_file="$2"
    local output_name="$3"
    local format="$4"
    
    echo "Calling broad peaks for $output_name..."
    echo "  Treatment: $treatment_files"
    echo "  Control: $control_file"
    echo "  Format: $format"
    
    macs2 callpeak \
        -t $treatment_files \
        -c $control_file \
        -f $format \
        -g hs \
        -n $output_name \
        --outdir ${OUTPUT_DIR} \
        -q 0.05 \
        --broad \
        --broad-cutoff 0.1 \
        --keep-dup all \
        2> ${OUTPUT_DIR}/${output_name}_macs2_broad.log
    
    if [ $? -eq 0 ]; then
        echo "  SUCCESS: Broad peak calling completed for $output_name"
        # Count peaks
        if [ -f "${OUTPUT_DIR}/${output_name}_peaks.broadPeak" ]; then
            peak_count=$(wc -l < "${OUTPUT_DIR}/${output_name}_peaks.broadPeak")
            echo "  Found $peak_count broad peaks"
        fi
        # Also report gapped peaks if generated
        if [ -f "${OUTPUT_DIR}/${output_name}_peaks.gappedPeak" ]; then
            gapped_count=$(wc -l < "${OUTPUT_DIR}/${output_name}_peaks.gappedPeak")
            echo "  Found $gapped_count gapped peaks"
        fi
    else
        echo "  ERROR: Broad peak calling failed for $output_name"
        echo "  Check log file: ${OUTPUT_DIR}/${output_name}_macs2_broad.log"
        tail -20 ${OUTPUT_DIR}/${output_name}_macs2_broad.log
    fi
    echo ""
}

# Broad peak calling for combined replicates (primary analysis for overlaps)
echo "Processing combined replicates for broad peaks..."
echo "=================================================="

# TES samples
run_macs2_broad \
    "${BAM_DIR}/TES-1_filtered.bam ${BAM_DIR}/TES-2_filtered.bam ${BAM_DIR}/TES-3_filtered.bam" \
    "${BAM_DIR}/IggMs_filtered.bam" \
    "TES" \
    "$FORMAT"

# TESmut samples
run_macs2_broad \
    "${BAM_DIR}/TESmut-1_filtered.bam ${BAM_DIR}/TESmut-2_filtered.bam ${BAM_DIR}/TESmut-3_filtered.bam" \
    "${BAM_DIR}/IggMs_filtered.bam" \
    "TESmut" \
    "$FORMAT"

# TEAD1 samples
run_macs2_broad \
    "${BAM_DIR}/TEAD1-1_filtered.bam ${BAM_DIR}/TEAD1-2_filtered.bam ${BAM_DIR}/TEAD1-3_filtered.bam" \
    "${BAM_DIR}/IggRb_filtered.bam" \
    "TEAD1" \
    "$FORMAT"

echo ""
echo "Processing individual samples for broad peaks..."
echo "================================================="

# Individual TES samples
for i in 1 2 3; do
    run_macs2_broad \
        "${BAM_DIR}/TES-${i}_filtered.bam" \
        "${BAM_DIR}/IggMs_filtered.bam" \
        "TES-${i}" \
        "$FORMAT"
done

# Individual TESmut samples
for i in 1 2 3; do
    run_macs2_broad \
        "${BAM_DIR}/TESmut-${i}_filtered.bam" \
        "${BAM_DIR}/IggMs_filtered.bam" \
        "TESmut-${i}" \
        "$FORMAT"
done

# Individual TEAD1 samples
for i in 1 2 3; do
    run_macs2_broad \
        "${BAM_DIR}/TEAD1-${i}_filtered.bam" \
        "${BAM_DIR}/IggRb_filtered.bam" \
        "TEAD1-${i}" \
        "$FORMAT"
done

echo ""
echo "Step 3: Summary of broad peak results"
echo "------------------------------------------"

# Summary report
echo "Broad peak files generated:"
ls -lh ${OUTPUT_DIR}/*_peaks.broadPeak 2>/dev/null | awk '{print "  "$9": "$5" ("$1")"}'

echo ""
echo "Broad peak counts:"
for peak_file in ${OUTPUT_DIR}/*_peaks.broadPeak; do
    if [ -f "$peak_file" ]; then
        sample=$(basename $peak_file _peaks.broadPeak)
        count=$(wc -l < $peak_file)
        printf "  %-20s: %6d broad peaks\n" "$sample" "$count"
    fi
done

echo ""
echo "Gapped peak counts (if available):"
for peak_file in ${OUTPUT_DIR}/*_peaks.gappedPeak; do
    if [ -f "$peak_file" ]; then
        sample=$(basename $peak_file _peaks.gappedPeak)
        count=$(wc -l < $peak_file)
        printf "  %-20s: %6d gapped peaks\n" "$sample" "$count"
    fi
done

echo ""
echo "Log files:"
ls -lh ${OUTPUT_DIR}/*.log 2>/dev/null | awk '{print "  "$9}'

echo ""
echo "Comparison with narrow peaks:"
echo "-----------------------------"
NARROW_DIR="${BASE_DIR}/results/05_peaks_narrow"
if [ -d "$NARROW_DIR" ]; then
    echo "Narrow vs Broad peak comparison:"
    for condition in TES TESmut TEAD1; do
        narrow_file="${NARROW_DIR}/${condition}_peaks.narrowPeak"
        broad_file="${OUTPUT_DIR}/${condition}_peaks.broadPeak"
        
        if [ -f "$narrow_file" ] && [ -f "$broad_file" ]; then
            narrow_count=$(wc -l < "$narrow_file")
            broad_count=$(wc -l < "$broad_file")
            printf "  %-8s: %6d narrow -> %6d broad peaks\n" "$condition" "$narrow_count" "$broad_count"
        fi
    done
fi

echo ""
echo "=========================================="
echo "Broad peak calling analysis completed!"
echo "These peaks are better suited for overlap analysis"
echo "as they capture broader binding domains."
echo "Date: $(date)"
echo "=========================================="

# Optional: Generate a comparative diagnostic plot if R is available
if command -v R &> /dev/null; then
    echo ""
    echo "Generating comparative diagnostic plots..."
    Rscript -e "
    library(ggplot2, quietly=TRUE)
    
    # Read broad peak counts
    broad_files <- list.files('${OUTPUT_DIR}', pattern='_peaks.broadPeak$', full.names=TRUE)
    narrow_files <- list.files('${NARROW_DIR}', pattern='_peaks.narrowPeak$', full.names=TRUE)
    
    if(length(broad_files) > 0) {
        broad_counts <- sapply(broad_files, function(x) length(readLines(x)))
        names(broad_counts) <- gsub('_peaks.broadPeak', '', basename(names(broad_counts)))
        
        narrow_counts <- rep(0, length(broad_counts))
        names(narrow_counts) <- names(broad_counts)
        
        if(length(narrow_files) > 0) {
            narrow_temp <- sapply(narrow_files, function(x) length(readLines(x)))
            narrow_names <- gsub('_peaks.narrowPeak', '', basename(names(narrow_temp)))
            narrow_counts[narrow_names] <- narrow_temp
        }
        
        # Create comparison plot
        pdf('${OUTPUT_DIR}/broad_vs_narrow_peaks_comparison.pdf', width=12, height=6)
        
        # Prepare data for plotting
        plot_data <- data.frame(
            Sample = rep(names(broad_counts), 2),
            Count = c(narrow_counts, broad_counts),
            Type = rep(c('Narrow', 'Broad'), each=length(broad_counts))
        )
        
        # Create barplot
        p <- ggplot(plot_data, aes(x=Sample, y=Count, fill=Type)) +
            geom_bar(stat='identity', position='dodge') +
            theme(axis.text.x = element_text(angle=45, hjust=1)) +
            labs(title='Peak Count Comparison: Narrow vs Broad Peaks',
                 subtitle='Broad peaks should show better overlap in downstream analysis',
                 x='Sample', y='Number of Peaks') +
            scale_fill_manual(values=c('Narrow'='lightblue', 'Broad'='darkblue'))
        
        print(p)
        dev.off()
        cat('Comparative diagnostic plot saved to ${OUTPUT_DIR}/broad_vs_narrow_peaks_comparison.pdf\n')
    }
    " 2>/dev/null || echo "R plotting skipped (optional)"
fi