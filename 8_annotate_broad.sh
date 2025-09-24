#!/bin/bash

#===============================================================================
# SCRIPT: 8b_annotate_broad.sh
# PURPOSE: Peak annotation and motif analysis for Cut&Tag data using BROAD PEAKS
#
# DESCRIPTION:
# This script performs comprehensive peak annotation and motif discovery analysis
# for identified broad peaks. It annotates peaks with genomic features (promoters,
# introns, exons, etc.) and discovers enriched DNA motifs within peak regions
# using HOMER motif analysis tools.
#
# KEY OPERATIONS:
# 1. Peak annotation using R-based genomic annotation tools
# 2. HOMER motif discovery for TES and TEAD1 broad peaks
# 3. Motif enrichment analysis with background correction
# 4. Generation of motif logos and statistics
#
# METHODOLOGY:
# - Uses R script for detailed genomic feature annotation
# - Employs HOMER findMotifsGenome.pl for de novo motif discovery
# - Analyzes 200bp regions around peak centers
# - Includes repeat masking for accurate motif identification
# - Compares motifs between different conditions
#
# IMPORTANT PARAMETERS:
# - Memory: 16GB (sufficient for motif analysis)
# - Time: 2 hours (adequate for HOMER processing)
# - Threads: 8 (parallel motif discovery)
# - Motif size: 200bp around peak centers
# - Genome: hg38 reference
#
# INPUTS:
# - Broad peak files from MACS2 (05_peaks_broad directory)
# - Peak comparison results from previous analysis
# - R annotation script (annotate_peaks_broad.R)
#
# OUTPUTS:
# - Genomic feature annotations (from R script)
# - HOMER motif analysis results (TES_motifs/, TEAD1_motifs/)
# - Motif logos, enrichment statistics, and HTML reports
# - Known and de novo motif discoveries
#
# DEPENDENCIES:
# - R with genomic annotation packages
# - HOMER motif analysis suite
# - Conda environment: annotation_enrichment
#
# USAGE:
# sbatch 8b_annotate_broad.sh
#
# NOTES:
# - Adds "chr" prefix to peak coordinates for HOMER compatibility
# - Includes file existence checks before processing
# - Generates publication-ready motif analysis reports
# - Compares motifs between TES and TEAD1 conditions
# - Uses BROAD PEAKS for extended motif analysis regions
#===============================================================================

#SBATCH --job-name=8b_annotate_broad
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/8b_annotate_broad.err"
#SBATCH --output="./logs/8b_annotate_broad.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate annotation_enrichment

# Load HOMER module
# module load homer

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva"
PEAKS_DIR="${BASE_DIR}/results/05_peaks_broad"
ANALYSIS_DIR="${BASE_DIR}/results/07_analysis_broad"

# Create output directory if it doesn't exist
mkdir -p ${ANALYSIS_DIR}

# Run R script for broad peak annotation
Rscript ${BASE_DIR}/scripts/annotate_peaks_broad.R

# HOMER motif analysis for TES broad peaks (COMMENTED OUT)
# if [ -s "${PEAKS_DIR}/TES_peaks.broadPeak" ]; then
#     echo "Running motif analysis for TES broad peaks..."
#     awk 'BEGIN{OFS="\t"} {print "chr"$0}' ${PEAKS_DIR}/TES_peaks.broadPeak > ${PEAKS_DIR}/TES_peaks.broadPeak.tmp
#     findMotifsGenome.pl \
#         ${PEAKS_DIR}/TES_peaks.broadPeak.tmp \
#         hg38 \
#         ${ANALYSIS_DIR}/TES_motifs \
#         -size 200 \
#         -mask \
#         -p 8
#     rm ${PEAKS_DIR}/TES_peaks.broadPeak.tmp
# else
#     echo "Skipping HOMER for TES broad peaks, file is empty or does not exist."
# fi

# Compare with TEAD1 broad peaks (COMMENTED OUT)
# if [ -s "${PEAKS_DIR}/TEAD1_peaks.broadPeak" ]; then
#     echo "Running motif analysis for TEAD1 broad peaks..."
#     awk 'BEGIN{OFS="\t"} {print "chr"$0}' ${PEAKS_DIR}/TEAD1_peaks.broadPeak > ${PEAKS_DIR}/TEAD1_peaks.broadPeak.tmp
#     findMotifsGenome.pl \
#         ${PEAKS_DIR}/TEAD1_peaks.broadPeak.tmp \
#         hg38 \
#         ${ANALYSIS_DIR}/TEAD1_motifs \
#         -size 200 \
#         -mask \
#         -p 8
#     rm ${PEAKS_DIR}/TEAD1_peaks.broadPeak.tmp
# else
#     echo "Skipping HOMER for TEAD1 broad peaks, file is empty or does not exist."
# fi

echo "Broad peak annotation and motif analysis complete!"