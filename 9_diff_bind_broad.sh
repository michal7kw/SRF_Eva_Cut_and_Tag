#!/bin/bash

#===============================================================================
# SCRIPT: 9b_diff_bind_broad.sh
# PURPOSE: Differential binding analysis for Cut&Tag data using DiffBind with BROAD PEAKS
#
# DESCRIPTION:
# This script performs statistical differential binding analysis between
# experimental conditions using the DiffBind R package with broad peaks. It identifies
# regions with significantly different binding patterns between TES, TESmut,
# and TEAD1 conditions, providing statistical validation of binding differences.
#
# KEY OPERATIONS:
# 1. Sample correlation analysis and quality assessment
# 2. Peak consensus set generation across replicates
# 3. Read count normalization and statistical modeling
# 4. Differential binding site identification
# 5. Generation of MA plots, PCA plots, and heatmaps
# 6. Export of significantly different binding sites
#
# METHODOLOGY:
# - Uses DiffBind R package for comprehensive analysis
# - Employs DESeq2 or edgeR for statistical testing
# - Performs TMM normalization for count data
# - Applies FDR correction for multiple testing
# - Generates consensus peak sets from replicates
#
# IMPORTANT PARAMETERS:
# - Memory: 32GB (required for large matrix operations)
# - Time: 4 hours (sufficient for statistical analysis)
# - Threads: 16 (parallel processing for count operations)
# - FDR threshold: typically 0.05 (set in R script)
# - Fold change threshold: typically 2-fold (set in R script)
#
# INPUTS:
# - Sample sheet with BAM files and broad peak files (diffbind_samplesheet_broad.csv)
# - Filtered BAM files from alignment pipeline
# - Broad peak files from MACS2 analysis
# - Configuration parameters in R script
#
# OUTPUTS:
# - Differential binding results (CSV files)
# - Quality control plots (PCA, correlation heatmaps)
# - MA plots showing fold changes and significance
# - Venn diagrams of overlapping differential sites
# - Normalized count matrices
#
# DEPENDENCIES:
# - R with DiffBind, DESeq2, and visualization packages
# - Conda environment: diffbind_analysis
# - Sample sheet configuration file for broad peaks
#
# USAGE:
# sbatch 9b_diff_bind_broad.sh
#
# NOTES:
# - Requires properly formatted sample sheet for broad peaks
# - Generates publication-ready statistical plots
# - Provides rigorous statistical validation of binding differences
# - Results inform biological interpretation of Cut&Tag data
# - Uses BROAD PEAKS for enhanced sensitivity in differential analysis
#===============================================================================

#SBATCH --job-name=9_diff_bind_broad
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/9_diff_bind_broad.err"
#SBATCH --output="./logs/9_diff_bind_broad.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate diffbind_analysis

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva"

# Run DiffBind analysis for broad peaks
Rscript ${BASE_DIR}/scripts/diffbind_analysis_broad.R

echo "Differential binding analysis with broad peaks complete!"