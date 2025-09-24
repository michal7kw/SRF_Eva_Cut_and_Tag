#!/bin/bash

#===============================================================================
# SCRIPT: 9_diff_bind.sh
# PURPOSE: Differential binding analysis for Cut&Tag data using DiffBind
#
# DESCRIPTION:
# This script performs statistical differential binding analysis between
# experimental conditions using the DiffBind R package. It identifies
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
# - Sample sheet with BAM files and peak files (diffbind_samplesheet.csv)
# - Filtered BAM files from alignment pipeline
# - Peak files from MACS2 analysis
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
# - Sample sheet configuration file
#
# USAGE:
# sbatch 9_diff_bind.sh
#
# NOTES:
# - Requires properly formatted sample sheet
# - Generates publication-ready statistical plots
# - Provides rigorous statistical validation of binding differences
# - Results inform biological interpretation of Cut&Tag data
#===============================================================================

#SBATCH --job-name=9_diffbind
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/9_diffbind.err"
#SBATCH --output="./logs/9_diffbind.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate diffbind_analysis

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva"

# Run DiffBind analysis
Rscript ${BASE_DIR}/scripts/diffbind_analysis.R

echo "Differential binding analysis complete!"