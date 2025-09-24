#!/bin/bash

#===============================================================================
# SCRIPT: run_additional_analyses.sh
# PURPOSE: Execute additional Cut&Tag quality control and analysis scripts
#
# DESCRIPTION:
# This script runs supplementary analyses that are critical for Cut&Tag
# data quality assessment and biological interpretation. It includes:
# 1. Fragment size analysis
# 2. Replicate reproducibility (IDR) analysis
# 3. Enhanced GO enrichment analysis
#===============================================================================

#SBATCH --job-name=additional_analyses
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/additional_analyses.err"
#SBATCH --output="./logs/additional_analyses.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate annotation_enrichment

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva"
cd ${BASE_DIR}

echo "=========================================="
echo "Running Additional Cut&Tag Analyses"
echo "Date: $(date)"
echo "=========================================="

# 1. Fragment size analysis
echo ""
echo "Step 1: Fragment size analysis"
echo "--------------------------------"
if [ -f "scripts/fragment_analysis.R" ]; then
    echo "Running fragment size analysis..."
    Rscript scripts/fragment_analysis.R
    if [ $? -eq 0 ]; then
        echo "Fragment size analysis completed successfully"
    else
        echo "Warning: Fragment size analysis encountered errors"
    fi
else
    echo "Warning: fragment_analysis.R script not found"
fi

# 2. Replicate reproducibility analysis
echo ""
echo "Step 2: Replicate reproducibility (IDR) analysis"
echo "------------------------------------------------"
if [ -f "scripts/idr_analysis.R" ]; then
    echo "Running IDR analysis..."
    Rscript scripts/idr_analysis.R
    if [ $? -eq 0 ]; then
        echo "IDR analysis completed successfully"
    else
        echo "Warning: IDR analysis encountered errors"
    fi
else
    echo "Warning: idr_analysis.R script not found"
fi

# 3. Enhanced GO enrichment analysis
echo ""
echo "Step 3: Enhanced GO enrichment analysis"
echo "--------------------------------------"
if [ -f "scripts/annotate_peaks.R" ]; then
    echo "Re-running peak annotation with enhanced GO analysis..."
    Rscript scripts/annotate_peaks.R
    if [ $? -eq 0 ]; then
        echo "Enhanced GO analysis completed successfully"
    else
        echo "Warning: Enhanced GO analysis encountered errors"
    fi
else
    echo "Warning: annotate_peaks.R script not found"
fi

echo ""
echo "=========================================="
echo "Additional analyses completed!"
echo "Date: $(date)"
echo "=========================================="

# Summary of outputs
echo ""
echo "Generated outputs:"
echo "------------------"
echo "Fragment analysis:"
ls -la results/fragment_analysis/ 2>/dev/null | tail -n +2 | awk '{print "  "$9}'

echo ""
echo "IDR analysis:"
ls -la results/idr_analysis/ 2>/dev/null | tail -n +2 | awk '{print "  "$9}'

echo ""
echo "Enhanced GO analysis:"
ls -la results/07_analysis_narrow/*GO* 2>/dev/null | awk '{print "  "$9}'