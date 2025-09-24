#!/bin/bash

#===============================================================================
# SCRIPT: master_script_narrow.sh
# PURPOSE: Master script for complete Cut&Tag analysis pipeline using NARROW PEAKS
#
# DESCRIPTION:
# This master script orchestrates the entire Cut&Tag analysis pipeline using
# narrow peak calling methodology. It submits all pipeline steps with proper
# job dependencies to ensure sequential execution and optimal resource utilization.
# The pipeline includes quality control, alignment, filtering, narrow peak calling,
# visualization, annotation, differential binding analysis, and reporting.
#
# PIPELINE WORKFLOW (Narrow Track Only):
# 1. Narrow peak calling (MACS2)
# 2. Peak comparison and visualization
# 3. Peak annotation and motif analysis
# 4. Differential binding analysis
# 5. Replicate combination and consensus
#
# KEY FEATURES:
# - Uses SLURM job dependencies for optimal scheduling
# - Narrow peak calling for precise binding site identification
# - Comprehensive quality control and validation
# - Publication-ready visualizations and reports
# - Automatic error handling and job monitoring
#
# OUTPUTS:
# - All results stored in results/ with _narrow suffix
# - Complete analysis pipeline from raw data to final report
# - Both PDF and PNG format plots for maximum compatibility
#
# USAGE:
# ./master_script_narrow.sh
#
# NOTES:
# - Run AFTER master_script_broad.sh has completed setup and shared components
# - This script only processes narrow peak-specific analysis steps
# - Assumes BigWig tracks and directory structure already exist
# - Use narrow peaks for precise transcription factor binding analysis
#===============================================================================

echo "=========================================="
echo "Cut&Tag Analysis Pipeline - NARROW PEAKS"
echo "Starting narrow peak-specific analysis"
echo "Date: $(date)"
echo "=========================================="

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva"
cd $BASE_DIR

# Ensure logs directory exists (may already exist from broad pipeline)
mkdir -p logs

echo "Step 1: Narrow peak calling"
JOB1=$(sbatch 5_peak_calling_narrow.sh | awk '{print $4}')
echo "Narrow peak calling job: $JOB1"

echo "Step 2: Peak comparison and visualization"
JOB2=$(sbatch --dependency=afterok:$JOB1 7_compare_narrow.sh | awk '{print $4}')
echo "Narrow peak analysis job: $JOB2"

echo "Step 3: Peak annotation and motif analysis"
JOB3=$(sbatch --dependency=afterok:$JOB1 8_annotate_narrow.sh | awk '{print $4}')
echo "Narrow peak annotation job: $JOB3"

echo "Step 4: Differential binding analysis"
JOB4=$(sbatch --dependency=afterok:$JOB1 9_diff_bind_narrow.sh | awk '{print $4}')
echo "Narrow peak DiffBind job: $JOB4"

echo "Step 5: Replicate combination and consensus"
JOB5=$(sbatch --dependency=afterok:$JOB1 11_combine_replicates_narrow.sh | awk '{print $4}')
echo "Narrow peak replicate combination job: $JOB5"

echo ""
echo "=========================================="
echo "NARROW PEAK PIPELINE SUBMITTED SUCCESSFULLY"
echo "=========================================="
echo "Pipeline type: NARROW PEAKS (Complementary Track)"
echo "Total jobs submitted: 5"
echo "Final job ID: $JOB5"
echo ""
echo "Job dependencies configured for optimal execution:"
echo "Peak Calling → [Analysis + Annotation + DiffBind + Consensus]"
echo ""
echo "Output directories:"
echo "- Narrow peaks: results/05_peaks_narrow/"
echo "- Analysis: results/07_analysis_narrow/"
echo "- Consensus: results/11_combined_replicates_narrow/"
echo ""
echo "Note: This is the complementary narrow peak track."
echo "Shared components (BigWig, SES overlap, reports) handled by broad pipeline."
echo ""
echo "Monitor progress with: squeue -u $USER"
echo "Check completion with: sacct -j $JOB5"
echo "=========================================="