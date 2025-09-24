#!/bin/bash

#===============================================================================
# SCRIPT: master_script_broad.sh
# PURPOSE: Master script for complete Cut&Tag analysis pipeline using BROAD PEAKS
#
# DESCRIPTION:
# This master script orchestrates the entire Cut&Tag analysis pipeline using
# broad peak calling methodology. It submits all pipeline steps with proper
# job dependencies to ensure sequential execution and optimal resource utilization.
# The pipeline includes quality control, alignment, filtering, broad peak calling,
# visualization, annotation, differential binding analysis, and reporting.
#
# PIPELINE WORKFLOW:
# 1. Setup and initialization
# 2. Quality control (FastQC)
# 3. Adapter trimming
# 4. Read alignment
# 5. BAM filtering and deduplication
# 6. Broad peak calling (MACS2)
# 7. BigWig track generation
# 8. Peak comparison and visualization
# 9. Peak annotation and motif analysis
# 10. Differential binding analysis
# 11. Replicate combination and consensus
# 12. Pairwise overlap analysis
# 13. SES overlap analysis
# 14. Summary report generation
#
# KEY FEATURES:
# - Uses SLURM job dependencies for optimal scheduling
# - Broad peak calling for enhanced overlap detection and sensitivity
# - Comprehensive quality control and validation
# - Publication-ready visualizations and reports
# - Automatic error handling and job monitoring
#
# OUTPUTS:
# - All results stored in results/ with _broad suffix
# - Complete analysis pipeline from raw data to final report
# - Both PDF and PNG format plots for maximum compatibility
#
# USAGE:
# ./master_script_broad.sh
#
# NOTES:
# - Ensure all input data is properly configured
# - Monitor job progress through SLURM queue
# - Final job completion indicates successful pipeline execution
# - Use broad peaks for improved overlap analysis and chromatin domain detection
#===============================================================================

echo "=========================================="
echo "Cut&Tag Analysis Pipeline - BROAD PEAKS"
echo "Starting complete analysis workflow"
echo "Date: $(date)"
echo "=========================================="

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva"
cd $BASE_DIR

# Ensure logs directory exists
mkdir -p logs

echo "Step 1: Setting up environment and directories"
JOB0=$(sbatch setup.sh | awk '{print $4}')
echo "Setup job: $JOB0"

# echo "Step 2: Quality control analysis"
# JOB1=$(sbatch --dependency=afterok:$JOB0 1_qc.sh | awk '{print $4}')
# echo "FastQC job: $JOB1"

# echo "Step 3: Adapter trimming"
# JOB2=$(sbatch --dependency=afterok:$JOB1 2_adapter_trim.sh | awk '{print $4}')
# echo "Trimming job: $JOB2"

# echo "Step 4: Read alignment"
# JOB3=$(sbatch --dependency=afterok:$JOB2 3_align.sh | awk '{print $4}')
# echo "Alignment job: $JOB3"

# echo "Step 5: BAM filtering and deduplication"
# JOB4=$(sbatch --dependency=afterok:$JOB3 4_filter.sh | awk '{print $4}')
# echo "Filtering job: $JOB4"

echo "Step 6: Broad peak calling"
JOB5=$(sbatch --dependency=afterok:$JOB0 5_peak_calling_broad.sh | awk '{print $4}')
echo "Broad peak calling job: $JOB5"

echo "Step 7: BigWig track generation"
JOB6=$(sbatch --dependency=afterok:$JOB0 6_bigwig.sh | awk '{print $4}')
echo "BigWig job: $JOB6"

echo "Step 8: Peak comparison and visualization"
JOB7=$(sbatch --dependency=afterok:$JOB5:$JOB6 7_compare_broad.sh | awk '{print $4}')
echo "Broad peak analysis job: $JOB7"

echo "Step 9: Peak annotation and motif analysis"
JOB8=$(sbatch --dependency=afterok:$JOB5 8_annotate_broad.sh | awk '{print $4}')
echo "Broad peak annotation job: $JOB8"

echo "Step 10: Differential binding analysis"
JOB9=$(sbatch --dependency=afterok:$JOB5:$JOB6 9_diff_bind_broad.sh | awk '{print $4}')
echo "Broad peak DiffBind job: $JOB9"

echo "Step 11: Replicate combination and consensus"
JOB10=$(sbatch --dependency=afterok:$JOB5:$JOB6 11_combine_replicates_broad.sh | awk '{print $4}')
echo "Broad peak replicate combination job: $JOB10"

echo "Step 12: Pairwise overlap analysis"
JOB11=$(sbatch --dependency=afterok:$JOB10 12_pairwise_overlap.sh | awk '{print $4}')
echo "Pairwise overlap job: $JOB11"

echo "Step 13: Combine BigWig tracks"
JOB12=$(sbatch --dependency=afterok:$JOB6 13_combine_bigwig.sh | awk '{print $4}')
echo "Combine BigWig job: $JOB12"

echo "Step 14: SES overlap analysis"
JOB13=$(sbatch --dependency=afterok:$JOB5 13_ses_overlap.sh | awk '{print $4}')
echo "SES overlap job: $JOB13"

echo "Step 15: SES overlap plots"
JOB14=$(sbatch --dependency=afterok:$JOB13 14_ses_overlap_plots.sh | awk '{print $4}')
echo "SES overlap plots job: $JOB14"

echo "Step 16: Generate comprehensive summary"
JOB15=$(sbatch --dependency=afterok:$JOB7:$JOB8:$JOB9:$JOB11:$JOB12:$JOB13:$JOB14 10_summary.sh | awk '{print $4}')
echo "Summary report job: $JOB15"

echo "Step 17: Generate MultiQC report"
JOB16=$(sbatch --dependency=afterok:$JOB15 multiqc.sh | awk '{print $4}')
echo "MultiQC job: $JOB16"

echo ""
echo "=========================================="
echo "BROAD PEAK PIPELINE SUBMITTED SUCCESSFULLY"
echo "=========================================="
echo "Pipeline type: BROAD PEAKS"
echo "Total jobs submitted: 18"
echo "Final job ID: $JOB16"
echo ""
echo "Job dependencies configured for optimal execution:"
echo "QC → Trim → Align → Filter → [Peak Calling + BigWig] → Analysis"
echo ""
echo "Output directories:"
echo "- Broad peaks: results/05_peaks_broad/"
echo "- Analysis: results/07_analysis_broad/"
echo "- Consensus: results/11_combined_replicates_broad/"
echo ""
echo "Key advantages of broad peak analysis:"
echo "- Enhanced sensitivity for overlap detection"
echo "- Better capture of extended chromatin domains"
echo "- Improved consensus peak calling"
echo "- Optimal for Cut&Tag chromatin accessibility analysis"
echo ""
echo "Monitor progress with: squeue -u $USER"
echo "Check completion with: sacct -j $JOB16"
echo "=========================================="