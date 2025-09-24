#!/bin/bash
# Save as: run_full_pipeline.sh

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva"
cd $BASE_DIR

# Submit jobs with dependencies
JOB0=$(sbatch scripts/0_setup.sh | awk '{print $4}')
echo "Setup job: $JOB0"

JOB1=$(sbatch --dependency=afterok:$JOB0 scripts/1_fastqc.sh | awk '{print $4}')
echo "FastQC job: $JOB1"

JOB2=$(sbatch --dependency=afterok:$JOB1 scripts/2_trim.sh | awk '{print $4}')
echo "Trimming job: $JOB2"

JOB3=$(sbatch --dependency=afterok:$JOB2 scripts/3_align.sh | awk '{print $4}')
echo "Alignment job: $JOB3"

JOB4=$(sbatch --dependency=afterok:$JOB3 scripts/4_filter.sh | awk '{print $4}')
echo "Filtering job: $JOB4"

JOB5=$(sbatch --dependency=afterok:$JOB4 scripts/5_peaks.sh | awk '{print $4}')
echo "Peak calling job: $JOB5"

JOB6=$(sbatch --dependency=afterok:$JOB4 scripts/6_bigwig.sh | awk '{print $4}')
echo "BigWig job: $JOB6"

JOB7=$(sbatch --dependency=afterok:$JOB5:$JOB6 scripts/7_analysis.sh | awk '{print $4}')
echo "Analysis job: $JOB7"

JOB8=$(sbatch --dependency=afterok:$JOB7 scripts/8_annotate.sh | awk '{print $4}')
echo "Annotation job: $JOB8"

JOB9=$(sbatch --dependency=afterok:$JOB5 scripts/9_diffbind.sh | awk '{print $4}')
echo "DiffBind job: $JOB9"

JOB10=$(sbatch --dependency=afterok:$JOB8:$JOB9 scripts/10_report.sh | awk '{print $4}')
echo "Report job: $JOB10"

echo "Pipeline submitted! Final job ID: $JOB10"