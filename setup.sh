#!/bin/bash
#SBATCH --job-name=setup_env
#SBATCH --account=kubacki.michal
#SBATCH --mem=4GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/setup.err"
#SBATCH --output="./logs/setup.out"

# Create directory structure for dual-track pipeline
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG"
cd $BASE_DIR

mkdir -p scripts logs config
mkdir -p results/{01_fastqc,02_trimmed,03_aligned,04_filtered}
mkdir -p results/{05_peaks_narrow,05_peaks_broad}
mkdir -p results/06_bigwig
mkdir -p results/{07_analysis_narrow,07_analysis_broad}
mkdir -p results/{11_combined_replicates_narrow,11_combined_replicates_broad}
mkdir -p results/{12_pairwise_overlap,13_ses_overlap}
mkdir -p results/{multiqc,report}

# Create sample configuration file
cat > config/samples.txt << EOF
TES-1
TES-2
TES-3
TESmut-1
TESmut-2
TESmut-3
TEAD1-1
TEAD1-2
TEAD1-3
IggMs
IggRb
EOF

# Create group configuration
cat > config/groups.txt << EOF
TES:TES-1,TES-2,TES-3:IggMs
TESmut:TESmut-1,TESmut-2,TESmut-3:IggMs
TEAD1:TEAD1-1,TEAD1-2,TEAD1-3:IggRb
EOF

echo "Dual-track Cut&Tag pipeline setup complete!"
echo ""
echo "Created directory structure:"
echo "- Common stages: 01_fastqc, 02_trimmed, 03_aligned, 04_filtered, 06_bigwig"
echo "- Narrow track: 05_peaks_narrow, 07_analysis_narrow, 11_combined_replicates_narrow"
echo "- Broad track: 05_peaks_broad, 07_analysis_broad, 11_combined_replicates_broad"
echo "- Comparative: 12_pairwise_overlap, 13_ses_overlap"
echo "- Reports: multiqc, report"
echo ""
echo "Sample configuration files created in config/"
echo "Ready for dual-track analysis pipeline execution!"