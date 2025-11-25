#!/bin/bash

#SBATCH --job-name=15_cell_death_proliferation_analysis
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/15_cell_death_proliferation_analysis.err"
#SBATCH --output="./logs/15_cell_death_proliferation_analysis.out"

# SCRIPT: run_cell_death_analysis.sh
# PURPOSE: Execute cell death and proliferation pathway analysis
#
# DESCRIPTION:
# This script runs specialized pathway analysis focused on cell death and
# proliferation mechanisms. This is particularly relevant for TES/TEAD1
# research in glioblastoma, where these transcription factors play crucial
# roles in cell survival, apoptosis, and proliferation control.
#
# KEY ANALYSES:
# 1. Focused GO enrichment for cell death and proliferation
# 2. KEGG pathway analysis for cancer-related pathways
# 3. Reactome pathway analysis
# 4. Comparative analysis between TES, TESmut, and TEAD1
# 5. Gene overlap and functional similarity analysis


# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate annotation_enrichment

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG"
cd ${BASE_DIR}

echo "=========================================="
echo "Cell Death & Proliferation Analysis"
echo "Focus: TES/TEAD1 in Glioblastoma"
echo "Date: $(date)"
echo "=========================================="
echo ""
echo "NOTE: This analysis uses annotated peaks from step 8."
echo "Step 8 now automatically uses high-quality consensus peaks"
echo "when available, ensuring the most reproducible results."
echo ""

# Check if required input files exist
echo "Checking input files..."
echo "----------------------"

ANNOTATION_DIR="results/07_analysis_narrow"
required_files=(
    "${ANNOTATION_DIR}/TES_peaks_annotated.csv"
    "${ANNOTATION_DIR}/TESmut_peaks_annotated.csv"
    "${ANNOTATION_DIR}/TEAD1_peaks_annotated.csv"
)

missing_files=0
for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        echo "✓ Found: $file"
    else
        echo "✗ Missing: $file"
        missing_files=$((missing_files + 1))
    fi
done

if [ $missing_files -gt 0 ]; then
    echo ""
    echo "Warning: $missing_files required files are missing"
    echo "Please run the peak annotation script first:"
    echo "sbatch 8_annotate_narrow.sh"
    echo ""
    echo "Attempting to run analysis with available files..."
fi

# Run the cell death and proliferation analysis
echo ""
echo "Running cell death and proliferation analysis..."
echo "------------------------------------------------"

if [ -f "scripts/cell_death_proliferation_analysis.R" ]; then
    echo "Executing specialized pathway analysis..."

    # Set R library path and run analysis
    export R_LIBS_USER="/beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/annotation_enrichment/lib/R/library"

    Rscript scripts/cell_death_proliferation_analysis.R

    if [ $? -eq 0 ]; then
        echo "✓ Cell death and proliferation analysis completed successfully"
    else
        echo "✗ Analysis encountered errors - check log files"
        exit 1
    fi
else
    echo "✗ Error: cell_death_proliferation_analysis.R script not found"
    exit 1
fi

echo ""
echo "=========================================="
echo "Analysis Summary"
echo "=========================================="

# Check output directory
OUTPUT_DIR="results/cell_death_proliferation"
if [ -d "$OUTPUT_DIR" ]; then
    echo ""
    echo "Generated files:"
    echo "---------------"

    # Count different types of outputs
    pdf_count=$(find "$OUTPUT_DIR" -name "*.pdf" | wc -l)
    png_count=$(find "$OUTPUT_DIR" -name "*.png" | wc -l)
    csv_count=$(find "$OUTPUT_DIR" -name "*.csv" | wc -l)

    echo "PDF plots: $pdf_count"
    echo "PNG plots: $png_count"
    echo "CSV results: $csv_count"

    echo ""
    echo "Key output files:"
    echo "----------------"

    # List important files
    if [ -f "$OUTPUT_DIR/analysis_summary.csv" ]; then
        echo "✓ Analysis summary: $OUTPUT_DIR/analysis_summary.csv"
    fi

    if [ -f "$OUTPUT_DIR/gene_overlap_stats.csv" ]; then
        echo "✓ Gene overlap statistics: $OUTPUT_DIR/gene_overlap_stats.csv"
    fi

    # List GO results
    go_files=$(find "$OUTPUT_DIR" -name "*_GO_results.csv" | wc -l)
    if [ $go_files -gt 0 ]; then
        echo "✓ GO enrichment results: $go_files files"
    fi

    # List KEGG results
    kegg_files=$(find "$OUTPUT_DIR" -name "*_KEGG_*results.csv" | wc -l)
    if [ $kegg_files -gt 0 ]; then
        echo "✓ KEGG pathway results: $kegg_files files"
    fi

    # List Reactome results
    reactome_files=$(find "$OUTPUT_DIR" -name "*_Reactome_*results.csv" | wc -l)
    if [ $reactome_files -gt 0 ]; then
        echo "✓ Reactome pathway results: $reactome_files files"
    fi

    echo ""
    echo "Biological insights to examine:"
    echo "------------------------------"
    echo "• Cell death pathways regulated by TES vs TEAD1"
    echo "• Proliferation mechanisms specific to each factor"
    echo "• Apoptosis regulation differences between TES and TESmut"
    echo "• Cancer-related pathway enrichments"
    echo "• Gene overlap patterns between conditions"

else
    echo "✗ Output directory not created - analysis may have failed"
fi

echo ""
echo "=========================================="
echo "Analysis completed: $(date)"
echo "Results location: $OUTPUT_DIR"
echo "=========================================="

# Quick summary if analysis summary exists
if [ -f "$OUTPUT_DIR/analysis_summary.csv" ]; then
    echo ""
    echo "Quick Results Summary:"
    echo "---------------------"
    # Display the summary table (skip header)
    tail -n +2 "$OUTPUT_DIR/analysis_summary.csv" | while IFS=, read -r sample analysis terms; do
        printf "%-10s %-20s %s significant terms\n" "$sample" "$analysis" "$terms"
    done
fi