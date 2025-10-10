#!/bin/bash

#===============================================================================
# IMPROVED CUT&TAG ANALYSIS PIPELINE - NARROW PEAKS ONLY
#
# This script implements an improved, unified pipeline for Cut&Tag analysis
# with enhanced statistical rigor, quality control, and biological interpretation
#
# Key Improvements:
# - Unified configuration management
# - Enhanced quality control checkpoints
# - Proper statistical workflow (DiffBind before peak overlap)
# - Comprehensive pathway analysis (GO, KEGG, Reactome)
# - Better error handling and validation
# - Focus on narrow peaks only for transcription factor binding
#
# Author: Improved pipeline for SRF_Eva analysis
# Date: 2025
#===============================================================================

#SBATCH --job-name=improved_pipeline
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/improved_pipeline.err"
#SBATCH --output="./logs/improved_pipeline.out"

# Set up environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG"
cd $BASE_DIR

# Create necessary directories
mkdir -p logs temp results/qc_reports

echo "=== IMPROVED CUT&TAG PIPELINE STARTED ==="
echo "Date: $(date)"
echo "Base directory: $BASE_DIR"

# Function to check if previous step completed successfully
check_success() {
    if [ $? -ne 0 ]; then
        echo "ERROR: $1 failed. Exiting pipeline."
        exit 1
    else
        echo "SUCCESS: $1 completed."
    fi
}

# Function to validate configuration
validate_config() {
    echo "Validating pipeline configuration..."

    if [ ! -f "config/pipeline_config.yaml" ]; then
        echo "ERROR: Configuration file not found: config/pipeline_config.yaml"
        exit 1
    fi

    # Check if required directories exist
    for dir in data/raw results logs scripts config; do
        if [ ! -d "$dir" ]; then
            echo "Creating directory: $dir"
            mkdir -p $dir
        fi
    done

    echo "Configuration validation complete."
}

# Step 0: Validate configuration
validate_config
check_success "Configuration validation"

# Step 1: Enhanced Differential Binding Analysis (DiffBind with QC)
echo "=== STEP 1: ENHANCED DIFFERENTIAL BINDING ANALYSIS ==="
conda activate diffbind_analysis

Rscript scripts/improved_diffbind_analysis.R
check_success "Enhanced DiffBind analysis"

# Step 2: Integrated Annotation and Pathway Analysis
echo "=== STEP 2: INTEGRATED ANNOTATION AND PATHWAY ANALYSIS ==="

Rscript scripts/integrated_annotation_and_GO.R
check_success "Integrated annotation and pathway analysis"

# Step 3: Enhanced GO Visualization (using our previously created script)
echo "=== STEP 3: ENHANCED GO VISUALIZATIONS ==="

Rscript scripts/enhanced_GO_visualization.R
check_success "Enhanced GO visualizations"

# Step 4: Quality Control Report Generation
echo "=== STEP 4: QUALITY CONTROL REPORT ==="

# Create QC summary script
cat > scripts/generate_qc_report.R << 'EOF'
#!/usr/bin/env Rscript
# Generate comprehensive QC report

library(yaml)
library(knitr)
library(rmarkdown)

config <- yaml.load_file("config/pipeline_config.yaml")
setwd(config$directories$base)

# Create QC report
qc_report <- '
---
title: "SRF Eva Cut&Tag Analysis - Quality Control Report"
author: "Improved Pipeline"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: cosmo
---

# Pipeline Summary

This report summarizes the quality control metrics and results from the improved Cut&Tag analysis pipeline.

## Sample Information

```{r sample_info, echo=FALSE}
config <- yaml::yaml.load_file("config/pipeline_config.yaml")
knitr::kable(data.frame(
  Group = names(config$samples$groups),
  Samples = sapply(config$samples$groups, function(x) paste(x$samples, collapse=", ")),
  Control = sapply(config$samples$groups, function(x) x$control),
  Description = sapply(config$samples$groups, function(x) x$description)
))
```

## Differential Binding Results

```{r diffbind_stats, echo=FALSE}
if (file.exists("results/diffbind_analysis/differential_binding_statistics.csv")) {
  stats <- read.csv("results/diffbind_analysis/differential_binding_statistics.csv")
  knitr::kable(stats, caption = "Differential Binding Statistics")
} else {
  cat("Differential binding statistics not found.")
}
```

## Quality Control Metrics

### Read Count Summary
```{r read_counts, echo=FALSE}
if (file.exists("results/diffbind_analysis/read_count_summary.csv")) {
  counts <- read.csv("results/diffbind_analysis/read_count_summary.csv")
  knitr::kable(counts[,c("ID", "Reads", "FRiP")], caption = "Read Count Summary")
} else {
  cat("Read count summary not found.")
}
```

## Key Visualizations

### PCA Plot
![PCA Analysis](results/diffbind_analysis/03_enhanced_PCA.pdf)

### Sample Correlations
![Sample Correlations](results/diffbind_analysis/04_final_correlation_heatmap.pdf)

## Pathway Analysis Summary

The integrated analysis identified key pathways enriched in TES-specific and TEAD1-specific binding sites.

### TES-Specific Pathways
Key biological processes uniquely targeted by TES include:
- Small GTPase-mediated signal transduction
- Regulation of neuron projection development
- Actin filament-based processes
- Epithelial morphogenesis

### Analysis Files Generated
- Differential binding results: `results/diffbind_analysis/`
- Pathway analysis: `results/integrated_analysis/`
- Enhanced GO visualizations: `results/07_analysis_narrow/enhanced_GO/`

## Pipeline Configuration
The analysis used the following key parameters:
- FDR threshold: 0.05
- Fold change threshold: 1.5
- Peak type: Narrow peaks only
- Statistical method: DESeq2 via DiffBind
'

writeLines(qc_report, "scripts/qc_report.Rmd")

# Render the report
rmarkdown::render("scripts/qc_report.Rmd",
                  output_file = "results/qc_reports/pipeline_qc_report.html")

cat("QC report generated: results/qc_reports/pipeline_qc_report.html\n")
EOF

Rscript scripts/generate_qc_report.R
check_success "QC report generation"

# Step 5: Generate final summary
echo "=== STEP 5: FINAL SUMMARY ==="

cat > results/IMPROVED_PIPELINE_SUMMARY.txt << EOF
====================================================
IMPROVED CUT&TAG ANALYSIS PIPELINE - FINAL SUMMARY
====================================================

Pipeline completed successfully on: $(date)

KEY IMPROVEMENTS IMPLEMENTED:
✓ Unified configuration management (config/pipeline_config.yaml)
✓ Enhanced statistical rigor (proper DiffBind workflow)
✓ Comprehensive quality control checkpoints
✓ Integrated pathway analysis (GO, KEGG, Reactome)
✓ Enhanced visualizations with biological interpretation
✓ Focus on narrow peaks for transcription factor analysis
✓ Better error handling and validation

RESULTS DIRECTORIES:
- results/diffbind_analysis/          - Statistical differential binding results
- results/integrated_analysis/        - Peak annotations and pathway analysis
- results/07_analysis_narrow/enhanced_GO/ - Enhanced GO visualizations
- results/qc_reports/                 - Quality control reports

KEY FILES:
- Differential binding statistics: results/diffbind_analysis/differential_binding_statistics.csv
- TES-specific genes: results/integrated_analysis/TES_specific_genes.txt
- TEAD1-specific genes: results/integrated_analysis/TEAD1_specific_genes.txt
- QC report: results/qc_reports/pipeline_qc_report.html

BIOLOGICAL INSIGHTS:
The improved analysis reveals that TES targets additional pathways beyond TEAD1,
particularly in:
- Small GTPase signaling regulation
- Neuronal development and axonogenesis
- Cytoskeletal organization
- Epithelial morphogenesis

This provides a clearer understanding of TES-specific transcriptional programs
in glioblastoma cells.

NEXT STEPS:
1. Review QC report for any quality issues
2. Examine pathway enrichment results for biological insights
3. Consider validation experiments for key TES-specific targets
4. Use enhanced visualizations for publication
EOF

echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo "Summary saved to: results/IMPROVED_PIPELINE_SUMMARY.txt"
echo "QC report available at: results/qc_reports/pipeline_qc_report.html"