#!/bin/bash
#SBATCH --job-name=promoter_heatmaps
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/promoter_heatmaps.out
#SBATCH --error=logs/promoter_heatmaps.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

echo "=========================================="
echo "PROMOTER HEATMAP GENERATION"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo ""
echo "Task: Visualize TES/TEAD1 binding at UP vs DOWN gene promoters"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG

# Create output directory
mkdir -p results/08_promoter_heatmaps

# Activate deepTools environment
echo "Loading deepTools environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate tg  # deepTools is installed in tg environment

# =============================================================================
# STEP 1: PREPARE GENE LISTS FROM INTEGRATIVE ANALYSIS
# =============================================================================

echo ""
echo "=== STEP 1: Preparing Gene Lists ==="
echo ""

# First, run R script to extract gene coordinates
cat > /tmp/prepare_gene_lists_$SLURM_JOB_ID.R << 'EOF'
library(dplyr)
library(GenomicRanges)
library(rtracklayer)

# Load GTF annotation
gtf <- rtracklayer::import("/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf")

# Convert to data frame and filter for genes only
genes_df <- as.data.frame(gtf) %>%
  filter(type == "gene") %>%
  select(seqnames, start, end, strand, gene_id, gene_name)

# Clean gene IDs (remove version)
genes_df$gene_id_clean <- gsub("\\..*", "", genes_df$gene_id)

# Load integrative analysis results
tes_direct <- read.csv("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/output/results/direct_targets/TES_direct_targets_all_genes.csv")
tead1_direct <- read.csv("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/output/results/direct_targets/TEAD1_direct_targets_all_genes.csv")

# Split by direction
tes_up <- tes_direct[tes_direct$log2FoldChange > 0, ]
tes_down <- tes_direct[tes_direct$log2FoldChange < 0, ]
tead1_up <- tead1_direct[tead1_direct$log2FoldChange > 0, ]
tead1_down <- tead1_direct[tead1_direct$log2FoldChange < 0, ]

cat(sprintf("TES UP: %d genes, TES DOWN: %d genes\n", nrow(tes_up), nrow(tes_down)))
cat(sprintf("TEAD1 UP: %d genes, TEAD1 DOWN: %d genes\n", nrow(tead1_up), nrow(tead1_down)))

# Function to create BED file
create_bed <- function(gene_data, genes_annotation, output_file) {
  # Merge with annotation
  merged <- merge(gene_data, genes_annotation,
                  by.x = "ensembl_id", by.y = "gene_id_clean",
                  all.x = TRUE)

  # Remove genes without coordinates
  merged <- merged[!is.na(merged$start), ]

  # Define promoter regions (TSS +/- 2kb)
  merged$promoter_start <- ifelse(merged$strand == "+",
                                  merged$start - 2000,
                                  merged$end - 2000)
  merged$promoter_end <- ifelse(merged$strand == "+",
                                merged$start + 2000,
                                merged$end + 2000)

  # Ensure valid coordinates
  merged$promoter_start[merged$promoter_start < 0] <- 0

  # Create BED format (chr, start, end, name, score, strand)
  bed <- data.frame(
    chr = merged$seqnames,
    start = merged$promoter_start,
    end = merged$promoter_end,
    name = merged$gene_name,
    score = abs(merged$log2FoldChange) * 100,  # Scale for visibility
    strand = merged$strand
  )

  # Sort by chromosome and start
  bed <- bed[order(bed$chr, bed$start), ]

  # Write BED file
  write.table(bed, output_file,
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

  cat(sprintf("Created %s: %d regions\n", output_file, nrow(bed)))
  return(nrow(bed))
}

# Create BED files for each category
outdir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/08_promoter_heatmaps"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

n_tes_up <- create_bed(tes_up, genes_df, file.path(outdir, "TES_UP_promoters.bed"))
n_tes_down <- create_bed(tes_down, genes_df, file.path(outdir, "TES_DOWN_promoters.bed"))
n_tead1_up <- create_bed(tead1_up, genes_df, file.path(outdir, "TEAD1_UP_promoters.bed"))
n_tead1_down <- create_bed(tead1_down, genes_df, file.path(outdir, "TEAD1_DOWN_promoters.bed"))

cat("\nBED file creation complete!\n")
EOF

# Temporarily switch to R environment for gene list preparation
echo "Switching to R environment for gene list preparation..."
conda activate r_chipseq_env

# Run R script
Rscript /tmp/prepare_gene_lists_$SLURM_JOB_ID.R

if [ $? -ne 0 ]; then
    echo "ERROR: Gene list preparation failed!"
    exit 1
fi

# Switch back to deepTools environment (tg)
echo "Switching back to deepTools environment..."
conda activate tg

# =============================================================================
# STEP 2: COMPUTE COVERAGE MATRICES WITH DEEPTOOLS
# =============================================================================

echo ""
echo "=== STEP 2: Computing Coverage Matrices ==="
echo ""

OUTDIR="results/08_promoter_heatmaps"
BIGWIG_DIR="results/06_bigwig"

# Check if BigWig files exist
if [ ! -d "$BIGWIG_DIR" ]; then
    echo "ERROR: BigWig directory not found: $BIGWIG_DIR"
    exit 1
fi

# Find BigWig files (pattern: TES-1_CPM.bw, TES-2_CPM.bw, etc.)
TES_BIGWIGS=$(ls ${BIGWIG_DIR}/TES-[0-9]*_CPM.bw 2>/dev/null | tr '\n' ' ')
TEAD1_BIGWIGS=$(ls ${BIGWIG_DIR}/TEAD1-[0-9]*_CPM.bw 2>/dev/null | tr '\n' ' ')

if [ -z "$TES_BIGWIGS" ] || [ -z "$TEAD1_BIGWIGS" ]; then
    echo "ERROR: BigWig files not found!"
    echo "TES BigWigs: $TES_BIGWIGS"
    echo "TEAD1 BigWigs: $TEAD1_BIGWIGS"
    exit 1
fi

echo "TES BigWigs: $TES_BIGWIGS"
echo "TEAD1 BigWigs: $TEAD1_BIGWIGS"
echo ""

# Combine BED files for single analysis
cat ${OUTDIR}/TES_UP_promoters.bed > ${OUTDIR}/TES_promoters_combined.bed
cat ${OUTDIR}/TES_DOWN_promoters.bed >> ${OUTDIR}/TES_promoters_combined.bed

cat ${OUTDIR}/TEAD1_UP_promoters.bed > ${OUTDIR}/TEAD1_promoters_combined.bed
cat ${OUTDIR}/TEAD1_DOWN_promoters.bed >> ${OUTDIR}/TEAD1_promoters_combined.bed

# Compute matrix for TES targets
echo "Computing matrix for TES targets..."
computeMatrix reference-point \
    --referencePoint TSS \
    -b 10000 -a 10000 \
    --binSize 10 \
    -R ${OUTDIR}/TES_UP_promoters.bed ${OUTDIR}/TES_DOWN_promoters.bed \
    -S $TES_BIGWIGS $TEAD1_BIGWIGS \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/TES_targets_matrix.gz \
    --outFileSortedRegions ${OUTDIR}/TES_targets_sorted.bed \
    -p 16 \
    2>&1 | tee ${OUTDIR}/computeMatrix_TES.log

if [ $? -ne 0 ]; then
    echo "WARNING: TES matrix computation had issues, continuing..."
fi

# Compute matrix for TEAD1 targets
echo "Computing matrix for TEAD1 targets..."
computeMatrix reference-point \
    --referencePoint TSS \
    -b 10000 -a 10000 \
    --binSize 10 \
    -R ${OUTDIR}/TEAD1_UP_promoters.bed ${OUTDIR}/TEAD1_DOWN_promoters.bed \
    -S $TES_BIGWIGS $TEAD1_BIGWIGS \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/TEAD1_targets_matrix.gz \
    --outFileSortedRegions ${OUTDIR}/TEAD1_targets_sorted.bed \
    -p 16 \
    2>&1 | tee ${OUTDIR}/computeMatrix_TEAD1.log

if [ $? -ne 0 ]; then
    echo "WARNING: TEAD1 matrix computation had issues, continuing..."
fi

# =============================================================================
# STEP 3: GENERATE HEATMAPS
# =============================================================================

echo ""
echo "=== STEP 3: Generating Heatmaps ==="
echo ""

# Heatmap for TES targets
if [ -f "${OUTDIR}/TES_targets_matrix.gz" ]; then
    echo "Creating TES targets heatmap..."
    plotHeatmap -m ${OUTDIR}/TES_targets_matrix.gz \
        -out ${OUTDIR}/TES_targets_heatmap.pdf \
        --colorMap RdYlBu_r \
        --whatToShow 'heatmap and colorbar' \
        --zMin 0 --zMax 10 \
        --heatmapHeight 15 \
        --refPointLabel "TSS" \
        --regionsLabel "TES UP" "TES DOWN" \
        --samplesLabel "TES-1" "TES-2" "TES-3" "TEAD1-1" "TEAD1-2" "TEAD1-3" \
        --plotTitle "TES/TEAD1 Binding at TES Target Promoters (UP vs DOWN)" \
        2>&1 | tee ${OUTDIR}/plotHeatmap_TES.log

    # Profile plot for TES targets
    echo "Creating TES targets profile plot..."
    plotProfile -m ${OUTDIR}/TES_targets_matrix.gz \
        -out ${OUTDIR}/TES_targets_profile.pdf \
        --perGroup \
        --colors "#E31A1C" "#FB9A99" "#FDBF6F" "#1F78B4" "#A6CEE3" "#B2DF8A" \
        --refPointLabel "TSS" \
        --regionsLabel "TES UP" "TES DOWN" \
        --plotTitle "Average TES/TEAD1 Binding Profile at TES Targets" \
        2>&1 | tee ${OUTDIR}/plotProfile_TES.log
fi

# Heatmap for TEAD1 targets
if [ -f "${OUTDIR}/TEAD1_targets_matrix.gz" ]; then
    echo "Creating TEAD1 targets heatmap..."
    plotHeatmap -m ${OUTDIR}/TEAD1_targets_matrix.gz \
        -out ${OUTDIR}/TEAD1_targets_heatmap.pdf \
        --colorMap RdYlBu_r \
        --whatToShow 'heatmap and colorbar' \
        --zMin 0 --zMax 10 \
        --heatmapHeight 15 \
        --refPointLabel "TSS" \
        --regionsLabel "TEAD1 UP" "TEAD1 DOWN" \
        --samplesLabel "TES-1" "TES-2" "TES-3" "TEAD1-1" "TEAD1-2" "TEAD1-3" \
        --plotTitle "TES/TEAD1 Binding at TEAD1 Target Promoters (UP vs DOWN)" \
        2>&1 | tee ${OUTDIR}/plotHeatmap_TEAD1.log

    # Profile plot for TEAD1 targets
    echo "Creating TEAD1 targets profile plot..."
    plotProfile -m ${OUTDIR}/TEAD1_targets_matrix.gz \
        -out ${OUTDIR}/TEAD1_targets_profile.pdf \
        --perGroup \
        --colors "#E31A1C" "#FB9A99" "#FDBF6F" "#1F78B4" "#A6CEE3" "#B2DF8A" \
        --refPointLabel "TSS" \
        --regionsLabel "TEAD1 UP" "TEAD1 DOWN" \
        --plotTitle "Average TES/TEAD1 Binding Profile at TEAD1 Targets" \
        2>&1 | tee ${OUTDIR}/plotProfile_TEAD1.log
fi

# =============================================================================
# COMPLETION
# =============================================================================

echo ""
echo "=========================================="
echo "PROMOTER HEATMAP GENERATION COMPLETE"
echo "=========================================="
echo "Completed: $(date)"
echo ""
echo "Output directory: results/08_promoter_heatmaps/"
echo ""
echo "Generated files:"
ls -lh ${OUTDIR}/*.pdf 2>/dev/null || echo "  (Check for PDF files manually)"
echo ""
echo "Key outputs:"
echo "  - TES_targets_heatmap.pdf"
echo "  - TES_targets_profile.pdf"
echo "  - TEAD1_targets_heatmap.pdf"
echo "  - TEAD1_targets_profile.pdf"
