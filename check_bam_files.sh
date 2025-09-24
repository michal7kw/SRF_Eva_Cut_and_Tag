#!/bin/bash

# Quick diagnostic script for BAM files
BAM_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva/results/04_filtered"

echo "BAM File Diagnostic Report"
echo "=========================="
echo "Directory: $BAM_DIR"
echo ""

# Check if directory exists
if [ ! -d "$BAM_DIR" ]; then
    echo "ERROR: Directory does not exist!"
    exit 1
fi

# Summary table header
printf "%-25s %-10s %-12s %-12s %-10s\n" "File" "Size" "Total Reads" "Paired %" "Status"
printf "%-25s %-10s %-12s %-12s %-10s\n" "----" "----" "-----------" "--------" "------"

# Check each BAM file
for bam in ${BAM_DIR}/*.bam; do
    if [ -f "$bam" ]; then
        filename=$(basename $bam)
        filesize=$(ls -lh $bam | awk '{print $5}')
        
        # Count reads
        total_reads=$(samtools view -c $bam 2>/dev/null)
        paired_reads=$(samtools view -f 1 -c $bam 2>/dev/null)
        
        # Calculate percentage
        if [ "$total_reads" -gt 0 ]; then
            paired_pct=$((paired_reads * 100 / total_reads))
            if [ "$paired_pct" -gt 50 ]; then
                status="PAIRED"
            else
                status="SINGLE"
            fi
        else
            paired_pct=0
            status="EMPTY!"
        fi
        
        printf "%-25s %-10s %-12d %-12s %-10s\n" \
            "$filename" "$filesize" "$total_reads" "${paired_pct}%" "$status"
    fi
done

echo ""
echo "Detailed Check for First File:"
echo "-------------------------------"

# Detailed check on first treatment file
FIRST_BAM="${BAM_DIR}/TEAD1-1_filtered.bam"
if [ -f "$FIRST_BAM" ]; then
    echo "File: $(basename $FIRST_BAM)"
    echo ""
    echo "Flagstat output:"
    samtools flagstat $FIRST_BAM
    echo ""
    echo "First 3 reads:"
    samtools view $FIRST_BAM | head -3 | cut -f1-6
fi