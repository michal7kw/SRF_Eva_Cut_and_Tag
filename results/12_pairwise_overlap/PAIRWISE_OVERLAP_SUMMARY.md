# Pairwise Binding Site Overlap Analysis Summary

**Generated on:** Wed Sep 17 05:02:34 PM CEST 2025  
**Peak Type Used:** consensus_broad peaks (CONSENSUS_BROAD)  
**Peak Files Source:** results/11_combined_replicates_broad/peaks

## Overview

This analysis compares binding sites between TES, TESmut, and TEAD1 conditions using pairwise overlap analysis with consensus_broad peaks. Note: Consider using broad or consensus peaks for better overlap sensitivity.

## Peak Counts by Condition

- **TES**: 27605 peaks
- **TESmut**: 19864 peaks
- **TEAD1**: 18404 peaks

## Pairwise Overlap Results

### TES vs TESmut
- TES peaks overlapping with TESmut: **11597**
- TESmut peaks overlapping with TES: **11597**
- TES-specific peaks (not in TESmut): **16087**
- TESmut-specific peaks (not in TES): **8361**

### TES vs TEAD1
- TES peaks overlapping with TEAD1: **8885**
- TEAD1 peaks overlapping with TES: **8885**
- TES-specific peaks (not in TEAD1): **18735**
- TEAD1-specific peaks (not in TES): **9649**

### TESmut vs TEAD1
- TESmut peaks overlapping with TEAD1: **10562**
- TEAD1 peaks overlapping with TESmut: **10562**
- TESmut-specific peaks (not in TEAD1): **9303**
- TEAD1-specific peaks (not in TESmut): **8020**

## Files Generated

### Overlap Files
- 
: Peak regions that overlap between conditions
- 
: Condition-specific peaks with no overlap

### Visualization Files
- 
: Venn diagrams showing overlap relationships
- 
: Additional visualization plots

### Statistics Files
- 
: Detailed numerical overlap statistics
- 
: This comprehensive summary

### Analysis Files
- 
: Signal matrices for overlap regions (if BigWig files available)

## Interpretation Notes

1. **TES vs TESmut**: Comparison shows binding specificity of TES construct
2. **TES vs TEAD1**: Reveals similarity between TES and endogenous TEAD1 binding
3. **TESmut vs TEAD1**: Control comparison to validate TESmut background levels

## Next Steps

1. Review Venn diagrams for visual overlap assessment
2. Analyze signal intensity at overlap regions using provided matrices
3. Perform functional annotation of condition-specific peak sets
4. Consider motif analysis in overlap vs unique regions


## Peak Type Notes

- **Current Analysis**: Used consensus_broad peaks from results/11_combined_replicates_broad/peaks
- **Recommendation**: Consider running broad peak analysis (5b_broad_peak_calling.sh) for improved overlap sensitivity
- **Alternative Available**: Run 5b_broad_peak_calling.sh to generate broad peaks

