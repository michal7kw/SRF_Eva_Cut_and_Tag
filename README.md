# SRF_Eva Cut&Tag Analysis Pipeline

**Current Version:** 2.0 (Dec 2024)
**Maintainer:** kubacki.michal@hsr.it

## 🚀 Quick Start

### Run Complete Pipeline
The master script handles everything from FastQ to final diffbind analysis.
```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG
./master_rerun_all.sh
```

### Common Maintenance Commands
```bash
# Check job status
squeue -u $USER

# Monitor logs
tail -f logs/*.out
```

---

## 📂 Pipeline Structure

Scripts are in the main directory and run sequentially.

| Step | Script Name | Purpose | Output Location |
|:---|:---|:---|:---|
| **01** | `1_qc.sh` | Quality control of raw reads | `results/01_fastqc/` |
| **02** | `2_adapter_trim.sh` | Adapter trimming | `results/02_trimmed/` |
| **03** | `3_align.sh` | Bowtie2 alignment | `results/03_aligned/` |
| **04** | `4_filter.sh` | Deduplication & QC filtering | `results/04_filtered/` |
| **05** | `5_peak_calling_narrow.sh` | MACS2 Narrow Peak Calling | `results/05_peaks_narrow/` |
| **06** | `6_bigwig.sh` | Generate coverage tracks | `results/06_bigwig/` |
| **07** | `7_compare_narrow.sh` | Comparison logic setup | (Intermediate) |
| **08** | `8_annotate_narrow.sh` | Peak Annotation | `results/07_analysis_narrow/` |
| **09** | `9_diff_bind_narrow.sh` | Differential Binding (DiffBind) | `results/07_analysis_narrow/` |
| **11** | `11_combine_replicates...` | Merge replicate peaks | `results/11_combined_replicates.../` |
| **12** | `12_pairwise_overlap...` | Peak overlap stats | `results/12_pairwise_overlap.../` |
| **13** | `13_ses_overlap.sh` | Super-Enhancer integration | `results/13_ses_overlap/` |
| **16** | `16_run_additional...` | IDR & Fragment size | `results/idr_analysis/` |

---

## 🛠️ Configuration & Customization

### Key Analysis Groups
- **Control**: `Nestin_Ctrl`, `Emx1_Ctrl` (IgG controls usually handled securely)
- **Experimental**: `TES`, `TEAD1`, `YAP1`

### Troubleshooting
-   **Missing Dependencies**: Ensure `conda activate cutntag` is active.
-   **Job Failures**: Check `logs/slurm-[jobid].err`.
-   **DiffBind Errors**: Check `results/07_analysis_narrow/DiffBind_*.pdf` for QC failure clues.

---

## 📊 Key Results Checklist
1.  **DiffBind Report**: `results/07_analysis_narrow/DiffBind_report.csv`
2.  **QC Report**: `results/multiqc/multiqc_report.html`
3.  **Peak Overlaps**: `results/12_pairwise_overlap_narrow/pairwise_overlap_heatmap.pdf`
4.  **SE Integration**: `results/13_ses_overlap/plots/`
