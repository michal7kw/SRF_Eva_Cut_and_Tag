library(VennDiagram)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Peak counts from shell script
tes_total=17796
tesmut_total=13871
tead1_total=11163

# Direct overlap counts (peaks with any overlap)
tes_overlap_tesmut=6577
tesmut_overlap_tes=6540
tes_overlap_tead1=4263
tead1_overlap_tes=4220
tesmut_overlap_tead1=6248
tead1_overlap_tesmut=6192

# For Venn diagrams, use overlap counts (peaks that have any overlap)
tes_tesmut_overlap=6577
tes_tead1_overlap=4263
tesmut_tead1_overlap=6248

tes_unique_tesmut=11219
tesmut_unique_tes=7331
tes_unique_tead1=13533
tead1_unique_tes=6943
tesmut_unique_tead1=7623
tead1_unique_tesmut=4971

# Summit-based overlaps
summit_tes_tesmut=6381
summit_tes_tead1=4024
summit_tesmut_tead1=5997

# Reciprocal overlaps (>=50%)
reciprocal_tes_tesmut=5037
reciprocal_tes_tead1=2792
reciprocal_tesmut_tead1=4456

# Jaccard indices
jaccard_tes_tesmut=0.262
jaccard_tes_tead1=0.173
jaccard_tesmut_tead1=0.333

cat(paste("Narrow Peak counts:\n"))
cat(paste("TES:", tes_total, "\n"))
cat(paste("TESmut:", tesmut_total, "\n"))
cat(paste("TEAD1:", tead1_total, "\n\n"))

cat(paste("Direct Overlaps:\n"))
cat(paste("TES overlapping TESmut:", tes_overlap_tesmut, "\n"))
cat(paste("TES overlapping TEAD1:", tes_overlap_tead1, "\n"))
cat(paste("TESmut overlapping TEAD1:", tesmut_overlap_tead1, "\n\n"))

# Create Venn diagrams for all pairwise comparisons
output_dir <- "results/12_pairwise_overlap_narrow/venn"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. TES vs TESmut
cat("Creating TES vs TESmut Venn diagram...\n")
png(file.path(output_dir, "TES_vs_TESmut_narrow_venn.png"),
    width = 3000, height = 3000, res = 300)
venn.plot1 <- draw.pairwise.venn(
    area1 = tes_total,
    area2 = tesmut_total,
    cross.area = tes_tesmut_overlap,
    category = c("TES", "TESmut"),
    col = "transparent",
    fill = c("lightblue", "lightcoral"),
    alpha = 0.6,
    cex = 1.5,
    fontfamily = "serif",
    main = "TES vs TESmut Narrow Peak Overlap",
    main.cex = 1.8
)
grid.draw(venn.plot1)
dev.off()

# 2. TES vs TEAD1
cat("Creating TES vs TEAD1 Venn diagram...\n")
png(file.path(output_dir, "TES_vs_TEAD1_narrow_venn.png"),
    width = 3000, height = 3000, res = 300)
venn.plot2 <- draw.pairwise.venn(
    area1 = tes_total,
    area2 = tead1_total,
    cross.area = tes_tead1_overlap,
    category = c("TES", "TEAD1"),
    col = "transparent",
    fill = c("lightblue", "lightgreen"),
    alpha = 0.6,
    cex = 1.5,
    fontfamily = "serif",
    main = "TES vs TEAD1 Narrow Peak Overlap",
    main.cex = 1.8
)
grid.draw(venn.plot2)
dev.off()

# 3. TESmut vs TEAD1
cat("Creating TESmut vs TEAD1 Venn diagram...\n")
png(file.path(output_dir, "TESmut_vs_TEAD1_narrow_venn.png"),
    width = 3000, height = 3000, res = 300)
venn.plot3 <- draw.pairwise.venn(
    area1 = tesmut_total,
    area2 = tead1_total,
    cross.area = tesmut_tead1_overlap,
    category = c("TESmut", "TEAD1"),
    col = "transparent",
    fill = c("lightcoral", "lightgreen"),
    alpha = 0.6,
    cex = 1.5,
    fontfamily = "serif",
    main = "TESmut vs TEAD1 Narrow Peak Overlap",
    main.cex = 1.8
)
grid.draw(venn.plot3)
dev.off()

# Create comprehensive summary statistics table with multiple overlap methods
stats_table <- data.frame(
    Comparison = c("TES vs TESmut", "TES vs TEAD1", "TESmut vs TEAD1"),
    Direct_Overlap = c(tes_tesmut_overlap, tes_tead1_overlap, tesmut_tead1_overlap),
    Summit_Overlap = c(summit_tes_tesmut, summit_tes_tead1, summit_tesmut_tead1),
    Reciprocal_50pct = c(reciprocal_tes_tesmut, reciprocal_tes_tead1, reciprocal_tesmut_tead1),
    First_Only = c(tes_unique_tesmut, tes_unique_tead1, tesmut_unique_tead1),
    Second_Only = c(tesmut_unique_tes, tead1_unique_tes, tead1_unique_tesmut),
    First_Total = c(tes_total, tes_total, tesmut_total),
    Second_Total = c(tesmut_total, tead1_total, tead1_total),
    Overlap_Pct_First = round(c(tes_tesmut_overlap/tes_total*100,
                                tes_tead1_overlap/tes_total*100,
                                tesmut_tead1_overlap/tesmut_total*100), 1),
    Overlap_Pct_Second = round(c(tes_tesmut_overlap/tesmut_total*100,
                                 tes_tead1_overlap/tead1_total*100,
                                 tesmut_tead1_overlap/tead1_total*100), 1),
    Jaccard_Index = c(jaccard_tes_tesmut, jaccard_tes_tead1, jaccard_tesmut_tead1)
)

# Write statistics to file
write.table(stats_table, file.path(output_dir, "narrow_peak_overlap_statistics.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)

# Create bar plot comparing overlap methods
overlap_comparison <- data.frame(
    Comparison = rep(c("TES vs TESmut", "TES vs TEAD1", "TESmut vs TEAD1"), 3),
    Method = rep(c("Direct Overlap", "Summit-based (200bp)", "Reciprocal (50%)"), each = 3),
    Count = c(tes_tesmut_overlap, tes_tead1_overlap, tesmut_tead1_overlap,
              summit_tes_tesmut, summit_tes_tead1, summit_tesmut_tead1,
              reciprocal_tes_tesmut, reciprocal_tes_tead1, reciprocal_tesmut_tead1)
)

p_methods <- ggplot(overlap_comparison, aes(x = Comparison, y = Count, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Comparison of Overlap Methods for Narrow Peaks",
         y = "Number of Overlapping Peaks", x = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set2")

ggsave(file.path(output_dir, "overlap_methods_comparison.pdf"), p_methods, width = 10, height = 6)
ggsave(file.path(output_dir, "overlap_methods_comparison.png"), p_methods, width = 10, height = 6, dpi = 300)

cat("Narrow peak Venn diagram analysis completed.\n")
cat("\nSummary Statistics:\n")
print(stats_table)

cat("\nInterpretation Guide:\n")
cat("- Direct Overlap: Any genomic overlap (most permissive)\n")
cat("- Summit-based: Peak centers within 200bp (appropriate for TF binding)\n")
cat("- Reciprocal 50%: At least 50% mutual overlap (most stringent)\n")
cat("- Jaccard Index: Overlap similarity score (0-1, higher = more similar)\n")
