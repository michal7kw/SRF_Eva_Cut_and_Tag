library(VennDiagram)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Pass counts from shell script, ensuring they are not empty
tes_total=18404
tesmut_total=18404
tead1_total=18404

tes_tesmut_intersection=10562
tes_tead1_intersection=10562
tesmut_tead1_intersection=10562

tes_unique_tesmut=9303
tesmut_unique_tes=8020
tes_unique_tead1=9303
tead1_unique_tes=8020
tesmut_unique_tead1=9303
tead1_unique_tesmut=8020


cat(paste("Peak counts:\n"))
cat(paste("TES:", tes_total, "\n"))
cat(paste("TESmut:", tesmut_total, "\n"))
cat(paste("TEAD1:", tead1_total, "\n\n"))

# Create Venn diagrams for all pairwise comparisons
output_dir <- "results/12_pairwise_overlap/venn"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. TES vs TESmut
cat("Creating TES vs TESmut Venn diagram...\n")
png(file.path(output_dir, "TES_vs_TESmut_consensus_broad_venn.png"),
    width = 3000, height = 3000, res = 300) # Set resolution and dimensions
venn.plot1 <- draw.pairwise.venn(
    area1 = tes_total,
    area2 = tesmut_total,
    cross.area = tes_tesmut_intersection,
    category = c("TES", "TESmut"),
    col = "transparent",
    fill = c("lightblue", "lightcoral"),
    alpha = 0.6,
    cex = 1.5,
    fontfamily = "serif",
    main = "TES vs TESmut consensus_broad Peak Overlap",
    main.cex = 1.8
)
grid.draw(venn.plot1)
dev.off()

# 2. TES vs TEAD1
cat("Creating TES vs TEAD1 Venn diagram...\n")
png(file.path(output_dir, "TES_vs_TEAD1_consensus_broad_venn.png"),
    width = 3000, height = 3000, res = 300) # Set resolution and dimensions
venn.plot2 <- draw.pairwise.venn(
    area1 = tes_total,
    area2 = tead1_total,
    cross.area = tes_tead1_intersection,
    category = c("TES", "TEAD1"),
    col = "transparent",
    fill = c("lightblue", "lightgreen"),
    alpha = 0.6,
    cex = 1.5,
    fontfamily = "serif",
    main = "TES vs TEAD1 consensus_broad Peak Overlap",
    main.cex = 1.8
)
grid.draw(venn.plot2)
dev.off()

# 3. TESmut vs TEAD1
cat("Creating TESmut vs TEAD1 Venn diagram...\n")
png(file.path(output_dir, "TESmut_vs_TEAD1_consensus_broad_venn.png"),
    width = 3000, height = 3000, res = 300) # Set resolution and dimensions
venn.plot3 <- draw.pairwise.venn(
    area1 = tesmut_total,
    area2 = tead1_total,
    cross.area = tesmut_tead1_intersection,
    category = c("TESmut", "TEAD1"),
    col = "transparent",
    fill = c("lightcoral", "lightgreen"),
    alpha = 0.6,
    cex = 1.5,
    fontfamily = "serif",
    main = "TESmut vs TEAD1 consensus_broad Peak Overlap",
    main.cex = 1.8
)
grid.draw(venn.plot3)
dev.off()

# Calculate detailed overlap statistics
cat("Calculating overlap statistics...\n")

# Create summary statistics table
stats_table <- data.frame(
    Comparison = c("TES vs TESmut", "TES vs TEAD1", "TESmut vs TEAD1"),
    Overlap = c(tes_tesmut_intersection, tes_tead1_intersection, tesmut_tead1_intersection),
    First_Only = c(tes_unique_tesmut, tes_unique_tead1, tesmut_unique_tead1),
    Second_Only = c(tesmut_unique_tes, tead1_unique_tes, tead1_unique_tesmut),
    First_Total = c(tes_total, tes_total, tesmut_total),
    Second_Total = c(tesmut_total, tead1_total, tead1_total),
    Overlap_Percent_First = round(c(tes_tesmut_intersection/tes_total*100,
                                   tes_tead1_intersection/tes_total*100,
                                   tesmut_tead1_intersection/tesmut_total*100), 1),
    Overlap_Percent_Second = round(c(tes_tesmut_intersection/tesmut_total*100,
                                    tes_tead1_intersection/tead1_total*100,
                                    tesmut_tead1_intersection/tead1_total*100), 1)
)

# Write statistics to file
write.table(stats_table, file.path(output_dir, "consensus_broad_peak_overlap_statistics.txt"), 
            sep="\t", row.names=FALSE, quote=FALSE)

cat("Venn diagram analysis completed.\n")
cat("Files created in:", file.path(output_dir), "\n")

# Print summary to console
print(stats_table)
