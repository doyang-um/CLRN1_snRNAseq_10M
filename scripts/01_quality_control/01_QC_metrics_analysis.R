# ===================================================================
# Quality Control (QC) Metrics Analysis - snRNA-seq
# ===================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)

# Load your Seurat object if needed
# retina <- readRDS("CLRN1_Retina_CellType_Annotated_Filtered.rds")

# ===================================================================
# Extract QC Metrics
# ===================================================================

# Get metadata with QC metrics
qc_data <- retina@meta.data

# Add sample labels for clarity
qc_data$Sample <- ifelse(qc_data$condition == "WT", "WT", "CLRN1 KO")

cat("=== QC Metrics Summary ===\n")
cat("\nTotal nuclei analyzed:", nrow(qc_data), "\n")
cat("WT nuclei:", sum(qc_data$condition == "WT"), "\n")
cat("KO nuclei:", sum(qc_data$condition == "CLRN1_KO"), "\n")

# Summary statistics
qc_summary <- qc_data %>%
  group_by(Sample) %>%
  summarise(
    n_nuclei = n(),
    mean_nCount = mean(nCount_RNA),
    median_nCount = median(nCount_RNA),
    mean_nFeature = mean(nFeature_RNA),
    median_nFeature = median(nFeature_RNA),
    mean_percent_mt = mean(percent.mt),
    median_percent_mt = median(percent.mt),
    .groups = 'drop'
  )

print(qc_summary)
write.csv(qc_summary, "QC_Summary_Statistics.csv", row.names = FALSE)

# ===================================================================
# 1. UMI Counts per Nucleus (nCount_RNA)
# ===================================================================

# Violin plot - UMI counts
p1 <- ggplot(qc_data, aes(x = Sample, y = nCount_RNA, fill = Sample)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_y_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "UMI Counts per Nucleus",
    x = "Sample",
    y = "UMI Counts (log10 scale)"
  ) +
  scale_fill_manual(values = c("WT" = "#619CFF", "CLRN1 KO" = "#F8766D")) +
  stat_summary(fun = median, geom = "point", size = 3, color = "black")

ggsave("QC_UMI_Counts_Violin.pdf", p1, width = 6, height = 6)

# Histogram - UMI distribution
p2 <- ggplot(qc_data, aes(x = nCount_RNA, fill = Sample)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_x_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    title = "Distribution of UMI Counts",
    x = "UMI Counts (log10 scale)",
    y = "Number of Nuclei"
  ) +
  scale_fill_manual(values = c("WT" = "#619CFF", "CLRN1 KO" = "#F8766D"))

ggsave("QC_UMI_Counts_Histogram.pdf", p2, width = 8, height = 6)

# ===================================================================
# 2. Gene Counts per Nucleus (nFeature_RNA)
# ===================================================================

# Violin plot - Gene counts
p3 <- ggplot(qc_data, aes(x = Sample, y = nFeature_RNA, fill = Sample)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_y_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Genes Detected per Nucleus",
    x = "Sample",
    y = "Number of Genes (log10 scale)"
  ) +
  scale_fill_manual(values = c("WT" = "#619CFF", "CLRN1 KO" = "#F8766D")) +
  stat_summary(fun = median, geom = "point", size = 3, color = "black")

ggsave("QC_Gene_Counts_Violin.pdf", p3, width = 6, height = 6)

# Histogram - Gene distribution
p4 <- ggplot(qc_data, aes(x = nFeature_RNA, fill = Sample)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_x_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    title = "Distribution of Genes per Nucleus",
    x = "Number of Genes (log10 scale)",
    y = "Number of Nuclei"
  ) +
  scale_fill_manual(values = c("WT" = "#619CFF", "CLRN1 KO" = "#F8766D"))

ggsave("QC_Gene_Counts_Histogram.pdf", p4, width = 8, height = 6)

# ===================================================================
# 3. Mitochondrial Percentage
# ===================================================================

# Violin plot - Mitochondrial percentage
p5 <- ggplot(qc_data, aes(x = Sample, y = percent.mt, fill = Sample)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Mitochondrial Gene Percentage",
    x = "Sample",
    y = "% Mitochondrial Genes"
  ) +
  scale_fill_manual(values = c("WT" = "#619CFF", "CLRN1 KO" = "#F8766D")) +
  stat_summary(fun = median, geom = "point", size = 3, color = "black") +
  geom_hline(yintercept = 15, linetype = "dashed", color = "red", size = 0.5)

ggsave("QC_Mitochondrial_Percent_Violin.pdf", p5, width = 6, height = 6)

# Histogram - Mitochondrial distribution
p6 <- ggplot(qc_data, aes(x = percent.mt, fill = Sample)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    title = "Distribution of Mitochondrial Percentage",
    x = "% Mitochondrial Genes",
    y = "Number of Nuclei"
  ) +
  scale_fill_manual(values = c("WT" = "#619CFF", "CLRN1 KO" = "#F8766D")) +
  geom_vline(xintercept = 15, linetype = "dashed", color = "red", size = 0.5)

ggsave("QC_Mitochondrial_Percent_Histogram.pdf", p6, width = 8, height = 6)

# ===================================================================
# 4. Correlation Plots
# ===================================================================

# UMI vs Gene count (colored by sample)
p7 <- ggplot(qc_data, aes(x = nCount_RNA, y = nFeature_RNA, color = Sample)) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    title = "UMI Counts vs Genes Detected",
    x = "UMI Counts (log10 scale)",
    y = "Number of Genes (log10 scale)"
  ) +
  scale_color_manual(values = c("WT" = "#619CFF", "CLRN1 KO" = "#F8766D")) +
  geom_smooth(method = "lm", se = FALSE)

ggsave("QC_UMI_vs_Genes_Scatter.pdf", p7, width = 8, height = 6)

# Mitochondrial % vs UMI count
p8 <- ggplot(qc_data, aes(x = nCount_RNA, y = percent.mt, color = Sample)) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_x_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    title = "UMI Counts vs Mitochondrial Percentage",
    x = "UMI Counts (log10 scale)",
    y = "% Mitochondrial Genes"
  ) +
  scale_color_manual(values = c("WT" = "#619CFF", "CLRN1 KO" = "#F8766D")) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "red", size = 0.5)

ggsave("QC_UMI_vs_Mito_Scatter.pdf", p8, width = 8, height = 6)

# ===================================================================
# 5. Combined QC Panel (Publication Figure)
# ===================================================================

# Create a comprehensive panel
combined_qc <- (p1 | p3 | p5) / (p2 | p4 | p6)

ggsave("QC_Combined_Panel.pdf", combined_qc, width = 18, height = 10)

# Alternative layout - 2x2 grid with key metrics
key_qc <- (p1 | p3) / (p5 | p7)

ggsave("QC_Key_Metrics_Panel.pdf", key_qc, width = 12, height = 10)

# ===================================================================
# 6. QC Metrics by Cell Type
# ===================================================================

# Prepare data
qc_celltype <- qc_data %>%
  filter(!is.na(cell_type))

# UMI counts by cell type
p9 <- ggplot(qc_celltype, aes(x = cell_type, y = nCount_RNA, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  scale_y_log10() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    title = "UMI Counts by Cell Type",
    x = "Cell Type",
    y = "UMI Counts (log10 scale)"
  ) +
  scale_fill_manual(values = c("WT" = "#619CFF", "CLRN1 KO" = "#F8766D"))

ggsave("QC_UMI_by_CellType.pdf", p9, width = 12, height = 6)

# Genes by cell type
p10 <- ggplot(qc_celltype, aes(x = cell_type, y = nFeature_RNA, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  scale_y_log10() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    title = "Genes Detected by Cell Type",
    x = "Cell Type",
    y = "Number of Genes (log10 scale)"
  ) +
  scale_fill_manual(values = c("WT" = "#619CFF", "CLRN1 KO" = "#F8766D"))

ggsave("QC_Genes_by_CellType.pdf", p10, width = 12, height = 6)

# Mitochondrial % by cell type
p11 <- ggplot(qc_celltype, aes(x = cell_type, y = percent.mt, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    title = "Mitochondrial % by Cell Type",
    x = "Cell Type",
    y = "% Mitochondrial Genes"
  ) +
  scale_fill_manual(values = c("WT" = "#619CFF", "CLRN1 KO" = "#F8766D")) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "red", size = 0.5)

ggsave("QC_Mito_by_CellType.pdf", p11, width = 12, height = 6)

# ===================================================================
# 7. Detailed Statistics Table
# ===================================================================

# Overall statistics
overall_stats <- data.frame(
  Metric = c(
    "Total Nuclei",
    "WT Nuclei", 
    "KO Nuclei",
    "Mean UMI (WT)",
    "Mean UMI (KO)",
    "Median UMI (WT)",
    "Median UMI (KO)",
    "Mean Genes (WT)",
    "Mean Genes (KO)",
    "Median Genes (WT)",
    "Median Genes (KO)",
    "Mean % MT (WT)",
    "Mean % MT (KO)",
    "Median % MT (WT)",
    "Median % MT (KO)"
  ),
  Value = c(
    nrow(qc_data),
    sum(qc_data$condition == "WT"),
    sum(qc_data$condition == "CLRN1_KO"),
    round(mean(qc_data$nCount_RNA[qc_data$condition == "WT"]), 1),
    round(mean(qc_data$nCount_RNA[qc_data$condition == "CLRN1_KO"]), 1),
    round(median(qc_data$nCount_RNA[qc_data$condition == "WT"]), 1),
    round(median(qc_data$nCount_RNA[qc_data$condition == "CLRN1_KO"]), 1),
    round(mean(qc_data$nFeature_RNA[qc_data$condition == "WT"]), 1),
    round(mean(qc_data$nFeature_RNA[qc_data$condition == "CLRN1_KO"]), 1),
    round(median(qc_data$nFeature_RNA[qc_data$condition == "WT"]), 1),
    round(median(qc_data$nFeature_RNA[qc_data$condition == "CLRN1_KO"]), 1),
    round(mean(qc_data$percent.mt[qc_data$condition == "WT"]), 2),
    round(mean(qc_data$percent.mt[qc_data$condition == "CLRN1_KO"]), 2),
    round(median(qc_data$percent.mt[qc_data$condition == "WT"]), 2),
    round(median(qc_data$percent.mt[qc_data$condition == "CLRN1_KO"]), 2)
  )
)

write.csv(overall_stats, "QC_Overall_Statistics.csv", row.names = FALSE)

# Per cell type statistics
celltype_stats <- qc_celltype %>%
  group_by(cell_type, Sample) %>%
  summarise(
    n_nuclei = n(),
    mean_UMI = round(mean(nCount_RNA), 1),
    median_UMI = round(median(nCount_RNA), 1),
    mean_genes = round(mean(nFeature_RNA), 1),
    median_genes = round(median(nFeature_RNA), 1),
    mean_mito = round(mean(percent.mt), 2),
    median_mito = round(median(percent.mt), 2),
    .groups = 'drop'
  )

write.csv(celltype_stats, "QC_CellType_Statistics.csv", row.names = FALSE)

# ===================================================================
# 8. Summary Report
# ===================================================================

cat("\n===============================================\n")
cat("QC METRICS SUMMARY REPORT\n")
cat("===============================================\n\n")

cat("Dataset Overview:\n")
cat("  Total nuclei:", nrow(qc_data), "\n")
cat("  WT nuclei:", sum(qc_data$condition == "WT"), "\n")
cat("  CLRN1 KO nuclei:", sum(qc_data$condition == "CLRN1_KO"), "\n\n")

cat("UMI Counts (per nucleus):\n")
cat("  WT - Mean:", round(mean(qc_data$nCount_RNA[qc_data$condition == "WT"]), 1), 
    "  Median:", round(median(qc_data$nCount_RNA[qc_data$condition == "WT"]), 1), "\n")
cat("  KO - Mean:", round(mean(qc_data$nCount_RNA[qc_data$condition == "CLRN1_KO"]), 1),
    "  Median:", round(median(qc_data$nCount_RNA[qc_data$condition == "CLRN1_KO"]), 1), "\n\n")

cat("Genes Detected (per nucleus):\n")
cat("  WT - Mean:", round(mean(qc_data$nFeature_RNA[qc_data$condition == "WT"]), 1),
    "  Median:", round(median(qc_data$nFeature_RNA[qc_data$condition == "WT"]), 1), "\n")
cat("  KO - Mean:", round(mean(qc_data$nFeature_RNA[qc_data$condition == "CLRN1_KO"]), 1),
    "  Median:", round(median(qc_data$nFeature_RNA[qc_data$condition == "CLRN1_KO"]), 1), "\n\n")

cat("Mitochondrial Percentage:\n")
cat("  WT - Mean:", round(mean(qc_data$percent.mt[qc_data$condition == "WT"]), 2), "%",
    "  Median:", round(median(qc_data$percent.mt[qc_data$condition == "WT"]), 2), "%\n")
cat("  KO - Mean:", round(mean(qc_data$percent.mt[qc_data$condition == "CLRN1_KO"]), 2), "%",
    "  Median:", round(median(qc_data$percent.mt[qc_data$condition == "CLRN1_KO"]), 2), "%\n\n")

cat("Quality Assessment:\n")
high_mito <- sum(qc_data$percent.mt > 15)
cat("  Nuclei with >15% mitochondrial reads:", high_mito, 
    "(", round(high_mito/nrow(qc_data)*100, 1), "%)\n")

cat("\n===============================================\n")
cat("OUTPUT FILES GENERATED:\n")
cat("===============================================\n")
cat("Individual Metrics:\n")
cat("  - QC_UMI_Counts_Violin.pdf\n")
cat("  - QC_Gene_Counts_Violin.pdf\n")
cat("  - QC_Mitochondrial_Percent_Violin.pdf\n")
cat("  - QC_UMI_vs_Genes_Scatter.pdf\n")
cat("  - QC_UMI_vs_Mito_Scatter.pdf\n\n")
cat("Combined Panels:\n")
cat("  - QC_Combined_Panel.pdf (comprehensive)\n")
cat("  - QC_Key_Metrics_Panel.pdf (publication-ready)\n\n")
cat("Cell Type Analysis:\n")
cat("  - QC_UMI_by_CellType.pdf\n")
cat("  - QC_Genes_by_CellType.pdf\n")
cat("  - QC_Mito_by_CellType.pdf\n\n")
cat("Data Tables:\n")
cat("  - QC_Summary_Statistics.csv\n")
cat("  - QC_Overall_Statistics.csv\n")
cat("  - QC_CellType_Statistics.csv\n\n")

cat("Analysis complete!\n")