#!/usr/bin/env Rscript
# =============================================================================
# Supplementary Figure 5: Quality Control and Sample Metadata
# - Panel A: QC metrics by genotype (violin plots)
# - Panel B: UMAP by donor ID
# - Panel C: Cell counts per donor
# - Panel D: Cell type distribution table
# - Panel E: Marker gene expression dot plot
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(gridExtra)
library(RColorBrewer)

# =============================================================================
# CONFIGURATION
# =============================================================================

setwd("/home/doyang/turbo/CLRN1 WT VS KO 10M SnRNAseq/")
output_dir <- "Figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load Seurat object with donor IDs
cat("Loading Seurat object...\n")
retina <- readRDS("CLRN1_Retina_with_DonorIDs.rds")

# Filter out doublets and unassigned
cat("Filtering doublets and unassigned cells...\n")
retina_clean <- subset(retina, subset = donor_id %in% c("donor0", "donor1", "donor2"))

cat("Cells after filtering:", ncol(retina_clean), "\n")

# =============================================================================
# PANEL A: QC METRICS BY GENOTYPE
# =============================================================================

cat("\n=== Creating Panel A: QC Metrics ===\n")

# Prepare data
qc_data <- retina_clean@meta.data %>%
  select(condition, nFeature_RNA, nCount_RNA, percent.mt)

# Create individual violin plots
p_nfeature <- ggplot(qc_data, aes(x = condition, y = nFeature_RNA, fill = condition)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("CLRN1_KO" = "#E64B35", "WT" = "#4DBBD5")) +
  labs(title = "nFeature_RNA", x = "Identity", y = "") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 10))

p_ncount <- ggplot(qc_data, aes(x = condition, y = nCount_RNA, fill = condition)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("CLRN1_KO" = "#E64B35", "WT" = "#4DBBD5")) +
  labs(title = "nCount_RNA", x = "Identity", y = "") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 10))

p_percent_mt <- ggplot(qc_data, aes(x = condition, y = percent.mt, fill = condition)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("CLRN1_KO" = "#E64B35", "WT" = "#4DBBD5")) +
  labs(title = "percent.mt", x = "Identity", y = "") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 10))

panel_a <- p_nfeature | p_ncount | p_percent_mt

cat("✓ Panel A created\n")

# =============================================================================
# PANEL B: UMAP BY DONOR ID
# =============================================================================

cat("\n=== Creating Panel B: UMAP by Donor ===\n")

# Define donor colors
donor_colors <- c(
  "donor0" = "#E64B35",
  "donor1" = "#4DBBD5", 
  "donor2" = "#00A087",
  "doublet" = "#F39B7F",
  "unassigned" = "#B0B0B0"
)

panel_b <- DimPlot(retina_clean, 
                   group.by = "donor_id",
                   cols = donor_colors,
                   pt.size = 0.5) +
  labs(title = "") +
  theme(legend.position = "right",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9))

cat("✓ Panel B created\n")

# =============================================================================
# PANEL C: CELL COUNTS PER DONOR
# =============================================================================

cat("\n=== Creating Panel C: Cell Counts per Donor ===\n")

# Count cells per donor and genotype
donor_counts <- retina_clean@meta.data %>%
  group_by(donor_id, condition) %>%
  summarise(count = n(), .groups = 'drop')

# Add donor labels
donor_counts$donor_label <- paste0(donor_counts$donor_id, "\n(n=", donor_counts$count, ")")

panel_c <- ggplot(donor_counts, aes(x = donor_id, y = count, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = count), position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("CLRN1_KO" = "#E64B35", "WT" = "#4DBBD5"),
                    name = "Genotype",
                    labels = c("CLRN1_KO" = "KO", "WT" = "WT")) +
  labs(title = "Cell Counts per Donor",
       x = "Donor",
       y = "Number of Cells") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top")

cat("✓ Panel C created\n")

# =============================================================================
# PANEL D: CELL TYPE DISTRIBUTION TABLE
# =============================================================================

cat("\n=== Creating Panel D: Cell Type Distribution Table ===\n")

# Create cell type by sample table
sample_celltype <- retina_clean@meta.data %>%
  mutate(sample = paste0(condition, "_", donor_id)) %>%
  group_by(sample, cell_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = cell_type, values_from = count, values_fill = 0)

# Add totals
sample_celltype$Total <- rowSums(sample_celltype[, -1])

# Create a table plot
table_grob <- tableGrob(sample_celltype, rows = NULL, 
                        theme = ttheme_default(base_size = 8,
                                               core = list(fg_params = list(hjust = 1, x = 0.95)),
                                               colhead = list(fg_params = list(fontface = "bold"))))

panel_d <- ggplot() + 
  annotation_custom(table_grob) +
  theme_void()

cat("✓ Panel D created\n")

# =============================================================================
# PANEL E: MARKER GENE DOT PLOT
# =============================================================================

cat("\n=== Creating Panel E: Marker Gene Expression ===\n")

# Define marker genes for each cell type
marker_genes <- c(
  # Rods
  "RHO", "NR2E3", "ROM1",
  # Cones  
  "ARR3", "OPN1SW", "OPN1MW",
  # Müller glia
  "RLBP1", "GLUL", "SLC1A3",
  # Bipolar
  "VSX2", "CABP5", "TRPM1",
  # Amacrine
  "GAD1", "SLC6A9", "CHAT",
  # RGC
  "RBPMS", "THY1", "NEFL",
  # Horizontal
  "ONECUT1", "LHX1", "PROX1"
)

# Check which markers are present
markers_present <- marker_genes[marker_genes %in% rownames(retina_clean)]
cat("Using", length(markers_present), "of", length(marker_genes), "marker genes\n")

# Create dot plot
panel_e <- DotPlot(retina_clean, 
                   features = markers_present,
                   group.by = "cell_type",
                   cols = c("lightgrey", "blue"),
                   dot.scale = 8) +
  coord_flip() +
  labs(x = "Features", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.position = "right",
        legend.box = "vertical")

cat("✓ Panel E created\n")

# =============================================================================
# COMBINE ALL PANELS
# =============================================================================

cat("\n=== Combining all panels ===\n")

# Create layout
# Top row: Panel A (QC violins) and Panel B (UMAP)
top_row <- panel_a | panel_b

# Middle row: Panel C (donor counts) and Panel D (table)
middle_row <- panel_c | panel_d

# Bottom row: Panel E (dot plot) - full width
bottom_row <- panel_e

# Combine all
supp_fig5 <- top_row / middle_row / bottom_row +
  plot_layout(heights = c(1, 1, 1.2)) +
  plot_annotation(
    title = "Supplementary Figure 5: Quality Control and Sample Metadata",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

cat("✓ Figure combined\n")

# =============================================================================
# SAVE FIGURE
# =============================================================================

cat("\n=== Saving figure ===\n")

ggsave(file.path(output_dir, "SuppFig5_QC_Metadata_Complete.pdf"),
       supp_fig5, width = 14, height = 12, device = cairo_pdf)

ggsave(file.path(output_dir, "SuppFig5_QC_Metadata_Complete.png"),
       supp_fig5, width = 14, height = 12, dpi = 300)

cat("✓ Figure saved\n")

# Save individual panels
ggsave(file.path(output_dir, "SuppFig5_PanelA_QC.pdf"), 
       panel_a, width = 8, height = 3, device = cairo_pdf)
ggsave(file.path(output_dir, "SuppFig5_PanelB_UMAP_Donor.pdf"), 
       panel_b, width = 5, height = 4, device = cairo_pdf)
ggsave(file.path(output_dir, "SuppFig5_PanelC_DonorCounts.pdf"), 
       panel_c, width = 5, height = 4, device = cairo_pdf)
ggsave(file.path(output_dir, "SuppFig5_PanelE_Markers.pdf"), 
       panel_e, width = 8, height = 6, device = cairo_pdf)

cat("✓ Individual panels saved\n")

# =============================================================================
# SAVE METADATA STATISTICS
# =============================================================================

cat("\n=== Saving metadata statistics ===\n")

# QC statistics
qc_stats <- retina_clean@meta.data %>%
  group_by(condition) %>%
  summarise(
    n_cells = n(),
    mean_nFeature = mean(nFeature_RNA),
    median_nFeature = median(nFeature_RNA),
    mean_nCount = mean(nCount_RNA),
    median_nCount = median(nCount_RNA),
    mean_percent_mt = mean(percent.mt),
    median_percent_mt = median(percent.mt)
  )

write.csv(qc_stats, file.path(output_dir, "SuppFig5_QC_Statistics.csv"), 
          row.names = FALSE)

# Cell type counts
celltype_stats <- retina_clean@meta.data %>%
  group_by(cell_type, condition) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = condition, values_from = count, values_fill = 0)

write.csv(celltype_stats, file.path(output_dir, "SuppFig5_CellType_Counts.csv"), 
          row.names = FALSE)

# Donor distribution
write.csv(sample_celltype, file.path(output_dir, "SuppFig5_Donor_Distribution.csv"), 
          row.names = FALSE)

cat("✓ Statistics saved\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n", rep("=", 72), "\n", sep="")
cat("SUPPLEMENTARY FIGURE 5 GENERATION COMPLETE\n")
cat(rep("=", 72), "\n\n", sep="")

cat("OUTPUT FILES:\n")
cat("1. SuppFig5_QC_Metadata_Complete.pdf/png - Full figure\n")
cat("2. SuppFig5_PanelA_QC.pdf - QC violin plots\n")
cat("3. SuppFig5_PanelB_UMAP_Donor.pdf - UMAP by donor\n")
cat("4. SuppFig5_PanelC_DonorCounts.pdf - Cell counts\n")
cat("5. SuppFig5_PanelE_Markers.pdf - Marker expression\n")
cat("6. SuppFig5_*.csv - Statistics tables\n\n")

cat("SUMMARY STATISTICS:\n")
print(qc_stats)
cat("\n")

cat("CELL TYPE DISTRIBUTION:\n")
print(celltype_stats)
cat("\n")

cat(rep("=", 72), "\n", sep="")
cat("READY FOR MANUSCRIPT!\n")
cat(rep("=", 72), "\n")
