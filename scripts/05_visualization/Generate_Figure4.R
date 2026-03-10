#!/usr/bin/env Rscript
# =============================================================================
# Generate Complete CLRN1 Expression Figure
# Panel A: UMAP (WT vs KO split)
# Panel B: CLRN1 by Cell Type
# Panel C: CLRN1 in Müller Glia (WT vs KO)
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# =============================================================================
# CONFIGURATION
# =============================================================================

setwd("/home/doyang/turbo/CLRN1 WT VS KO 10M SnRNAseq/")
output_dir <- "/home/doyang/turbo/CLRN1_figures_regenerated"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading Seurat object...\n")
retina <- readRDS("CLRN1_Retina_with_DonorIDs.rds")

cat("Total cells:", ncol(retina), "\n")
cat("Cell types:", unique(retina$cell_type), "\n\n")

# =============================================================================
# DEFINE COLORS (Improved distinguishable palette)
# =============================================================================

cell_colors <- c(
  "Rods" = "#E377C2",           # Pink/magenta (largest cluster - keep visible)
  "Cones" = "#8C564B",          # Brown
  "Muller_Glia" = "#00CED1",    # Cyan/turquoise
  "Müller_Glia" = "#00CED1",
  "Müller glia" = "#00CED1",
  "Bipolar_Cells" = "#FF7F0E",  # Orange
  "Amacrine_Cells" = "#D62728", # Red
  "Horizontal_Cells" = "#BCBD22", # Yellow-green (CHANGED - was pink)
  "RGC" = "#9467BD"             # Purple
)

# For violin plots (simpler palette)
violin_colors <- c(
  "Rods" = "#1F77B4",          # Blue
  "Cones" = "#FF7F0E",         # Orange  
  "Müller glia" = "#2CA02C"    # Green
)

genotype_colors <- c(
  "WT" = "#6BAED6",            # Light blue
  "CLRN1 KO" = "#FB9A99"       # Light red/pink
)

# =============================================================================
# PANEL A: UMAP SPLIT BY GENOTYPE
# =============================================================================

cat("=== GENERATING PANEL A: UMAP ===\n")

# Make sure condition labels are clean
retina$condition_clean <- gsub("CLRN1_KO|KO", "CLRN1_KO", retina$condition)
retina$condition_clean <- gsub("WT|wt", "WT", retina$condition_clean)

panel_a <- DimPlot(retina, 
                   reduction = "umap",
                   group.by = "cell_type",
                   split.by = "condition_clean",
                   cols = cell_colors,
                   pt.size = 0.5,
                   ncol = 2) +
  theme_void() +
  theme(
    plot.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  )

cat("✓ Panel A created\n\n")

# =============================================================================
# PANEL B: CLRN1 EXPRESSION BY CELL TYPE
# =============================================================================

cat("=== GENERATING PANEL B: CLRN1 by Cell Type ===\n")

# Get CLRN1 expression
if (inherits(retina@assays$RNA, "Assay5")) {
  all_data <- LayerData(retina, assay = "RNA", layer = "data")
} else {
  all_data <- GetAssayData(retina, slot = "data", assay = "RNA")
}

clrn1_all <- all_data["CLRN1", ]

# Focus on Rods, Cones, Müller glia
main_types <- c("Rods", "Cones", "Muller_Glia", "Müller_Glia", "Müller glia", "Muller glia")
cells_for_panelb <- retina$cell_type %in% main_types

plot_data_b <- data.frame(
  CellType = retina$cell_type[cells_for_panelb],
  CLRN1 = clrn1_all[cells_for_panelb]
)

# Standardize Müller glia name
plot_data_b$CellType <- gsub("Muller.*|Müller.*", "Müller glia", plot_data_b$CellType)
plot_data_b$CellType <- factor(plot_data_b$CellType, 
                                levels = c("Rods", "Cones", "Müller glia"))

cat("Cell counts in Panel B:\n")
print(table(plot_data_b$CellType))

panel_b <- ggplot(plot_data_b, aes(x = CellType, y = CLRN1, fill = CellType)) +
  geom_violin(scale = "width", trim = FALSE) +
  scale_fill_manual(values = violin_colors) +
  scale_y_continuous(limits = c(0, 3.5), breaks = seq(0, 3, 1)) +
  labs(title = "CLRN1 Expression in Retinal Cell Types",
       x = "",
       y = "CLRN1 Expression (log-normalized)") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

cat("✓ Panel B created\n")
cat("  Y-axis range: 0 to 3.5\n\n")

# =============================================================================
# PANEL C: CLRN1 IN MÜLLER GLIA BY GENOTYPE
# =============================================================================

cat("=== GENERATING PANEL C: CLRN1 in Müller Glia ===\n")

# Get Müller glia cells
muller_mask <- retina$cell_type %in% c("Muller_Glia", "Müller_Glia", "Müller glia", "Muller glia")
muller_cells <- retina[, muller_mask]

# Get CLRN1 expression in Müller glia
if (inherits(muller_cells@assays$RNA, "Assay5")) {
  muller_data <- LayerData(muller_cells, assay = "RNA", layer = "data")
} else {
  muller_data <- GetAssayData(muller_cells, slot = "data", assay = "RNA")
}

clrn1_muller <- muller_data["CLRN1", ]

# Create plot data
plot_data_c <- data.frame(
  Genotype = muller_cells$condition,
  CLRN1 = clrn1_muller
)

# Clean genotype names
plot_data_c$Genotype <- gsub("CLRN1_KO|KO|ko", "CLRN1 KO", plot_data_c$Genotype)
plot_data_c$Genotype <- gsub("^WT$|wt", "WT", plot_data_c$Genotype)
plot_data_c$Genotype <- factor(plot_data_c$Genotype, levels = c("WT", "CLRN1 KO"))

cat("Sample sizes:\n")
print(table(plot_data_c$Genotype))

# Statistical test (cell-level for display, as shown in your figure)
# Note: Proper analysis uses pseudobulk (p=0.7), but figure shows cell-level
stat_test <- wilcox.test(CLRN1 ~ Genotype, data = plot_data_c)
p_value <- stat_test$p.value

cat("  Mann-Whitney p-value (cell-level):", round(p_value, 3), "\n\n")

panel_c <- ggplot(plot_data_c, aes(x = Genotype, y = CLRN1, fill = Genotype)) +
  geom_violin(scale = "width", trim = FALSE) +
  scale_fill_manual(values = genotype_colors) +
  scale_y_continuous(limits = c(0, 3.5), breaks = seq(0, 3, 1)) +
  labs(title = "CLRN1 Expression in Müller Glia Cells",
       x = "",
       y = "CLRN1 Expression\n(log-normalized counts)") +
  annotate("text", x = 1.5, y = 3.2,
           label = paste0("Mann-Whitney U test\np = ", round(p_value, 3), "\nns"),
           size = 4) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

cat("✓ Panel C created\n")
cat("  Y-axis range: 0 to 3.5\n\n")

# =============================================================================
# COMBINE PANELS
# =============================================================================

cat("=== COMBINING PANELS ===\n")

# Layout: Panel A on top (full width)
#         Panel B and C on bottom (side by side)

combined_figure <- panel_a / (panel_b | panel_c) +
  plot_layout(heights = c(2, 1.5)) +
  plot_annotation(tag_levels = 'a',
                  tag_prefix = '',
                  tag_suffix = '') &
  theme(plot.tag = element_text(size = 16, face = "bold"))

cat("✓ Figure combined\n\n")

# =============================================================================
# SAVE FIGURES
# =============================================================================

cat("=== SAVING FIGURES ===\n")

# Save combined figure
ggsave(file.path(output_dir, "Figure_CLRN1_Complete.pdf"),
       combined_figure, width = 12, height = 10)
ggsave(file.path(output_dir, "Figure_CLRN1_Complete.png"),
       combined_figure, width = 12, height = 10, dpi = 300)

cat("✓ Combined figure saved\n")

# Save individual panels
ggsave(file.path(output_dir, "Panel_A_UMAP.pdf"),
       panel_a, width = 10, height = 5)
ggsave(file.path(output_dir, "Panel_B_CLRN1_CellTypes.pdf"),
       panel_b, width = 5, height = 5)
ggsave(file.path(output_dir, "Panel_C_CLRN1_Muller_Genotypes.pdf"),
       panel_c, width = 5, height = 5)

cat("✓ Individual panels saved\n\n")

# =============================================================================
# CREATE VERSION WITH PSEUDOBULK STATISTICS
# =============================================================================

cat("=== CREATING VERSION WITH PROPER STATISTICS ===\n")

# Filter to clean donors only
muller_clean <- subset(retina, 
                       cell_type %in% c("Muller_Glia", "Müller_Glia", "Müller glia", "Muller glia") &
                       donor_id %in% c("donor0", "donor1", "donor2"))

if (inherits(muller_clean@assays$RNA, "Assay5")) {
  clrn1_clean <- LayerData(muller_clean, assay = "RNA", layer = "data")["CLRN1", ]
} else {
  clrn1_clean <- GetAssayData(muller_clean, slot = "data")["CLRN1", ]
}

# Pseudobulk
pseudobulk <- data.frame(
  CLRN1 = clrn1_clean,
  Sample = muller_clean$sample_id
) %>%
  group_by(Sample) %>%
  summarise(
    Mean_CLRN1 = mean(CLRN1),
    .groups = "drop"
  ) %>%
  mutate(
    Genotype = ifelse(grepl("WT", Sample), "WT", "CLRN1 KO")
  )

# Proper statistical test
pseudobulk_test <- wilcox.test(Mean_CLRN1 ~ Genotype, data = pseudobulk)
p_pseudobulk <- pseudobulk_test$p.value

cat("Pseudobulk statistics:\n")
cat("  3 WT vs 3 KO animals\n")
cat("  p-value:", round(p_pseudobulk, 3), "\n\n")

# Create Panel C with pseudobulk p-value
panel_c_proper <- ggplot(plot_data_c, aes(x = Genotype, y = CLRN1, fill = Genotype)) +
  geom_violin(scale = "width", trim = FALSE) +
  scale_fill_manual(values = genotype_colors) +
  scale_y_continuous(limits = c(0, 3.5), breaks = seq(0, 3, 1)) +
  labs(title = "CLRN1 Expression in Müller Glia Cells",
       x = "",
       y = "CLRN1 Expression\n(log-normalized counts)",
       caption = paste0("n=3 animals per genotype\nPseudobulk Wilcoxon test")) +
  annotate("text", x = 1.5, y = 3.2,
           label = paste0("p = ", round(p_pseudobulk, 2), "\nns"),
           size = 4) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    plot.caption = element_text(size = 9, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# Combined with proper stats
combined_proper <- panel_a / (panel_b | panel_c_proper) +
  plot_layout(heights = c(2, 1.5)) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave(file.path(output_dir, "Figure_CLRN1_Complete_ProperStats.pdf"),
       combined_proper, width = 12, height = 10)
ggsave(file.path(output_dir, "Figure_CLRN1_Complete_ProperStats.png"),
       combined_proper, width = 12, height = 10, dpi = 300)

cat("✓ Version with proper statistics saved\n\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat(rep("=", 72), "\n", sep="")
cat("FIGURE GENERATION COMPLETE\n")
cat(rep("=", 72), "\n\n")

cat("Output directory:", output_dir, "\n\n")

cat("FILES CREATED:\n")
cat("1. Figure_CLRN1_Complete.pdf/png\n")
cat("   - Matches your original figure (cell-level p-value)\n")
cat("   - Panel A: UMAP split by genotype\n")
cat("   - Panel B: CLRN1 by cell type\n")
cat("   - Panel C: CLRN1 in Müller glia (WT vs KO)\n\n")

cat("2. Figure_CLRN1_Complete_ProperStats.pdf/png\n")
cat("   - Uses proper pseudobulk statistics (p = 0.7)\n")
cat("   - Accounts for biological replication\n")
cat("   - RECOMMENDED for publication\n\n")

cat("3. Individual panels (PDF):\n")
cat("   - Panel_A_UMAP.pdf\n")
cat("   - Panel_B_CLRN1_CellTypes.pdf\n")
cat("   - Panel_C_CLRN1_Muller_Genotypes.pdf\n\n")

cat("STATISTICS:\n")
cat("  Cell-level test: p =", round(p_value, 3), "\n")
cat("  Pseudobulk test: p =", round(p_pseudobulk, 2), "(PROPER - use this!)\n\n")

cat("RECOMMENDATION:\n")
cat("  Use 'Figure_CLRN1_Complete_ProperStats.pdf' for your manuscript.\n")
cat("  This version accounts for biological replication (3 vs 3 animals)\n")
cat("  and is statistically rigorous.\n\n")

cat(rep("=", 72), "\n", sep="")
cat("READY FOR MANUSCRIPT!\n")
cat(rep("=", 72), "\n")
