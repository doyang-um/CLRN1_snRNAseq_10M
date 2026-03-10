#!/usr/bin/env Rscript
# =============================================================================
# Generate Figure 5: Differential Expression Volcano Plots
# - Volcano plots for each cell type (Rods, Cones, Müller Glia, etc.)
# - Bar chart showing number of DEGs per cell type
# =============================================================================

library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(tidyr)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Set paths
base_dir <- "/home/doyang/turbo/CLRN1_snRNAseq_GitHub"
results_dir <- file.path(base_dir, "results/DESeq2")
output_dir <- "/home/doyang/turbo/CLRN1_figures_regenerated"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define cell types (order matters for plotting)
# These match your actual file names!
cell_types <- c("Rods", "Cones", "Muller", 
                "Amacrine", "RGC", "Bipolar", 
                "Horizontal")

# Display names
cell_type_labels <- c(
  "Rods" = "Rods",
  "Cones" = "Cones",
  "Muller" = "Müller Glia",
  "Amacrine" = "Amacrine",
  "RGC" = "RGC",
  "Bipolar" = "Bipolar Cells",
  "Horizontal" = "Horizontal Cells"
)

# Colors for volcano plots
volcano_colors <- c(
  "Down in KO" = "#4DBBD5",    # Blue
  "Up in KO" = "#E64B35",       # Red
  "Not Significant" = "#B0B0B0" # Gray
)

# =============================================================================
# FUNCTION: CREATE VOLCANO PLOT
# =============================================================================

create_volcano_plot <- function(de_results, cell_type_name, n_genes_to_label = 10) {
  
  # Calculate -log10(p-value)
  de_results$neg_log10_p <- -log10(de_results$padj)
  
  # Replace infinite values
  max_p <- max(de_results$neg_log10_p[!is.infinite(de_results$neg_log10_p)], na.rm = TRUE)
  de_results$neg_log10_p[is.infinite(de_results$neg_log10_p)] <- max_p + 2
  
  # Define significance
  de_results$significance <- "Not Significant"
  de_results$significance[de_results$padj < 0.05 & de_results$log2FoldChange < 0] <- "Down in KO"
  de_results$significance[de_results$padj < 0.05 & de_results$log2FoldChange > 0] <- "Up in KO"
  
  de_results$significance <- factor(de_results$significance, 
                                     levels = c("Down in KO", "Not Significant", "Up in KO"))
  
  # Count DEGs
  n_degs <- sum(de_results$padj < 0.05, na.rm = TRUE)
  
  # Select top genes to label (by significance)
  sig_genes <- de_results %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    head(n_genes_to_label)
  
  # Create plot
  p <- ggplot(de_results, aes(x = log2FoldChange, y = neg_log10_p, color = significance)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = volcano_colors) +
    
    # Add threshold lines
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    
    # Add gene labels
    geom_text_repel(
      data = sig_genes,
      aes(label = gene),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      color = "black",
      segment.size = 0.2
    ) +
    
    # Labels and theme
    labs(
      title = paste0(cell_type_labels[cell_type_name], " (", n_degs, " DEGs)"),
      x = "Log2 Fold Change (KO/WT)",
      y = "-log10(p-value)"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  return(list(plot = p, n_degs = n_degs))
}

# =============================================================================
# LOAD DATA AND CREATE VOLCANO PLOTS
# =============================================================================

cat("=== LOADING DESEQ2 RESULTS ===\n")

volcano_plots <- list()
deg_counts <- data.frame(
  CellType = character(),
  N_DEGs = numeric(),
  stringsAsFactors = FALSE
)

for (ct in cell_types) {
  
  # Construct file name (adjust based on your actual file names)
  file_name <- paste0(ct, "_Pseudobulk_DESeq2_Results.csv")
  file_path <- file.path(results_dir, file_name)
  
  # Check if file exists
  if (!file.exists(file_path)) {
    cat("  WARNING: File not found:", file_name, "\n")
    cat("  Trying alternative names...\n")
    
    # Try alternative naming
    alt_names <- c(
      paste0(cell_type_labels[ct], "_Pseudobulk_DESeq2_Results.csv"),
      paste0(gsub("_", " ", ct), "_Pseudobulk_DESeq2_Results.csv")
    )
    
    found <- FALSE
    for (alt_name in alt_names) {
      alt_path <- file.path(results_dir, alt_name)
      if (file.exists(alt_path)) {
        file_path <- alt_path
        found <- TRUE
        cat("  Found:", alt_name, "\n")
        break
      }
    }
    
    if (!found) {
      cat("  SKIPPING:", ct, "\n\n")
      next
    }
  }
  
  # Load data
  cat("  Loading:", basename(file_path), "\n")
  de_data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Check if gene column exists
  if (!"gene" %in% colnames(de_data)) {
    if ("Gene" %in% colnames(de_data)) {
      de_data$gene <- de_data$Gene
    } else if (rownames(de_data)[1] != "1") {
      de_data$gene <- rownames(de_data)
    } else {
      de_data$gene <- de_data[, 1]
    }
  }
  
  # Check for required columns
  required_cols <- c("log2FoldChange", "padj")
  if (!all(required_cols %in% colnames(de_data))) {
    cat("  ERROR: Missing required columns in", basename(file_path), "\n")
    cat("  Available columns:", paste(colnames(de_data), collapse = ", "), "\n\n")
    next
  }
  
  # Remove NAs
  de_data <- de_data %>%
    filter(!is.na(padj), !is.na(log2FoldChange))
  
  # Create volcano plot
  result <- create_volcano_plot(de_data, ct)
  volcano_plots[[ct]] <- result$plot
  
  # Store DEG count
  deg_counts <- rbind(deg_counts, data.frame(
    CellType = cell_type_labels[ct],
    N_DEGs = result$n_degs
  ))
  
  cat("  DEGs found:", result$n_degs, "\n\n")
}

# =============================================================================
# CREATE DEG SUMMARY BAR CHART
# =============================================================================

cat("=== CREATING DEG SUMMARY BAR CHART ===\n")

# Order by number of DEGs
deg_counts$CellType <- factor(deg_counts$CellType, 
                               levels = deg_counts$CellType[order(deg_counts$N_DEGs, decreasing = TRUE)])

# Create bar chart
deg_barplot <- ggplot(deg_counts, aes(x = CellType, y = N_DEGs, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = N_DEGs), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c(
    "Rods" = "#4DBBD5",
    "Amacrine" = "#B565A7",
    "Cones" = "#E64B35",
    "Horizontal Cells" = "#E89AC7",
    "RGC" = "#8B4513",
    "Müller Glia" = "#00A087",
    "Bipolar Cells" = "#3C5488"
  )) +
  labs(
    title = "Differentially Expressed Genes per Cell Type",
    subtitle = "Pseudobulk DESeq2 (padj < 0.05)",
    x = "",
    y = "Number of DEGs"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  ylim(0, max(deg_counts$N_DEGs) * 1.15)

cat("✓ Bar chart created\n\n")

# =============================================================================
# COMBINE ALL PLOTS
# =============================================================================

cat("=== COMBINING PLOTS ===\n")

# Arrange volcano plots in 3x3 grid with bar chart in bottom right
combined_figure <- (
  volcano_plots[["Rods"]] | volcano_plots[["Cones"]] | volcano_plots[["Muller"]]
) / (
  volcano_plots[["Amacrine"]] | volcano_plots[["RGC"]] | volcano_plots[["Bipolar"]]
) / (
  volcano_plots[["Horizontal"]] | deg_barplot | plot_spacer()
) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    title = "Differential Gene Expression Analysis: CLRN1 KO vs WT",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

cat("✓ Figure combined\n\n")

# =============================================================================
# SAVE FIGURE
# =============================================================================

cat("=== SAVING FIGURE ===\n")

# Save combined figure
ggsave(file.path(output_dir, "Figure5_Volcano_Plots_Complete.pdf"),
       combined_figure, width = 14, height = 12)
ggsave(file.path(output_dir, "Figure5_Volcano_Plots_Complete.png"),
       combined_figure, width = 14, height = 12, dpi = 300)

cat("✓ Combined figure saved\n")

# Save individual volcano plots
for (ct in names(volcano_plots)) {
  file_name <- paste0("Volcano_", ct, ".pdf")
  ggsave(file.path(output_dir, file_name),
         volcano_plots[[ct]], width = 5, height = 5)
}

cat("✓ Individual volcano plots saved\n")

# Save bar chart
ggsave(file.path(output_dir, "DEG_Summary_BarChart.pdf"),
       deg_barplot, width = 6, height = 5)

cat("✓ Bar chart saved\n\n")

# =============================================================================
# SAVE DEG SUMMARY TABLE
# =============================================================================

write.csv(deg_counts, 
          file.path(output_dir, "DEG_Summary_Table.csv"),
          row.names = FALSE)

cat("✓ Summary table saved\n\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat(rep("=", 72), "\n", sep="")
cat("FIGURE 5 GENERATION COMPLETE\n")
cat(rep("=", 72), "\n\n")

cat("Output directory:", output_dir, "\n\n")

cat("FILES CREATED:\n")
cat("1. Figure5_Volcano_Plots_Complete.pdf/png - Main figure\n")
cat("2. Volcano_[CellType].pdf - Individual volcano plots\n")
cat("3. DEG_Summary_BarChart.pdf - Bar chart only\n")
cat("4. DEG_Summary_Table.csv - DEG counts per cell type\n\n")

cat("DEG SUMMARY:\n")
print(deg_counts)
cat("\n")

cat("Total DEGs across all cell types:", sum(deg_counts$N_DEGs), "\n\n")

cat(rep("=", 72), "\n", sep="")
cat("READY FOR MANUSCRIPT!\n")
cat(rep("=", 72), "\n")
