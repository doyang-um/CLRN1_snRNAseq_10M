# ===================================================================
# Volcano Plots - Müller Glia and Cones with Specific Gene Highlighting
# ===================================================================

library(ggplot2)
library(ggrepel)
library(tidyverse)

# ===================================================================
# Load DE Results
# ===================================================================

# Load DE results for Müller Glia and Cones
muller_de <- read.csv("DE_Genes_Muller_Glia_KO_vs_WT.csv")
cones_de <- read.csv("DE_Genes_Cones_KO_vs_WT.csv")

cat("=== DE Results Loaded ===\n")
cat("Müller Glia DE genes:", nrow(muller_de), "\n")
cat("Cones DE genes:", nrow(cones_de), "\n\n")

# ===================================================================
# Define Genes of Interest
# ===================================================================

# Genes to highlight
downregulated_genes <- c("CTNNA2", "TENM2")
upregulated_genes <- c("CRACR2A", "DLG4")
all_genes_of_interest <- c(downregulated_genes, upregulated_genes)

cat("Genes to highlight:\n")
cat("  Downregulated:", paste(downregulated_genes, collapse = ", "), "\n")
cat("  Upregulated:", paste(upregulated_genes, collapse = ", "), "\n\n")

# ===================================================================
# Prepare Data for Volcano Plot
# ===================================================================

prepare_volcano_data <- function(de_data, cell_type) {
  
  # Add gene column if not present
  if (!"gene" %in% colnames(de_data)) {
    de_data$gene <- de_data$Gene
  }
  
  # Calculate -log10(p-value)
  de_data$neg_log10_p <- -log10(de_data$p_val_adj)
  
  # Replace infinite values
  de_data$neg_log10_p[is.infinite(de_data$neg_log10_p)] <- max(de_data$neg_log10_p[!is.infinite(de_data$neg_log10_p)]) + 10
  
  # Create significance categories
  de_data$significance <- "Not Significant"
  de_data$significance[de_data$p_val_adj < 0.05 & de_data$avg_log2FC > 0.5] <- "Upregulated"
  de_data$significance[de_data$p_val_adj < 0.05 & de_data$avg_log2FC < -0.5] <- "Downregulated"
  
  # Mark genes of interest
  de_data$gene_of_interest <- "Other"
  de_data$gene_of_interest[de_data$gene %in% downregulated_genes] <- "Downregulated (Highlighted)"
  de_data$gene_of_interest[de_data$gene %in% upregulated_genes] <- "Upregulated (Highlighted)"
  
  # Label for plotting (only genes of interest)
  de_data$label <- ""
  de_data$label[de_data$gene %in% all_genes_of_interest] <- de_data$gene[de_data$gene %in% all_genes_of_interest]
  
  # Check which genes were found
  found_genes <- de_data$gene[de_data$gene %in% all_genes_of_interest]
  missing_genes <- all_genes_of_interest[!all_genes_of_interest %in% de_data$gene]
  
  cat(cell_type, ":\n")
  cat("  Found genes:", paste(found_genes, collapse = ", "), "\n")
  if (length(missing_genes) > 0) {
    cat("  Missing genes:", paste(missing_genes, collapse = ", "), "\n")
  }
  cat("\n")
  
  return(de_data)
}

# Prepare data
muller_volcano <- prepare_volcano_data(muller_de, "Müller Glia")
cones_volcano <- prepare_volcano_data(cones_de, "Cones")

# ===================================================================
# Function to Create Volcano Plot
# ===================================================================

create_volcano_plot <- function(volcano_data, cell_type, filename) {
  
  # Define colors
  colors <- c(
    "Not Significant" = "grey70",
    "Upregulated" = "#FC4E07",
    "Downregulated" = "#00AFBB",
    "Upregulated (Highlighted)" = "#E41A1C",
    "Downregulated (Highlighted)" = "#377EB8"
  )
  
  # Create base plot
  p <- ggplot(volcano_data, aes(x = avg_log2FC, y = neg_log10_p)) +
    # Plot all points first (grey and colored)
    geom_point(data = subset(volcano_data, significance == "Not Significant"),
               aes(color = significance), alpha = 0.3, size = 1) +
    geom_point(data = subset(volcano_data, significance == "Upregulated"),
               aes(color = significance), alpha = 0.5, size = 1.5) +
    geom_point(data = subset(volcano_data, significance == "Downregulated"),
               aes(color = significance), alpha = 0.5, size = 1.5) +
    # Highlight genes of interest with larger, more prominent points
    geom_point(data = subset(volcano_data, gene_of_interest == "Upregulated (Highlighted)"),
               aes(color = gene_of_interest), size = 4, shape = 19) +
    geom_point(data = subset(volcano_data, gene_of_interest == "Downregulated (Highlighted)"),
               aes(color = gene_of_interest), size = 4, shape = 19) +
    # Add gene labels
    geom_text_repel(aes(label = label),
                    size = 5,
                    fontface = "bold.italic",
                    box.padding = 0.5,
                    point.padding = 0.5,
                    segment.color = "black",
                    segment.size = 0.5,
                    max.overlaps = 20,
                    min.segment.length = 0) +
    # Add threshold lines
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", size = 0.5) +
    # Colors
    scale_color_manual(values = colors,
                       name = "Category",
                       breaks = c("Downregulated (Highlighted)", 
                                  "Upregulated (Highlighted)",
                                  "Downregulated", 
                                  "Upregulated", 
                                  "Not Significant"),
                       labels = c("Downregulated (Key)", 
                                  "Upregulated (Key)",
                                  "Downregulated", 
                                  "Upregulated", 
                                  "Not Significant")) +
    # Theme and labels
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      panel.grid.major = element_line(color = "grey95"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste(cell_type, "- Differential Expression"),
      subtitle = "CLRN1 KO vs WT",
      x = "Log2 Fold Change",
      y = "-Log10 (Adjusted P-value)"
    ) +
    # Add count annotations
    annotate("text", x = min(volcano_data$avg_log2FC) * 0.8, 
             y = max(volcano_data$neg_log10_p) * 0.95,
             label = paste0("Downregulated: ", 
                            sum(volcano_data$significance == "Downregulated")),
             hjust = 0, size = 4, fontface = "bold", color = "#00AFBB") +
    annotate("text", x = max(volcano_data$avg_log2FC) * 0.5, 
             y = max(volcano_data$neg_log10_p) * 0.95,
             label = paste0("Upregulated: ", 
                            sum(volcano_data$significance == "Upregulated")),
             hjust = 0, size = 4, fontface = "bold", color = "#FC4E07")
  
  # Save plot
  ggsave(filename, p, width = 10, height = 8)
  
  return(p)
}

# ===================================================================
# Create Volcano Plots
# ===================================================================

cat("Creating volcano plots...\n\n")

# Müller Glia volcano plot
p_muller <- create_volcano_plot(muller_volcano, "Müller Glia", 
                                "Volcano_Muller_Glia_Highlighted.pdf")

# Cones volcano plot
p_cones <- create_volcano_plot(cones_volcano, "Cones",
                               "Volcano_Cones_Highlighted.pdf")

# ===================================================================
# Create Side-by-Side Comparison
# ===================================================================

library(patchwork)

# Create simplified versions for side-by-side
create_simplified_volcano <- function(volcano_data, cell_type) {
  
  colors <- c(
    "Not Significant" = "grey70",
    "Upregulated" = "#FC4E07",
    "Downregulated" = "#00AFBB",
    "Upregulated (Highlighted)" = "#E41A1C",
    "Downregulated (Highlighted)" = "#377EB8"
  )
  
  p <- ggplot(volcano_data, aes(x = avg_log2FC, y = neg_log10_p)) +
    geom_point(data = subset(volcano_data, significance == "Not Significant"),
               aes(color = significance), alpha = 0.3, size = 0.8) +
    geom_point(data = subset(volcano_data, significance == "Upregulated"),
               aes(color = significance), alpha = 0.5, size = 1) +
    geom_point(data = subset(volcano_data, significance == "Downregulated"),
               aes(color = significance), alpha = 0.5, size = 1) +
    geom_point(data = subset(volcano_data, gene_of_interest == "Upregulated (Highlighted)"),
               aes(color = gene_of_interest), size = 3, shape = 19) +
    geom_point(data = subset(volcano_data, gene_of_interest == "Downregulated (Highlighted)"),
               aes(color = gene_of_interest), size = 3, shape = 19) +
    geom_text_repel(aes(label = label),
                    size = 4,
                    fontface = "bold.italic",
                    box.padding = 0.3,
                    max.overlaps = 20) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.4) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", size = 0.4) +
    scale_color_manual(values = colors) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.position = "none",
      panel.grid.major = element_line(color = "grey95")
    ) +
    labs(
      title = cell_type,
      x = "Log2 Fold Change",
      y = "-Log10 (Adj. P-value)"
    )
  
  return(p)
}

p_muller_simple <- create_simplified_volcano(muller_volcano, "Müller Glia")
p_cones_simple <- create_simplified_volcano(cones_volcano, "Cones")

combined_volcano <- p_muller_simple | p_cones_simple

ggsave("Volcano_Muller_Cones_Combined.pdf", combined_volcano, width = 16, height = 7)

# ===================================================================
# Extract Statistics for Highlighted Genes
# ===================================================================

extract_gene_stats <- function(volcano_data, cell_type) {
  
  gene_stats <- volcano_data %>%
    filter(gene %in% all_genes_of_interest) %>%
    select(gene, avg_log2FC, p_val_adj, pct.1, pct.2) %>%
    mutate(
      Cell_Type = cell_type,
      Direction = ifelse(avg_log2FC > 0, "Upregulated", "Downregulated")
    ) %>%
    arrange(avg_log2FC)
  
  return(gene_stats)
}

muller_stats <- extract_gene_stats(muller_volcano, "Müller Glia")
cones_stats <- extract_gene_stats(cones_volcano, "Cones")

all_stats <- rbind(muller_stats, cones_stats)

write.csv(all_stats, "Highlighted_Genes_Statistics.csv", row.names = FALSE)

cat("=== Statistics for Highlighted Genes ===\n")
print(all_stats)
cat("\n")

# ===================================================================
# Create Summary Table
# ===================================================================

summary_table <- data.frame(
  Cell_Type = character(),
  Total_DE_Genes = integer(),
  Upregulated = integer(),
  Downregulated = integer(),
  Highlighted_Found = character(),
  stringsAsFactors = FALSE
)

# Müller Glia
muller_summary <- data.frame(
  Cell_Type = "Müller Glia",
  Total_DE_Genes = sum(muller_volcano$p_val_adj < 0.05 & abs(muller_volcano$avg_log2FC) > 0.5),
  Upregulated = sum(muller_volcano$significance == "Upregulated"),
  Downregulated = sum(muller_volcano$significance == "Downregulated"),
  Highlighted_Found = paste(muller_volcano$gene[muller_volcano$gene %in% all_genes_of_interest], collapse = ", ")
)

# Cones
cones_summary <- data.frame(
  Cell_Type = "Cones",
  Total_DE_Genes = sum(cones_volcano$p_val_adj < 0.05 & abs(cones_volcano$avg_log2FC) > 0.5),
  Upregulated = sum(cones_volcano$significance == "Upregulated"),
  Downregulated = sum(cones_volcano$significance == "Downregulated"),
  Highlighted_Found = paste(cones_volcano$gene[cones_volcano$gene %in% all_genes_of_interest], collapse = ", ")
)

summary_table <- rbind(muller_summary, cones_summary)

write.csv(summary_table, "Volcano_Plot_Summary.csv", row.names = FALSE)

cat("=== Volcano Plot Summary ===\n")
print(summary_table)
cat("\n")

# ===================================================================
# Summary Report
# ===================================================================

cat("===============================================\n")
cat("VOLCANO PLOT ANALYSIS COMPLETE\n")
cat("===============================================\n\n")

cat("Highlighted Genes:\n")
cat("  Downregulated: ", paste(downregulated_genes, collapse = ", "), "\n")
cat("  Upregulated: ", paste(upregulated_genes, collapse = ", "), "\n\n")

cat("OUTPUT FILES GENERATED:\n")
cat("  - Volcano_Muller_Glia_Highlighted.pdf\n")
cat("  - Volcano_Cones_Highlighted.pdf\n")
cat("  - Volcano_Muller_Cones_Combined.pdf (side-by-side)\n")
cat("  - Highlighted_Genes_Statistics.csv\n")
cat("  - Volcano_Plot_Summary.csv\n\n")

cat("Key Findings:\n")
for (i in 1:nrow(all_stats)) {
  gene <- all_stats$gene[i]
  ct <- all_stats$Cell_Type[i]
  fc <- round(all_stats$avg_log2FC[i], 2)
  pval <- format.pval(all_stats$p_val_adj[i], digits = 2)
  direction <- all_stats$Direction[i]
  
  cat("  ", gene, "in", ct, ":", direction, "(Log2FC =", fc, ", p =", pval, ")\n")
}

cat("\nAnalysis complete!\n")