# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Set working directory
setwd("/home/doyang/turbo/CLRN1 WT VS KO 10M SnRNAseq/")

# Load the Seurat object with donor IDs
retina <- readRDS("CLRN1_Retina_with_DonorIDs.rds")

# Check it loaded
retina

# Define a nice color palette for cell types
cell_colors <- c(
  "Rods" = "#1f77b4",
  "Cones" = "#ff7f0e", 
  "Muller_Glia" = "#2ca02c",
  "Bipolar_Cells" = "#d62728",
  "Amacrine_Cells" = "#9467bd",
  "RGC" = "#8c564b",
  "Horizontal_Cells" = "#e377c2"
)

# ============================================
# FIGURE 1: Cell Type Identification
# ============================================

# Panel A: UMAP colored by cell type
fig1a <- DimPlot(retina, reduction = "umap", group.by = "cell_type", 
                 cols = cell_colors, pt.size = 0.5) +
  ggtitle("A. Cell Type Identification") +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14))

# Panel B: UMAP split by condition
fig1b <- DimPlot(retina, reduction = "umap", group.by = "cell_type",
                 split.by = "condition", cols = cell_colors, pt.size = 0.3) +
  ggtitle("B. WT vs CLRN1 KO") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

# Display
fig1a
fig1b


# Create a figures directory
dir.create("Figures", showWarnings = FALSE)

# Save Figure 1A
ggsave("Figures/Fig1A_UMAP_CellTypes.pdf", fig1a, width = 8, height = 6)
ggsave("Figures/Fig1A_UMAP_CellTypes.png", fig1a, width = 8, height = 6, dpi = 300)

# Save Figure 1B
ggsave("Figures/Fig1B_UMAP_WT_vs_KO.pdf", fig1b, width = 12, height = 5)
ggsave("Figures/Fig1B_UMAP_WT_vs_KO.png", fig1b, width = 12, height = 5, dpi = 300)

cat("Figures saved to:", getwd(), "/Figures/\n")
list.files("Figures")

# Panel C: Dot plot of canonical markers
markers <- c("RHO", "NR2E3",           # Rods
             "ARR3", "OPN1SW",          # Cones
             "RLBP1", "GLUL",           # Müller Glia
             "VSX2", "GRM6",            # Bipolar
             "GAD1", "SLC6A9",          # Amacrine
             "RBPMS", "THY1",           # RGC
             "ONECUT1", "LHX1")         # Horizontal

fig1c <- DotPlot(retina, features = markers, group.by = "cell_type") +
  RotatedAxis() +
  ggtitle("C. Canonical Marker Expression") +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 9))

fig1c
ggsave("Figures/Fig1C_DotPlot_Markers.pdf", fig1c, width = 10, height = 5)

# Panel D: Cell type proportions
cell_props <- retina@meta.data %>%
  group_by(condition, cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(proportion = count / sum(count) * 100)

fig1d <- ggplot(cell_props, aes(x = condition, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cell_colors) +
  labs(x = "", y = "Proportion (%)", fill = "Cell Type",
       title = "D. Cell Type Composition") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "right")

fig1d
ggsave("Figures/Fig1D_CellType_Proportions.pdf", fig1d, width = 6, height = 5)

# Combine all Figure 1 panels
fig1_combined <- (fig1a | fig1d) / (fig1b) / (fig1c)
ggsave("Figures/Figure1_Complete.pdf", fig1_combined, width = 14, height = 14)

cat("Figure 1 panels saved!\n")
list.files("Figures")

# Load DESeq2 results if not already loaded
rods_results <- read.csv("Rods_Pseudobulk_DESeq2_Results.csv")
cones_results <- read.csv("Cones_Pseudobulk_DESeq2_Results.csv")
muller_results <- read.csv("Muller_Pseudobulk_DESeq2_Results.csv")
amacrine_results <- read.csv("Amacrine_Pseudobulk_DESeq2_Results.csv")

# Function to create volcano plot
make_volcano <- function(results, title, genes_to_label = NULL, padj_cutoff = 0.05, fc_cutoff = 1) {
  
  # Add significance categories
  results$significance <- "Not Significant"
  results$significance[results$padj < padj_cutoff & results$log2FoldChange > fc_cutoff] <- "Up in KO"
  results$significance[results$padj < padj_cutoff & results$log2FoldChange < -fc_cutoff] <- "Down in KO"
  
  # Create plot
  p <- ggplot(results, aes(x = log2FoldChange, y = -log10(pvalue), color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up in KO" = "#d62728", 
                                  "Down in KO" = "#1f77b4", 
                                  "Not Significant" = "grey70")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey40") +
    labs(title = title, x = "Log2 Fold Change (KO/WT)", y = "-Log10(p-value)") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position = "bottom")
  
  # Add gene labels if specified
  if (!is.null(genes_to_label)) {
    label_data <- results[results$gene %in% genes_to_label, ]
    p <- p + ggrepel::geom_text_repel(data = label_data, 
                                      aes(label = gene),
                                      size = 3, 
                                      max.overlaps = 20,
                                      color = "black")
  }
  
  return(p)
}

# Load ggrepel for labels
library(ggrepel)

# Panel A: Rods volcano plot
top_rods <- c("CWF19L2", "GPC6", "CFAP54", "COL8A1", "ELMOD1", "SPON1")
fig2a <- make_volcano(rods_results, "A. Rods (232 DEGs)", genes_to_label = top_rods)
fig2a

# Panel B: Cones volcano plot
top_cones <- c("TENM2", "THSD7B", "SLC24A2", "PLCB1", "SAG", "CFAP54")
fig2b <- make_volcano(cones_results, "B. Cones (68 DEGs)", genes_to_label = top_cones)
fig2b

# Panel C: Müller Glia volcano plot (use p-value since few pass padj)
top_muller <- c("CTNNA2", "DPP6", "CRACR2A")
fig2c <- make_volcano(muller_results, "C. Müller Glia (1 DEG)", genes_to_label = top_muller, padj_cutoff = 0.05)
fig2c

# Panel D: Summary bar plot of DEG counts
deg_summary <- data.frame(
  Cell_Type = c("Rods", "Cones", "Amacrine", "Müller Glia"),
  DEGs = c(232, 68, 200, 1)
)
deg_summary$Cell_Type <- factor(deg_summary$Cell_Type, levels = deg_summary$Cell_Type)

fig2d <- ggplot(deg_summary, aes(x = Cell_Type, y = DEGs, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = DEGs), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Rods" = "#1f77b4", "Cones" = "#ff7f0e", 
                               "Amacrine" = "#9467bd", "Müller Glia" = "#2ca02c")) +
  labs(title = "D. DEGs per Cell Type (padj < 0.05)", x = "", y = "Number of DEGs") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "none") +
  ylim(0, 260)

fig2d

# Save individual panels
ggsave("Figures/Fig2A_Volcano_Rods.pdf", fig2a, width = 6, height = 5)
ggsave("Figures/Fig2B_Volcano_Cones.pdf", fig2b, width = 6, height = 5)
ggsave("Figures/Fig2C_Volcano_Muller.pdf", fig2c, width = 6, height = 5)
ggsave("Figures/Fig2D_DEG_Summary.pdf", fig2d, width = 5, height = 5)

# Combine Figure 2
fig2_combined <- (fig2a | fig2b) / (fig2c | fig2d)
ggsave("Figures/Figure2_Complete.pdf", fig2_combined, width = 12, height = 10)
ggsave("Figures/Figure2_Complete.png", fig2_combined, width = 12, height = 10, dpi = 300)

cat("Figure 2 saved!\n")
list.files("Figures")


# Load GO results
load("GO_enrichment_results.RData")

# Function to create GO dot plot
make_go_dotplot <- function(go_result, title, n_terms = 10) {
  
  if(is.null(go_result) || nrow(as.data.frame(go_result)) == 0) {
    return(NULL)
  }
  
  df <- as.data.frame(go_result)[1:min(n_terms, nrow(as.data.frame(go_result))), ]
  
  # Parse GeneRatio
  df$GeneRatio_numeric <- sapply(df$GeneRatio, function(x) {
    parts <- as.numeric(strsplit(x, "/")[[1]])
    return(parts[1]/parts[2])
  })
  
  # Truncate long descriptions
  df$Description <- ifelse(nchar(df$Description) > 40, 
                           paste0(substr(df$Description, 1, 37), "..."),
                           df$Description)
  
  # Reorder by p.adjust
  df$Description <- factor(df$Description, levels = rev(df$Description))
  
  p <- ggplot(df, aes(x = GeneRatio_numeric, y = Description, 
                      size = Count, color = -log10(p.adjust))) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red", name = "-Log10(padj)") +
    scale_size_continuous(range = c(3, 8), name = "Gene Count") +
    labs(title = title, x = "Gene Ratio", y = "") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 11),
          axis.text.y = element_text(size = 9))
  
  return(p)
}

# Panel A: Rods Upregulated (ECM/fibrosis)
fig3a <- make_go_dotplot(go_rods_up, "A. Rods - Upregulated in KO", n_terms = 10)
fig3a

# Panel B: Cones Downregulated (synaptic/phototransduction)
fig3b <- make_go_dotplot(go_cones_down, "B. Cones - Downregulated in KO", n_terms = 10)
fig3b

# Panel C: Amacrine Upregulated (ER stress/synaptic vesicle)
fig3c <- make_go_dotplot(go_amacrine_up, "C. Amacrine - Upregulated in KO", n_terms = 10)
fig3c

# Panel D: Müller Glia Downregulated (glutamate signaling)
fig3d <- make_go_dotplot(go_muller_down, "D. Müller Glia - Downregulated in KO", n_terms = 10)
fig3d

# Save individual panels
ggsave("Figures/Fig3A_GO_Rods_Up.pdf", fig3a, width = 7, height = 5)
ggsave("Figures/Fig3B_GO_Cones_Down.pdf", fig3b, width = 7, height = 5)
ggsave("Figures/Fig3C_GO_Amacrine_Up.pdf", fig3c, width = 7, height = 5)
ggsave("Figures/Fig3D_GO_Muller_Down.pdf", fig3d, width = 7, height = 5)

# Combine Figure 3
fig3_combined <- (fig3a | fig3b) / (fig3c | fig3d)
ggsave("Figures/Figure3_Complete.pdf", fig3_combined, width = 14, height = 10)
ggsave("Figures/Figure3_Complete.png", fig3_combined, width = 14, height = 10, dpi = 300)

cat("Figure 3 saved!\n")

# ============================================
# FIGURE 4: Key Gene Expression Patterns
# ============================================

# Panel A: Feature plot - CWF19L2 (top Rod DEG)
fig4a <- FeaturePlot(retina, features = "CWF19L2", split.by = "condition", 
                     pt.size = 0.5, order = TRUE) &
  scale_color_gradientn(colors = c("grey90", "blue", "darkblue")) &
  theme(plot.title = element_text(face = "bold", size = 12))

fig4a

# Panel B: Feature plot - TENM2 (top Cone DEG)
fig4b <- FeaturePlot(retina, features = "TENM2", split.by = "condition",
                     pt.size = 0.5, order = TRUE) &
  scale_color_gradientn(colors = c("grey90", "orange", "red")) &
  theme(plot.title = element_text(face = "bold", size = 12))

fig4b

# Panel C: Feature plot - CTNNA2 (top Müller DEG)
fig4c <- FeaturePlot(retina, features = "CTNNA2", split.by = "condition",
                     pt.size = 0.5, order = TRUE) &
  scale_color_gradientn(colors = c("grey90", "green", "darkgreen")) &
  theme(plot.title = element_text(face = "bold", size = 12))

fig4c

# Panel D: Violin plots for top DEGs by cell type
# Subset to relevant cell types for cleaner visualization
key_genes <- c("CWF19L2", "TENM2", "CTNNA2")

fig4d <- VlnPlot(retina, features = key_genes, group.by = "cell_type", 
                 split.by = "condition", pt.size = 0, 
                 cols = c("WT" = "#4DAF4A", "CLRN1_KO" = "#E41A1C")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig4d

# Save individual panels
ggsave("Figures/Fig4A_FeaturePlot_CWF19L2.pdf", fig4a, width = 8, height = 4)
ggsave("Figures/Fig4B_FeaturePlot_TENM2.pdf", fig4b, width = 8, height = 4)
ggsave("Figures/Fig4C_FeaturePlot_CTNNA2.pdf", fig4c, width = 8, height = 4)
ggsave("Figures/Fig4D_ViolinPlot_KeyGenes.pdf", fig4d, width = 12, height = 8)

# Alternative: Individual violin plots for each gene with better clarity
fig4d_cwf <- VlnPlot(retina, features = "CWF19L2", group.by = "cell_type", 
                     split.by = "condition", pt.size = 0,
                     cols = c("WT" = "#4DAF4A", "CLRN1_KO" = "#E41A1C")) +
  ggtitle("CWF19L2 Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

fig4d_tenm2 <- VlnPlot(retina, features = "TENM2", group.by = "cell_type", 
                       split.by = "condition", pt.size = 0,
                       cols = c("WT" = "#4DAF4A", "CLRN1_KO" = "#E41A1C")) +
  ggtitle("TENM2 Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

fig4d_ctnna2 <- VlnPlot(retina, features = "CTNNA2", group.by = "cell_type", 
                        split.by = "condition", pt.size = 0,
                        cols = c("WT" = "#4DAF4A", "CLRN1_KO" = "#E41A1C")) +
  ggtitle("CTNNA2 Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

# Save individual violin plots
ggsave("Figures/Fig4D1_Violin_CWF19L2.pdf", fig4d_cwf, width = 8, height = 4)
ggsave("Figures/Fig4D2_Violin_TENM2.pdf", fig4d_tenm2, width = 8, height = 4)
ggsave("Figures/Fig4D3_Violin_CTNNA2.pdf", fig4d_ctnna2, width = 8, height = 4)

# Combine Figure 4
fig4_combined <- (fig4a) / (fig4b) / (fig4c)
ggsave("Figures/Figure4_FeaturePlots.pdf", fig4_combined, width = 8, height = 12)
ggsave("Figures/Figure4_FeaturePlots.png", fig4_combined, width = 8, height = 12, dpi = 300)

# Alternative combined with violins
fig4_violin_combined <- fig4d_cwf / fig4d_tenm2 / fig4d_ctnna2
ggsave("Figures/Figure4_ViolinPlots.pdf", fig4_violin_combined, width = 10, height = 10)
ggsave("Figures/Figure4_ViolinPlots.png", fig4_violin_combined, width = 10, height = 10, dpi = 300)

cat("Figure 4 saved!\n")
list.files("Figures")
