# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

# Set working directory
setwd("/home/doyang/turbo/CLRN1 WT VS KO 10M SnRNAseq/")

# Load Seurat object
retina <- readRDS("CLRN1_Retina_with_DonorIDs.rds")

# Define HSP90 and chaperone genes
hsp90_genes <- c("HSP90AB1", "HSP90AA1", "HSP90B1")
chaperone_genes <- c("HSP90AB1", "HSP90AA1", "HSP90B1", "HSPA5", "HSPD1", "CLU", "CANX", "PDIA3")

# ============================================
# Panel C: Heatmap of log2FC and significance across all cell types
# ============================================

# Build data from DESeq2 results
amacrine_results <- read.csv("Amacrine_Pseudobulk_DESeq2_Results.csv")
rgc_res <- read.csv("RGC_Pseudobulk_DESeq2_Results.csv")
horiz_res <- read.csv("Horizontal_Pseudobulk_DESeq2_Results.csv")
bipolar_res <- read.csv("Bipolar_Pseudobulk_DESeq2_Results.csv")
rods_results <- read.csv("Rods_Pseudobulk_DESeq2_Results.csv")
cones_results <- read.csv("Cones_Pseudobulk_DESeq2_Results.csv")
muller_results <- read.csv("Muller_Pseudobulk_DESeq2_Results.csv")

extract_chaperone <- function(res, cell_type, genes) {
  df <- res[res$gene %in% genes, c("gene", "log2FoldChange", "padj")]
  df$cell_type <- cell_type
  return(df)
}

heatmap_data <- rbind(
  extract_chaperone(amacrine_results, "Amacrine", chaperone_genes),
  extract_chaperone(rgc_res, "RGC", chaperone_genes),
  extract_chaperone(horiz_res, "Horizontal", chaperone_genes),
  extract_chaperone(bipolar_res, "Bipolar", chaperone_genes),
  extract_chaperone(rods_results, "Rods", chaperone_genes),
  extract_chaperone(cones_results, "Cones", chaperone_genes),
  extract_chaperone(muller_results, "Müller Glia", chaperone_genes)
)

# Add significance labels
heatmap_data$sig <- ""
heatmap_data$sig[heatmap_data$padj < 0.05 & !is.na(heatmap_data$padj)] <- "*"
heatmap_data$sig[heatmap_data$padj < 0.01 & !is.na(heatmap_data$padj)] <- "**"
heatmap_data$sig[heatmap_data$padj < 0.001 & !is.na(heatmap_data$padj)] <- "***"

# Order cell types: inner retina first, then photoreceptors, then glia
heatmap_data$cell_type <- factor(heatmap_data$cell_type, 
                                 levels = c("Amacrine", "RGC", "Horizontal", "Bipolar",
                                            "Rods", "Cones", "Müller Glia"))

heatmap_data$gene <- factor(heatmap_data$gene, 
                            levels = c("HSP90AB1", "HSP90AA1", "HSP90B1", 
                                       "HSPA5", "HSPD1", "CLU", "CANX", "PDIA3"))

fig_heatmap <- ggplot(heatmap_data, aes(x = gene, y = cell_type, fill = log2FoldChange)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig), size = 5, vjust = 0.5) +
  geom_text(aes(label = round(log2FoldChange, 2)), size = 2.5, vjust = 1.8, color = "grey30") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "Log2FC\n(KO/WT)", limits = c(-0.5, 2)) +
  labs(title = "Chaperone / UPR Gene Expression Changes",
       subtitle = "Pseudobulk DESeq2  |  * p.adj<0.05  ** p.adj<0.01  *** p.adj<0.001",
       x = "", y = "") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank()) +
  # Add a divider between inner retina and photoreceptors
  annotate("segment", x = 0.5, xend = 8.5, y = 3.5, yend = 3.5, 
           linetype = "dashed", color = "grey50") +
  annotate("text", x = 8.8, y = 2, label = "Inner\nRetina", size = 2.5, hjust = 0, color = "grey40") +
  annotate("text", x = 8.8, y = 5.5, label = "Outer Retina\n& Glia", size = 2.5, hjust = 0, color = "grey40")

fig_heatmap
ggsave("Figures/Fig_Chaperone_Heatmap.pdf", fig_heatmap, width = 8, height = 5)
ggsave("Figures/Fig_Chaperone_Heatmap.png", fig_heatmap, width = 8, height = 5, dpi = 300)


cat("\n✅ All HSP90/chaperone figures saved!\n")
list.files("Figures")