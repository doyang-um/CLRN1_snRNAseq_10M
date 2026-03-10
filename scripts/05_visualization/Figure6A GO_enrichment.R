# ============================================
# COMPREHENSIVE GO ENRICHMENT HEATMAP
# All cell types in one figure
# ============================================

library(ggplot2)
library(dplyr)

# Build a combined pathway summary for all cell types
# Select key representative pathways
pathway_summary <- data.frame(
  Cell_Type = c(
    # Rods Up
    rep("Rods ↑", 3),
    # Cones Down
    rep("Cones ↓", 3),
    # Amacrine Up
    rep("Amacrine ↑", 3),
    # RGC Up
    rep("RGC ↑", 3),
    # Horizontal Up
    rep("Horizontal ↑", 3),
    # Müller Down
    rep("Müller ↓", 3)
  ),
  Pathway = c(
    # Rods Up
    "ECM organization", "Connective tissue dev.", "Cartilage development",
    # Cones Down
    "Synapse organization", "Visual perception", "cAMP/PKA signaling",
    # Amacrine Up
    "Protein folding / ER stress", "Glycolytic process", "Synaptic vesicle cycle",
    # RGC Up
    "Aerobic respiration", "Presynaptic endocytosis", "Chaperone assembly",
    # Horizontal Up
    "Response to heat", "Ca²⁺ channel regulation", "Lactate metabolism",
    # Müller Down
    "Glutamate receptor signaling", "Receptor localization to synapse", "GABAergic transmission"
  ),
  neg_log10_padj = c(
    # Rods Up
    -log10(2.44e-06), -log10(3.67e-05), -log10(7.64e-06),
    # Cones Down
    -log10(0.006), -log10(0.026), -log10(0.008),
    # Amacrine Up
    -log10(5.17e-05), -log10(1.43e-04), -log10(3.07e-04),
    # RGC Up
    -log10(7.13e-05), -log10(7.13e-05), -log10(6.15e-04),
    # Horizontal Up
    -log10(0.023), -log10(0.023), -log10(0.023),
    # Müller Down
    -log10(0.031), -log10(0.031), -log10(0.034)
  ),
  Category = c(
    rep("ECM / Fibrosis", 3),
    rep("Synaptic / Visual", 3),
    rep("Stress / Metabolic", 3),
    rep("Stress / Metabolic", 3),
    rep("Stress / Metabolic", 3),
    rep("Glial Support", 3)
  )
)

# Order cell types
pathway_summary$Cell_Type <- factor(pathway_summary$Cell_Type, 
                                     levels = c("Müller ↓", "Rods ↑", "Cones ↓", 
                                                "Amacrine ↑", "RGC ↑", "Horizontal ↑"))

# Create bubble plot
fig_go_summary <- ggplot(pathway_summary, aes(x = Cell_Type, y = Pathway, 
                                               size = neg_log10_padj, color = Category)) +
  geom_point() +
  scale_size_continuous(range = c(4, 12), name = "-Log10(padj)") +
  scale_color_manual(values = c("ECM / Fibrosis" = "#d62728",
                                 "Synaptic / Visual" = "#1f77b4",
                                 "Stress / Metabolic" = "#ff7f0e",
                                 "Glial Support" = "#2ca02c")) +
  labs(title = "GO Enrichment Summary Across All Retinal Cell Types",
       subtitle = "CLRN1 KO vs WT", x = "", y = "") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 10),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey90"))

fig_go_summary
ggsave("Figures/Fig_GO_Summary_AllCellTypes.pdf", fig_go_summary, width = 10, height = 8)
ggsave("Figures/Fig_GO_Summary_AllCellTypes.png", fig_go_summary, width = 10, height = 8, dpi = 300)

# ============================================
# UPDATED COMPLETE SUPPLEMENTARY TABLE
# ============================================

complete_go_table <- data.frame(
  Cell_Type = c("Rods", "Rods", "Cones", "Amacrine", "RGC", "Horizontal", "Horizontal", "Müller Glia"),
  Direction = c("Up in KO", "Down in KO", "Down in KO", "Up in KO", "Up in KO", "Up in KO", "Down in KO", "Down in KO"),
  n_DEGs = c(119, 113, 50, 131, 27, 23, 14, 59),
  Top_GO_Term = c("ECM organization", "No significant terms", "Synapse organization",
                   "Protein folding", "Presynaptic endocytosis", "Response to heat",
                   "Smooth muscle relaxation", "Glutamate receptor signaling"),
  GO_padj = c("2.4e-06", "NA", "0.006", "5.2e-05", "7.1e-05", "0.023", "0.014", "0.031"),
  Biological_Theme = c("Fibrotic remodeling", "–", "Synaptic/visual dysfunction",
                        "ER stress & metabolic shift", "Metabolic stress & synaptic", 
                        "Stress response", "–", "Impaired glial support")
)

print(complete_go_table)
write.csv(complete_go_table, "Complete_GO_Summary_AllCellTypes.csv", row.names = FALSE)

cat("\n✅ All GO summary figures and tables saved!\n")
list.files("Figures")