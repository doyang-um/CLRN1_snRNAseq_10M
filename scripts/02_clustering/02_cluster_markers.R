# ===================================================================
# Cluster Identification Markers - Publication Figures
# ===================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)

# Load your Seurat object if needed
# retina <- readRDS("CLRN1_Retina_CellType_Annotated_Filtered.rds")

# ===================================================================
# Define Key Retinal Cell Type Markers
# ===================================================================

# Select the most established, canonical markers for each cell type
retinal_markers <- list(
  Rods = c("RHO", "NRL", "SAG", "PDE6A"),
  Cones = c("ARR3", "OPN1SW", "OPN1MW", "GNAT2"),
  Bipolar_Cells = c("VSX2", "VSX1", "CABP5", "GRM6"),
  Horizontal_Cells = c("ONECUT1", "ONECUT2", "LHX1", "PROX1"),
  Amacrine_Cells = c("GAD1", "GAD2", "SLC6A9", "PAX6"),
  RGC = c("RBPMS", "POU4F2", "SNCG", "THY1"),
  Muller_Glia = c("RLBP1", "SLC1A3", "GLUL", "VIM")
)

# Check which markers are actually in the dataset
all_genes <- rownames(retina)

available_markers <- list()
for (cell_type in names(retinal_markers)) {
  markers <- retinal_markers[[cell_type]]
  found <- markers[markers %in% all_genes]
  if (length(found) > 0) {
    available_markers[[cell_type]] <- found
  }
}

cat("=== Available Markers by Cell Type ===\n")
for (cell_type in names(available_markers)) {
  cat(cell_type, ":", paste(available_markers[[cell_type]], collapse = ", "), "\n")
}

# Create a flat list of all available markers
all_markers_vec <- unique(unlist(available_markers))

cat("\nTotal markers to plot:", length(all_markers_vec), "\n")

# ===================================================================
# 1. Comprehensive Dot Plot - All Cell Types and Markers
# ===================================================================

# Set order of cell types for better visualization
cell_type_order <- c("Rods", "Cones", "Bipolar_Cells", "Horizontal_Cells", 
                     "Amacrine_Cells", "RGC", "Muller_Glia")

# Filter to only include cell types that exist
cell_type_order <- cell_type_order[cell_type_order %in% unique(retina$cell_type)]

# Set factor levels
retina$cell_type <- factor(retina$cell_type, levels = cell_type_order)

# Create comprehensive dot plot
p_dot_all <- DotPlot(retina, 
                     features = all_markers_vec,
                     group.by = "cell_type",
                     dot.scale = 8) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Retinal Cell Type Marker Expression",
    x = "Marker Genes",
    y = "Cell Type"
  ) +
  scale_color_gradient(low = "lightgrey", high = "darkblue")

ggsave("Cluster_Markers_DotPlot_All.pdf", p_dot_all, width = 14, height = 8)

# ===================================================================
# 2. Dot Plot by Cell Type Category (Organized)
# ===================================================================

# Create a more organized dot plot with markers grouped by cell type
marker_order <- c()
for (ct in cell_type_order) {
  if (ct %in% names(available_markers)) {
    marker_order <- c(marker_order, available_markers[[ct]])
  }
}

marker_order <- unique(marker_order)

p_dot_organized <- DotPlot(retina,
                           features = marker_order,
                           group.by = "cell_type",
                           dot.scale = 8) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  labs(
    title = "Canonical Retinal Cell Type Markers",
    x = "Marker Genes",
    y = "Cell Type"
  ) +
  scale_color_gradientn(colors = c("lightgrey", "yellow", "red", "darkred"),
                        name = "Average\nExpression")

ggsave("Cluster_Markers_DotPlot_Organized.pdf", p_dot_organized, width = 16, height = 8)

# ===================================================================
# 3. Feature Plots - Key Markers for Each Cell Type
# ===================================================================

# Select top 2 markers per cell type for feature plots
feature_plot_markers <- c()
for (ct in names(available_markers)) {
  markers <- available_markers[[ct]][1:min(2, length(available_markers[[ct]]))]
  feature_plot_markers <- c(feature_plot_markers, markers)
}

feature_plot_markers <- unique(feature_plot_markers)

cat("\nCreating feature plots for", length(feature_plot_markers), "markers...\n")

# Create feature plots in batches of 4
n_batches <- ceiling(length(feature_plot_markers) / 4)

for (i in 1:n_batches) {
  start_idx <- (i - 1) * 4 + 1
  end_idx <- min(i * 4, length(feature_plot_markers))
  markers_batch <- feature_plot_markers[start_idx:end_idx]
  
  p_feat <- FeaturePlot(retina, 
                        features = markers_batch,
                        ncol = 2,
                        pt.size = 0.3) &
    theme(
      plot.title = element_text(face = "italic", size = 12),
      legend.position = "right"
    )
  
  ggsave(paste0("Cluster_Markers_FeaturePlot_Batch", i, ".pdf"), 
         p_feat, width = 12, height = 12)
}

# ===================================================================
# 4. Cell Type-Specific Feature Plot Panels
# ===================================================================

# Create individual panels for each major cell type
major_celltypes <- c("Rods", "Cones", "Muller_Glia", "RGC", "Bipolar_Cells")

for (ct in major_celltypes) {
  if (ct %in% names(available_markers) && length(available_markers[[ct]]) >= 2) {
    
    markers <- available_markers[[ct]][1:min(4, length(available_markers[[ct]]))]
    
    p_ct <- FeaturePlot(retina,
                        features = markers,
                        ncol = 2,
                        pt.size = 0.4) &
      theme(
        plot.title = element_text(face = "italic", size = 12),
        legend.position = "right"
      )
    
    ggsave(paste0("Markers_FeaturePlot_", ct, ".pdf"), 
           p_ct, width = 12, height = 10)
  }
}

# ===================================================================
# 5. Combined UMAP with Cell Type Labels + Marker Feature Plots
# ===================================================================

# UMAP with cell type labels
p_umap <- DimPlot(retina, 
                  group.by = "cell_type",
                  label = TRUE,
                  repel = TRUE,
                  pt.size = 0.5,
                  label.size = 5) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom"
  ) +
  labs(title = "Retinal Cell Types") +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))

# Select one key marker per cell type
key_markers_per_type <- c()
for (ct in names(available_markers)) {
  key_markers_per_type <- c(key_markers_per_type, available_markers[[ct]][1])
}
key_markers_per_type <- unique(key_markers_per_type)

# Create feature plots for key markers
p_features_key <- FeaturePlot(retina,
                              features = key_markers_per_type[1:min(6, length(key_markers_per_type))],
                              ncol = 3,
                              pt.size = 0.3) &
  theme(
    plot.title = element_text(face = "italic", size = 11),
    legend.position = "right"
  )

# Combine UMAP and feature plots
combined_panel <- p_umap / p_features_key + 
  plot_layout(heights = c(1, 2))

ggsave("Cluster_Identification_Combined_Panel.pdf", 
       combined_panel, width = 16, height = 14)

# ===================================================================
# 6. Heatmap of Marker Expression
# ===================================================================

library(pheatmap)

# Get average expression per cell type
Idents(retina) <- "cell_type"
avg_exp <- AverageExpression(retina, 
                             features = all_markers_vec,
                             group.by = "cell_type")

# Extract the RNA assay data
exp_matrix <- as.matrix(avg_exp$RNA)

# Scale the data
exp_matrix_scaled <- t(scale(t(exp_matrix)))

# Create annotation for markers (which cell type they represent)
marker_annotation <- data.frame(
  Marker_For = character(nrow(exp_matrix_scaled)),
  row.names = rownames(exp_matrix_scaled)
)

for (gene in rownames(exp_matrix_scaled)) {
  for (ct in names(available_markers)) {
    if (gene %in% available_markers[[ct]]) {
      marker_annotation[gene, "Marker_For"] <- ct
      break
    }
  }
}

# Create color palette
ann_colors <- list(
  Marker_For = c(
    Rods = "#E41A1C",
    Cones = "#377EB8",
    Bipolar_Cells = "#4DAF4A",
    Horizontal_Cells = "#984EA3",
    Amacrine_Cells = "#FF7F00",
    RGC = "#FFFF33",
    Muller_Glia = "#A65628"
  )
)

# Create heatmap
pdf("Cluster_Markers_Heatmap.pdf", width = 10, height = 12)
pheatmap(exp_matrix_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_row = marker_annotation,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Retinal Cell Type Marker Expression",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 9,
         fontsize_col = 10,
         angle_col = 45)
dev.off()

# ===================================================================
# 7. Violin Plots for Key Markers
# ===================================================================

# Select top marker for each cell type
top_markers <- c()
for (ct in names(available_markers)) {
  top_markers <- c(top_markers, available_markers[[ct]][1])
}
top_markers <- unique(top_markers)

# Create violin plots
p_violin <- VlnPlot(retina,
                    features = top_markers[1:min(8, length(top_markers))],
                    group.by = "cell_type",
                    pt.size = 0,
                    ncol = 4) &
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(face = "italic", size = 11)
  )

ggsave("Cluster_Markers_ViolinPlot.pdf", p_violin, width = 16, height = 8)

# ===================================================================
# 8. Ridge Plot (Distribution of Expression)
# ===================================================================

library(ggridges)

# Select a few key markers for ridge plot
ridge_markers <- top_markers[1:min(4, length(top_markers))]

for (marker in ridge_markers) {
  
  plot_data <- FetchData(retina, vars = c(marker, "cell_type"))
  colnames(plot_data)[1] <- "expression"
  
  p_ridge <- ggplot(plot_data, aes(x = expression, y = cell_type, fill = cell_type)) +
    geom_density_ridges(alpha = 0.7) +
    theme_ridges() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14)
    ) +
    labs(
      title = paste(marker, "Expression Distribution"),
      x = "Expression Level",
      y = "Cell Type"
    )
  
  ggsave(paste0("Marker_RidgePlot_", marker, ".pdf"), 
         p_ridge, width = 10, height = 6)
}

# ===================================================================
# 9. Summary Table of Markers
# ===================================================================

# Create a summary table
marker_summary <- data.frame(
  Cell_Type = character(),
  Markers_Used = character(),
  N_Cells = integer(),
  stringsAsFactors = FALSE
)

for (ct in names(available_markers)) {
  n_cells <- sum(retina$cell_type == ct)
  markers_str <- paste(available_markers[[ct]], collapse = ", ")
  
  marker_summary <- rbind(marker_summary, data.frame(
    Cell_Type = ct,
    Markers_Used = markers_str,
    N_Cells = n_cells
  ))
}

write.csv(marker_summary, "Cluster_Identification_Markers_Summary.csv", row.names = FALSE)

print("\n=== Cluster Identification Markers Summary ===")
print(marker_summary)

# ===================================================================
# Summary Report
# ===================================================================

cat("\n===============================================\n")
cat("CLUSTER IDENTIFICATION FIGURES GENERATED\n")
cat("===============================================\n\n")

cat("Dot Plots:\n")
cat("  - Cluster_Markers_DotPlot_All.pdf (all markers)\n")
cat("  - Cluster_Markers_DotPlot_Organized.pdf (organized by cell type)\n\n")

cat("Feature Plots:\n")
cat("  - Cluster_Markers_FeaturePlot_Batch[1-N].pdf (all markers)\n")
cat("  - Markers_FeaturePlot_[CellType].pdf (per cell type)\n\n")

cat("Combined Panels:\n")
cat("  - Cluster_Identification_Combined_Panel.pdf (UMAP + markers)\n\n")

cat("Additional Visualizations:\n")
cat("  - Cluster_Markers_Heatmap.pdf (expression heatmap)\n")
cat("  - Cluster_Markers_ViolinPlot.pdf (distribution)\n")
cat("  - Marker_RidgePlot_[Gene].pdf (density plots)\n\n")

cat("Summary Data:\n")
cat("  - Cluster_Identification_Markers_Summary.csv\n\n")

cat("Analysis complete!\n")
cat("Recommended for publication: Cluster_Markers_DotPlot_Organized.pdf\n")