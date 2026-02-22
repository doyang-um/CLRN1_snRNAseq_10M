# ===================================================================
# Retinal Cell Type Identification Using Known Markers
# ===================================================================

library(Seurat)
library(tidyverse)
library(patchwork)

# Load your Seurat object (if not already loaded)
# retina <- readRDS("CLRN1_Retina_Seurat_Object.rds")

# ===================================================================
# Define Known Retinal Cell Type Markers
# ===================================================================

# Comprehensive list of retinal cell type markers
# Note: Gene names might vary slightly in rabbit - we'll check what's available

retinal_markers <- list(
  Rods = c("RHO", "NRL", "NR2E3", "SAG", "PDE6A", "PDE6B", "GNAT1", "CNGA1", "ROM1"),
  Cones = c("ARR3", "OPN1SW", "OPN1MW", "OPN1LW", "GNAT2", "PDE6C", "PDE6H", "GNGT2"),
  Photoreceptors_General = c("RCVRN", "RECOVERIN", "CRX", "OTX2"),
  Bipolar_Cells = c("VSX1", "VSX2", "CABP5", "TRPM1", "GRM6", "PRKCA", "SCGN", "ISL1"),
  Horizontal_Cells = c("ONECUT1", "ONECUT2", "LHX1", "CALB1", "PROX1"),
  Amacrine_Cells = c("GAD1", "GAD2", "SLC6A9", "PAX6", "TFAP2A", "TFAP2B", "CHAT"),
  Retinal_Ganglion_Cells = c("RBPMS", "RBPMS2", "POU4F1", "POU4F2", "SNCG", "THY1", "NEFL", "NEFM"),
  Muller_Glia = c("RLBP1", "SLC1A3", "GLUL", "VIM", "APOE", "CLU", "CRALBP"),
  Microglia = c("C1QA", "C1QB", "C1QC", "CX3CR1", "TMEM119", "P2RY12", "AIF1", "ITGAM", "CD68"),
  Astrocytes = c("GFAP", "AQP4", "S100B", "ALDH1L1"),
  RPE = c("RPE65", "BEST1", "MERTK", "RDH5"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "FLT1"),
  Pericytes = c("PDGFRB", "ACTA2", "RGS5", "CSPG4")
)

# ===================================================================
# Check Which Markers Are Present in Your Data
# ===================================================================

all_genes <- rownames(retina)

# Function to find available markers
check_markers <- function(marker_list, available_genes) {
  lapply(marker_list, function(markers) {
    found <- markers[markers %in% available_genes]
    missing <- markers[!markers %in% available_genes]
    list(found = found, missing = missing)
  })
}

marker_check <- check_markers(retinal_markers, all_genes)

# Print what markers we found
cat("\n=== Marker Availability Check ===\n")
for (cell_type in names(marker_check)) {
  cat("\n", cell_type, ":\n", sep = "")
  cat("  Found (", length(marker_check[[cell_type]]$found), "): ", 
      paste(marker_check[[cell_type]]$found, collapse = ", "), "\n", sep = "")
  if (length(marker_check[[cell_type]]$missing) > 0) {
    cat("  Missing (", length(marker_check[[cell_type]]$missing), "): ", 
        paste(marker_check[[cell_type]]$missing, collapse = ", "), "\n", sep = "")
  }
}

# ===================================================================
# Create a Cleaned Marker List (Only Available Genes)
# ===================================================================

available_markers <- list(
  Rods = marker_check$Rods$found,
  Cones = marker_check$Cones$found,
  Photoreceptors_General = marker_check$Photoreceptors_General$found,
  Bipolar_Cells = marker_check$Bipolar_Cells$found,
  Horizontal_Cells = marker_check$Horizontal_Cells$found,
  Amacrine_Cells = marker_check$Amacrine_Cells$found,
  Retinal_Ganglion_Cells = marker_check$Retinal_Ganglion_Cells$found,
  Muller_Glia = marker_check$Muller_Glia$found,
  Microglia = marker_check$Microglia$found
)

# Remove empty categories
available_markers <- available_markers[sapply(available_markers, length) > 0]

# ===================================================================
# Visualize Markers Across Clusters
# ===================================================================

# Dot plot for all available markers
all_markers_vec <- unique(unlist(available_markers))

if (length(all_markers_vec) > 0) {
  p_dot_all <- DotPlot(retina, features = all_markers_vec) + 
    RotatedAxis() +
    ggtitle("Expression of Retinal Cell Type Markers Across Clusters")
  
  ggsave("Cell_Type_Markers_DotPlot_All.pdf", p_dot_all, width = 16, height = 8)
  
  # Create separate dot plots for each major cell type
  for (cell_type in names(available_markers)) {
    if (length(available_markers[[cell_type]]) > 0) {
      p <- DotPlot(retina, features = available_markers[[cell_type]]) + 
        RotatedAxis() +
        ggtitle(paste(cell_type, "Markers"))
      
      filename <- paste0("Markers_", gsub(" ", "_", cell_type), ".pdf")
      ggsave(filename, p, width = 10, height = 6)
    }
  }
}

# ===================================================================
# Feature Plots for Key Markers
# ===================================================================

# Select top markers for each cell type (if available)
key_markers <- c()
for (cell_type in names(available_markers)) {
  if (length(available_markers[[cell_type]]) > 0) {
    key_markers <- c(key_markers, available_markers[[cell_type]][1:min(2, length(available_markers[[cell_type]]))])
  }
}

if (length(key_markers) > 0) {
  # Create feature plots in batches of 4
  n_batches <- ceiling(length(key_markers) / 4)
  
  for (i in 1:n_batches) {
    start_idx <- (i - 1) * 4 + 1
    end_idx <- min(i * 4, length(key_markers))
    markers_batch <- key_markers[start_idx:end_idx]
    
    p_features <- FeaturePlot(retina, features = markers_batch, ncol = 2)
    ggsave(paste0("FeaturePlot_Markers_Batch", i, ".pdf"), p_features, width = 12, height = 12)
  }
}

# ===================================================================
# Calculate Module Scores for Each Cell Type
# ===================================================================

# Add module scores for each cell type
for (cell_type in names(available_markers)) {
  if (length(available_markers[[cell_type]]) > 0) {
    module_name <- paste0(gsub(" ", "_", cell_type), "_Score")
    retina <- AddModuleScore(
      object = retina,
      features = list(available_markers[[cell_type]]),
      name = module_name
    )
  }
}

# Get the actual column names (Seurat adds "1" to the end)
score_columns <- grep("_Score1$", colnames(retina@meta.data), value = TRUE)

# Feature plots of module scores
if (length(score_columns) > 0) {
  n_batches <- ceiling(length(score_columns) / 4)
  
  for (i in 1:n_batches) {
    start_idx <- (i - 1) * 4 + 1
    end_idx <- min(i * 4, length(score_columns))
    scores_batch <- score_columns[start_idx:end_idx]
    
    p_scores <- FeaturePlot(retina, features = scores_batch, ncol = 2)
    ggsave(paste0("ModuleScores_Batch", i, ".pdf"), p_scores, width = 12, height = 12)
  }
}

# ===================================================================
# Calculate Average Expression and Scores Per Cluster
# ===================================================================

# Get average expression of markers per cluster
avg_exp <- AverageExpression(retina, features = all_markers_vec, group.by = "seurat_clusters")

# Create a heatmap of average expression
library(pheatmap)

if (length(all_markers_vec) > 0) {
  exp_matrix <- as.matrix(avg_exp$RNA)
  
  # Scale the expression values
  exp_matrix_scaled <- t(scale(t(exp_matrix)))
  
  pdf("Cell_Type_Markers_Heatmap.pdf", width = 12, height = 10)
  pheatmap(exp_matrix_scaled,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           main = "Average Expression of Cell Type Markers per Cluster",
           color = colorRampPalette(c("blue", "white", "red"))(100))
  dev.off()
}

# Calculate average module scores per cluster
if (length(score_columns) > 0) {
  cluster_scores <- retina@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(across(all_of(score_columns), mean, .names = "{col}"))
  
  # Save the scores
  write.csv(cluster_scores, "Cluster_CellType_Scores.csv", row.names = FALSE)
  
  # Create a heatmap of module scores
  score_matrix <- as.matrix(cluster_scores[, -1])
  rownames(score_matrix) <- cluster_scores$seurat_clusters
  colnames(score_matrix) <- gsub("_Score1", "", colnames(score_matrix))
  
  pdf("Cell_Type_Scores_Heatmap.pdf", width = 10, height = 8)
  pheatmap(t(score_matrix),
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           main = "Cell Type Module Scores per Cluster",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           display_numbers = TRUE,
           number_format = "%.2f")
  dev.off()
  
  # Print the top scoring cell type for each cluster
  cat("\n=== Predicted Cell Types per Cluster ===\n")
  for (i in 1:nrow(score_matrix)) {
    cluster_num <- cluster_scores$seurat_clusters[i]
    scores <- score_matrix[i, ]
    top_type <- names(scores)[which.max(scores)]
    top_score <- max(scores)
    
    cat(sprintf("Cluster %s: %s (score: %.3f)\n", 
                cluster_num, 
                gsub("_", " ", top_type), 
                top_score))
  }
}

# ===================================================================
# Compare Your Cluster Markers with Known Cell Types
# ===================================================================

# Load your cluster markers
if (file.exists("Top5_Markers_per_Cluster.csv")) {
  top_markers <- read.csv("Top5_Markers_per_Cluster.csv")
  
  cat("\n=== Top Markers Per Cluster ===\n")
  for (cluster_id in unique(top_markers$cluster)) {
    cluster_markers <- top_markers %>% filter(cluster == cluster_id)
    cat("\nCluster", cluster_id, ":\n")
    cat("  Top genes:", paste(cluster_markers$gene[1:min(5, nrow(cluster_markers))], collapse = ", "), "\n")
  }
}

# ===================================================================
# Save Updated Seurat Object with Module Scores
# ===================================================================

saveRDS(retina, "CLRN1_Retina_Annotated.rds")

cat("\n=== Analysis Complete ===\n")
cat("Check the following files:\n")
cat("  - Cell_Type_Markers_DotPlot_All.pdf\n")
cat("  - Cell_Type_Scores_Heatmap.pdf\n")
cat("  - Cluster_CellType_Scores.csv\n")
cat("  - Individual marker plots for each cell type\n")