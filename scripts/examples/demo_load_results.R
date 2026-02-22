#!/usr/bin/env Rscript
#------------------------------------------------------------------------------
# Demo: Explore Pre-computed DESeq2 Results
# 
# This script demonstrates how to load and analyze the pre-computed results.
# No large data files or package installation needed!
#------------------------------------------------------------------------------

# Simple base R - no packages needed!
cat("\n====================================\n")
cat("CLRN1 snRNA-seq Results Explorer\n")
cat("====================================\n\n")

# 1. Müller Glia Results
cat("1. Müller Glia (Primary CLRN1 Expression Site)\n")
cat("   Loading results...\n")
muller <- read.csv("results/DESeq2/Muller_Pseudobulk_DESeq2_Results.csv")

sig_muller <- muller[!is.na(muller$padj) & muller$padj < 0.05, ]
cat("   Total genes tested:", nrow(muller), "\n")
cat("   Significant DEGs (padj<0.05):", nrow(sig_muller), "\n\n")

if (nrow(sig_muller) > 0) {
  cat("   Significant genes:\n")
  for (i in 1:nrow(sig_muller)) {
    cat("   -", sig_muller$gene[i], 
        "| log2FC:", round(sig_muller$log2FoldChange[i], 2),

  "| log2FC:", round(sig_muller$log2FoldChange[i], 2),
        "| padj:", format(sig_muller$padj[i], scientific=TRUE, digits=3), "\n")
  }
}
# 2. Rod Photoreceptors
cat("\n2. Rod Photoreceptors (Secondary Degeneration)\n")
rods <- read.csv("results/DESeq2/Rods_Pseudobulk_DESeq2_Results.csv")
sig_rods <- rods[!is.na(rods$padj) & rods$padj < 0.05, ]

cat("   Total genes tested:", nrow(rods), "\n")
cat("   Significant DEGs:", nrow(sig_rods), "\n")
cat("   - Upregulated:", sum(sig_rods$log2FoldChange > 0), "\n")
cat("   - Downregulated:", sum(sig_rods$log2FoldChange < 0), "\n")

# Top 5 downregulated
down_rods <- sig_rods[sig_rods$log2FoldChange < 0, ]
down_rods <- down_rods[order(down_rods$log2FoldChange), ]
cat("\n   Top 5 most downregulated:\n")
for (i in 1:min(5, nrow(down_rods))) {
  cat("   -", down_rods$gene[i], 
      "| log2FC:", round(down_rods$log2FoldChange[i], 2),
      "| padj:", format(down_rods$padj[i], scientific=TRUE, digits=3), "\n")
}

# 3. Cone Photoreceptors
cat("\n3. Cone Photoreceptors\n")
cones <- read.csv("results/DESeq2/Cones_Pseudobulk_DESeq2_Results.csv")
sig_cones <- cones[!is.na(cones$padj) & cones$padj < 0.05, ]

cat("   Significant DEGs:", nrow(sig_cones), "\n")

# Check for TENM2
tenm2 <- cones[cones$gene == "TENM2", ]
if (nrow(tenm2) > 0) {
  cat("\n   Key finding - TENM2 (synaptic organizer):\n")
  cat("   - log2FC:", round(tenm2$log2FoldChange, 2), "(most severe change observed)\n")
  cat("   - padj:", format(tenm2$padj, scientific=TRUE, digits=3), "\n")
}

# 4. Amacrine Cells (Compensation)
cat("\n4. Amacrine Cells (Inner Retinal Compensation)\n")
amacrine <- read.csv("results/DESeq2/Amacrine_Pseudobulk_DESeq2_Results.csv")
sig_amacrine <- amacrine[!is.na(amacrine$padj) & amacrine$padj < 0.05, ]

cat("   Significant DEGs:", nrow(sig_amacrine), "\n")
cat("   - Upregulated:", sum(sig_amacrine$log2FoldChange > 0), "\n")
cat("   - Downregulated:", sum(sig_amacrine$log2FoldChange < 0), "\n")

# 5. Summary Table
cat("\n5. Summary Across All Cell Types\n")
cat("   =====================================\n")
cat("   Cell Type       | DEGs (padj<0.05)\n")
cat("   =====================================\n")

cell_types <- c("Muller", "Rods", "Cones", "Amacrine", "Horizontal", "Bipolar", "RGC")
for (ct in cell_types) {
  file <- paste0("results/DESeq2/", ct, "_Pseudobulk_DESeq2_Results.csv")
  if (file.exists(file)) {
    data <- read.csv(file)
    n_sig <- sum(!is.na(data$padj) & data$padj < 0.05)
    cat(sprintf("   %-15s | %d\n", ct, n_sig))
  }
}
cat("   =====================================\n")

cat("\n====================================\n")
cat("✓ Demo completed successfully!\n")
cat("====================================\n\n")

cat("Next steps:\n")
cat("- Explore individual CSV files in results/DESeq2/\n")
cat("- See full analysis pipeline in scripts/\n")
cat("- Download raw data from GEO: GSE[XXXXX]\n\n")
