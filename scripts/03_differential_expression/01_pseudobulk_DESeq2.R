# Load DESeq2
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
library(DESeq2)

# Function to run pseudobulk DESeq2 for a cell type
run_pseudobulk_deseq2 <- function(seurat_obj, cell_type_name) {
  
  # Subset to cell type
  cells <- subset(seurat_obj, cell_type == cell_type_name)
  
  # Check sample sizes
  sample_counts <- table(cells$sample_id)
  print(paste("Running DESeq2 for:", cell_type_name))
  print(sample_counts)
  
  # Aggregate counts by sample
  counts_matrix <- AggregateExpression(cells, 
                                       group.by = "sample_id",
                                       assays = "RNA",
                                       slot = "counts")$RNA
  
  # Create sample metadata
  sample_info <- data.frame(
    sample_id = colnames(counts_matrix),
    condition = ifelse(grepl("WT", colnames(counts_matrix)), "WT", "KO")
  )
  rownames(sample_info) <- sample_info$sample_id
  sample_info$condition <- factor(sample_info$condition, levels = c("WT", "KO"))
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                colData = sample_info,
                                design = ~ condition)
  
  # Filter low counts
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get results
  res <- results(dds, contrast = c("condition", "KO", "WT"))
  res <- as.data.frame(res)
  res$gene <- rownames(res)
  res <- res[order(res$padj), ]
  
  return(list(dds = dds, results = res))
}

# Run for Rods (largest cell population)
rods_deseq <- run_pseudobulk_deseq2(retina_clean, "Rods")
head(rods_deseq$results, 20)

# Run for all major cell types
muller_deseq <- run_pseudobulk_deseq2(retina_clean, "Muller_Glia")
cones_deseq <- run_pseudobulk_deseq2(retina_clean, "Cones")
amacrine_deseq <- run_pseudobulk_deseq2(retina_clean, "Amacrine_Cells")

# Summary of significant DEGs (padj < 0.05)
cat("\n=== Summary of DEGs (padj < 0.05) ===\n")
cat("Rods:", sum(rods_deseq$results$padj < 0.05, na.rm = TRUE), "DEGs\n")
cat("Müller Glia:", sum(muller_deseq$results$padj < 0.05, na.rm = TRUE), "DEGs\n")
cat("Cones:", sum(cones_deseq$results$padj < 0.05, na.rm = TRUE), "DEGs\n")
cat("Amacrine:", sum(amacrine_deseq$results$padj < 0.05, na.rm = TRUE), "DEGs\n")

# View top Müller glia DEGs
cat("\n=== Top Müller Glia DEGs ===\n")
head(muller_deseq$results, 10)

# View top Cones DEGs
cat("\n=== Top Cones DEGs ===\n")
head(cones_deseq$results, 10)

# Save all DESeq2 results
write.csv(rods_deseq$results, "Rods_Pseudobulk_DESeq2_Results.csv", row.names = FALSE)
write.csv(muller_deseq$results, "Muller_Pseudobulk_DESeq2_Results.csv", row.names = FALSE)
write.csv(cones_deseq$results, "Cones_Pseudobulk_DESeq2_Results.csv", row.names = FALSE)
write.csv(amacrine_deseq$results, "Amacrine_Pseudobulk_DESeq2_Results.csv", row.names = FALSE)

# Compare: Check if original Müller DEGs are in pseudobulk results
original_muller_degs <- c("CRACR2A", "CTNNA2", "DLG4", "CYP2J2", "RSU1", "CFL2", "PIKFYVE", "SERINC1", "PSMD3")
cat("\n=== Original Müller DEGs in Pseudobulk Results ===\n")
muller_deseq$results[muller_deseq$results$gene %in% original_muller_degs, c("gene", "log2FoldChange", "pvalue", "padj")]

# Compare: Check TENM2 in cones
cat("\n=== TENM2 in Cones (your gene of interest) ===\n")
cones_deseq$results[cones_deseq$results$gene == "TENM2", ]

# Compare: Check CWF19L2 in Rods
cat("\n=== CWF19L2 in Rods ===\n")
rods_deseq$results[rods_deseq$results$gene == "CWF19L2", ]

cat("\nAll results saved!\n")


# Save all DESeq2 results
# Compare: Check if original Müller DEGs are in pseudobulk results
original_muller_degs <- c("CRACR2A", "CTNNA2", "DLG4", "CYP2J2", "RSU1", "CFL2", "PIKFYVE", "SERINC1", "PSMD3")
cat("\n=== Original Müller DEGs in Pseudobulk Results ===\n")
muller_deseq$results[muller_deseq$results$gene %in% original_muller_degs, c("gene", "log2FoldChange", "pvalue", "padj")]

# Compare: Check TENM2 in cones
cat("\n=== TENM2 in Cones ===\n")
cones_deseq$results[cones_deseq$results$gene == "TENM2", ]

# Compare: Check CWF19L2 in Rods
cat("\n=== CWF19L2 in Rods ===\n")
rods_deseq$results[rods_deseq$results$gene == "CWF19L2", ]

# Create final summary table
final_summary <- data.frame(
  Cell_Type = c("Rods", "Cones", "Muller_Glia", "Amacrine"),
  Total_DEGs_padj05 = c(
    sum(rods_deseq$results$padj < 0.05, na.rm = TRUE),
    sum(cones_deseq$results$padj < 0.05, na.rm = TRUE),
    sum(muller_deseq$results$padj < 0.05, na.rm = TRUE),
    sum(amacrine_deseq$results$padj < 0.05, na.rm = TRUE)
  ),
  Top_Gene = c("CWF19L2", "TENM2", "CTNNA2", 
               amacrine_deseq$results$gene[1]),
  Top_Gene_log2FC = c(
    rods_deseq$results[rods_deseq$results$gene == "CWF19L2", "log2FoldChange"],
    cones_deseq$results[cones_deseq$results$gene == "TENM2", "log2FoldChange"],
    muller_deseq$results[muller_deseq$results$gene == "CTNNA2", "log2FoldChange"],
    amacrine_deseq$results$log2FoldChange[1]
  ),
  Top_Gene_padj = c(
    rods_deseq$results[rods_deseq$results$gene == "CWF19L2", "padj"],
    cones_deseq$results[cones_deseq$results$gene == "TENM2", "padj"],
    muller_deseq$results[muller_deseq$results$gene == "CTNNA2", "padj"],
    amacrine_deseq$results$padj[1]
  )
)

print(final_summary)

# Save
write.csv(final_summary, "Pseudobulk_DESeq2_Summary.csv", row.names = FALSE)

# Save the updated Seurat object
saveRDS(retina_clean, "CLRN1_Retina_Demultiplexed_Clean.rds")

cat("\n✅ Analysis complete! Files saved:\n")
cat("- Rods_Pseudobulk_DESeq2_Results.csv\n")
cat("- Cones_Pseudobulk_DESeq2_Results.csv\n")
cat("- Muller_Pseudobulk_DESeq2_Results.csv\n")
cat("- Amacrine_Pseudobulk_DESeq2_Results.csv\n")
cat("- Pseudobulk_DESeq2_Summary.csv\n")
cat("- CLRN1_Retina_Demultiplexed_Clean.rds\n")
if (!requireNamespace("DESeq2", quietly = TRUE)) {
    BiocManager::install("DESeq2")
}
# Load required packages
# Install these first if needed (see repository README)
library(DESeq2)
library(tidyverse)

