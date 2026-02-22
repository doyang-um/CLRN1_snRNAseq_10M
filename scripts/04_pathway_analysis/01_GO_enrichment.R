# ===================================================================
# FINAL Comprehensive GO Enrichment + GSEA Analysis
# Maximum sensitivity for Rods, Cones, Muller Glia
# ===================================================================

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)
library(tidyverse)

# Create output directory
output_dir <- "GO_Final_Analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

cat("===================================================================\n")
cat("COMPREHENSIVE GO ENRICHMENT + GSEA ANALYSIS\n")
cat("===================================================================\n\n")

# ===================================================================
# PARAMETERS
# ===================================================================

CELL_TYPES <- c("Rods", "Cones", "Muller_Glia")
PADJ_THRESHOLD <- 0.05
LOGFC_THRESHOLD <- 0.25

# Very relaxed GO parameters
GO_PVALUE_CUTOFF <- 0.2  # Very permissive
GO_QVALUE_CUTOFF <- 0.5  # Very permissive

# GSEA parameters
GSEA_PVAL_CUTOFF <- 0.05

cat("Standard GO Enrichment Parameters:\n")
cat("  GO p-value cutoff:", GO_PVALUE_CUTOFF, "\n")
cat("  GO q-value cutoff:", GO_QVALUE_CUTOFF, "\n\n")

# ===================================================================
# LOAD DATA AND PREPARE GENE LISTS
# ===================================================================

cat("Loading DE results...\n")

de_data <- list()
for (cell_type in CELL_TYPES) {
  file <- paste0("CellType_DE_Analysis_Results/DE_Genes_", cell_type, "_KO_vs_WT.csv")
  if (file.exists(file)) {
    de_data[[cell_type]] <- read.csv(file)
    cat("  ✓", cell_type, "\n")
  }
}

cat("\n")

# ===================================================================
# FUNCTION: Convert to Entrez
# ===================================================================

convert_to_entrez <- function(gene_symbols) {
  gene_symbols <- unique(gene_symbols[!is.na(gene_symbols)])
  entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols,
                       column = "ENTREZID", keytype = "SYMBOL",
                       multiVals = "first")
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  return(entrez_ids)
}

# ===================================================================
# STANDARD GO ENRICHMENT (Over-representation Analysis)
# ===================================================================

cat("===================================================================\n")
cat("PART 1: STANDARD GO OVER-REPRESENTATION ANALYSIS\n")
cat("===================================================================\n\n")

go_summary <- data.frame()

for (cell_type in names(de_data)) {
  
  cat("───────────────────────────────────────────────────────────────\n")
  cat(cell_type, "\n")
  cat("───────────────────────────────────────────────────────────────\n")
  
  de <- de_data[[cell_type]]
  
  # Upregulated
  up_genes <- de %>% filter(p_val_adj < PADJ_THRESHOLD, avg_log2FC > LOGFC_THRESHOLD)
  up_entrez <- convert_to_entrez(up_genes$gene)
  
  # Downregulated  
  down_genes <- de %>% filter(p_val_adj < PADJ_THRESHOLD, avg_log2FC < -LOGFC_THRESHOLD)
  down_entrez <- convert_to_entrez(down_genes$gene)
  
  # Background
  all_entrez <- convert_to_entrez(de$gene)
  
  cat("  Up:", length(up_entrez), "genes | Down:", length(down_entrez), "genes\n\n")
  
  # Analyze UP genes
  if (length(up_entrez) >= 10) {
    cat("  Analyzing UPREGULATED genes...\n")
    
    for (ont in c("BP", "CC", "MF")) {
      go_res <- enrichGO(gene = up_entrez, universe = all_entrez,
                         OrgDb = org.Hs.eg.db, ont = ont,
                         pAdjustMethod = "BH",
                         pvalueCutoff = GO_PVALUE_CUTOFF,
                         qvalueCutoff = GO_QVALUE_CUTOFF,
                         readable = TRUE)
      
      if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
        n_terms <- nrow(as.data.frame(go_res))
        cat("    ✓", ont, ":", n_terms, "terms\n")
        
        # Save
        write.csv(as.data.frame(go_res),
                  file.path(output_dir, paste0(cell_type, "_UP_GO_", ont, ".csv")),
                  row.names = FALSE)
        
        # Plot for BP
        if (ont == "BP") {
          go_df <- as.data.frame(go_res) %>%
            arrange(p.adjust) %>%
            head(20) %>%
            mutate(Description = str_trunc(Description, 60))
          
          p <- ggplot(go_df, aes(x = reorder(Description, -log10(p.adjust)),
                                 y = -log10(p.adjust))) +
            geom_bar(stat = "identity", fill = "#E64B35", alpha = 0.8) +
            coord_flip() + theme_minimal() +
            labs(title = paste(cell_type, "- Upregulated GO BP"),
                 x = NULL, y = "-log10(Adj P)")
          
          ggsave(file.path(output_dir, paste0(cell_type, "_UP_BP_Barplot.pdf")),
                 p, width = 12, height = 8)
        }
        
        go_summary <- rbind(go_summary, data.frame(
          Cell_Type = cell_type, Direction = "UP", Ontology = ont,
          N_Terms = n_terms
        ))
      } else {
        cat("    ✗", ont, ": No terms\n")
      }
    }
    cat("\n")
  }
  
  # Analyze DOWN genes
  if (length(down_entrez) >= 10) {
    cat("  Analyzing DOWNREGULATED genes...\n")
    
    for (ont in c("BP", "CC", "MF")) {
      go_res <- enrichGO(gene = down_entrez, universe = all_entrez,
                         OrgDb = org.Hs.eg.db, ont = ont,
                         pAdjustMethod = "BH",
                         pvalueCutoff = GO_PVALUE_CUTOFF,
                         qvalueCutoff = GO_QVALUE_CUTOFF,
                         readable = TRUE)
      
      if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
        n_terms <- nrow(as.data.frame(go_res))
        cat("    ✓", ont, ":", n_terms, "terms\n")
        
        write.csv(as.data.frame(go_res),
                  file.path(output_dir, paste0(cell_type, "_DOWN_GO_", ont, ".csv")),
                  row.names = FALSE)
        
        if (ont == "BP") {
          go_df <- as.data.frame(go_res) %>%
            arrange(p.adjust) %>%
            head(20) %>%
            mutate(Description = str_trunc(Description, 60))
          
          p <- ggplot(go_df, aes(x = reorder(Description, -log10(p.adjust)),
                                 y = -log10(p.adjust))) +
            geom_bar(stat = "identity", fill = "#4DBBD5", alpha = 0.8) +
            coord_flip() + theme_minimal() +
            labs(title = paste(cell_type, "- Downregulated GO BP"),
                 x = NULL, y = "-log10(Adj P)")
          
          ggsave(file.path(output_dir, paste0(cell_type, "_DOWN_BP_Barplot.pdf")),
                 p, width = 12, height = 8)
        }
        
        go_summary <- rbind(go_summary, data.frame(
          Cell_Type = cell_type, Direction = "DOWN", Ontology = ont,
          N_Terms = n_terms
        ))
      } else {
        cat("    ✗", ont, ": No terms\n")
      }
    }
    cat("\n")
  }
}

# Save GO summary
if (nrow(go_summary) > 0) {
  write.csv(go_summary, file.path(output_dir, "GO_Summary.csv"), row.names = FALSE)
}

# ===================================================================
# PART 2: GSEA (Gene Set Enrichment Analysis)
# ===================================================================

cat("\n===================================================================\n")
cat("PART 2: GENE SET ENRICHMENT ANALYSIS (GSEA)\n")
cat("===================================================================\n\n")

cat("GSEA uses ALL genes (not just significant ones)\n")
cat("Genes are ranked by: sign(log2FC) * (-log10(pvalue))\n\n")

gsea_summary <- data.frame()

for (cell_type in names(de_data)) {
  
  cat("───────────────────────────────────────────────────────────────\n")
  cat(cell_type, "- GSEA\n")
  cat("───────────────────────────────────────────────────────────────\n")
  
  de <- de_data[[cell_type]]
  
  # Create ranked gene list
  # Ranking metric: sign(log2FC) * (-log10(pvalue))
  de$rank_metric <- sign(de$avg_log2FC) * (-log10(de$p_val + 1e-300))
  
  # Convert to Entrez and create named vector
  gene_list <- de$rank_metric
  names(gene_list) <- de$gene
  
  # Convert names to Entrez
  entrez_mapping <- mapIds(org.Hs.eg.db,
                           keys = names(gene_list),
                           column = "ENTREZID",
                           keytype = "SYMBOL",
                           multiVals = "first")
  
  # Keep only mapped genes
  gene_list_entrez <- gene_list[!is.na(entrez_mapping)]
  names(gene_list_entrez) <- entrez_mapping[!is.na(entrez_mapping)]
  
  # Remove duplicates (keep highest ranked)
  gene_list_entrez <- gene_list_entrez[!duplicated(names(gene_list_entrez))]
  
  # Sort by rank
  gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
  
  cat("  Ranked gene list:", length(gene_list_entrez), "genes\n")
  cat("  Range:", round(min(gene_list_entrez), 2), "to",
      round(max(gene_list_entrez), 2), "\n\n")
  
  # Run GSEA for GO BP
  cat("  Running GSEA for GO Biological Process...\n")
  
  gsea_go <- gseGO(geneList = gene_list_entrez,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pvalueCutoff = GSEA_PVAL_CUTOFF,
                   pAdjustMethod = "BH",
                   verbose = FALSE)
  
  if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
    
    n_terms <- nrow(as.data.frame(gsea_go))
    cat("  ✓ Found", n_terms, "enriched GO BP terms\n\n")
    
    # Save results
    write.csv(as.data.frame(gsea_go),
              file.path(output_dir, paste0(cell_type, "_GSEA_GO_BP.csv")),
              row.names = FALSE)
    
    # Show top enriched pathways
    gsea_df <- as.data.frame(gsea_go) %>%
      arrange(pvalue) %>%
      head(10)
    
    cat("  Top 10 enriched pathways:\n")
    for (i in 1:nrow(gsea_df)) {
      cat("    ", i, ". ", gsea_df$Description[i], "\n", sep = "")
      cat("       NES =", round(gsea_df$NES[i], 2),
          "| p =", format(gsea_df$pvalue[i], scientific = TRUE, digits = 2), "\n")
    }
    cat("\n")
    
    # Create enrichment plot for top pathway
    if (nrow(gsea_df) > 0) {
      top_pathway_id <- gsea_df$ID[1]
      
      p_gsea <- gseaplot2(gsea_go,
                          geneSetID = top_pathway_id,
                          title = paste(cell_type, "-", gsea_df$Description[1]),
                          pvalue_table = TRUE)
      
      ggsave(file.path(output_dir, paste0(cell_type, "_GSEA_Top_Pathway.pdf")),
             p_gsea, width = 10, height = 8)
    }
    
    # Create dot plot
    p_dot <- dotplot(gsea_go, showCategory = 20, font.size = 10) +
      ggtitle(paste(cell_type, "- GSEA GO BP"))
    
    ggsave(file.path(output_dir, paste0(cell_type, "_GSEA_Dotplot.pdf")),
           p_dot, width = 12, height = 10)
    
    gsea_summary <- rbind(gsea_summary, data.frame(
      Cell_Type = cell_type,
      N_Pathways = n_terms,
      Top_Pathway = gsea_df$Description[1],
      Top_NES = round(gsea_df$NES[1], 2),
      Top_Pvalue = gsea_df$pvalue[1]
    ))
    
  } else {
    cat("  ✗ No enriched pathways found\n\n")
  }
}

# Save GSEA summary
if (nrow(gsea_summary) > 0) {
  write.csv(gsea_summary, file.path(output_dir, "GSEA_Summary.csv"),
            row.names = FALSE)
}

# ===================================================================
# FINAL SUMMARY
# ===================================================================

cat("===================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("===================================================================\n\n")

cat("Output directory:", output_dir, "\n\n")

if (nrow(go_summary) > 0) {
  cat("Standard GO Enrichment Results:\n")
  print(go_summary)
  cat("\n")
} else {
  cat("Standard GO Enrichment: No results\n\n")
}

if (nrow(gsea_summary) > 0) {
  cat("GSEA Results:\n")
  print(gsea_summary)
  cat("\n")
} else {
  cat("GSEA: No results\n\n")
}

cat("===================================================================\n")
cat("FILES GENERATED:\n")
cat("===================================================================\n")
cat("  Standard GO:\n")
cat("    - [CellType]_[UP/DOWN]_GO_[BP/CC/MF].csv\n")
cat("    - [CellType]_[UP/DOWN]_BP_Barplot.pdf\n")
cat("    - GO_Summary.csv\n\n")
cat("  GSEA:\n")
cat("    - [CellType]_GSEA_GO_BP.csv\n")
cat("    - [CellType]_GSEA_Dotplot.pdf\n")
cat("    - [CellType]_GSEA_Top_Pathway.pdf\n")
cat("    - GSEA_Summary.csv\n")
cat("===================================================================\n")

# Check the CC terms for Rods
rods_cc <- read.csv("GO_Final_Analysis/Rods_UP_GO_CC.csv")
cat("Rods Upregulated - Top 20 Cellular Component terms:\n")
print(rods_cc[1:20, c("Description", "p.adjust", "Count", "GeneRatio")])

# Check the MF terms
rods_mf <- read.csv("GO_Final_Analysis/Rods_UP_GO_MF.csv")
cat("\nRods Upregulated - Molecular Function terms:\n")
print(rods_mf[, c("Description", "p.adjust", "Count", "GeneRatio")])