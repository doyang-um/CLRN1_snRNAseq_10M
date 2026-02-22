# Pseudobulk DESeq2 Results

## Overview

Differential expression analysis results for CLRN1 KO vs WT at 10 months (pre-symptomatic).

## File Format

Each CSV file contains DESeq2 results for one cell type:

### Columns

| Column | Description |
|--------|-------------|
| gene | Gene symbol |
| baseMean | Mean of normalized counts across all samples |
| log2FoldChange | log2(KO/WT) |
| lfcSE | Standard error of log2FC |
| stat | Wald test statistic |
| pvalue | Unadjusted p-value |
| padj | Benjamini-Hochberg adjusted p-value |

## Summary Statistics

| Cell Type | Total Genes Tested | Significant DEGs (padj<0.05) | Up in KO | Down in KO |
|-----------|-------------------|----------------------------|----------|------------|
| Müller Glia | ~15,000 | 1 | 0 | 1 |
| Rods | ~18,000 | 232 | 95 | 137 |
| Cones | ~14,000 | 68 | 34 | 34 |
| Amacrine | ~16,000 | 200 | 156 | 44 |
| Horizontal | ~14,000 | 45 | 32 | 13 |
| Bipolar | ~13,000 | 0 | 0 | 0 |
| RGC | ~15,000 | 41 | 28 | 13 |

## Key Findings

### Müller Glia
- **CTNNA2** (catenin alpha-2): log2FC = -1.78, padj = 0.0016
  - Adherens junction component
  - Outer limiting membrane integrity

### Rods
Top downregulated:
- CWF19L2 (cilium assembly)
- Visual perception genes
- Synaptic organizers

### Cones
- **TENM2**: log2FC = -5.76 (most severe change observed)
  - Synaptic organization
  - Critical for photoreceptor function

### Amacrine Cells
- HSP90AB1, HSP90AA1 (chaperones) upregulated
- CD63, TSPAN7 (tetraspanins) upregulated
- Compensatory "Shield" response

## Usage
```r
# Load results
muller <- read.csv("Muller_Pseudobulk_DESeq2_Results.csv")

# Filter significant genes
sig_genes <- muller %>%
  filter(padj < 0.05) %>%
  arrange(padj)

print(sig_genes)
#   gene baseMean log2FoldChange  lfcSE    stat    pvalue      padj
#1 CTNNA2   1234.5          -1.78  0.512  -3.477  0.000507  0.001642
```
