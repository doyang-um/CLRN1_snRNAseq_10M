# CLRN1 Knockout Rabbit Retina snRNA-seq Analysis

[![DOI](badge-will-be-added)]()

## Overview

Complete analysis code and results for:

**"Single-nucleus transcriptomics reveals Müller glia-driven non-cell-autonomous mechanisms in CLRN1-associated retinal degeneration"**

David O. Yang, et al.  
*Journal of Clinical Investigation* (2026)

## Experimental Design

- **Model**: CLRN1 knockout rabbit (23-bp frameshift deletion)
- **Timepoint**: 10 months (pre-symptomatic, before OCT/ERG changes)
- **Samples**: n=3 wild-type, n=3 knockout
- **Platform**: 10x Genomics Chromium Single Cell 3' v3.1
- **Tissue**: Retina visual streak region (2-3mm ventral to medullary rays)
- **Sequencing**: Illumina NovaSeq X Plus
- **Total nuclei analyzed**: 6,542 post-QC

## Key Findings

### The "Anchor-Shield" Model

1. **CLRN1 is exclusively expressed in Müller glia** (not photoreceptors)
   - Challenges assumption of photoreceptor-autonomous disease

2. **Primary defect: "Anchor" failure in Müller glia**
   - CTNNA2 downregulation (log2FC = -1.78, padj = 0.0016)
   - Outer limiting membrane adherens junction disruption

3. **Secondary photoreceptor vulnerability**
   - Rods: 232 DEGs (synaptic loss: TENM2, CWF19L2)
   - Cones: 68 DEGs (severe TENM2 downregulation, log2FC = -5.76)

4. **Inner retinal "Shield" compensation**
   - Amacrine: 200 DEGs (HSP90AB1, CD63, TSPAN7 upregulation)
   - Horizontal: 45 DEGs
   - RGC: 41 DEGs
   - Bipolar: 0 DEGs (protected)

### Cell Type Distribution

| Cell Type | Nuclei Count | DEGs (padj<0.05) | Direction |
|-----------|--------------|------------------|-----------|
| Müller Glia | ~800 | 1 | ↓ CTNNA2 |
| Rods | ~3,500 | 232 | Mixed |
| Cones | ~600 | 68 | Mostly ↓ |
| Amacrine | ~900 | 200 | Mostly ↑ |
| Horizontal | ~300 | 45 | Mixed |
| Bipolar | ~300 | 0 | - |
| RGC | ~150 | 41 | Mixed |

## Repository Structure
```
CLRN1_snRNAseq_GitHub/
├── README.md
├── LICENSE
├── results/                              # Pre-computed results
│   ├── DESeq2/                          # Main differential expression
│   │   ├── Muller_Pseudobulk_DESeq2_Results.csv
│   │   ├── Rods_Pseudobulk_DESeq2_Results.csv
│   │   ├── Cones_Pseudobulk_DESeq2_Results.csv
│   │   ├── Amacrine_Pseudobulk_DESeq2_Results.csv
│   │   ├── Horizontal_Pseudobulk_DESeq2_Results.csv
│   │   ├── Bipolar_Pseudobulk_DESeq2_Results.csv
│   │   └── RGC_Pseudobulk_DESeq2_Results.csv
│   ├── GO_enrichment/                   # Pathway analysis
│   ├── QC_metrics/                      # Quality control stats
│   └── supplementary/                   # Additional tables
├── scripts/                              # Analysis pipeline
│   ├── 01_quality_control/
│   ├── 02_clustering/
│   ├── 03_differential_expression/
│   ├── 04_pathway_analysis/
│   └── 05_visualization/
├── data/
│   └── sample_metadata.csv
└── docs/
    └── analysis_pipeline.md
```

## Data Availability

### Raw Data

- **GEO Accession**: GSE[XXXXX]
- **URL**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE[XXXXX]

Raw FASTQ files, processed count matrices, and Seurat objects available on GEO.

### Results Files

Pre-computed differential expression and pathway analysis results are included in this repository under `results/`.

## Software Requirements

### Environment
- R ≥ 4.2.0
- RStudio (recommended)

### R Packages
```r
# CRAN packages
install.packages(c(
  "Seurat",        # ≥ 4.3.0
  "tidyverse",
  "ggplot2",
  "pheatmap",
  "viridis",
  "patchwork",
  "cowplot"
))

# Bioconductor packages
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c(
  "DESeq2",           # ≥ 1.38.0
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
))
```

## Installation
```bash
# Clone repository
git clone https://github.com/davidoyang/CLRN1_snRNAseq_GitHub.git
cd CLRN1_snRNAseq_GitHub

# Install R packages (if needed)
Rscript -e "source('install_packages.R')"
```

## Usage

### Using Pre-computed Results

The easiest way to explore findings:
```r
# Load main DESeq2 results
muller <- read.csv("results/DESeq2/Muller_Pseudobulk_DESeq2_Results.csv")
rods <- read.csv("results/DESeq2/Rods_Pseudobulk_DESeq2_Results.csv")
cones <- read.csv("results/DESeq2/Cones_Pseudobulk_DESeq2_Results.csv")

# Filter significant genes
muller_sig <- muller %>% filter(padj < 0.05)
# Result: CTNNA2 only

# Top rod DEGs
rods %>% 
  filter(padj < 0.05) %>% 
  arrange(padj) %>% 
  head(10)
```

### Reproducing Full Analysis

Download data from GEO first, then:
```bash
# 1. Quality control
Rscript scripts/01_quality_control/01_QC_metrics_analysis.R

# 2. Cell type identification
Rscript scripts/02_clustering/01_cell_type_identification.R

# 3. Differential expression
Rscript scripts/03_differential_expression/01_pseudobulk_DESeq2.R

# 4. Pathway analysis
Rscript scripts/04_pathway_analysis/01_GO_enrichment.R

# 5. Generate figures
Rscript scripts/05_visualization/01_generate_figures.R
```

## Key Analysis Details

### Pseudobulk Differential Expression

- **Method**: DESeq2 with pseudobulk aggregation
- **Aggregation**: Sum counts by sample within each cell type
- **Design**: ~ genotype
- **Filtering**: Genes with ≥10 counts in ≥3 samples
- **Multiple testing**: Benjamini-Hochberg FDR correction
- **Threshold**: padj < 0.05

### Quality Control Filters

- Genes detected per nucleus: >500
- Total counts per nucleus: 500-50,000
- Mitochondrial content: <5%
- Doublet removal: SNP-based with Vireo

### Cell Type Annotation

Based on canonical markers:
- **Müller glia**: GLUL, RLBP1, SLC1A3
- **Rods**: RHO, NRL, SAG
- **Cones**: ARR3, OPN1SW, GNAT2
- **Amacrine**: GAD1, SLC6A9, CHAT
- **Horizontal**: LHX1, ONECUT1
- **Bipolar**: VSX2, CABP5, GRM6
- **RGC**: RBPMS, POU4F1, NEFL

## Citation
```bibtex
@article{Yang2026,
  title={Single-nucleus transcriptomics reveals Müller glia-driven non-cell-autonomous mechanisms in CLRN1-associated retinal degeneration},
  author={Yang, David O and [...]},
  journal={Journal of Clinical Investigation},
  year={2026},
  doi={[will be added]}
}
```

## Contact

- **David O. Yang, M.D., Ph.D.**
- Email: doyang@umich.edu
- Institution: University of Michigan
- **Issues**: https://github.com/davidoyang/CLRN1_snRNAseq_GitHub/issues

## License

MIT License - see LICENSE file

## Acknowledgments

- University of Michigan Advanced Genomics Core
- Great Lakes High Performance Computing
- Funding: [Add grant numbers]

## Related Resources

- Manuscript preprint: [Link when available]
- Interactive browser: [If applicable]
- Longitudinal phenotyping data: [Separate repository if applicable]
