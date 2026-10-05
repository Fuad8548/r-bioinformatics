# r-bioinformatics

[![R Version](https://img.shields.io/badge/R-4.3%2B-blue.svg)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.18%2B-green.svg)](https://bioconductor.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A structured, self-paced curriculum and hands-on code repository for mastering R and Bioconductor in computational biology, statistical genomics, single-cell analysis, and epigenomics.

---

## 📌 Overview

This repository provides modular, reproducible workflows for genomic analysis using R. It covers fundamental data wrangling, differential gene expression, functional enrichment, single-cell transcriptomics, and publication-ready visualization.

---

## 🛠️ Prerequisites & Installation

* **R (>= 4.3.0)** and **RStudio** or **VS Code** with R extension.
* **renv** for reproducible project dependencies.

### Quick Start

```bash
# Clone the repository
git clone [https://github.com/FUAD8548/r-bioinformatics.git](https://github.com/FUAD8548/r-bioinformatics.git)
cd r-bioinformatics
```

## Run automated workspace initializer
Rscript setup.R

## Repository Structure

r-bioinformatics/
├── README.md
├── r-bioinformatics.Rproj
├── renv.lock
├── setup.R                   # Automated environment & directory generator
├── data/
│   ├── raw/                  # Read-only input datasets
│   └── processed/            # RDS, SummarizedExperiment, Seurat objects
├── R/
│   ├── utils.R               # Reusable helper functions
│   └── plotting_theme.R      # Custom publication-ready ggplot2 themes
├── modules/                  # 10 core technical learning modules
└── outputs/
    ├── figures/              # Exported plots (PDF, PNG)
    └── tables/               # Differential expression summaries, CSVs

## 🎓 Learning Modules
Module 01: Tidyverse & Base R Foundations — Data wrangling with dplyr/tidyr and plotting with ggplot2.

Module 02: Bioconductor & S4 Architecture — S4 classes, Biostrings, and SummarizedExperiment.

Module 03: Genomic Ranges & Annotations — Genomic intervals with GenomicRanges, IRanges, and rtracklayer.

Module 04: NGS Data Handling & Alignment — Reading SAM/BAM files with Rsamtools and GenomicAlignments.

Module 05: Bulk RNA-Seq & Differential Expression — Count normalization and contrast testing with DESeq2 and edgeR.

Module 06: Functional Enrichment Analysis — ORA and GSEA pathway analysis using clusterProfiler and pathview.

Module 07: Single-Cell RNA-Seq (scRNA-Seq) — Quality control, UMAP clustering, and marker analysis with Seurat.

Module 08: Epigenomics & Peak Analysis — Differential binding with DiffBind and peak annotation via ChIPseeker.

Module 09: Microbiome & Metagenomics — Taxonomic profiling and diversity analysis using phyloseq.

Module 10: Capstone Project & Reporting — End-to-end integrated analysis with Quarto (.qmd) and Shiny.


## 🔗 Related Track
🐍 Python Track: git clone [https://github.com/FUAD8548/r-bioinformatics.git](https://github.com/FUAD8548/r-bioinformatics.git) — Sequence alignment algorithms, genomic parsing, and structural biology in Python.


---

### `setup.R`

```r
# Automated Workspace Initializer for r-bioinformatics

cat("=== Initializing 'r-bioinformatics' Repository Workspace ===\n\n")

# 1. Directory Tree Definition
directories <- c(
  "data/raw",
  "data/processed",
  "R",
  "outputs/figures",
  "outputs/tables",
  "modules/01_tidyverse_foundations",
  "modules/02_bioconductor_s4",
  "modules/03_genomic_ranges",
  "modules/04_ngs_processing",
  "modules/05_rnaseq_deseq2",
  "modules/06_enrichment_analysis",
  "modules/07_scrnaseq_seurat",
  "modules/08_chipseq_epigenomics",
  "modules/09_microbiome_phyloseq",
  "modules/10_capstone_project"
)

# 2. Folder Creation and Tracking setup
for (dir in directories) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat(sprintf("[CREATED] %s\n", dir))
  } else {
    cat(sprintf("[EXISTS]  %s\n", dir))
  }
  
  # Add .gitkeep so empty directories stay tracked in Git
  gitkeep_path <- file.path(dir, ".gitkeep")
  if (!file.exists(gitkeep_path)) {
    file.create(gitkeep_path)
  }
}
```

### A concise reference summary of the key concepts learned across Modules 01 through 06

|          Module           |                                Core Mental Model                                 |                                               Key S4 classes & Functions                                                |                                     Practical Output/ Utility                                      |
| :-----------------------: | :------------------------------------------------------------------------------: | :---------------------------------------------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------: |
|  01. Tidyverse & Base R   |            Tidy data wrangling and grammar-of-graphics visualization             | dplyr (filter, mutate, group_by, summarize), tidyr (pivot_longer, pivot_wider), ggplot2 (geom_point, geom_boxplot, aes) |           Cleaned expression matrices and metadata data frames; publication-ready plots.           |
|   02. Bioconductor & S4   | Slot-based object encapsulation separating count assays from row/column metadata | SummarizedExperiment, Biostrings (DNAStringSet), slot access via @ or accessor methods (assay(), colData(), rowData())  |      Standardized omics data containers combining expression values with sample annotations.       |
|    03. Genomic Ranges     |            1-based genomic interval math and coordinate intersections            |                      GRanges, IRanges, rtracklayer::import(), findOverlaps(), reduce(), nearest()                       |      Annotation parsing (GTF/GFF/BED) and overlap queries between peaks, genes, and variants.      |
| 04. NGS Data & Alignment  |        Memory-efficient streaming and indexing of aligned sequence reads         |                  Rsamtools (BamFile, ScanBamParam), GenomicAlignments (readGAlignments()), coverage()                   | Efficient extraction of read alignments, mapping quality filtering, and genomic coverage profiles. |
|    05. Bulk RNA-Seq DE    |      Negative Binomial modeling of raw read counts with variance shrinkage       |                         DESeq2 (DESeqDataSetFromMatrix, DESeq(), results(), lfcShrink()), edgeR                         |  Normalized counts, $log_2\text{FC}$ estimates, adjusted $p$-values (FDR), and volcano/MA plots.   |
| 06. Functional Enrichment |       Mapping DE genes or ranked metrics to biological pathway ontologies        |                       clusterProfiler (enrichGO, enrichKEGG, GSEA), msigdbr, enrichplot, pathview                       |       ORA $p$-values, GSEA Normalized Enrichment Scores (NES), network plots, and KEGG maps.       |


### Core Data Flow Across Modules 01–06

```text
Raw Reads (.bam) ─────────► Genomic Coordinates (.gtf / GRanges)
        │                                     │
   [Module 04]                           [Module 03]
        │                                     │
        └──────────────────┬──────────────────┘
                           ▼
            SummarizedExperiment Container
            ├── Assays: Raw Count Matrix
            └── Metadata: Sample Phenotypes
                       [Module 02]
                           │
                           ▼
          Differential Expression Analysis
          ├── Size Factor Normalization
          └── Negative Binomial Wald Test
                       [Module 05]
                           │
                           ▼
          Significantly DE Genes & Rank Metric
                           │
                           ▼
            Pathway & Functional Enrichment
            ├── Over-Representation Analysis (ORA)
            └── Gene Set Enrichment Analysis (GSEA)
                       [Module 06]
                           │
                           ▼
       Tidy Data Wrangling & Publication Figures
                       [Module 01]
```

















































