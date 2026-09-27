# 1. Synthetic Dataset Creation (RNA-Seq Metadata & Expression)

```r
# Sample metadata frame
metadata <- tibble(
  sample_id = paste0("SMP_", 1:12),
  condition = rep(c("Control", "Treated"), each = 6),
  cell_line = rep(c("Hela", "HEK293"), times = 6),
  batch = sample(c("Batch_A", "Batch_B"), 12, replace = TRUE),
  rin_score = round(runif(12, min = 6.8, max = 9.9), 1) # RNA Integrity Number
)
```

## Explanation:
The code uses the `tibble()` function (a modern, user-friendly version of a standard R `data.frame`) to build 5 specific columns for 12 hypothetical samples: 
- `sample_id = paste0("SMP_", 1:12)`
  - What it does: Generates 12 unique sample identifiers labeled `SMP_1`, `SMP_2`, ..., up to `SMP_12`
  - Why it matters: This acts as the "Primary Key" or unique anchor used to link these sample traits directly to their genetic expression data later on.
- `condition = rep(c("Control", "Treated"), each = 6)`
  - What it does: Creates a vector that repeats "Control" 6 times, followed by "Treated" 6 times.
  - Why it matters: This establishes our primary biological comparison group (the independent variable).
- ``cell_line = rep(c("Hela", "HEK293"), times = 6)``
  - What it does: Alternates between two widely used cell lines (`Hela`, `HEK293`, `Hela`, `HEK293`...) for all 12 rows.
  - Why it matters: This introduces a biological covariate. Cancers behave differently depending on the cell type background, so tracking this allows us to control for baseline differences between HeLa (cervical cancer origin) and HEK293 (human embryonic kidney origin).
- `batch = sample(c("Batch_A", "Batch_B"), 12, replace = TRUE)`
  - What it does: Randomly assigns each sample to either `Batch_A` or `Batch_B`.
  - Why it matters: It simulates "batch effects"—technical variations that happen when samples are processed on different days, by different lab technicians, or using different reagent kits.
- `rin_score = round(runif(12, min = 6.8, max = 9.9), 1)`
  - What it does: Uses `runif()` to generate 12 random numbers between 6.8 and 9.9, rounding them to one decimal place.
  - Why it matters: RIN stands for RNA Integrity Number. It ranges from 1 (completely degraded RNA) to 10 (perfectly intact RNA). A score above 7.0 is typically required for high-quality sequencing. This column acts as a Quality Control (QC) metric.


```r
# Gene expression data in wide matrix format (log2 TPM)
expression_wide <- tibble(
  sample_id = paste0("SMP_", 1:12),
  TP53      = rnorm(12, mean = 8.5, sd = 0.6),
  BRCA1     = rnorm(12, mean = ifelse(metadata$condition == "Treated", 11.2, 7.8), sd = 0.8),
  EGFR      = rnorm(12, mean = ifelse(metadata$condition == "Treated", 5.4, 9.1), sd = 0.7),
  GAPDH     = rnorm(12, mean = 13.0, sd = 0.2) # Housekeeping control
)
```

## Explanation:
This codebase simulates a synthetic RNA-Seq gene expression dataset for 12 samples across 4 specific genes. It uses statistical sampling to model realistic biological patterns, such as baseline housekeepers and treatment-responsive biomarkers.
- `expression_wide <- tibble(...)`
  This creates a data frame structured in a "wide" format, where each row represents an individual biological sample, and columns contain the expression measurements for specific genes.
- `sample_id = paste0("SMP_", 1:12)`
  Generates a column of unique sample names matching our metadata exactly (`SMP_1` through `SMP_12`).
- `TP53 = rnorm(12, mean = 8.5, sd = 0.6)`
  Simulates expression for the `TP53` tumor suppressor gene. The `rnorm()` function draws 12 random numbers from a normal (Gaussian) distribution. In this case, TP53 has a steady, moderate expression baseline (average of 8.5) that does not change between Control and Treated samples.
- `BRCA1 = rnorm(12, mean = ifelse(metadata$condition == "Treated", 11.2, 7.8), sd = 0.8)`
Simulates **up-regulation** for the DNA repair gene `BRCA1`.
  - The `ifelse()` statement checks our `metadata` data frame.
  - If a sample is "**Treated**", its average expression jumps to a high level (11.2).
  - If it is a "Control", it stays lower (7.8)
- `EGFR = rnorm(12, mean = ifelse(metadata$condition == "Treated", 5.4, 9.1), sd = 0.7)`
Simulates **down-regulation** for the growth receptor gene `EGFR`. Under the "**Treated**" condition, its expression drops to an average of 5.4, compared to 9.1 in the "**Control**" group.
- `GAPDH = rnorm(12, mean = 13.0, sd = 0.2)`
Simulates a standard **housekeeping gene (GAPDH)**. Housekeeping genes maintain vital cellular structures and exhibit high (`mean = 13.0`), stable expression with incredibly low variance (`sd = 0.2`) regardless of experimental conditions. They serve as excellent normalization controls.

## Biological Significance of the Unit: log2 TPM
The comment mentions the data is in `log2 TPM` (Transcripts Per Million)
- **TPM** normalizes sequencing data by adjusting for both the length of a gene and the total sequencing depth of that run, allowing us to directly compare different samples.
- `log2` transformation is a standard bioinformatics practice. Raw sequencing reads scale exponentially; transforming them to a $\[\log _{2}\]$ scale squashes outliers, stabilizes data variance, and makes a $\(1\text{-unit change}\)$ mathematically equal to a 2-fold difference in gene expression. 
