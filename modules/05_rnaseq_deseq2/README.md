# Module 05: Bulk RNA-Seq & Differential Expression with DESeq2

## Description: End-to-end differential gene expression workflow using DESeq2, covering data pre-filtering, design formulas, dispersion fitting, hypothesis testing, and Log2 Fold Change (LFC) shrinkage.

## The Core Goal
We have a table of raw read counts for thousands of genes across multiple samples (e.g., 3 Control vs. 3 Treated). We want to determine which genes are significantly higher or lower in expression due to treatment, rather than random biological or technical noise.

## The 5-Step Engine of DESeq2
1. **Why Raw Integer Counts? (The Count Model)**
RNA-Seq data consists of discrete read counts, not continuous values.
   - **Why not simple $t$-tests?** Standard $t$-tests assume data follows a normal (bell-curve) distribution. Count data is heavily skewed, discrete, and bound by zero.
   - **The Negative Binomial Distribution**: DESeq2 models counts using a Negative Binomial distribution because it accounts for two types of variance:
        $\text{Variance} = \text{Mean} + (\text{Dispersion} \times \text{Mean}^2)$
      - Shot noise (Poisson randomness from sequencing).
      - Biological variance (differences between actual biological replicates).

2. **Median-of-Ratios Normalization** (`sizeFactors`)
Sample A might have 20 million reads total, while Sample B got 40 million. We can't compare raw counts directly.
   - Why NOT CPM / Total Reads? If one single gene (like albumin) makes up 30% of Sample A's reads, it artificially inflates total counts and makes every other gene look downregulated.
   - How DESeq2 Fixes It: It computes a pseudo-reference sample using the geometric mean across all samples for each gene. Each sample's size factor is calculated as the median ratio of its counts to the reference. This renders normalization immune to a few hyper-expressed outlier genes.

3. Dispersion Shrinkage (Information Borrowing)
With only 3 replicates per group, estimating the true variance for 20,000 genes individually is unreliable—some genes will randomly show zero variance, producing false positives.
   - DESeq2's Solution: It assumes genes with similar expression levels share similar dispersion. It fits a trend curve across all genes and "shrinks" individual gene variance estimates toward the global curve.

```text
Gene Dispersion ──► Shrink toward Global Trend Line ──► Stable Variance Estimate
```

4. **Differential Expression & Multiple Testing (padj)**
DESeq2 fits a Generalized Linear Model (GLM) for each gene and performs a Wald Test to see if $\text{Log}_2(\text{Fold Change}) \neq 0$.
   - **The False Discovery Rate (FDR)**: Testing 20,000 genes at a standard $p < 0.05$ significance level means we'll get 1,000 false positives by pure chance.
   - **Adjusted $p$-value (`padj`)**: DESeq2 applies the Benjamini-Hochberg correction to control the False Discovery Rate. A `padj < 0.05` means that no more than 5% of the genes called "significant" are expected to be false positives.

5. **Log2 Fold Change (LFC) Shrinkage**
If a gene goes from 1 read in Control to 4 reads in Treated, that looks like a $\text{Log}_2(\text{Fold Change}) = 2$ (a 4-fold jump). But 1 read vs 4 reads is mostly random noise.
   - `lfcShrink()` pulls noisy, low-count fold changes back toward zero while leaving high-count, high-confidence genes untouched. This is critical for accurate gene ranking and volcano plots.

## How the Code Translates to This Logic

|           Function            |                                            What it actually does                                             |
| :---------------------------: | :----------------------------------------------------------------------------------------------------------: |
|   DESeqDataSetFromMatrix()    |                Locks the raw count matrix and sample metadata into a synchronized S4 object.                 |
| dds$condition <- relevel(...) |    Tells R which group is the baseline (Control) so positive fold-changes mean "upregulated in Treated".     |
|       dds <- DESeq(dds)       | Executes steps 2, 3, and 4 automatically: estimates size factors, fits dispersions, and runs Wald GLM tests. |
|          lfcShrink()          |                       Executes step 5: dampens noisy fold changes for low-count genes.                       |


### loading core libaries
```r
suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(tidyverse)
})

set.seed(42)
```


# 1. Synthetic Dataset Generation (Unnormalized Counts)

```r
n_genes <- 1000
n_samples <- 6

# Generate baseline Poisson counts across genes
counts_matrix <- matrix(
  rpois(n_genes * n_samples, lambda = 100),
  nrow = n_genes,
  ncol = n_samples,
  dimnames = list(
    paste0("GENE_", sprintf("%04d", 1:n_genes)),
    paste0("SAMPLE_", 1:n_samples)
  )
)

# Inject true differential expression into the first 100 genes for Treated group
counts_matrix[1:100, 4:6] <- as.integer(counts_matrix[1:100, 4:6] * runif(100, 2.5, 5.0))

# Sample metadata
col_data <- DataFrame(
  condition = factor(rep(c("Control", "Treated"), each = 3)),
  batch     = factor(rep(c("B1", "B2", "B1"), times = 2)),
  row.names = colnames(counts_matrix)
)
```

### Explanations
In this step, we are creating a mock RNA-Seq experiment from scratch. Generating synthetic data is a standard way to test differential expression workflows because we know exactly which genes are altered. This lets we confirm if `DESeq2` can accurately find them later.

1. **Setting Up the Matrix Dimensions**
```r
n_genes <- 1000
n_samples <- 6
```
- What it means: We are simulating a mini-genome containing 1,000 genes across 6 distinct samples.

2. **Simulating Background Noise (The Poisson Matrix)**

   - `rpois(..., lambda = 100)`: In RNA-Seq, read counts are discrete integers (we can't have half a sequence read). The **Poisson distribution** (`rpois`) is used here to generate random whole numbers that fluctuate naturally around an average baseline value of **100 counts per gene**.
   - `dimnames = list(...)`: This clean code labels our rows sequentially from `GENE_0001` to `GENE_1000` and columns from `SAMPLE_1` to `SAMPLE_6`.

3. **Injecting the Biological Ground Truth (True DE Genes `counts_matrix[1:100, 4:6]`)**
   - **What it does**: This is where we rig the experiment. We select the first **100 genes** across the last 3 samples (`4:6`, which will be our `Treated` group).
   - **The Math**: We multiply their baseline counts by a random scaling factor between 2.5x and 5.0x (`runif(100, 2.5, 5.0)`).
   - **The Biological Result**: We have just engineered 100 up-regulated genes in the treated group. The remaining 900 genes (`GENE_0101` to `GENE_1000`) remain unchanged between groups, acting as our true biological negatives.

4. **Creating the Experimental Metadata Table** (`col_data`)
Every differential expression model requires a design table explaining what our samples actually represent. `DataFrame` is an S4 metadata container from `SummarizedExperiment`. If we print `col_data`, it will map out our experimental design precisely:
   - `condition`: Classifies samples 1–3 as `Control` and samples 4–6 as `Treated`.
   - `batch`: Injects a common real-world complication—batch effects. It simulates that the samples were processed across two separate experimental runs (`B1` and `B2`), which creates artificial technical variance that `DESeq2` will need to correct for.


# 2. Constructing DESeqDataSet & Pre-Filtering

```r
# Build DESeqDataSet directly from count matrix and metadata
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData   = col_data,
  design    = ~ condition
)

# Pre-filter low count genes (at least 10 reads across all samples)
keep_genes <- rowSums(counts(dds)) >= 10
dds <- dds[keep_genes, ]

cat("--- Pre-Filtered DESeqDataSet Dimensions ---\n")
print(dim(dds))

# Explicitly lock the reference level for statistical contrasts
dds$condition <- relevel(dds$condition, ref = "Control")
```

### Explanations
In this step, we are preparing our synthetic matrix for statistical analysis by wrapping it in a specialized object, removing low-quality genes, and locking in our control group.

1. **Building the** `DESeqDataSet` (`DESeqDataSetFromMatrix`)
   - **What it does**: This links our raw data (`counts_matrix`) with our sample data (`col_data`).
   - `design = ~ condition`: This is the core formula telling DESeq2 exactly what we want to measure. By setting it to `~ condition`, we instruct the model to calculate how gene expression changes across our experimental states (`Control` vs. `Treated`), while ignoring the `batch` column for now.

2. **Pre-Filtering Low-Count Genes**
```r
keep_genes <- rowSums(counts(dds)) >= 10
dds <- dds[keep_genes, ]
```
- **Why this is done**: In real-world data, many genes have zero or very few sequence reads because they aren't expressed in our tissue. Keeping them wastes computing power and hurts our statistical sensitivity.
- **The Result (`1000 6`)**: Our output dimension remains exactly 1,000 genes across 6 samples. Because our synthetic baseline was generated with a high average count (`lambda = 100` from Section 1), every single gene easily passed the 10-read threshold, so zero genes were discarded.

3. **Setting the Reference Level (relevel)**
```r
dds$condition <- relevel(dds$condition, ref = "Control")
```
- **Why this is critical**: By default, R sorts factors alphabetically. Without this line, `DESeq2` would alphabetically choose Control as the baseline anyway, but if our groups were named `A_Treated` and `B_Control`, R would accidentally treat the treated group as the baseline!
- **What it fixes**: This explicitly locks `Control` as the comparison anchor. Consequently, any downstream fold-changes will be calculated logically as:
      $\text{Fold Change} = \frac{\text{Treated}}{\text{Control}}$

A positive log-fold change will mean a gene went up during treatment.


# 3. Running DESeq2 Analysis Pipeline

```r
# Executing size factor estimation, dispersion estimation, and GLM fitting
dds <- DESeq(dds)

# Inspect size factors generated by median-of-ratios normalization
cat("\nEstimated Size Factors per Sample:\n")
print(sizeFactors(dds))
```

## Output: 
```bash
SAMPLE_1  SAMPLE_2  SAMPLE_3  SAMPLE_4  SAMPLE_5  SAMPLE_6 
0.9896984 0.9923298 0.9912324 1.0106375 1.0183030 1.0151069
```

## Explanation
In this step, we execute `DESeq(dds)`, which is the core processing engine of the entire workflow. Instead of making we run multiple individual math equations, this master function wraps three major statistical steps into a single command to find true differential expression.

1. **Estimating Size Factors**
   - The Problem: Samples never have the exact same number of total sequencing reads because of random differences during library preparation. If Sample 4 has twice as many total reads as Sample 1 simply because it sat in the sequencer longer, its raw numbers will look inflated.
   - The Solution: DESeq2 uses a robust method called median-of-ratios normalization. It creates a pseudo-reference genome across all our samples, calculates ratio scores for every gene, and calculates a scaling baseline factor (Size Factor) for each column.

2. **Estimating Dispersions**
   - The Problem: In sequencing data, genes with higher expression naturally have higher, more unpredictable variance. Traditional statistical tests (like standard t-tests) fall apart here.
   - The Solution: DESeq2 measures dispersion—a specialized parameter that quantifies how much a gene’s expression fluctuates across biological replicates relative to its average expression. It uses a method called "shrinkage" to pool information across all 1,000 genes, sharing statistical strength to make variance estimates highly accurate even with small sample sizes.

3. GLM Fitting & Wald Test 
- The Problem: Our counts are discrete whole numbers, meaning standard linear curves don't fit them properly.- The Solution: It fits a **Generalized Linear Model (GLM)** following a Negative Binomial distribution for every gene. Once the curve is shaped based on our design (~ condition), it runs a Wald test to calculate a p-value for each gene, asking: "*Is the change in expression between Control and Treated statistically greater than zero?*" 

### Output
**What our Size Factors output tells us:**
  - Values near `1.0`: Because we generated our synthetic matrix using a perfectly uniform baseline (`lambda = 100` for every single sample), the overall sequencing depth across all 6 samples is almost perfectly equal.
  - How DESeq2 uses them: When DESeq2 performs its downstream calculations, it divides the raw count values of `SAMPLE_1` by `0.9896984` to slightly downscale them, and multiplies/divides SAMPLE_6 counts to scale them up. This puts all 6 samples on a perfectly level playing field before statistical testing.


# 4. Extracting Results & Applying LFC Shrinkage

```r
# Extract standard Wald test results for Treated vs Control at FDR threshold 0.05
res_raw <- results(dds, contrast = c("condition", "Treated", "Control"), alpha = 0.05)

cat("\n--- Raw Differential Expression Summary ---\n")
summary(res_raw)

# Apply Log2 Fold Change Shrinkage (normal / apeglm) to penalize noisy low-count genes
res_shrunk <- lfcShrink(
  dds = dds,
  contrast = c("condition", "Treated", "Control"),
  type = "normal"
)

# Convert results into a tidy data frame
de_results_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column(var = "gene_id") %>%
  drop_na(padj) %>%
  arrange(padj)

cat("\nTop 5 Differentially Expressed Genes:\n")
print(head(de_results_df, 5))
```

### Output
```bash
gene_id     baseMean log2FoldChange     lfcSE     stat       pvalue
1 GENE_0004 316.7264       2.542961 0.1415551 17.93193 6.642728e-72
2 GENE_0031 304.1964       2.373154 0.1400972 16.91322 3.595022e-64
3 GENE_0067 289.5774       2.294850 0.1400201 16.36512 3.393426e-60
4 GENE_0024 303.5361       2.299423 0.1418843 16.18375 6.566794e-59
5 GENE_0003 276.9302       2.158899 0.1388479 15.52820 2.235580e-54
          padj
1 6.642728e-69
2 1.797511e-61
3 1.131142e-57
4 1.641699e-56
5 4.471159e-52
```

### Explanation
In this step, we are extracting the statistical results of our differential expression model, applying an essential noise-reduction technique (shrinkage), and ranking our genes to find the top hits.
1. **Extracting the Comparison** (`results`)
```r
res_raw <- results(dds, contrast = c("condition", "Treated", "Control"), alpha = 0.05)
```

  - `contrast`: This explicitly tells the function to calculate the log2 fold change as $\(\log_2(\text{Treated} / \text{Control})\)$.
  - `alpha = 0.05`: This sets our False Discovery Rate (FDR) target threshold to 5%. DESeq2 uses this value behind the scenes to optimize independent filtering of genes that have too few counts to ever achieve statistical significance. 

2. **Why Log2 Fold Change Shrinkage is Necessary**(`lfcShrink`)
```r
res_shrunk <- lfcShrink(dds = dds, contrast = c(...), type = "normal")
```

   - **The Problem**: In sequencing datasets, genes with very low or noisy read counts can randomly show massive fold changes just by chance (e.g., a gene going from 1 read to 5 reads looks like a massive 5x increase, but it is actually just background noise). If we sort our results purely by raw fold change, these noisy genes will crowd the top of our list.
   - The Solution: `lfcShrink` applies a Bayesian information-sharing framework. It pulls down ("shrinks") the fold-change values of highly variable, low-count genes toward zero, while leaving reliable high-count genes completely untouched. This step is critical for downstream visualizations like Volcano plots. 

3. **Tidying up the Data Frame**
```r
de_results_df <- as.data.frame(res_shrunk) %>% ... %>% drop_na(padj) %>% arrange(padj)
```
We convert the specialized S4 object into a standard tidyverse `data.frame`, strip out genes that couldn't be statistically tested (which get assigned an `NA` in adjusted p-values), and sort them from the most significant to least significant using `padj`.

4. **Decoding our Top 5 Hits (The Output)**
Looking at the `gene_id` names in our output table: `GENE_0004`, `GENE_0031`, `GENE_0067`, `GENE_0024`, `GENE_0003`. Our simulation successfully worked! All five of these top-ranked genes fall cleanly between 1 and 100—the exact window where we manually injected artificial differential expression in Section 1.
Let's break down the metrics for our top row (`GENE_0004`):

   - `baseMean (316.72)`(**The Baseline Filter**): The average normalized read count for this gene across all 6 samples combined. It is a high number, which gives the statistical test high confidence.
	It acts as a safety filter. If a gene has a huge fold change but a baseMean near 0, it means it was just random noise in a single sample. DESeq2 uses baseMean to weed out low-expression junk genes before computing statistics.

   - `log2FoldChange (2.54)`: A value of `2.54` translates to a real-world expression multiplier of $\(2^{2.54} \approx 5.8\)$. This tells us that `GENE_0004` is heavily up-regulated in our Treated group compared to our Control group. 
   - 
   - `lfcSE (0.14)`: The standard error of the log2 fold change. At `0.14`, it is tiny, indicating that the expression change was highly consistent across our replicates.
     - Why it's necessary: It acts as the penalty metric. If our 3 Treated samples have wildly different numbers, the lfcSE will be high. This is the exact number used during LFC Shrinkage (lfcShrink) to push shaky, unpredictable fold changes back toward zero so they don't corrupt our plot.
   - `stat` (**The Wald Statistic**): This is our raw mathematical score. It is calculated simply as: $\(\text{stat} = \frac{\text{log2FoldChange}}{\text{lfcSE}}\)$.
     - Why it's necessary: It scales the fold change against its uncertainty. A high `stat` score means the gene's change is huge and remarkably consistent across all samples. This score is the raw mathematical input used to compute the next step: the p-value.

   - `pvalue (6.64e-72)`: The raw probability that is the direct result of the Wald statistic. It tells we how likely it is to see this change completely by random chance.

   - `padj (6.64e-69)`: The adjusted p-value (using the **Benjamini-Hochberg** correction). Because we are testing 1,000 genes at once, we run into the multiple testing problem. The `padj` controls our false discovery rate. A value this low (6.64 × 10⁻⁶⁹) means this gene is a highly dependable, true biological discovery.


# Explanation for the volcano plot
1. **The Threshold Cutoff Lines**
   - The Horizontal Line (`-log10(0.05`): This line rests low on our Y-axis at the value `1.30`. Any gene point floating above this horizontal boundary is statistically significant (`padj < 0.05`).
   - **The Vertical Lines** (`-1` and `1`): These segment our chart into 3 vertical zones. Genes must cross outside these lines to hit a biological fold-change threshold of doubling or halving expression levels.
   - **Interpretation**: Our top 5 genes (`GENE_0004`, `GENE_0031`, etc.) will sit high above the horizontal line and far to the right of the `1` vertical line, landing deeply in the upper-right "Up-regulated" quadrant.

**Why we used `padj` instead of `pvalue`**?

- `pvalue` (Raw P-value): This evaluates one single gene in isolation. A p-value of 0.05 means there is a 5% chance that this gene's difference is a random fluke.
- `padj` (Adjusted P-value / FDR): This corrects for the False Discovery Rate when testing thousands of genes at the exact same time.

Imagine we put 20,000 people in a room and ask them all to flip a coin 10 times. By pure random chance, a handful of those 20,000 people will flip 10 heads in a row. In isolation, their raw pvalue looks amazing. But in reality, it was just a statistical certainty because we ran the test so many times.

When analyzing RNA-Seq, we are testing roughly 20,000 genes simultaneously. If we use a raw `pvalue < 0.05`, about **1,000 genes will look significant purely by random chance (false positives)**.
The `padj` column applies a mathematical correction (like the Benjamini-Hochberg method) to penalize the raw p-values based on how many tests we ran. It guarantees that if we select genes where `padj < 0.05`, only 5% of our entire final list will be false positives.
Summary: We always use `padj` for our Y-axis and our filtering. Using raw `pvalue` will fill our results with **false-positive noise**.




























