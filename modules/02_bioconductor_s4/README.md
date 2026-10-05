# Module 02: Bioconductor & S4 Architecture
## Description: Understanding S4 object design, sequence manipulation with Biostrings, and multi-assay data management with SummarizedExperiment.

Bioconductor is an open-source, community-driven ecosystem based on the R statistical programming language that provides software tools for the analysis and comprehension of high-throughput genomic and biological data. 
It provides three main types of packages:
  - **Software**: Analytical tools for processing data like RNA-Seq, single-cell genomics, proteomics, and CRISPR screens.
  - **AnnotationData**: Curated databases that map raw genetic indicators to usable biological information (e.g., matching a probe ID to a human gene name).
  - **ExperimentData**: Real-world biological datasets used for teaching, testing, and benchmarking

S4 is a formal, strict object-oriented programming (OOP) system within R that Bioconductor uses as its primary infrastructure to ensure that complex biological datasets remain organized, valid, and interoperable across different packages. S4 enforces strict rules to prevent "silent errors" (errors that ruin calculations without stopping the code)
The three core components of S4 are:
  - **Classes**: Formal definitions that describe the structure components (called slots), and data validation rules for objects.
  - **Generic function**: Special functions that act as a polite interface, deciding which specific task to run based on the type of data passed to them.
  - **Methods**: The actual blocks of code that executes a specific operation for a defined class under the guidance of a generic function. 

## Common S4 objects:
- `GRanges`: Explicitly manages genomic intervals and chromosomal coordinates.
- `Biostrings`: Efficiently stores and manipulates massive DNA, RNA, or amino acid sequences.
- `SummarizedExperiment`: A powerful matrix container that ties together expression data (like counts) alongside feature and patient metadata.

## 1. Why S4? (The Bioconductor "Safety Net")
Standard R uses S3 classes or simple lists. A list in R lets we store anything, anywhere. That flexibility is dangerous in computational biology: if a user renames a column or deletes sample IDs, a pipeline running overnight will silently break or produce garbage results.
Bioconductor built the S4 Class System to enforce strict structural contracts:

- **Slots**: S4 objects store specific data types in fixed internal compartments called "slots" (e.g., `@assays`, `@colData`).
- **Validation**: We cannot inject mismatched data into an S4 object. If we try to create a dataset where the sample metadata doesn't match the expression matrix columns, R throws an immediate error.

## 2. `SummarizedExperiment`: The 3-Way Synchronization Lock
This is the single most important object in all of Bioconductor. Almost every RNA-Seq (`DESeq2`), single-cell (`Seurat`/`SingleCellExperiment`), or epigenomics pipeline relies on it.

Instead of keeping expression matrices, sample tables, and gene tables in three separate files—risking indexing mismatches—`SummarizedExperiment` locks them together in a 3D grid:

```text
colData(se)
            [Sample Metadata: Treatment, Batch, Age]
                         │
                         ▼
                     ┌───────┐
                     │Samples│ (Columns: 1 to N)
                     └───────┘
                     ┌───────┐
rowData(se) ──────►  │       │
[Gene Metadata:      │Assays │ (Matrix: Counts, TPM)
 Symbol, Chr, GC]    │       │
                     └───────┘
                     ▲
                     │
                     │ Rows: 1 to M (Genes / Features)
```

When we subset a `SummarizedExperiment`, R automatically updates all three dimensions simultaneously:

```r
# Filter for genes on Chromosome 1 AND only Treated samples:
se_sub <- se[rowData(se)$chromosome == "chr1", colData(se)$condition == "Treated"]
```

In standard R, we would have to carefully slice three separate matrices and metadata tables by hand. Here, R slices the expression matrix, trims the sample table (`colData`), and trims the gene table (`rowData`) in a single line without breaking row/column alignment.

## 3. The Accessor Rule: Getters vs. Direct Slots (@)
In the script, we noticed functions like `colData(se)` and `assay(se, "counts")` instead of `se@colData`.
The Rule: We should never use `@` directly in our scripts; because in S4, direct slot access (`se@colData`) bypasses object validation. Bioconductor developers write getter/setter functions (`colData()`, `assay()`, `rowData()`) so that if the internal software changes in a future package update, our code won't break. Getters serve as the safe public API.

## 4. `Biostrings`: Why Plain R Strings ("`ATGC`") Fail
Why not just use standard R vectors like `c("ATGC", "GCTA")` for sequences?
  1. **Memory Overhead**: Standard R string vectors carry massive internal overhead. Storing millions of genomic reads as standard R strings will exhaust our RAM almost immediately.
  2. **Binary Encoding**: A `DNAStringSet` encodes nucleotides using binary representations (2 bits per base instead of a full byte per character), dramatically shrinking memory footprint.
  3. **Genomic Intelligence**: Plain strings don't know biological rules. A `DNAStringSet` enables instant C-accelerated genomic operations across millions of reads:

```r
reverseComplement(dna_seqs) # Instant C-level operation
letterFrequency(dna_seqs, letters = "GC") # C-accelerated GC calculation
```

## 1. Biostrings & Sequence Manipulation (S4 Classes)
```r
dna_seqs <- DNAStringSet(c(
  Read_1 = "ATGCGATCGATCGATCGATCG",
  Read_2 = "GCTAGCTAGCTAGCTA",
  Read_3 = "NNNATGCGATCGATCG",
  Read_4 = "ATGCGATCGATCGAAA"
))

cat("--- Biostrings DNAStringSet ---\n")
print(dna_seqs)

cat("\nSequence Lengths:\n")
print(width(dna_seqs))

cat("\nGC Content Calculation:\n")
gc_content <- letterFrequency(dna_seqs, letters = "GC", as.prob = TRUE)
print(gc_content)

cat("\nReverse Complement of Sequences:\n")
rev_comp <- reverseComplement(dna_seqs)
print(rev_comp)
```

## Explanations:
**Part 1: Creating the S4 Data Container**
- `DNAStringSet(...)` converts a standard R character vector into a highly memory-efficient **S4 collection of DNA sequences**.
- **S4 Validation**: This class automatically enforces strict biological rules. It reads our sequences and validates that they only contain standardized genetic characters (`A`, `T`, `C`, `G`, and `N` for unknown bases). If we accidentally included an invalid letter like `X`, the S4 system would throw an immediate error.

## Part 2: Sequence Operations & Calculations
1. **Sequence Lengths**
```r
print(width(dna_seqs))
```
- `width()` is an S4 accessor method that instantly extracts the exact length (number of base pairs) of each sequence.
- **Output**: It will return a vector of integers corresponding to the lengths: 21 for Read_1, and 16 for Read_2, Read_3, and Read_4.

2. **GC Content Calculation**
```r
gc_content <- letterFrequency(dna_seqs, letters = "GC", as.prob = TRUE)
```
- `letterFrequency(...)` scans each sequence and counts the occurrences of specific nucleotides.
- `letters = "GC"` targets Guanine (G) and Cytosine (C).
- `as.prob = TRUE` returns the value as a percentage/proportion (between 0 and 1) instead of raw counts.
- Output:
  - **Read_1**: 0.476 (~47.6% GC content)
  - **Read_2**: 0.500 (50% GC content)
  - **Read_3 & Read_4**: 0.375 (37.5% GC content. Note: For Read_3, the "N" bases count toward the total length, which lowers its overall GC proportion).

3. **Reverse Complement**
```r
rev_comp <- reverseComplement(dna_seqs)
```

- `reverseComplement()` simulates a critical biological process. It flips each sequence backwards (3' to 5' direction reversed to 5' to 3') and replaces each base with its complementary matching base pair (A $\[\leftrightarrow \]$ T, C $\[\leftrightarrow \]$ G). The S4 class natively knows that an unknown base N pairs with an `N`.
- Output Example: `Read_4` ("`ATGCGATCGATCGAAA`") is flipped and transformed into "`TTTCGATCGATCGCAT`".


## 2. Construction of SummarizedExperiment Objects
# 0. Setting the Seed
```r
set.seed(123)
```
- Ensures reproducibility. Because this script generates random data (using Poisson and uniform distributions below), setting a seed guarantees we get the exact same numbers every time we run it.

1. **The Primary Count Matrix (The Matrix)**
```r
n_genes <- 100
n_samples <- 6
counts_matrix <- matrix(
  rpois(n_genes * n_samples, lambda = 50),
  nrow = n_genes, ncol = n_samples,
  dimnames = list(...)
)
```
- **What it represents**: This simulates raw gene expression data (e.g., RNA-Seq read counts) for 100 genes across 6 samples.
- `rpois(..., lambda = 50)`: Generates random numbers using a Poisson distribution centered around 50, mimicking how DNA/RNA reads are naturally distributed.
- `dimnames`: Sets the row names as `GENE_001` to `GENE_100` and column names as `SAMPLE_1` to `SAMPLE_6`.


2. **Sample Metadata (colData)**
```r
col_data <- DataFrame(
  condition = factor(rep(c("Control", "Treated"), each = 3)),
  batch     = factor(rep(c("B1", "B2", "B1"), times = 2)),
  lib_size  = colSums(counts_matrix),
  row.names = colnames(counts_matrix)
)
```
- **What it represents**: Information about our samples/patients (columns of the matrix).
- `DataFrame()`: An S4 specific version of R's standard data frame, optimized for Bioconductor.
- **Variables**: It assigns 3 samples to a "Control" group and 3 to a "Treated" group, tracks experimental "batch" effects, and calculates the total sequencing depth (`lib_size`) for each sample.
- **Crucial Rule**: The `row.names` of `colData` must match the column names of our count matrix.

3. **Feature Annotation Metadata (rowData)**
```r
row_data <- DataFrame(
  gene_symbol = paste0("Gene", 1:n_genes),
  chromosome  = sample(c("chr1", "chr2", "chrX"), n_genes, replace = TRUE),
  gc_pct      = runif(n_genes, 0.35, 0.65),
  row.names   = rownames(counts_matrix)
)
```

- **What it represents**: Information about our genes (rows of the matrix).
- **Variables**: It assigns human-readable symbols, randomly distributes the genes across three chromosomes (chr1, chr2, chrX), and simulates a GC content percentage between 35% and 65%.
- **Crucial Rule**: The `row.names` of `rowData` must match the row names of our count matrix.

4. **Assembling the `SummarizedExperiment`**

```r
se <- SummarizedExperiment(
  assays  = list(counts = counts_matrix, logcounts = log2(counts_matrix + 1)),
  colData = col_data,
  rowData = row_data
)
```
- This builds the final S4 object.
- `assays`: Notice that it takes a list of matrices. It stores the raw `counts` matrix and computes a normalized `logcounts` matrix simultaneously. Both matrices must share the exact same dimensions.
- **S4 Integrity Check**: When we run this, the S4 system automatically checks that our row names and column names perfectly align. If they don't match, it halts and throws an error to protect our data integrity.

**The Output Structure**
When we `print(se)`, R shows a clean S4 summary detailing:
- **class**: SummarizedExperiment
- **dim**: 100 rows (genes), 6 columns (samples)
- **metadata**: Any experimental details attached
- **assays**: `counts`, `logcounts`
- **rownames / colnames**: The identifiers used
- **rowData names**: `gene_symbol`, `chromosome`, `gc_pct`
- **colData names**: `condition`, `batch`, `lib_size`


## 3. Interacting with S4 Slots via Getter/Setter Accessors
This code block demonstrates how to interact with, extract data from, and subset an assembled `SummarizedExperiment` object using formal S4 accessor methods.

1. Extracting Matrix Data (`assays`)
```r
raw_counts <- assay(se, "counts")
log_counts <- assay(se, "logcounts")
```
- **What it does**: Instead of directly digging into the internal code of the object, we use the standard getter function `assay()`. It looks inside the `se` container and pulls out specific numeric tables by their names ("`counts`" or "`logcounts`").
- **Why it matters**: This ensures us safely grab the raw or normalized data as standard R matrices ready for calculation, without accidentally corrupting the master `se` object.

2. **Accessing Metadata Slots (`colData` and `rowData`)**
```r
sample_info <- colData(se)
gene_info   <- rowData(se)
```
- `colData(se)`: Extracts the entire sample metadata table (conditions, batches, library sizes) representing the columns of our experiment.
- `rowData(se)`: Extracts the entire feature annotation table (chromosomes, gene symbols, GC percentage) representing the rows of our experiment.

3. **Coordinated S4 Subsetting (The Matrix Filter)**
```r
se_filtered <- se[rowData(se)$chromosome == "chr1", colData(se)$condition == "Treated"]
```
This single line showcases the true power of Bioconductor's S4 design. It filters the entire experiment down using standard 2-dimensional matrix coordinates: `object[rows, columns]`.

- **The Rows Rule** (`rowData(se)$chromosome == "chr1"`): Finds all genes located on chromosome 1.
- **The Columns Rule** (`colData(se)$condition == "Treated"`): Finds all experimental samples belonging to the "Treated" group.

**The S4 Magic**: When we run this query, R doesn't just subset the metadata. It **simultaneously cuts down** the counts matrix, the logcounts matrix, the sample rows, and the gene rows all at once. Everything stays perfectly in sync automatically.

4. **Dimensions Output**
```r
cat("\nFiltered SummarizedExperiment Dimensions (chr1 & Treated only):\n")
print(dim(se_filtered))
```
- **What it shows**: This prints out the size of our newly sliced `se_filtered` object.
- **The Math**: Because we filtered for "Treated" samples, the columns will instantly drop from 6 down to 3. The number of rows will drop from 100 down to whatever number of genes were randomly assigned to "chr1" during the generation step.

## Core Takeaway
1. Compress sequences with `Biostrings` instead of standard R text.
2. Bundle matrices + metadata into a `SummarizedExperiment` so we never lose track of sample-to-gene mappings.
3. Interact via getters (`assay()`, `colData()`, `rowData()`) to ensure our code remains robust and maintainable.









































