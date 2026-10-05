suppressPackageStartupMessages({
    library(Biostrings)
    library(SummarizedExperiment)
    library(S4Vectors)
})

# 1. Create a DNAStringSet object containing synthetic genomic reads
dna_seqs <- DNAStringSet(c(
    Read_1 = "ATGCGATCGATCGATCGATCG",
    Read_2 = "GCTAGCTAGCTAGCTA",
    Read_3 = "NNNATGCGATCGATCG",
    Read_4 = "ATGCGATCGATCGAAA"
))

cat("--- Biostrings DNAStringSet ---\n")
print(dna_seqs)

# 2. S4 Accessors and Sequence Operations
# Sequence Operations & Calculations
print(width(dna_seqs))

# GC Content Calculation
gc_content <- letterFrequency(dna_seqs, letters = "GC", as.prob = TRUE)
print(gc_content)

# Reverse Complement
rev_comp <- reverseComplement(dna_seqs)
print(rev_comp)


# 2. Construction of SummarizedExperiment Objects
# setting the seed to ensure reproducibility
set.seed(123)

# 2.1 Primary Count Matrix (Rows = Genes, Columns = Samples)
n_genes <- 100
n_samples <- 6

counts_matrix <- matrix(
    rpois(n_genes * n_samples, lambda = 50),
    nrow = n_genes,
    ncol = n_samples,
    dimnames = list(
        paste0("GENE_", sprintf("%03d", 1:n_genes)),
        paste0("SAMPLE_", 1:n_samples)
    )
)

# 2.2 Sample Metadata (colData)
col_data <- DataFrame(
    condition = factor(rep(c("Control", "Treated"), each = 3)),
    batch     = factor(rep(c("B1", "B2", "B1"), times = 2)),
    lib_size  = colSums(counts_matrix),
    row.names = colnames(counts_matrix)
)

# 2.3 Feature Annotation Metadata (rowData)
row_data <- DataFrame(
    gene_symbol = paste0("Gene", 1:n_genes),
    chromosome  = sample(c("chr1", "chr2", "chrX"), n_genes, replace = TRUE),
    gc_pct      = runif(n_genes, 0.35, 0.65),
    row.names   = rownames(counts_matrix)
)

# 2.4 Assemble the SummarizedExperiment Object
se <- SummarizedExperiment(
    assays  = list(counts = counts_matrix, logcounts = log2(counts_matrix + 1)),
    colData = col_data,
    rowData = row_data
)

cat("\n--- SummarizedExperiment Overview ---\n")
print(se)


# 3. Interacting with S4 Slots via Getter/Setter Accessors
# 3.1 Extracting Matrix Data (`assays`)
raw_counts <- assay(se, "counts")
log_counts <- assay(se, "logcounts")

# 3.2 Accessing Metadata Slots (`colData` and `rowData`)
sample_info <- colData(se)
gene_info <- rowData(se)

# 3.3 Coordinated S4 Subsetting (The Matrix Filter)
se_filtered <- se[rowData(se)$chromosome == "chr1", colData(se)$condition == "Treated"]

# 4.4 Dimensions Output
print(dim(se_filtered))
