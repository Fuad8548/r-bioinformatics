suppressPackageStartupMessages({
    library(Rsamtools)
    library(GenomicAlignments)
    library(GenomicRanges)
    library(SummarizedExperiment)
})

# ==============================================================================
# 1. Connecting to BAM Files & Inspecting Headers
# ==============================================================================

# Locate example BAM file bundled with Rsamtools
bam_path <- system.file("extdata", "ex1.bam", package = "Rsamtools")

# Create a BamFile S4 pointer (does NOT load full file into memory)
bam_file <- BamFile(bam_path)
open(bam_file)

# Extract BAM Header (chromosome lengths, aligner metadata)
header <- scanBamHeader(bam_file)
cat("--- Reference Chromosomes in BAM Header ---\n")
print(header$targets)

close(bam_file)


# ==============================================================================
# 2. Targeted Import with ScanBamParam
# ==============================================================================

# Define a target region of interest (chr1:1 to 5000)
target_region <- GRanges("seq1", IRanges(start = 1, end = 5000))

# Define filtering parameters:
# - High mapping quality (MAPQ >= 20)
# - Exclude unmapped reads and PCR duplicates
# - Restrict to target_region
# - Keep specific fields (cigar, mapq, seq)
param <- ScanBamParam(
    which = target_region,
    what = c("rname", "strand", "pos", "cigar", "mapq"),
    flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE),
    mapqFilter = 20
)

# Read filtered alignments as a GAlignments S4 object
alignments <- readGAlignments(bam_path, param = param)

cat("\n--- Imported GAlignments Summary ---\n")
print(alignments)

# Inspect alignment metadata
cat("\nFirst 5 CIGAR Strings:\n")
print(head(cigar(alignments), 5))

cat("\nWidths of Aligned Reads on Genome:\n")
print(head(width(alignments), 5))


# ==============================================================================
# 3. Handling Spliced Reads & Exon Coordinates
# ==============================================================================

# Convert alignments into GRanges (expands CIGAR N operations into skipped ranges)
read_ranges <- granges(alignments)

cat("\n--- Converted Genomic Ranges from Alignments ---\n")
print(head(read_ranges, 3))


# ==============================================================================
# 4. Read Quantification over Features (summarizeOverlaps)
# ==============================================================================

# Create synthetic gene features (exons grouped by gene in a GRangesList)
gene_exons <- GRangesList(
    Gene_A = GRanges("seq1", IRanges(start = c(100, 800), end = c(500, 1200))),
    Gene_B = GRanges("seq1", IRanges(start = c(1500, 2200), end = c(1900, 2800)))
)

# Quantify read overlaps using 'Union' mode
# Union: A read counts toward a gene if it overlaps any exon of that gene
se_counts <- summarizeOverlaps(
    features = gene_exons,
    reads = bam_path,
    mode = "Union",
    singleEnd = TRUE,
    param = param
)

cat("\n--- Resulting SummarizedExperiment from Count Quantification ---\n")
print(se_counts)

cat("\nCounts Matrix:\n")
print(assays(se_counts)$counts)
