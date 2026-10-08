suppressPackageStartupMessages({
    library(GenomicRanges)
    library(IRanges)
    library(S4Vectors)
})

# ==============================================================================
# 1. Building GRanges Objects & Accessing Slots
# ==============================================================================

# Construct a GRanges object representing synthetic genes
genes_gr <- GRanges(
    seqnames = Rle(c("chr1", "chr1", "chr2")),
    ranges = IRanges(
        start = c(1000, 5000, 2000),
        end   = c(3000, 8000, 4000)
    ),
    strand = c("+", "-", "+"),
    gene_id = c("GENE_A", "GENE_B", "GENE_C"),
    score = c(85.2, 92.0, 45.1)
)

cat("--- GRanges Object Structure ---\n")
print(genes_gr)

# Accessing S4 Slots via Getter Functions (NEVER use @)
cat("\nExtracting Chromosomes (seqnames):\n")
print(seqnames(genes_gr))

cat("\nExtracting Coordinates (IRanges):\n")
print(ranges(genes_gr))

cat("\nExtracting Metadata DataFrame (mcols):\n")
print(mcols(genes_gr))


# ==============================================================================
# 2. Strand-Aware Range Transformations
# ==============================================================================

# 2.1 Extract Promoter Regions using flank()
# flank(x, width) gets adjacent regions. start = TRUE gets upstream region.
promoters_gr <- flank(genes_gr, width = 1000, start = TRUE)

cat("\n--- Promoters (1kb Upstream of TSS) ---\n")
# Notice how for GENE_B (minus strand), flank correctly extracts downstream coordinates
print(promoters_gr[, "gene_id"])


# 2.2 Intra-Range Operations (shift & promoters)
# promoters() directly gets [TSS - upstream, TSS + downstream]
promoters_explicit <- promoters(genes_gr, upstream = 2000, downstream = 200)

cat("\n--- Promoters (-2000bp to +200bp around TSS) ---\n")
print(promoters_explicit[, "gene_id"])


# ==============================================================================
# 3. Inter-Range Operations (reduce & disjoin)
# ==============================================================================

# Create overlapping ChIP-seq peaks
peaks_gr <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(
        start = c(1500, 2000, 7000),
        end   = c(2500, 3500, 9000)
    )
)

# reduce() merges overlapping ranges into single union intervals
merged_peaks <- reduce(peaks_gr)

cat("\n--- Merged Peaks (reduce) ---\n")
print(merged_peaks)

# disjoin() breaks overlapping ranges into distinct non-overlapping sub-intervals
disjoined_peaks <- disjoin(peaks_gr)

cat("\n--- Disjoined Intervals (disjoin) ---\n")
print(disjoined_peaks)


# ==============================================================================
# 4. Overlap & Intersection Analysis
# ==============================================================================

# Find which peaks overlap with our gene annotations
overlaps <- findOverlaps(query = peaks_gr, subject = genes_gr)

cat("\n--- Overlap Hits Object ---\n")
print(overlaps)

# Extract matching indices
queryIndices(overlaps) # Indices in peaks_gr
subjectIndices(overlaps) # Indices in genes_gr

# Filter peaks that directly overlap any gene using subsetByOverlaps
peaks_on_genes <- subsetByOverlaps(query = peaks_gr, subject = genes_gr)

cat("\n--- Peaks Overlapping Genes ---\n")
print(peaks_on_genes)
