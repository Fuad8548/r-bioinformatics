suppressPackageStartupMessages({
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(GenomicRanges)
    library(DiffBind)
    library(ggplot2)
})

set.seed(42)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# ==============================================================================
# PART 1: PEAK ANNOTATION WITH CHIPSEEKER
# ==============================================================================
cat("\n--- Part 1: Peak Annotation using ChIPseeker ---\n")

# Simulate a GRanges peak dataset (e.g., MACS2 output intervals)
chroms <- c("chr1", "chr2", "chr3")
starts <- c(1000000, 2500000, 5000000, 12000000, 18000000)
ends <- starts + sample(200:800, 5)

sim_peaks <- GRanges(
    seqnames = Rle(chroms, c(2, 2, 1)),
    ranges   = IRanges(start = starts, end = ends),
    strand   = Rle(c("*"))
)

cat("Simulated Peak Count:", length(sim_peaks), "\n")

# Annotate peaks relative to hg38 genomic features
peak_anno <- annotatePeak(
    peak         = sim_peaks,
    tssRegion    = c(-3000, 3000),
    TxDb         = txdb,
    annoDb       = "org.Hs.eg.db",
    verbose      = FALSE
)

# Convert annotation summary to a readable data frame
anno_df <- as.data.frame(peak_anno)
cat("\nAnnotated Peak Sample Output:\n")
print(head(anno_df[, c("seqnames", "start", "end", "annotation", "SYMBOL", "distanceToTSS")]))

# Diagnostic Visualizations
# 1. Barplot of genomic feature distributions (Promoter, Intron, Exon, etc.)
p_anno <- plotAnnoBar(peak_anno) +
    labs(title = "Genomic Distribution of Peaks")

# 2. Plot Distance to TSS
p_tss <- plotDistToTSS(peak_anno) +
    labs(title = "Peak Distribution Relative to TSS")

print(p_anno)
print(p_tss)

# 1. Save the genomic distribution barplot as a PDF
ggsave(
    filename = "genomic_distribution_barplot.pdf",
    plot     = p_anno,
    width    = 8,
    height   = 5
)

# 2. Save the distance to TSS distribution plot as a PDF
ggsave(
    filename = "distance_to_tss_plot.pdf",
    plot     = p_tss,
    width    = 8,
    height   = 5
)


# ==============================================================================
# PART 2: DIFFERENTIAL BINDING ANALYSIS WITH DIFFBIND
# ==============================================================================
cat("\n--- Part 2: Differential Binding Framework using DiffBind ---\n")

# Structure of a standard DiffBind Sample Sheet Data Frame
# Note: In real workflows, paths point to actual .bam and .bed files
samples_df <- data.frame(
    SampleID    = c("Control_R1", "Control_R2", "Treated_R1", "Treated_R2"),
    Tissue      = "HCT116",
    Factor      = "CTCF",
    Condition   = c("Control", "Control", "Treated", "Treated"),
    Replicate   = c(1, 2, 1, 2),
    bamReads    = c("bams/ctrl_1.bam", "bams/ctrl_2.bam", "bams/treat_1.bam", "bams/treat_2.bam"),
    Peaks       = c("peaks/ctrl_1.bed", "peaks/ctrl_2.bed", "peaks/treat_1.bed", "peaks/treat_2.bed"),
    PeakCaller  = "macs"
)

cat("DiffBind Sample Sheet Architecture:\n")
print(samples_df[, c("SampleID", "Condition", "Replicate", "Factor")])

cat("\nStandard DiffBind Execution Template (Pseudocode for Real BAM/BED Data):\n")
cat("
  # 1. Initialize DBA Object
  # db_obj <- dba(sampleSheet = samples_df)

  # 2. Count Reads in Consensus Peaks
  # db_counts <- dba.count(db_obj, bParallel = TRUE)

  # 3. Establish Normalization Strategy
  # db_norm <- dba.normalize(db_counts, normalize = DBA_NORM_LIB)

  # 4. Set Contrasts & Execute Differential Analysis
  # db_contrast <- dba.contrast(db_norm, categories = DBA_CONDITION)
  # db_analyzed <- dba.analyze(db_contrast, method = DBA_DESEQ2)

  # 5. Extract Differentially Bound Peaks (GRanges)
  # db_results <- dba.report(db_analyzed, th = 0.05)
")
