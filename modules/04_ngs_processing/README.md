# Module 04: NGS Data Handling & Alignment Parsing
## Description: Selective streaming of BAM alignments using Rsamtools, CIGAR parsing with GenomicAlignments, and feature counting.

After sequencing DNA or RNA, aligners (like STAR, Bowtie2, or BWA) map millions of short reads back to a reference genome, producing **SAM** (Sequence Alignment/Map) files or their compressed binary equivalents, **BAM** files.

In R, reading a 15 GB BAM file directly into memory will crash our session. Bioconductor solves this through random-access streaming and strand/CIGAR-aware S4 containers.

**Key Concepts to Internalize**
1. **SAM/BAM Files & Indexing** (`.bai`)
- **SAM**: Human-readable text format containing read ID, flag, chromosome, position, mapping quality (MAPQ), CIGAR string, and sequence.
- **BAM**: Binary, BGZF-compressed version of SAM.
- **BAI Index**: A spatial lookup table for BAM files. By pairing a `.bam` file with its `.bai` index, `Rsamtools` can perform random access—fetching reads from `chr1:1,000,000-1,050,000` without reading the rest of the file.

2. **Filtering on Import with** `ScanBamParam`
Instead of importing unmapped reads, low-quality reads, or PCR duplicates into R, we define a filter rule (`ScanBamParam`) at read-time:
- Filter by **Mapping Quality (MAPQ)**: E.g., `minq = 30` drops ambiguous multi-mapping reads.
- Filter by **SAM Flags**: E.g., `scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE)`.
- Filter by **Genomic Region**: Supply a `GRanges` target to stream only reads overlapping specific loci.

3. **CIGAR Strings & `GAlignments`**
RNA-Seq reads spanning introns contain split alignments represented by CIGAR strings (e.g., `50M1000N50M` = 50 matched bases, 1000 skipped bases [intron], 50 matched bases).
- `GenomicAlignments::readGAlignments()` parses CIGAR strings automatically, turning spliced alignments into intron-aware ranges.

4. **Feature Quantification with `summarizeOverlaps()`**
The function `summarizeOverlaps()` counts how many reads in a BAM file overlap each gene/exon in a `GRangesList`, returning a ready-to-use `SummarizedExperiment` for downstream differential expression (`DESeq2`).

**Core Script:**

### Loading required packages:
```r
suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicAlignments)
  library(GenomicRanges)
  library(SummarizedExperiment)
})
```


## 1. Connecting to BAM Files & Inspecting Headers

```r
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
```

### Explanations:
We are going to establish a connection to a BAM file and extracting its metadata without loading the heavy sequence alignment data into our computer's RAM.

1. **Locating the File:**
```r
bam_path <- system.file("extdata", "ex1.bam", package = "Rsamtools")
```
What it does: It finds the absolute file path to a tiny, sample BAM file (`ex1.bam`) that comes built-in with the `Rsamtools` package for testing purposes.

2. **Creating a Memory-Efficient Pointer**
```r
bam_file <- BamFile(bam_path)
open(bam_file)
```
   - What it does: Instead of reading the entire file into memory (which could easily crash our R session with real-world human genome data), `BamFile()` creates an S4 pointer object.
   - `open(bam_file)` establishes a streaming connection to the file on our disk. This allows R to jump directly to specific chromosomes or coordinates later on.

3. **Inspecting the Metadata (The Header)**
```r
header <- scanBamHeader(bam_file)
print(header$targets)
close(bam_file)
```
   - What it does: Every well-formed BAM file has a "Header" section at the very beginning. This header acts like a table of contents. `scanBamHeader()` reads this section and extracts `$targets`.

### Output:
```bash
seq1 seq2 
1575 1584
```
This tells we that our reference genome for this experiment consists of exactly two reference contigs/chromosomes:
   - `seq1`: has a total length of 1,575 base pairs.
   - `seq2`: has a total length of 1,584 base pairs.

Any alignment position inside this BAM file will fall within these coordinates. Once we have this information, close(bam_file) safely disconnects from the file.


# 2. Targeted Import with `ScanBamParam`

```r
# Define a target region of interest (chr1:1 to 5000)
target_region <- GRanges("seq1", IRanges(start = 1, end = 5000))

# Define filtering parameters:
# - High mapping quality (MAPQ >= 20)
# - Exclude unmapped reads and PCR duplicates
# - Restrict to target_region
# - Keep specific fields (cigar, mapq, seq)
param <- ScanBamParam(
  which   = target_region,
  what    = c("rname", "strand", "pos", "cigar", "mapq"),
  flag    = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE),
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
```

### Explanations:
In this step, we used `ScanBamParam` as a powerful camera lens to zoom in and import only a specific subset of sequencing reads that match our exact biological criteria.

1. **Breaking Down the Filter (`ScanBamParam`)**
Instead of swallowing the whole file, we instructed R to stream and extract reads matching four strict conditions:
   - `which = target_region:` Restricts the search area. We asked for `seq1:1-5000`. Since `seq1` is only 1,575 base pairs long (from Part 1), this captures every single read on `seq1`.
   - `flag = scanBamFlag(...)`: Quality control. It strips away noisy sequencing artifacts by throwing out unmapped reads and duplicates (identical reads created artificially during PCR amplification).
   - `mapqFilter = 20`: Confidence control. A mapping quality (MAPQ) score of 20 means there is a 99% probability that the read actually belongs to this genomic position. Low-confidence multi-mapping reads are ignored.
   - `what = c(...)`: Variable selection. We specified exactly which fields to load, which keeps memory footprint minimal.

### Output 1: The `GAlignments`
`readGAlignments()` returns a highly specialized Genomic Data Frame holding **1,476 matching reads**. The output table is split into two halves divided by a pipe character (`|`):
```bash
seqnames strand cigar qwidth start end width |  njunc | ...
[1] seq1    +    36M    36     1    36   36  |      0 | ...
```

**Left Side (The Global Coordinates)**
These are automatically calculated genomic boundaries managed by the `GenomicAlignments` package:
- `seqnames` / `strand`: Confirms that this read mapped perfectly to the forward (`+`) or reverse (`-`) strand of `seq1`.\
- `start` / `end`: The absolute, 1-based coordinates on the reference genome. Read `[1]` blankets base pairs 1 through 36. Read `[1476]` blankets base pairs 1,535 through 1,569.
- `qwidth` / `width`: `qwidth` is the length of the raw sequence query read itself. `width` is how much physical space it spans along the reference chromosome.

**Right Side (The Metadata Columns)**
Everything to the right of the `|` represents the raw raw data we specifically requested via the `what` argument.
  - `njunc`: Number of junctions (gaps or introns spanned by the read). They are all `0`, indicating continuous, unbroken genomic DNA or fully mapped exonic fragments.
  - `mapq`: The raw MAPQ scores. We can see values like `99` (extremely confident) and `63`, which easily cleared our minimum threshold of `20`.

### Output 2 & 3: CIGAR & Genomic Widths
```bash
Output 2: First 5 CIGAR Strings:
[1] "36M" "35M" "35M" "36M" "35M"
```
**CIGAR** stands for *Concise Idiosyncratic Gapped Alignment Report*. It uses a letter-coded short-hand to tell the story of how the read aligned:
  - `36M`: Means 36 Alignment Matches. The read matches the reference sequence consecutively without any gaps, insertions, or deletions.
  - `35M`: 35 continuous alignment matches.

```bash
Output 3: Widths of Aligned Reads on Genome:
[1] 36 35 35 36 35
```
This lists the physical genomic footprint of the first 5 reads. Because all 5 reads have a pure matching CIGAR string (`M` only) with zero insertions/deletions, the physical width on the genome exactly matches the lengths of the CIGAR strings.



# 3. Handling Spliced Reads & Exon Coordinates

```r
# Convert alignments into GRanges (expands CIGAR N operations into skipped ranges)
read_ranges <- granges(alignments)

cat("\n--- Converted Genomic Ranges from Alignments ---\n")
print(head(read_ranges, 3))
```

### Explanations:
In this step, we are executing a critical transition step in NGS analysis: converting read alignments into pure genomic coordinates (`GRanges`).

1. **What does `granges()` actually do?**
A `GAlignments` object tracks how a sequencer read is aligned to a reference (including CIGAR strings, mapping qualities, and flags). However, downstream tools—like peak callers, variant callers, or expression counters—only care about where the read is physically located on the genome.
`granges(alignments)` strips away the raw alignment metrics and simplifies the data down into standard genomic intervals (`seqnames`, `ranges`, and `strand`).

2. **The Power of `granges()` with "Spliced Reads"**
While the first three reads in our output are simple, continuous blocks (`1-36`, `3-37`, `5-39`), `granges()` does something incredibly intelligent under the hood when it encounters spliced reads (common in RNA-Seq data):
   - The CIGAR "N" Operation: If a read spans an intron, its CIGAR string might look like `10M100N25M` (10 base pairs match an exon, 100 base pairs are skipped as an intron, and 25 base pairs match the next exon).
   - The `granges()` Magic: When we run `granges()` on a spliced alignment, it automatically accounts for that `N` operation. It calculates the outermost genomic footprint (`start` to `end`), ensuring that the giant intron gap is recognized so that downstream overlap calculations with gene models remain perfectly accurate.



# 4. Read Quantification over Features (`summarizeOverlaps`)

```r
library(BiocParallel) # Load the Bioconductor parallel backend manager

# Create synthetic gene features (exons grouped by gene in a GRangesList)
gene_exons <- GRangesList(
  Gene_A = GRanges("seq1", IRanges(start = c(100, 800), end = c(500, 1200))),
  Gene_B = GRanges("seq1", IRanges(start = c(1500, 2200), end = c(1900, 2800)))
)

# Quantify read overlaps using 'Union' mode
# Union: A read counts toward a gene if it overlaps any exon of that gene
se_counts <- summarizeOverlaps(
  features = gene_exons,
  reads    = bam_path,
  mode     = "Union",
  singleEnd = TRUE,
  param    = param
)

cat("\n--- Resulting SummarizedExperiment from Count Quantification ---\n")
print(se_counts)

cat("\nCounts Matrix:\n")
print(assays(se_counts)$counts)
```

### Explanations
In this step, we performed **RNA-Seq read quantification**. We took our genomic sequence alignments (the BAM file) and mapped them against an annotation map (the exons) to count how many fragments belong to each gene.

1. **Defining the Genomic Annotation Structure ( `GRangesList`)**
- **What it is**: In real-world projects, this structure is usually loaded automatically from a **GTF/GFF annotation file** using packages like `GenomicFeatures`.
- **How it works**: It uses a `GRangesList` (a list of genomic ranges). This structure mimics real biology because genes are not solid blocks; they are split into exons.
  - `Gene_A` has two separate exons: one from base pair 100 to 500, and a second one from 800 to 1200.

2. **Quantifying Overlaps (`summarizeOverlaps`)**
This is the machine that matches our reads to our exons.
   - `mode = "Union"`: This defines the logical rule for counting. Under `Union`, if a sequencing read overlaps any part of any exon belonging to `Gene_A`, it scores a point for `Gene_A`. If a read splits across an intron but hits both exons, it still safely counts as exactly 1 hit for that gene.

### Output 
1. **The `RangedSummarizedExperiment` Container**
```bash
class: RangedSummarizedExperiment 
dim: 2 1 
assays(1): counts
rownames(2): Gene_A Gene_B
colnames(1): ex1.bam
```
Bioconductor packs the results into a `SummarizedExperiment` object. 
This is a clever "all-in-one" container widely used in downstream differential expression packages like **DESeq2** or **edgeR**. It acts like a 3D spreadsheet holding:
   - `dim: 2 1`: A matrix with 2 rows (our genes) and 1 column (our sample BAM file).
   - `assays(1): counts`: The actual mathematical matrix containing the raw expression numbers.
   - `rownames` & `colnames`: Labels mapping exactly to our features and our sample data files.

2. **The Counts Matrix Explained**
```bash
        ex1.bam
Gene_A     825
Gene_B      82
```

This is our final biological data payload:
- `Gene_A (825)`: Out of the 1,476 total filtered reads on `seq1` (from Part 2), 825 reads physically overlapped the exons of Gene_A. This indicates robust gene expression.
- `Gene_B (82)`: Even though Gene_B's coordinates (`1500` to `2800`) stretch way past the absolute end of `seq1` (`1575` bp), it managed to catch 82 reads. This means there is a pile of reads mapping right at the trailing tail edge of `seq1` (between 1500 and 1575) before the chromosome ends!









