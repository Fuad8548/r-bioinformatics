# Module 04: NGS Data Handling & Alignment Parsing
## Description: Selective streaming of BAM alignments using Rsamtools, CIGAR parsing with GenomicAlignments, and feature counting.

After sequencing DNA or RNA, aligners (like STAR, Bowtie2, or BWA) map millions of short reads back to a reference genome, producing **SAM** (Sequence Alignment/Map) files or their compressed binary equivalents, **BAM** files. In R, reading a 15 GB BAM file directly into memory will crash our session. Bioconductor solves this through random-access streaming and strand/CIGAR-aware S4 containers.

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

## Core Scripts

### Loading required packages:
```r
suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicAlignments)
  library(GenomicRanges)
  library(SummarizedExperiment)
})
```
**Explanations:**
|       Package        |                   Core Purpose                    |                                                                                          Function                                                                                           |
| :------------------: | :-----------------------------------------------: | :-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
|      Rsamtools       |         High-performance BAM file access          | Allows you to selectively stream and filter massive BAM files (e.g., loading only specific chromosomes or filtering out low-quality alignments) without overloading your computer's memory. |
|  GenomicAlignments   |     Reading and parsing aligned genomic data      |               Holds the functions needed to load reads into R (like readGAlignments) and parse CIGAR strings to accurately map where reads split, splice, or contain indels.                |
|    GenomicRanges     |  Representing and manipulating genomic intervals  |        Provides the foundational infrastructure (like GRanges objects) to define the boundaries of genomic features (genes, exons, or transcripts) and perform overlap calculations.        |
| SummarizedExperiment | Storing final feature-count matrices and metadata |  Acts as the ultimate container for your count data. It cleanly links your final matrix of read counts to both the feature descriptions (rows) and the sample/patient metadata (columns).   |


## 1. Connecting to BAM Files & Inspecting Headers
We are going to establish a connection to a BAM file and extracting its metadata without loading the heavy sequence alignment data into our computer's RAM.

1. **Locating the File:**
```r
bam_path <- system.file("extdata", "ex1.bam", package = "Rsamtools")
```
- It finds the absolute file path to a tiny, sample BAM file (`ex1.bam`) that comes built-in with the `Rsamtools` package for testing purposes.

1. **Creating a Memory-Efficient Pointer**
```r
bam_file <- BamFile(bam_path)
open(bam_file)
```
   - Instead of reading the entire file into memory (which could easily crash our R session with real-world human genome data), `BamFile()` creates an S4 pointer object.
   - `open(bam_file)` establishes a streaming connection to the file on our disk. This allows R to jump directly to specific chromosomes or coordinates later on.

1. **Inspecting the Metadata (The Header)**
```r
header <- scanBamHeader(bam_file)
close(bam_file)
```
   - Every well-formed BAM file has a "Header" section at the very beginning. This header acts like a table of contents. `scanBamHeader()` reads this section and extracts `$targets`.

### Output:
```bash
seq1 seq2 
1575 1584
```
This tells us that our reference genome for this experiment consists of exactly two reference contigs/chromosomes:
   - `seq1`: has a total length of 1,575 base pairs.
   - `seq2`: has a total length of 1,584 base pairs.

Any alignment position inside this BAM file will fall within these coordinates. Once we have this information, close(bam_file) safely disconnects from the file.


## 2. Targeted Import with `ScanBamParam`
In this step, we used `ScanBamParam` as a powerful camera lens to zoom in and import only a specific subset of sequencing reads that match our exact biological criteria.

1. **Breaking Down the Filter (`ScanBamParam`)**
Instead of swallowing the whole file, we instructed R to stream and extract reads matching four strict conditions:
   - `which = target_region:` Restricts the search area. We asked for `seq1:1-5000`. Since `seq1` is only 1,575 base pairs long (from Part 1), this captures every single read on `seq1`.
   - `flag = scanBamFlag(...)`: Quality control. It strips away noisy sequencing artifacts by throwing out unmapped reads and duplicates (identical reads created artificially during PCR amplification).
   - `mapqFilter = 20`: Confidence control. A mapping quality (MAPQ) score of 20 means there is a 99% probability that the read actually belongs to this genomic position. Low-confidence multi-mapping reads are ignored.
   - `what = c(...)`: Variable selection. We specified exactly which fields to load, which keeps memory footprint minimal.


## Output 1: The `GAlignments`

`readGAlignments()` returns a highly specialized Genomic Data Frame holding **1,476 matching reads**. The output table is split into two halves divided by a pipe character (`|`):
```bash
print(alignments)
seqnames strand    cigar    qwidth     start       end     width
seq1      +         36M        36         1        36        36
seq1      +         35M        35         3        37        35
seq1      +         35M        35         5        39        35
seq1      +         36M        36         6        41        36
seq1      +         35M        35         9        43        35

njunc |    rname   strand       pos       cigar      mapq
    0 |     seq1        +         1         36M        99
    0 |     seq1        +         3         35M        99
    0 |     seq1        +         5         35M        99
    0 |     seq1        +         6         36M        63
    0 |     seq1        +         9         35M        99
```


### Section 1: Genomic Coordinates & Core Alignment
The top block shows you where and how the reads mapped to the genome.
- `seqnames`: The chromosome or scaffold name. Here, your reads all mapped to seq1.
• `strand`: The DNA strand. + means the forward strand, and - means the reverse strand.
• `cigar`: The CIGAR string (Concise Idempotent Alignment Relation). It describes what the aligner did to match the read to the genome. For example, 36M means there are 36 consecutive Alignment Matches (or mismatches) between the read and the genome.
• `qwidth`: The actual length of the raw sequencing read itself (Query Width). For the first read, it's 36 base pairs long.
• `start` / `end`: The exact starting and ending coordinate positions on the reference genome (seq1) where this read maps. Notice how the first read covers positions 1 to 36.
• `width`: The total span of the reference genome covered by the read (end - start + 1).

### Section 2: Metadata & Junctions
The bottom block displays the specific data fields you requested using your ScanBamParam(what = ...) instruction, separated by a | line.
- `njunc`: The number of junctions (introns) crossed by the read. All values here are 0, meaning these reads are continuous and do not split across gaps/introns (common in DNA-seq or single-exon transcripts).
- `rname`: The reference sequence name (matches seqnames).
• pos: The leftmost alignment position on the reference (matches start).
• `mapq`: The Mapping Quality score. This tells you how confident the aligner is that the read belongs exactly there. A score of 99 is typically used by aligners to mean "highly confident / uniquely mapped," while 63 is also a very strong score (well above your filter of 20!).

**Summary:**
Your very first read is 36 base pairs long (`qwidth`), maps flawlessly to the forward strand (`+`) of `seq1` across positions 1 to 36 with 36 matches (`36M`), and has a perfect mapping confidence score of 99 (`mapq`).


## 3. Handling Spliced Reads & Exon Coordinates
In this step, we are executing a critical transition step in NGS analysis: converting read alignments into pure genomic coordinates (`GRanges`).

1. **What does `granges()` actually do?**
A `GAlignments` object tracks how a sequencer read is aligned to a reference (including CIGAR strings, mapping qualities, and flags). However, downstream tools—like peak callers, variant callers, or expression counters—only care about where the read is physically located on the genome.
`granges(alignments)` strips away the raw alignment metrics and simplifies the data down into standard genomic intervals (`seqnames`, `ranges`, and `strand`).

**Output:**
```bash
granges(alignments)

seqnames     ranges  strand
  <Rle>     <IRanges>  <Rle>
  seq1        1-36      +
  seq1        3-37      +
  seq1        5-39      +
```

1. **The Power of `granges()` with "Spliced Reads"**
While the first three reads in our output are simple, continuous blocks (`1-36`, `3-37`, `5-39`), `granges()` does something incredibly intelligent under the hood when it encounters spliced reads (common in RNA-Seq data):
   - The CIGAR "N" Operation: If a read spans an intron, its CIGAR string might look like `10M100N25M` (10 base pairs match an exon, 100 base pairs are skipped as an intron, and 25 base pairs match the next exon).
   - The `granges()` Magic: When we run `granges()` on a spliced alignment, it automatically accounts for that `N` operation. It calculates the outermost genomic footprint (`start` to `end`), ensuring that the giant intron gap is recognized so that downstream overlap calculations with gene models remain perfectly accurate.



# 4. Read Quantification over Features (`summarizeOverlaps`)
In this step, we performed **RNA-Seq read quantification**. We took our genomic sequence alignments (the BAM file) and mapped them against an annotation map (the exons) to count how many fragments belong to each gene.

1. **Defining the Genomic Annotation Structure ( `GRangesList`)**
- In real-world projects, this structure is usually loaded automatically from a **GTF/GFF annotation file** using packages like `GenomicFeatures`.
- It uses a `GRangesList` (a list of genomic ranges). This structure mimics real biology because genes are not solid blocks; they are split into exons.
  - `Gene_A` has two separate exons: one from base pair 100 to 500, and a second one from 800 to 1200.

2. **Quantifying Overlaps (`summarizeOverlaps`)**
This is the machine that matches our reads to our exons.
   - `mode = "Union"`: This defines the logical rule for counting. Under `Union`, if a sequencing read overlaps any part of any exon belonging to `Gene_A`, it scores a point for `Gene_A`. If a read splits across an intron but hits both exons, it still safely counts as exactly 1 hit for that gene.

### Output explanation
1. **The `RangedSummarizedExperiment` Container**
```bash
print(se_counts)

class: RangedSummarizedExperiment 
dim: 2 1 
assays(1): counts
rownames(2): Gene_A Gene_B
colnames(1): ex1.bam
```
Bioconductor packs the results into a `SummarizedExperiment` object. This is a clever "all-in-one" container widely used in downstream differential expression packages like **DESeq2** or **edgeR**. It acts like a 3D spreadsheet holding:
   - `dim: 2 1`: A matrix with 2 rows (our genes) and 1 column (our sample BAM file).
   - `assays(1): counts`: The actual mathematical matrix containing the raw expression numbers.
   - `rownames` & `colnames`: Labels mapping exactly to our features and our sample data files.


2. **The Counts Matrix**
```bash
print(assays(se_counts)$counts)

        ex1.bam
Gene_A     825
Gene_B      82
```

This is our final biological data payload:
```r
Gene_A = Exon 1 (100 to 500) and Exon 2 (800 to 1200)
Gene_B = Exon 1 (1500 to 1900) and Exon 2 (2200 to 2800)
```

- Gene B is looking for reads between positions 1500 to 2800. But our chromosome seq1 ends completely at position 1575!
- Gene A, however, sits perfectly in the middle of the action (positions 100–1200), right where the bulk of your sequencing reads were mapped. That is why Gene A racked up 825 overlaps, while Gene B only managed to catch 82 overlaps at the very edge of the chromosome.

**What this tells you biologically**
Gene A has roughly 10 times more reads mapping to it than Gene B. In a real experiment, this heavily implies that Gene A is much more actively expressed (turned on) in this sample than Gene B!







