# Module 08: Epigenomics & Peak Analysis with ChIPseeker & DiffBind
## Description: Peak annotation, TSS distribution profiling, and differential binding analysis framework for ChIP-seq / ATAC-seq data.

In epigenomics (ChIP-seq for transcription factors and histone modifications, or ATAC-seq for chromatin accessibility), peak callers like MACS2 produce genomic interval files (`.narrowPeak` or `.broadPeak`).

In R, epigenomics analysis centers around two primary goals:
1. **Functional Peak Annotation** (`ChIPseeker`): Mapping peak intervals to biological features (promoters, 5' UTRs, exons, introns, distal intergenic regions) and determining TSS (Transcription Start Site) distance distributions.

2. **Differential Affinity/Binding Analysis** (`DiffBind`): Quantifying read counts overlapping consensus peak sets across experimental conditions and conducting statistical differential testing (using `DESeq2` or `edgeR` backends).

3. **Peak Import & Genomic Range Construction:**
Imports BED/narrowPeak files into GRanges objects.

When we run a ChIP-seq or ATAC-seq experiment, the primary output after early processing (alignment and peak calling) is a list of genomic coordinates called "peaks." These peaks represent regions of the genome where a protein was bound (ChIP-seq) or where the chromatin is open/accessible (ATAC-seq).
However, raw genomic coordinates (e.g., `chr1: 1,234,567-1,235,890`) do not give you biological meaning by themselves. This framework allows bioinformaticians to answer three fundamental questions:

- **Peak Annotation**: Which genes are next to these coordinates? Are they falling inside an exon, an intron, or a promoter region? We do this to predict which genes are likely regulated by these DNA-binding proteins or open regions.
- **TSS Distribution Profiling**: Are these peaks localized around the Transcription Start Sites (TSS)? Transcription factors and transcriptional machinery heavily cluster right at the start of genes (promoters). Finding out if your data peaks near the TSS acts as a crucial biological quality control.
- **Differential Binding Analysis**: Does the protein bind more strongly or is the chromatin more open in Condition A vs. Condition B (e.g., Cancer vs. Normal cells)? This lets us pinpoint the exact genomic locations driving phenotypic changes or disease states.

1. **The Role of Each Package**
   - `ChIPseeker`: The core package used to map peak coordinates to genomic features and visualize how close peaks are to the TSS.
   - `TxDb.Hsapiens.UCSC.hg38.knownGene`: A specialized transcript database object containing the transcript coordinates for the human genome assembly hg38. It functions as the blueprint map for finding exons, introns, and promoters.
   - `org.Hs.eg.db`: The human genome wide annotation database. It maps abstract internal gene IDs (like Entrez IDs) to human-readable gene symbols (e.g., TP53).
   - `GenomicRanges`: The foundational Bioconductor infrastructure package used to store, manipulate, and compute overlaps on genomic intervals (chromosomes, starts, ends).
   - `DiffBind`: A statistical framework specialized in taking peak coordinates across multiple biological replicates and conditions to identify significantly changed, differentially bound regions.

## PART 1: Peak Annotation with** `ChIPseeker`
This section transforms abstract coordinates into concrete biological targets.
1. **Simulating Genomic Data** (`GRanges`)
```r
sim_peaks <- GRanges(
    seqnames = Rle(chroms, c(2, 2, 1)),
    ranges   = IRanges(start = starts, end = ends),
    strand   = Rle(c("*"))
)
```
- The code uses `GenomicRanges` to craft a mock dataset of 5 peaks spanning across Chromosomes 1, 2, and 3. `Rle` (Run-Length Encoding) is an optimization technique that stores repeating values efficiently. The strand is set to `*` because ChIP/ATAC-seq peak data usually doesn't have strand specificity (proteins bind DNA or open chromatin across both strands).
- Real-world equivalent: In a real pipeline, we won't build this manually. We will use `ChIPseeker::readPeakFile()` to import real `.bed` or `.narrowPeak` coordinate files produced by peak callers like MACS2.

2. **The Annotation Engine (`annotatePeak`)**
```r
peak_anno <- annotatePeak(
    peak         = sim_peaks,
    tssRegion    = c(-3000, 3000),
    TxDb         = txdb,
    annoDb       = "org.Hs.eg.db",
    verbose      = FALSE
)
```
- `tssRegion = c(-3000, 3000)`: Defines the boundaries of a Promoter region. Here, anything within 3,000 base pairs upstream or downstream of a gene's Transcription Start Site (TSS) is officially categorized as a promoter binding event.
- `TxDb` and `annoDb`: TxDb behaves like a map of coordinates (finding structural landmarks like exons, introns, and promoters). `annoDb` acts as a language translator, mapping machine IDs (Entrez ID) to human-readable gene identifiers (like `SYMBOL` or `GENENAME`).

3. **Output Tables and Plots**
```r
anno_df <- as.data.frame(peak_anno)
p_anno <- plotAnnoBar(peak_anno)
p_tss <- plotDistToTSS(peak_anno)
```

- `as.data.frame`: Converts the complex annotation object into a readable matrix showing exactly which genomic features (e.g., "Distal Intergenic", "Intron", "Promoter") each peak overlaps, alongside the closest gene symbol.
- `plotAnnoBar`: Produces a stacked bar chart illustrating percentage distributions. For instance, an ATAC-seq sample usually has peaks spread across promoters, introns, and intergenic regions, while a transcription factor ChIP-seq target often shows a massive skew towards the promoter region.
- `plotDistToTSS`: Generates a distribution summary illustrating the percentage of peaks falling within <1kb, 1-3kb, or 3-10kb from the nearest TSS.


## PART 2: Differential Binding Analysis with `DiffBind`
This part shifts from asking "where are my peaks?" to "how does my binding profile change across conditions?"

1. **The Sample Sheet Architecture**
```r
samples_df <- data.frame(
    SampleID    = c("Control_R1", "Control_R2", "Treated_R1", "Treated_R2"),
    Condition   = c("Control", "Control", "Treated", "Treated"),
    bamReads    = c("bams/ctrl_1.bam", ...),
    Peaks       = c("peaks/ctrl_1.bed", ...),
    ...
)
```
- `DiffBind` requires a strict tracking system called a sample sheet. This tracks every replicate, identifying which aligned sequencing reads file (`.bam`) pairs with which specific peak layout file (`.bed`).

2. **The Step-by-Step Pseudocode (The Real Matrix Execution)**
Since `DiffBind` requires physical files on our hard drive to run, our code blocks out the framework as a template. Here is how that framework functions mathematically behind the scenes:
- Step 1: `dba(sampleSheet = ...)`
Loads our tracking file. It evaluates our `.bed` sheets across all samples to construct a unified consensus peakset (a master list of all regions where binding was detected in at least a subset of samples).
- Step 2: `dba.count(...)`
The heavy lifting step. It opens up our large sequence alignment files (`.bam`) and explicitly counts the exact number of sequence reads that fall inside the boundaries of each consensus peak. It yields a raw count matrix (Rows = Consensus Peaks, Columns = Biological Samples).
- Step 3: `dba.normalize(...)`
Different sequencing runs generate different amounts of total data. This stage scales our raw counts against calculated library sizes to ensure we aren't mistaking a global sequencing boost for genuine biological variations.
- Step 4: `dba.contrast(...)` & `dba.analyze(...)`
Organizes our cohorts (e.g., Control vs Treated). It then hands the normalized count matrix over to reliable RNA-seq differential engines (like DESeq2 or EdgeR). It applies generalized linear models to test which peaks exhibit statistical differences between our groups.
- Step 5: `dba.report(..., th = 0.05)`
Filters the analysis matrix down to peaks possessing a False Discovery Rate (FDR) below 0.05. It outputs a final `GRanges` object of statistically significant, differentially bound genomic segments.

## Output
**Annotated Peak Sample Output:**

```bash
seqnames    start      end
1     chr1  1000000  1000760
2     chr1  2500000  2500520
3     chr2  5000000  5000352
4     chr2 12000000 12000273
5     chr3 18000000 18000427
                                           annotation    SYMBOL distanceToTSS
1                                    Promoter (<=1kb)      HES4             0
2                                    Promoter (1-2kb)     PLCH2          1353
3                                   Distal Intergenic LINC01249       -343502
4 Intron (ENST00000438292.5/100506457, intron 1 of 4) MIR3681HG         -6840
5    Intron (ENST00000624232.2/339862, intron 2 of 4)     BALR6        139078
```

## Explanation of the outputs
### Plot 1: Genomic Distribution of Peaks
This plot breaks down the structural regions of the genome where your peaks landed. Because our dataset contains exactly 5 peaks, each peak represents exactly 20% of the total distribution:
- **Promoter** ($\(\le \) 1kb)$ [20%]: Corresponds to Peak 1 (`HES4`), which sits precisely at the transcription start site (`distanceToTSS = 0`).
- **Promoter (1-2kb) [20%]**: Corresponds to Peak 2 (`PLCH2`), which is located `1,353 bp` away from the TSS.
- **1st Intron [20%]**: Corresponds to Peak 4 (`MIR3681HG`), which explicitly maps to intron 1 of 4.
- **Other Intron [20%]**: Corresponds to Peak 5 (`BALR6`), which maps to intron 2 of 4.
- **Distal Intergenic [20%]**: Corresponds to Peak 3 (`LINC01249`), which lands far outside any gene body in an intergenic desert (-343,502 bp away).


### Plot 2: Peak Distribution Relative to TSS
This plot maps out the spatial distance and direction (5' upstream vs. 3' downstream) of our peaks relative to the closest Transcription Start Site (TSS). The vertical line in the center represents the TSS (0 bp).
- **Center Green Band (0 bp - 1 kb) [20% total width]**: Occupies the center area from -10% to +10%. This is driven by Peak 1 (distanceToTSS = 0).
- **Right Teal Band (1 - 3 kb) [20% width]**: Located on the right side from 10% to 30%. This represents downstream binding, driven by Peak 2 (+1,353 bp).
- **Left Brown Band (5 - 10 kb) [20% width]**: Located on the left side from -10% to -30%. This represents upstream binding, driven by Peak 4 (-6,840 bp).
- **Outer Purple Bands (>100 kb) [40% total width]**: Distributed symmetrically at the flanks (-30% to -50% on the left, and 30% to 50% on the right). This represents our two extreme distal peaks: Peak 3 (-343,502 bp, upstream) and Peak 5 (+139,078 bp, downstream).


**DiffBind Sample Sheet Architecture:**
```bash
SampleID       Condition    Replicate Factor
1 Control_R1   Control         1   CTCF
2 Control_R2   Control         2   CTCF
3 Treated_R1   Treated         1   CTCF
4 Treated_R2   Treated         2   CTCF
```

**Biological Insight for our Pipeline (CTCF Factor)**
Our sample sheet indicates we are profiling CTCF (CCCTC-Binding Factor). CTCF is a well-known architectural protein responsible for mapping the 3D structure of the genome and organizing chromatin loops.
While many transcription factors cluster exclusively at the green promoter band (0-1kb), real CTCF datasets typically mirror this balanced distribution. We will see distinct peaks at promoters to regulate transcription, mixed with a massive abundance of peaks in introns and distal intergenic regions where CTCF functions as an insulator to loop DNA together.






















