# Module 03: Genomic Ranges & Annotations
## Description: Managing genomic intervals, strand-aware operations, range transformations, and overlapping features using GenomicRanges.

## Conceptual Blueprint: Why Genomic Ranges?

Genomic data is fundamentally spatial. Genes, exons, transcription factor binding sites (ChIP-seq peaks), and variants (SNPs) are all coordinates anchored to specific positions on chromosomes. 
In a standard R data frame or Python dictionary, checking if a 200 bp ChIP-seq peak overlaps with any of 30,000 promoter regions requires nested loops, taking $O(N \times M)$ time. Bioconductor’s `GenomicRanges` infrastructure uses optimized C-level interval trees (Interval Trees and NCList data structures) to reduce overlap searches across millions of features to $O(N \log M)$ milliseconds.

## Key Concepts to Internalize
1. **Coordinate Systems: 1-Based & Closed**
Unlike Python/BED files (0-based, half-open $[start, end)$), R and Bioconductor use 1-based, closed intervals $[start, end]$.
- An interval with `start = 100` and `end = 200` has a width of $200 - 100 + 1 = 101\text{ bp}$.

2. **IRanges vs. GRanges**
- **IRanges**: Pure integer ranges `[start, end]` without any genomic context (no chromosome, no strand).
- `GRanges`: Wraps an `IRanges` object with genomic awareness:

  - `seqnames`: Chromosome / Contig identifier (e.g., "`chr1`", "`chrX`").
  - `ranges`: An underlying `IRanges` object.
  - `strand`: Strand orientation (`"+"`, `"-"`, or `"*"` for unstranded/both).
  - `mcols`: Metadata columns (S4 DataFrame holding gene IDs, scores, annotations).

3. **Strand-Aware Transformations**
In genomic analysis, "upstream" depends entirely on the DNA strand:
   - On the `+` strand, the promoter is to the left of `start` ($start - width$).
   - On the `-` strand, transcription moves right-to-left, so the promoter is to the right of `end` ($end + width$).
`GRanges` methods automatically respect strand orientation when performing operations like `flank()`.

## loading packages
```r
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})
```

1. Buiding GRanges Objects & Accessing Slots
```r
# Construct a GRanges object representing synthetic genes
genes_gr <- GRanges(
  seqnames = Rle(c("chr1", "chr1", "chr2")),
  ranges   = IRanges(
    start = c(1000, 5000, 2000),
    end   = c(3000, 8000, 4000)
  ),
  strand   = c("+", "-", "+"),
  gene_id  = c("GENE_A", "GENE_B", "GENE_C"),
  score    = c(85.2, 92.0, 45.1)
)

cat("--- GRanges Object Structure ---\n")
print(genes_gr)
```

**Explanations:**
Here is what each argument inside `GRanges(...)` means:
- `seqnames = Rle(c("chr1", "chr1", "chr2"))`
  - **What it means**: This specifies which chromosome each gene is on.
  - **The Data**: Gene A is on Chromosome 1, Gene B is on Chromosome 1, and Gene C is on Chromosome 2.
  - Note: `Rle` stands for "Run-Length Encoding." It is just a highly efficient way for R to store repeated data (like "chr1" appearing multiple times) to save computer memory.

- `ranges = IRanges(start = ..., end = ...)`
  - **What it means**: This defines the exact boundaries (coordinates) of the genes on those chromosomes.
  - **The Data**:
    - Gene A spans from base pair 1,000 to 3,000.
    - Gene B spans from base pair 5,000 to 8,000.
    - Gene C spans from base pair 2,000 to 4,000.
  - strand = c("+", "-", "+")
    - **What it means**: DNA has two strands: a forward/plus (`+`) strand and a reverse/minus (`-`) strand. This tells R which direction the gene is facing.
  - `gene_id = ...` and `score = ...`
    - **What it means**: These are metadata columns. We can attach any extra information we want to our genomic ranges. Here, the code attaches custom names (`gene_id`) and a confidence or expression value (`score`) to each gene.

## What the Output Looks Like
When print(`genes_gr`) runs, R will print out a neat table that looks something like this:

```bash
GRanges object with 3 ranges and 2 metadata columns:
      seqnames    ranges strand |     gene_id     score
         <Rle> <IRanges>  <Rle> | <character> <numeric>
  [1]     chr1 1000-3000      + |      GENE_A      85.2
  [2]     chr1 5000-8000      - |      GENE_B      92.0
  [3]     chr2 2000-4000      + |      GENE_C      45.1
```

Notice the vertical bar (`|`). Everything to the left of the bar (`seqnames`, `ranges`, `strand`) is the core genomic info required for every `GRanges` object. Everything to the right (`gene_id`, `score`) is our custom metadata.

```r
# Accessing S4 Slots via Getter Functions (NEVER use @)
cat("\nExtracting Chromosomes (seqnames):\n")
print(seqnames(genes_gr))

cat("\nExtracting Coordinates (IRanges):\n")
print(ranges(genes_gr))

cat("\nExtracting Metadata DataFrame (mcols):\n")
print(mcols(genes_gr))
```

**Explanations:**
1. `seqnames(genes_gr)`
   - **What it does**: Extracts only the chromosome names from your data.
   - **What the output means**: It ignores the coordinates, strands, and metadata, giving us back just a list of the chromosomes: `chr1`, `chr1`, `chr2`.

Output: 
```bash
Lengths :       2    1
Values :     chr1 chr2
Levels(2) :  chr1 chr2
```

2. `ranges(genes_gr)`
   - **What it does**: Extracts only the start and end positions (the `IRanges` component), along with the width (length) of each region.
   - What the output means: R automatically calculates the width for you (End - Start + 1). You will see a clean summary of the boundaries:

Output:
```bash
          start       end     width
      <integer> <integer> <integer>
  [1]      1000      3000      2001
  [2]      5000      8000      3001
  [3]      2000      4000      2001
```

3. `mcols(genes_gr)`
   - **What it does**: Extracts only the custom metadata columns (everything to the right of the `|` vertical bar). `mcols` stands for "metadata columns".
   - **What the output means**: It strips away the genomic coordinates and leaves you with a standard data table containing just your custom variables: `gene_id` and `score`.

Output:
```bash
      gene_id     score
  <character> <numeric>
1      GENE_A      85.2
2      GENE_B      92.0
3      GENE_C      45.1
```

## 2. Strand-Aware Range Transformations

This section is all about finding the promoters of your genes. In biology, a promoter is a region of DNA located just before the start of a gene where transcription begins (the Transcription Start Site or TSS).
The term "Strand-Aware" is the most important concept here. Because DNA has two strands running in opposite directions, "before the gene" means different things depending on the strand:
- On the plus (`+`) strand, the gene goes left-to-right. The promoter sits to the left (lower coordinates).
- On the minus (`-`) strand, the gene goes right-to-left. The promoter sits to the right (higher coordinates).

### 2.1 Extract Promoter Regions using flank()

```r
# flank(x, width) gets adjacent regions. start = TRUE gets upstream region.
promoters_gr <- flank(genes_gr, width = 1000, start = TRUE)

cat("\n--- Promoters (1kb Upstream of TSS) ---\n")
# Notice how for GENE_B (minus strand), flank correctly extracts downstream coordinates
print(promoters_gr[, "gene_id"])
```

- The `flank()` function grabs a region adjacent to our existing range. By setting `start = TRUE`, we tell R to look upstream (before the TSS).
- **GENE_A** (`+` **strand, original: 1000 to 3000**): The TSS is at `1000`. Going 1,000 base pairs upstream (to the left) yields coordinates **0 to 999** (or 1 to 999 depending on formatting).
- **GENE_B** (`-` **strand, original: 5000 to 8000**): Because it is on the minus strand, the gene actually starts at `8000` and goes backward. Going 1,000 base pairs upstream (to the right) yields coordinates **8001 to 9000**.

**Output**:
```bash
      seqnames    ranges strand |     gene_id
         <Rle> <IRanges>  <Rle> | <character>
  [1]     chr1     0-999      + |      GENE_A
  [2]     chr1 8001-9000      - |      GENE_B
  [3]     chr2 1000-1999      + |      GENE_C
```


### 2.2 Intra-Range Operations (shift & promoters)

```r
# promoters() directly gets [TSS - upstream, TSS + downstream]
promoters_explicit <- promoters(genes_gr, upstream = 2000, downstream = 200)

cat("\n--- Promoters (-2000bp to +200bp around TSS) ---\n")
print(promoters_explicit[, "gene_id"])
```

**Finding Promoters using `promoters()`**

While `flank()` only looks strictly outside the gene boundaries, the built-in `promoters()` function allows you to create a window that captures data **both before and slightly inside the gene**.

```r
promoters_explicit <- promoters(genes_gr, upstream = 2000, downstream = 200)
```

This command tells R to capture **2000 base pairs upstream** (before the TSS) and **200 base pairs downstream** (inside the gene) to create a custom regulatory window around the TSS.
How it calculates the new coordinates:

- **GENE_A (`+` strand, TSS = 1000)**:
  - Go 2,000 bp left (upstream): `1000 - 2000 = -1000` (R defaults boundaries to 1 if they drop below 1).
  - Go 200 bp right (downstream inside the gene): `1000 + 200 = 1200`.
  - Final window: **1 to 1200**.
- **GENE_B (`-` strand, TSS = 8000):**
  - Go 2,000 bp right (upstream): `8000 + 2000 = 10000`.
  - Go 200 bp left (downstream inside the gene): `8000 - 200 = 7800`.
  - Final window: 7801 to 10000.

**Output:**
```bash
      seqnames     ranges strand |     gene_id
         <Rle>  <IRanges>  <Rle> | <character>
  [1]     chr1 -1000-1199      + |      GENE_A
  [2]     chr1 7801-10000      - |      GENE_B
  [3]     chr2     0-2199      + |      GENE_C
```

**Summary of Differences:**

- Use `flank()` when we want a region that stops exactly where the gene starts.
- Use `promoters()` when we want a window that spans across the exact start site of the gene.


## 3. Inter-Range Operations (reduce & disjoin)

```r
# Create overlapping ChIP-seq peaks
peaks_gr <- GRanges(
  seqnames = "chr1",
  ranges   = IRanges(
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
```

**Explanations:**

This section covers Inter-Range Operations, which look at how different intervals interact with each other. It creates a new set of data called `peaks_gr` to represent three **ChIP-seq peaks** (regions where proteins bind to DNA).
Notice that the first two peaks overlap:
  - Peak 1: 1500 to 2500
  - Peak 2: 2000 to 3500 (Starts before Peak 1 ends!)
  - Peak 3: 7000 to 9000 (Isolated)

`reduce()` and `disjoin()` are two different ways to clean up or flatten these overlapping regions.

1. `reduce()` — **The Merging Tool**
The `reduce()` function takes overlapping or touching intervals and **flattens them into a single, continuous interval**.
   - **What it does to our peaks**: It sees that `1500–2500` and `2000–3500` overlap. It merges them into one giant peak running from the very beginning of Peak 1 to the very end of Peak 2 (`1500–3500`). Peak 3 doesn't touch anything, so it stays exactly the same. 
   - **The Output (merged_peaks):**
 ```bash
 seqnames    ranges strand
     <Rle> <IRanges>  <Rle>
 [1]     chr1 1500-3500      *
 [2]     chr1 7000-9000      *
 ```
**When to use it**: When we just want to know the total footprint of protein binding and don't care about the individual peak boundaries anymore.


2. `disjoin()` — **The Slicing Tool**
The `disjoin()` function does the opposite. Instead of combining them, it **chops the overlapping regions into distinct, non-overlapping parts** based on where the boundaries cross.
- What it does to our peaks: It looks at the overlap between Peak 1 and Peak 2 and slices them into 3 clean, unique segments:
  1. The part unique to Peak 1: `1500 to 1999`
  2. The part where Peak 1 and Peak 2 overlap: `2000 to 2500`
  3. The part unique to Peak 2: `2501 to 3500`
    Again, Peak 3 (7000–9000) is left alone.
- **The Output** (`disjoined_peaks`):
```bash
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1 1500-1999      *
  [2]     chr1 2000-2500      *
  [3]     chr1 2501-3500      *
  [4]     chr1 7000-9000      *
```
**When to use it**: When you need to analyze the exact sections of DNA that are either shared between samples or strictly unique to one sample.


## 4. Overlap & Intersection Analysis
This section is the core of most genomic analyses: **cross-referencing two different datasets** to see where they intersect.
Here, we are comparing our **ChIP-seq peaks** (`peaks_gr`) against our **gene annotations** (`genes_gr`) to find out which peaks landed inside or near our genes.

```r
# Find which peaks overlap with our gene annotations
overlaps <- findOverlaps(query = peaks_gr, subject = genes_gr)

cat("\n--- Overlap Hits Object ---\n")
print(overlaps)

# Extract matching indices
queryHits(overlaps)   # Indices in peaks_gr
subjectHits(overlaps) # Indices in genes_gr

# Filter peaks that directly overlap any gene using subsetByOverlaps
peaks_on_genes <- subsetByOverlaps(x = peaks_gr, ranges = genes_gr)

cat("\n--- Peaks Overlapping Genes ---\n")
print(peaks_on_genes)
```

1. `findOverlaps()` — **The Matchmaker**
The `findOverlaps()` function looks at every range in the `query` and checks if it physically intersects with any range in the `subject`.

Instead of returning a new set of coordinates, it returns a specialized "Hits" object, which acts like a map of connections. It lists pairs of row numbers matching the query to the subject. If you print `overlaps`, it will show pairs like this:

```bash
Hits object with 3 hits and 0 metadata columns:
      queryHits subjectHits
      <integer>   <integer>
  [1]         1           1
  [2]         2           1
  [3]         3           2
```

- **Row 1**: Peak 1 (`1500-2500`) overlaps Gene A (`1000-3000`).
- **Row 2**: Peak 2 (`2000-3500`) also overlaps Gene A (`1000-3000`).
- (Peak 3 at `7000-9000` is on chr1, but Gene B is on the minus strand of `chr1`, and Gene C is on `chr2`—depending on settings, it won't match Gene B if strand strictness is turned on).

2. `queryHits()` & `subjectHits()` — **Extracting the Row Numbers**
These two functions allow you to pull those raw row numbers out of the Hits object so you can use them in standard R programming loop scripts or filters.

- `queryHits(overlaps)` gives us a simple vector of the row numbers from `peaks_gr` that found a match.
- `subjectHits(overlaps)` gives us the corresponding row numbers from `genes_gr` that were hit.

3. `subsetByOverlaps()` — **The Fast Filter**
If we don't care about the exact mapping pairs and **just want to keep the peaks that hit a gene**, we use `subsetByOverlaps()`.

   - **What it does**: It filters our original `peaks_gr` object. It throws away Peak 3 (because it didn't hit any genes) and keeps Peak 1 and Peak 2.
   - **The Output** (`peaks_on_genes`): A clean `GRanges` object containing only the peaks that successfully overlapped our genes dataset.

**Output:**
```bash
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1 1500-2500      *
  [2]     chr1 2000-3500      *
  [3]     chr1 7000-9000      *
```

We might remember that in our original data, **GENE_B** sits on `chr1` between `5000-8000`.
- Peak 3 is at `7000-9000`.
- Therefore, Peak 3 physically overlaps with **GENE_B** between coordinates `7000` and `8000`.

**Strand Column (*)**
the `strand` column for our peaks shows an asterisk (`*`). In the Bioconductor world, `*` means "unstranded" or "any strand".
Because our ChIP-seq peaks are unstranded, R completely ignored the fact that GENE_B is on the minus (`-`) strand. It only cared that the chromosome name (`chr1`) and the numerical coordinates overlapped.













