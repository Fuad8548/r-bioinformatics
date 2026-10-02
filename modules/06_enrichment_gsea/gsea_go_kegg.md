# Module 06 Functional Enrichment Analysis and GSEA

## Description: 
Differential expression analysis yields a list of thousands of p-values and fold-changes, but raw gene lists do not explain biological mechanisms. Functional Enrichment Analysis maps these gene lists onto curated databases (such as Gene Ontology, KEGG, and MSigDB) to identify active biological processes, metabolic pathways, and cellular structures.

## Conceptual Framework: ORA vs. GSEA
Functional analysis uses two primary statistical paradigms: **Over-Representation Analysis (ORA) and Gene Set Enrichment Analysis (GSEA)**.

```text
    ALL MEASURED GENES (Universe)
                 │
  ┌──────────────┴──────────────┐
  ▼                             ▼
Threshold-Based (ORA)         Threshold-Free (GSEA)
- Filter: padj < 0.05         - Rank ALL genes by metric 
- Discrete DEG list           - Continuous vector
- Fisher's Exact / Hypergeom  - Running Kolmogorov-Smirnov
```

1. Here we are going to discuss about GSEA:
GSEA evaluates changes across entire pathways without applying strict DEG significance cutoffs.
   - **Input**: A continuous vector containing **all measured genes**, ranked by a metric reflecting magnitude and direction of expression change. Common ranking metrics ($S_i$) include:
     - DESeq2 Wald Statistic (`stat` column)
     - Signed $-\log_{10}(p\text{-value})$: $\text{sign}(\log_2 \text{FoldChange}) \times -\log_{10}(p)$
   - **Algorithm Mechanics**:
     - **Rank Ordering**: Sort all $N$ genes from top upregulated to top downregulated: $L = \{g_1, g_2, \dots, g_N\}$.
     - **Running Enrichment Score (ES)**: Walk down $L$. Increase a cumulative running total when a gene is in Set $S$ (weighted by ranking metric); decrease it when a gene is absent.
     - **Peak Score (ES)**: The maximum deviation from zero corresponds to the Enrichment Score.
     - **Normalized Enrichment Score (NES)**: Adjusts ES for gene set size to allow direct comparison across pathways.
     - **Significance Evaluation**: Permutes sample labels (or gene labels) to calculate empirical false discovery rates ($q$-values).

## Core script
```r
# Description: Gene Set Enrichment Analysis (GSEA) using clusterProfiler, org.Hs.eg.db; gseGO() and gseKEGG() from clusterProfiler.

# 1. Synthetic Ranked Gene Vector Preparation

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(tidyverse)
})

set.seed(42)

# Fetch valid human Entrez IDs
all_entrez_ids <- keys(org.Hs.eg.db, keytype = "ENTREZID")[1:2500]

# Generate synthetic ranking metrics (e.g., DESeq2 Wald statistic or signed -log10 p-value)
stats_vector <- rnorm(length(all_entrez_ids), mean = 0, sd = 3)
names(stats_vector) <- all_entrez_ids

# Inject simulated pathway signal into top genes
stats_vector[1:200] <- stats_vector[1:200] + 4.5

# CRITICAL FOR GSEA: Vector must be named with Entrez IDs and sorted in DESCENDING order
ranked_genes <- sort(stats_vector, decreasing = TRUE)
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]

cat("Total ranked genes for GSEA:", length(ranked_genes), "\n")
cat("Top 3 ranked genes:\n")
print(head(ranked_genes, 3) %>% enframe(name = "Entrez_ID", value = "Metric_Score"))
```

### Explanation:
1. **Getting a List of Human Genes**
```r
all_entrez_ids <- keys(org.Hs.eg.db, keytype = "ENTREZID")[1:2500]
```
- What it does: It grabs the official ID numbers (Entrez IDs) for the first 2,500 human genes from a standard biological database (`org.Hs.eg.db`).

2. **Making Up Fake Test Results**
```r
stats_vector <- rnorm(length(all_entrez_ids), mean = 0, sd = 3)
names(stats_vector) <- all_entrez_ids
```
- It generates 2,500 random numbers using a normal distribution (bell curve) and assigns one number to each gene ID. In a real experiment, these numbers would represent how much a gene changed (e.g., if a drug turned a gene "up" or "down").

3. **Rigging the Data (Injecting a Signal)**
```r
stats_vector[1:200] <- stats_vector[1:200] + 4.5
```
- What it does: It artificially adds 4.5 to the scores of the first 200 genes. This guarantees that these specific 200 genes look heavily "activated" or significant.
- Why? When we run the actual GSEA software later, we want to make sure the software is working properly by seeing if it successfully detects this "rigged" group of genes.

4. **Sorting the Genes (The critical step for GSEA)**
```r
ranked_genes <- sort(stats_vector, decreasing = TRUE)
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]
```
- What it does: It sorts the entire list of genes from the highest score to the lowest score and removes any accidental duplicates.
- Why? GSEA software strictly requires the data to be ordered this way. It needs to look at the "winners" (highest scores) at the very top of the list.

5. **Printing the Results**
```r
cat("Total ranked genes for GSEA:", length(ranked_genes), "\n")
print(head(ranked_genes, 3) ... )
```
- What it does: It prints a summary to your screen, showing you the total number of genes (2,500) and a sneak peek at the top 3 highest-scoring genes in a neat little table.

### Output:
```bash
Entrez_ID     Metric_Score
  <chr>            <dbl>
1 142               12.6
2 18                11.4
3 2100              10.8
```
- **Entrez_ID** <chr>: This column lists the unique identification numbers for the human genes. The <chr> means "character" (text format).
  - For example, 142 is the official database ID for a real human gene(specifically, a gene named PARP1).
- **Metric_Score** <dbl>: This column shows how highly activated the gene is. The <dbl> stands for "double," which just means a number with decimals.
  - Gene 142 has the highest score of 12.6.


# 2. GSEA with Gene Ontology (GO)
```r
cat("\n--- Running GSEA for GO (Biological Process) ---\n")

gse_go_res <- gseGO(
  geneList      = ranked_genes,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",             # Options: "BP" (Biological Process), "MF", "CC", "ALL"
  keyType       = "ENTREZID",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  eps           = 1e-10,
  verbose       = FALSE
)

# Convert S4 gseaResult object to a data frame
gse_go_df <- as.data.frame(gse_go_res)

cat("Significantly Enriched GO BP Terms:", nrow(gse_go_df), "\n")
if (nrow(gse_go_df) > 0) {
  print(head(gse_go_df[, c("ID", "Description", "NES", "p.adjust")], 5))
}
```

**Explanation:**
This script takes your sorted list of genes and runs them through the `gseGO()` function (which comes from the `clusterProfiler` package) to find matching biological pathways.
1. `gseGO()` arguments:
- `geneList` = `ranked_genes`: This inputs our sorted, named vector of genes that we created in the first step.
- `OrgDb = org.Hs.eg.db`: This tells the function to use the Human genome database mapping tool to read our Entrez ID numbers.
- `ont = "BP"`: This limits the search strictly to Biological Processes (e.g., cell division, metabolism). We could alternatively look for "MF" (Molecular Functions) or "CC" (Cellular Components).
- `keyType = "ENTREZID"`: This informs the software that our gene list uses Entrez ID numbers (like 142) rather than gene symbols (like PARP1).
- `minGSSize = 15` and `maxGSSize = 500`: These filter out pathways that are either too small (fewer than 15 genes) or too broad/vague (more than 500 genes) to be biologically meaningful.
- `pvalueCutoff = 0.05`: This tells the code to only keep pathways with a statistically significant **p-value of less than 0.05**.
- `pAdjustMethod = "BH"`: This uses the Benjamini-Hochberg method to adjust the p-values. It controls for the "false discovery rate" because testing thousands of pathways at once can create accidental random matches.
- `eps = 1e-10`: This sets the boundary calculation accuracy for calculating p-values so the math doesn't stall out on highly significant pathways.
- `verbose = FALSE`: This tells R to run quietly without printing progress updates to the screen.

3. **Converting and Displaying the Data**
```r
gse_go_df <- as.data.frame(gse_go_res)
```
- What it does: The raw output of gseGO is a complex specialized object (an S4 object). This line flattens it into a standard Data Frame (spreadsheet format) so humans can easily read, filter, and print it.

### output
```bash
ID            Description       NES
GO:0034440    lipid oxidation   2.081343
GO:0009190    cyclic nucleotide biosynthetic process    2.039035
GO:0006164    purine nucleotide biosynthetic process    2.024051
GO:1903409    reactive oxygen species biosynthetic process    1.991242
GO:0009152    purine ribonucleotide biosynthetic process    1.975862
              p.adjust
GO:0034440    0.02906019
GO:0009190    0.03419405
GO:0006164    0.03169326
GO:1903409    0.03803520
GO:0009152    0.03169326
```

- `ID`: The unique identification code from the Gene Ontology (GO) database.
- `Description`: The actual biological name of the pathway.
- `NES` (**Normalized Enrichment Score**): This measures how heavily concentrated our top-ranked genes are within this pathway.• A positive NES (like 2.07) means the genes in this pathway were overwhelmingly found at the very top of our ranked list (the "winners").
- `p.adjust`: The adjusted p-value. It tells us the statistical confidence of the result, corrected for random chance. Biologists look for a value below 0.05 to consider a pathway "statistically significant."

# Visualization

```r
# 4A. Dotplot comparing Normalized Enrichment Scores (NES) across top GO terms
if (nrow(gse_go_df) > 0) {
  p_go_dot <- dotplot(gse_go_res, showCategory = 5, split = ".sign") +
    facet_grid(. ~ .sign) +
    ggtitle("GO BP GSEA: Activated vs Suppressed Pathways") +
    theme_minimal()
  
  print(p_go_dot)
  ggsave("gsea_go_score.pdf", plot = p_go_dot, width = 8, height = 6)
}
```

## Explanation:
**The GO Dotplot (p_go_dot)**
This section creates a Dotplot for our Gene Ontology results.
- `dotplot(gse_go_res, showCategory = 10, split = ".sign")`: This pulls the top 10 pathways and splits them based on their "sign" (whether they are Activated with a positive NES score, or Suppressed with a negative NES score).
- `facet_grid(. ~ .sign)`: This splits the visual window into two separate side-by-side panels: one panel showing activated pathways, and one showing suppressed pathways.
• What you will see: A chart with dots. The size of each dot represents how many genes overlap with that pathway, and the color intensity represents the statistical significance (p-value).

### What the plot shows:
1. **The X-Axis: Gene Ratio**
• What it means: The `GeneRatio` (tracked at the bottom from 0.35 to 0.50) tells us **what percentage of the core genes in that pathway** showed up in your highly ranked "rigged" list.
- **According to our output:**
	- `cyclic nucleotide biosynthetic process` and `reactive oxygen species`... are pulled all the way to the right at `0.50`. This means an impressive 50% of all the genes that control those pathways were sitting at the top of your ranked list.
	- lipid oxidation sits at `0.45` (45%).
	- The two purine pathways sit on the left near `0.31`.

2. **The Dot Color: Adjusted P-Value (`p.adjust`)**
- What it means: The color gradient measures statistical confidence. According to your legend, red/pink represents highly significant values (closer to `0.029`), while blue represents values approaching `0.038`.
- According to your output:
	- `lipid oxidation` has the lowest p-value (`0.02906`). On our plot, its dot is the most **vibrant red**, marking it as the most statistically reliable result.
	- `reactive oxygen species...`. has the highest p-value (`0.03803`). On our plot, it shifts noticeably to a cool **blue dot**.

3. **The Dot Size: Count**
• What it means: The scale on the right shows that the size of the circle is determined by the absolute number of genes from our data that belong to that specific pathway.
- According to our output:
	- Look at the `lipid oxidation` dot—it is the largest circle on the screen. This means it has the highest number of overlapping genes (around 16–18 genes) anchoring it down.
	- In contrast, the other dots are smaller, meaning fewer individual genes were required to calculate their enrichment.



# 3. GSEA with KEGG Pathways

```r
cat("\n--- Running GSEA for KEGG Pathways ---\n")

# Note: gseKEGG requires Entrez IDs by default for human ("hsa")
gse_kegg_res <- gseKEGG(
  geneList      = ranked_genes,
  organism      = "hsa",            # "hsa" for Homo sapiens, "mmu" for mouse, etc.
  keyType       = "ncbi-geneid",    # Matches Entrez IDs
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  eps           = 1e-10,
  verbose       = FALSE
)

# Convert to data frame
gse_kegg_df <- as.data.frame(gse_kegg_res)

cat("Significantly Enriched KEGG Pathways:", nrow(gse_kegg_df), "\n")
if (nrow(gse_kegg_df) > 0) {
  print(head(gse_kegg_df[, c("ID", "Description", "NES", "p.adjust")], 5))
}
```

### Output:
```bash
ID                        Description                       NES
hsa04923    Regulation of lipolysis in adipocytes           2.105008
hsa04611    Platelet activation                             2.082016
hsa00230    Purine metabolism                               2.054358
hsa01522    Endocrine resistance                            1.999560
hsa04213    Longevity regulating pathway - multiple species 1.946000
            p.adjust
hsa04923 0.002155450
hsa04611 0.002201219
hsa00230 0.002155450
hsa01522 0.008218669
hsa04213 0.015583606
```

**Script Explanations:**
This code takes your sorted gene list and looks for matching biochemical networks or molecular signaling pathways in the KEGG database.
- `gseKEGG(...)`: This function initiates the Gene Set Enrichment Analysis using the KEGG database instead of the GO database.
- `organism = "hsa"`: This specifies the organism. `hsa` is the official KEGG shortcode for **Homo sapiens** (human).
- `keyType = "ncbi-geneid"`: This tells the function that our gene list is named using official NCBI/Entrez ID numbers (like 142), which is the native format KEGG expects.
- `minGSSize = 15` and `maxGSSize = 500`: These filter out pathways that are either too small or too large to be helpful.
- `pvalueCutoff = 0.05`: This ensures the software only saves pathways with a statistically significant p-value.
- `pAdjustMethod = "BH"` and `eps = 1e-10`: This applies the **Benjamini-Hochberg** adjustment to control for false positives and sets mathematical precision boundaries.
- `gse_kegg_df <- as.data.frame(gse_kegg_res)`: This flattens the complex R object into a standard spreadsheet data frame.
- `if (nrow(gse_kegg_df) > 0) { ... }`: If R finds matching pathways, it extracts and displays only the top 5 rows (`head(..., 5)`) along with 4 specific columns `c(ID, Description, NES, and p.adjust)`.

**Output Explanation:**
- `ID`: The unique KEGG identification code. All human KEGG pathways start with hsa.
• `Description`: The specific biochemical process or disease network.
• `NES` (Normalized Enrichment Score): All values here are positive and high (ranging from 1.94 to 2.10). This means the genes inside these pathways were heavily crowded at the very top of our ranked list.
• `p.adjust`: All of these adjusted p-values are well below 0.05 (e.g., 0.002 is vastly lower than 0.05), proving these results are highly statistically valid.

**Breakdown of Top Pathways Found:**
1. `hsa04923` (Regulation of lipolysis in adipocytes): This pathway dictates how fat cells break down lipids. It has the highest score (NES = 2.10), meaning your highest-scoring genes heavily overlap with fat metabolism.


# Visualization
```r
# 4B. GSEA Running Score Plot for the top KEGG pathway
if (nrow(gse_kegg_df) > 0) {
  top_5_kegg_id <- gse_kegg_df$ID[1:5]
  
  p_kegg_run <- gseaplot2(
    gse_kegg_res,
    geneSetID = top_5_kegg_id,
    title     = "KEGG Pathway GSEA: Top 5 Activated Networks",
    pvalue_table  = TRUE
  )
  
  print(p_kegg_run)
  ggsave("gsea_kegg_score.pdf", plot = p_kegg_run, width = 8, height = 6)
}
```

**The KEGG Running Score Plot (`p_kegg_run`)**
This section creates the iconic "GSEA Plot" for our #1 top KEGG pathway (which was hsa04923: Regulation of lipolysis in adipocytes).
- `top_kegg_id <- gse_kegg_df$ID[1]`: This grabs the ID string (`"hsa04923"`) from the first row of our KEGG data frame automatically.
- `gseaplot2(...)`: This creates a three-tiered classic GSEA visualization:
	1. Top Panel: A green line showing the "Running Enrichment Score." Because our genes are highly activated, we will see this line spike upward drastically on the left side of the chart and then taper off.
	2. Middle Panel: A "barcode" plot. Every single vertical black line represents a gene from that specific lipid pathway. You will see a heavy cluster of black bars on the far-left side (the "winners" section).
	3. Bottom Panel: A ranking metric map showing your Wald statistics/scores plunging from high positive numbers down to negative numbers.
• ggsave(...): This automatically exports and saves that gorgeous KEGG running score chart as a clean, publication-ready PDF file named gsea_running_score.pdf in your current working folder.

### Explanation of the KEGG plot:
1. **Top Panel: The Running Enrichment Score (RES)**
   - What it shows: The line graph peaks near the top left. The software walks down our list of 2,500 genes one by one. Every time it hits a gene that belongs to the lipid pathway, the score goes up. When it hits a gene that doesn't, the score drops slightly.
   - Connection to our output: The peak height of this curve is used to calculate our NES (2.105008). Because we artificially added `4.5` to the top genes, the highest peak belongs to Regulation of lipolysis in adipocytes (`NES = 2.10`), while the lowest of the five will be Longevity regulating pathway (`NES = 1.94`).

2. **Middle Panel: The Barcode Hits**
   - What it shows: A series of black vertical lines (the "barcode"). Each line marks the exact location of a single gene from the Regulation of lipolysis in adipocytes pathway.
   - Connection to our output: We notice a dense cluster of black lines crowded on the far-left side of the graph. This visual clustering matches our highly significant `p.adjust` value of 0.002155. If the results were random, these black bars would be scattered evenly across the entire 0-to-2500 span.

3. **Bottom Panel: Ranked List Metric**
   - What it shows: This panel displays the raw scores we assigned to the 2,500 genes (plotted on the Y-axis from roughly 10 down to -10).
   - Connection to our output: The steep drop we see on the left shows the exact boundary where our top 200 "rigged" genes end and the rest of the randomly distributed genes begin.







































