# Code Functionality: Multi-Group Functional Enrichment Analysis

In high-throughput biological experiments (like RNA-seq), researchers rarely look at just a single treatment versus control. Instead, they track complex conditions such as **time-series progressions** (e.g., 6h vs 12h vs 24h) or **multiple parallel drug variants** (e.g., Treatment A vs B vs C).

Instead of running separate enrichment analyses for each group and trying to manually overlap Venn diagrams, this script uses `clusterProfiler::compareCluster()` to evaluate all groups side-by-side in a single calculation. It handles the statistical corrections uniformly and generates multi-panel comparative figures (dotplots and treeplots) directly suited for publication.

1. **Section 1: Synthetic Data Creation (Time-Series Simulation)**

```r
all_entrez <- keys(org.Hs.eg.db, keytype = "ENTREZID")[1:3000]

deg_time_6h  <- sample(all_entrez, 180)
deg_time_12h <- c(sample(deg_time_6h, 60), sample(all_entrez[301:1000], 140))
deg_time_24h <- c(sample(deg_time_12h, 40), sample(all_entrez[1001:2000], 160))
```
- Concept: Fetches 3,000 real human Entrez Gene IDs.
- Simulation Mechanics: Simulates a time-course experiment by sampling 180–200 genes per group. It purposefully programs explicit overlaps between early, middle, and late groups (e.g., 60 genes from 6h carry over into 12h) to mimic real cell-signaling dynamics.

2. **Section 2: List-Based Comparative ORA**
```r
ck_go <- compareCluster(
    geneClusters  = multi_group_genes,
    fun           = "enrichGO",
    ont           = "BP",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    readable      = TRUE
)
```
- `compareCluster()`: The engine of this module. It takes the named list of genes and loops over them using the function specified in `fun = "enrichGO"` (Gene Ontology).
- `ont = "BP"`: Targets the Biological Process branch of Gene Ontology.
- `readable = TRUE`: Automatically translates confusing Entrez vector numbers (e.g., 7157) into reader-friendly Gene Symbols (e.g., TP53) inside the final data tables.
- **Defensive Logic**: The script safely converts the result to a dataframe (`as.data.frame`) inside an if block to ensure that if an experiment returns zero significant pathways, the script skips the step gracefully instead of crashing. It saves the resulting table directly to `Comparative_GO_ORA_Results.csv`.

3. **Section 3: Formula Notation Interface (Tidy Data Layout)**
```r
df_long <- data.frame(Entrez = ..., Cluster = ...)
ck_formula <- compareCluster(Entrez ~ Cluster, data = df_long, ...)
```
- Purpose: Demonstrates advanced syntax capability. In production pipelines (like `DESeq2` or `edgeR` mixed with tidyverse), data rarely sits in neat lists; it is usually structured as a long table.
- `Entrez ~ Cluster`: Tells R: "*Analyze the Entrez genes grouped by the experimental label column named Cluster*." This allows us to skip complex data restructuring steps.

4. **Section 4: Multi-Group GSEA (Threshold-Free)**
```r
ck_gsea <- compareCluster(
    geneClusters = multi_ranked_list,
    fun          = "GSEA",
    TERM2GENE    = msig_h,
    pvalueCutoff = 0.1
)
```
- Purpose: Shifts from Over-Representation Analysis (which uses strict cutoffs like $\(p < 0.05\)$) to threshold-free Gene Set Enrichment Analysis (**GSEA**) across groups.
- Mechanics: It constructs ranked vectors for Treatment A, B, and C using a helper function, then runs multi-group GSEA against the MSigDB Hallmark database collection (msigdbr) to find moving biological pathways without using rigid filters. It saves the output to `Comparative_MSigDB_GSEA_Results.csv`.

5. **Section 5: Vector PDF Export and Data Visualizations**
This code block performs a Multi-Group Gene Set Enrichment Analysis (GSEA) across three distinct experimental conditions (`Treatment_A`, `Treatment_B`, and `Treatment_C`) using the **MSigDB Hallmark collection**.
Unlike the Over-Representation Analysis (ORA) code in the earlier sections—which uses arbitrary cutoffs (like $\(p < 0.05\)$) to evaluate isolated lists of genes—this module tracks **threshold-free expression continuums** for all conditions side by side.

- **The Helper Function & Target Database Mapping**
```r
set_ranked_vector_with_signal <- function(seed_val, signal_type) { ... }
```
- `rnorm(...)`: It initializes a clean vector of 2,000 human gene symbols, assigning each a random numeric differential expression value (resembling standard statistics like a DESeq2 Wald test or limma t-statistic).
• `msigdbr(...)`: Inside the function, the script programmatically pulls the MSigDB Hallmark collection ("H") for human cells. This acts as a biological mapping directory.

- **Simulating Unique Group-Specific Fingerprints**
To make sure the GSEA algorithm doesn't encounter mathematical noise and return zero results, the function uses conditional `if` branches to purposefully inject a strong biological signal (+5.0) into targeted sets of genes:
  - **Treatment A (Hypoxia Blueprint):**
    ```r
    stats[hypoxia_hits] <- stats[hypoxia_hits] + 5.0
    ```
    The script flags any gene that belongs to the official HALLMARK_HYPOXIA registry and pushes its score up. When this vector is sorted using sort(..., decreasing = TRUE), these Hypoxia genes are forced to cluster tightly at the very top of the list.
  - **Treatment B (Apoptosis Blueprint)**: It looks up the HALLMARK_APOPTOSIS registry and applies the (+5.0) shift exclusively to those genes.
  - **Treatment C (Glycolysis Tidy Blueprint)**: It uses clean tidyverse syntax (filter, pull) to find HALLMARK_GLYCOLYSIS genes, shifting them upward. 

3. **The `compareCluster` GSEA Engine**
```r
ck_gsea <- compareCluster(
    geneClusters = multi_ranked_list,
    fun          = "GSEA",
    TERM2GENE    = msig_h,
    pvalueCutoff = 0.05,
    verbose      = FALSE
)
```
This is the core calculation phase. The `compareCluster` wrapper changes its internal mathematics because we passed it `fun = "GSEA"`:
- It recognizes that our inputs are fully ranked numeric vectors rather than simple character lists of gene names.
- It systematically calculates a **running enrichment score (ES)** for every single pathway in the database against all three treatment lists.
- It corrects the calculated p-values globally across all groups simultaneously using the false discovery rate wrapper, applying a strict `pvalueCutoff = 0.05`.


## Results and Plots
1. **Comparative MSigDB Hallmark GSEA**
This chart displays the GSEA data. Because GSEA evaluates a sorted continuum of genes, the `GeneRatio` indicates how much of a pathway’s core membership shifted heavily to the top of that tr
- **Treatment_A (The Hypoxia Fingerprint)**: Our CSV shows a massive Normalized Enrichment Score for `HALLMARK_HYPOXIA` (**NES = 3.28**) with a highly significant adjusted `p-value` ($\(4.2 \times 10^{-09}\)$). In the plot, this is the massive, bright red circle at the top left, displaying a high `GeneRatio` (~0.85). It also captures a weaker cross-activation of `HALLMARK_GLYCOLYSIS` (NES = 2.06, blue dot) and `TGF_BETA_SIGNALING`.
- **Treatment_B (The Apoptosis Fingerprint)**: This group shifts `HALLMARK_APOPTOSIS` to the extreme top (**NES = 3.27**, large red dot) alongside `HALLMARK_P53_PATHWAY` (blue dot), representing classic cellular stress or programmed cell death responses.
- **Treatment_C (The Glycolysis Fingerprint)**: This treatment drives `HALLMARK_GLYCOLYSIS` to its highest level (**NES = 3.05**, large red dot), while showing a minor, secondary enrichment of `HYPOXIA`.

2. **Timecourse Dotplot — GO:BP ORA**
This plot shifts the focus to **Over-Representation Analysis (ORA)** across your timepoints. Unlike GSEA, the `Count` and dots represent the overlapping proportion of filtered target lists.
- **Time_6h (Early Responders)**: Our CSV logs `humoral immune response` as the top hit ($\(p.\text{adjust} = 0.002\)$, Count = 13 genes). The plot matches this as the purple circle at the very top left. It also shows early physiological changes through `positive regulation of cytosolic calcium ion concentration` ($\(p.\text{adjust} = 0.035\)$, Count = 8) and `amino acid metabolic process`.
- **Time_12h (Intermediate Immune Responders)**: The cellular focus moves deeply into adaptive immunity. The plot illustrates unique activations for `regulation of lymphocyte activation` (Count = 18), `regulation of T cell activation` (Count = 15), and `B cell activation` (Count = 12) that were completely absent at 6 hours.
- **Time_24h (Late Tissue Organizers)**: At 24 hours, the immune responses fade out, replaced by massive developmental and structural pathways. These processes have the **highest statistical significance** of the entire experiment ($\(p.\text{adjust} \approx 4.8 \times 10^{-06}\)$). Large, dark red circles confirm an intense genetic focus on `regulation of hormone levels` (Count = 20), `fatty acid metabolic process` (Count = 18), and structural tissue systems like `eye development` and `sensory system development` (Count = 15).

3. **Functional Treeplot — GO:BP ORA**
The treeplot collapses the text redundancy of our GO terms. It groups pathways into multi-colored blocks on the left based on shared gene components, using multi-colored pie charts to explain group behavior.
- **Cluster 5 (Pink Block - Top)**: Groups eye development and fatty acid metabolic processes. The pie circles are solid blue, indicating these mechanisms are driven exclusively by genes turning on at the `Time_24h` milestone.
- **Cluster 1 & 2 (Blue Block - Bottom)**: Groups regulation of T cell/lymphocyte activation. The pie circles are solid green, confirming this branch of lymphocyte proliferation triggers only within the `Time_12h` window.
• **Shared Signaling Triggers (Center Circles)**: Look closely at the split pie charts for `humoral immune response` and `positive regulation of cytosolic calcium ion concentration`. These circles contain distinct slices of Red (`Time_6h`), Green (`Time_12h`), and Blue (`Time_24h`). This provides biological proof of a sustained core workflow; while individual genes may change, the physical processes of calcium signaling and humoral defense remain continuously active throughout your entire experimental timeline.

















