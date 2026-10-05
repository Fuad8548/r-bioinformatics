# Advanced Pathway Visualization & Custom Lawets
# Description: Publication-ready pathway and network visualizations combining enrichplot (cnetplot, emapplot, treeplot, heatplot) and pathview.

### Scientific purpose of each package
1. `enrichplot` (Global & Structural Visualization)
When we run GSEA or ORA on thousands of genes, we often get 100+ "significant" pathways—many of which share 80% of the same genes (e.g., Cell Cycle, Mitotic Nuclear Division, DNA Replication).
   - **Reduces Term Redundancy** (`emapplot` / `treeplot`): Groups overlapping pathways into functional clusters so we can report 4 m ain biological themes instead of reading a list of 50 redundant terms.
   - **Shows Gene-Pathway Cross-Talk** (`cnetplot`): Reveals if a single key hub gene (e.g., TP53 or MYC) is driving the enrichment signal across multiple distinct pathways.
   - **Visualizes Directionality** (`heatplot`): Displays fold-changes for individual genes within each pathway side-by-side.

2. `pathview` (Mechanistic & Native Map Overlay)
`enrichplot` shows us that a pathway is enriched, but it doesn't show us how the cascade works. pathview downloads the official KEGG diagram for a pathway and paints the pathway boxes based on our experimental data.
   - **Positional Insight**: Distinguishes whether our up-regulated genes are upstream cell-surface receptors (membrane level) or downstream transcription factors (nucleus level).
   - **Multi-Layered Regulation**: Shows if an entire linear cascade is activated or if a specific feedback loop/inhibition block is occurring.

### Application Across Omics Disciplines
```r
                    High-Throughput Experiment 
               (RNA-seq / Proteomics / Metabolomics / GWAS)
                                   │
                                   ▼
                   Quantified Matrix (log2FC, p-vals)
                                   │
                                   ▼
             Functional Enrichment (clusterProfiler/ GSEA)
                                   │
         ┌─────────────────────────┴─────────────────────────┐
         ▼                                                   ▼
   enrichplot                                            pathview
 Global Network &                                 Mechanistic KEGG Maps
 Redundancy Reduction                             (Genes, Proteins, Metabolites)
```

1. **RNA-seq & Single-Cell RNA-seq** (Transcriptomics):
   - Usage: Map differential gene expression ($log_2\text{FC}$) to identify shifted cellular processes across cell types or treatment conditions.
2. **Proteomics** (Mass Spectrometry / Quantitative Proteomics):
   - Usage: Map protein abundance ratios ($log_2\text{FC}$) or spectral count changes.
   - Note: Convert UniProt/RefSeq IDs to Entrez IDs or Gene Symbols before feeding data into the pipeline.
3. **Metabolomics**:
   - Usage: `pathview` uniquely supports simultaneous dual-layer mapping! We can feed gene/protein expression via `gene.data` and metabolite concentrations via `cpd.data` (using KEGG Compound IDs like `C00031` for Glucose) to paint both genes and metabolites on the exact same KEGG pathway map.
4. **Genomics (GWAS, WES, CNV Analysis)**:
   - Usage: Identify functional pathways enriched for high-risk variants, mutated genes, or copy-number alterations (using $Z$-scores or mutation frequencies instead of fold-changes).


### Significance of Notable Parts of the Code
- **The Background Universe** (`universe = all_entrez`)
   - A list of all genes that were actually detected or capable of being measured in our experiment.
   - Significance: This is vital for mathematical accuracy. To know if 50 genes belonging to a "viral defense" pathway is a surprising finding, a hypergeometric statistical test needs to know if those 50 came out of a pool of 3,000 expressed genes or 30,000 total human genes. Omitting the true background universe leads to high rates of false positives.

- **Relaxed Statistical Cutoffs** (`pvalueCutoff = 1, qvalueCutoff = 1`)
  - Overriding the default filtering thresholds (p < 0.05).
  - Significance: In this specific script, because the data is randomly simulated, none of the pathways would naturally pass strict biological significance testing. Relaxing these parameters forces R to yield outputs anyway. In a real experiment, we would restore these to `0.05` to filter out background noise.

- **ID Translation** (`setReadable()`)
  - Mapping Entrez IDs (e.g., `4193`) to official Gene Symbols (e.g., MDM2).
  - Significance: Computational databases prefer numbers because they never change. Humans prefer symbols because they mean something. This function translates the machine-friendly data into human-friendly language right before plotting, preventing our final figures from being cluttered with unreadable digit blocks.

- **Semantic Similarity Calculations** (`pairwise_termsim()`)
  - Measuring how much overlap or meaning two biological terms share.
  - Significance: Databases like Gene Ontology are highly redundant. We might get five separate enriched terms that essentially mean "cell division." This function measures that redundancy so advanced plots (emapplot and treeplot) can cluster those sister terms together, simplifying our final interpretation.

### Advantages of This Pipeline
- **Standardized & Broadly Accepted**: `clusterProfiler` and `pathview` are gold-standard, peer-reviewed Bioconductor packages. Using them ensures our methodology aligns with current bioinformatics publication standards.
- **Diverse Visual Lawets**: Rather than relying entirely on basic bar charts, the pipeline offers multi-dimensional views. We can view interactions at the gene level (`cnetplot`), category level (`emapplot`), grouping level (`treeplot`), and inside cellular pathway maps (`pathview`).
- **Dynamic Data Integration**: The plots don't just show which pathways are active; they overlay our expression data (log₂FC values) directly onto the figures, showing the directional flow (upregulation vs. downregulation) of the biology.

### Limitations of This Approach
- Over-Representation Analysis (ORA) Blindspots: The pipeline relies on ORA, meaning we must pick an arbitrary threshold (like $\(\vert{}log_2FC\vert{} > 1.5\)$) to cut off our gene list. A gene with a value of 1.49 is completely thrown away, even though it might be biologically vital.
- Database Annotation Gaps: The analysis is only as good as the databases it queries. If a gene has not been thoroughly studied or annotated in the `org.Hs.eg.db` or KEGG repositories, it will be ignored, even if it is the driving force behind our experiment's phenotype.
- Static Pathway Representations: The `pathview` package draws data on top of static KEGG maps. It cannot capture real-time cellular dynamics, tissue-specific variations, or alternative protein splicing.


### Visualizations
1. **Gene-Concept Network** (`cnetplot`)
- The intricate web linking individual genes (peripheral dots) to their assigned KEGG pathway hubs (central colored circles).
- Key Observations:
  - Hub Centrality: The Pathways in cancer category sits at the center with a high `itemNum` (large circle size), meaning a vast proportion of our significant genes belong to it.
  - Gene Expression Overlay: The peripheral dots are colored on a scale from blue (downregulated) to red (upregulated). For instance, looking closely near Pathways in cancer, we can trace specific heavily upregulated genes (dark red) vs downregulated genes (dark blue) fueling that term.
- Interpretation: This lets us pinpoint pleiotropic genes—genes that connect to multiple disease pathways simultaneously (e.g., genes sitting on intersecting lines between Breast cancer and Pathways in cancer).

2. **Heatmap-like Plot** (`heatplot`)
- A grid matrix view of the `cnetplot` data. Pathways are on the y-axis, and individual Entrez gene numeric IDs are on the x-axis.
- Key Observations:
  - Sparsity vs Density: Pathways in cancer shows a dense, continuous sequence of vertical color bars, indicating massive gene coverage. Conversely, pathways like Arachidonic acid metabolism or Th17 cell differentiation are sparse, containing only a few active genes from our list.
  - Expression Patterns: The vertical bars change color dynamically based on fold change. We can easily scan to find clusters where a pathway is dominated by strongly upregulated (red) or downregulated (blue) genes.
- Interpretation: This is an excellent alternative to the cnetplot when we have too many genes, as it displays gene names and directional changes without turning into a tangled "hairball" network.

3. **GO Enrichment Map** (`emapplot`)
- A network showing how much our enriched GO terms overlap with one another. Lines (edges) connect terms that share a significant number of identical genes.
- Key Observations:
  - The Glucose/Inflammation Island: There is a highly interconnected cluster at the bottom right linking hexose transmembrane transport, regulation of D-glucose transmembrane transport, and inflammatory response to antigenic stimulus.
  - Orphan Terms: Categories like sulfur amino acid metabolic process and single fertilization float completely isolated in space.
- Interpretation: The connected cluster reveals that our simulated gene expression changes are driving a unified immuno-metabolic axis, whereas the orphan nodes represent completely unrelated, isolated biological events.

4. **Hierarchical Functional Treeplot** (`treeplot`)
- This plot groups redundant GO terms based on semantic similarity, creating clear functional blocks.
- Key Observations:
  - Cluster 1 (Green): Houses developmental and baseline structural terms like circulatory system development and cell differentiation involved in kidney development.
  - Cluster 2 (Blue): Groups metabolic and signaling pathways related to glucose control, specifically response to insulin and insulin-like growth factor receptor signaling.
  - Cluster 3 (Orange) & Cluster 4 (Purple): Isolate distinct modules for immune responses (inflammatory response to antigenic stimulus) and reproductive biology (single fertilization) respectively.
- Interpretation: Instead of reading 15 independent lines, the treeplot tells us that our gene list is split across four primary biological themes.


5. **Native KEGG Pathway Map** (`pathview`)
- What it shows: A high-resolution biological schematic of "Pathways in cancer" fetched directly from the KEGG database, with our simulated expression data overlaid on top.
- Key Observations:
  - Sub-Pathway Disruption: Looking at the rectangular protein boxes. They are split into colored segments (blue for down, red for up).
  - Apoptosis & Proliferation Blocks: In the center-right of the canvas, the core cancer hallmarks—"Evading apoptosis" and "Proliferation"—receive multiple incoming signals. Downregulated genes (blue boxes like BAD, BAX, or CASP9) and upregulated genes (red boxes) are scattered across the Jak-STAT, Wnt, and MAPK signaling cascades, illustrating exactly where the pathway's signal flow is altered.
- Interpretation: This is our most mechanistic figure. It moves away from abstract statistical charts to show exactly where cell signaling breaks down inside a physical cell diagram.



































