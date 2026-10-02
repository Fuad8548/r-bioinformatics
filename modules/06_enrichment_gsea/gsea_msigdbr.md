1. **Injecting Biological Signal & Sorting**
```r
stat_values[1:180] <- stat_values[1:180] + 4.5 
stat_values[181:300] <- stat_values[181:300] - 3.5 

ranked_genes <- sort(stat_values, decreasing = TRUE)
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]
```

- **Signal Injection**: By manually shifting the scores of genes 1–180 upwards and 181–300 downwards, the script creates a simulated signature (an up-regulated and a down-regulated module) for the GSEA algorithm to discover.
- `sort(..., decreasing = TRUE)`: Crucial. GSEA tracks a running sum calculation that walks down our gene list from top-activated to top-suppressed. If this vector is not strictly sorted in descending order, the algorithm will calculate incorrect metrics or fail.

2. **Database Queries & GSEA Execution**
```r
msig_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
    dplyr::select(gs_name, gene_symbol)
```

- `msigdbr()`: Programmatically pulls the **Molecular Signatures Database (MSigDB)** straight into R as a clean dataframe.
	- Category "H" pulls the **Hallmark Collection** (50 broad, well-defined biological states).
	- **Category "C2" / "CP:REACTOME"** pulls the Reactome Canonical Pathways (highly specific, physics-and-chemistry-based biochemical networks).
- `select(gs_name, gene_symbol)`: Trims the download into a clean 2-column key-value layout (`TERM2GENE`) required by `clusterProfiler`.


```r
gsea_reactome <- GSEA(
    geneList      = ranked_genes,
    TERM2GENE     = msig_reactome,
    minGSSize     = 10,
    maxGSSize     = 500,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    eps           = 1e-10,
    verbose       = FALSE
)
```

- `GSEA()`: Runs the core permutation math. It compares our `ranked_genes` against the custom target pathway lists (`TERM2GENE`).
• `minGSSize` / `maxGSSize`: Discards pathways that are too narrow (fewer than 10 matching genes) or too sweeping (more than 500 genes) to protect against mathematical bias.
• `pvalueCutoff = 0.05` & `pAdjustMethod = "BH"`: Applies a Benjamini-Hochberg False Discovery Rate correction to limit false positives while filtering out non-significant pathways.


3. **Defensive Cleanup & Parsing**
```r
if (!is.null(gsea_reactome) && nrow(as.data.frame(gsea_reactome)) > 0) {
    gsea_reactome@result$Description <- gsub("^REACTOME_", "", gsea_reactome@result$Description)
    res_reactome_df <- as.data.frame(gsea_reactome)
} else {
    res_reactome_df <- data.frame()
}
```
- Defensive Checking (if): If a small dataset yield zero significant hits, `clusterProfiler` returns a NULL object. Attempting to manipulate columns on an empty slot crashes a script. This block acts as a protective shield.
• `gsub("^REACTOME_", "", ...)`: Cleans up text clutter. It strips out the database prefix from the results object so that when we generate a plot, the labels print cleanly as "G_ALPHA_Z_SIGNALLING_EVENTS" instead of the verbose "REACTOME_G_ALPHA_Z_SIGNALLING_EVENTS"


## Output:
```bash
ID
REACTOME_G_ALPHA_Z_SIGNALLING_EVENTS                                                                                                
REACTOME_ADENYLATE_CYCLASE_INHIBITORY_PATHWAY                                                                                       
REACTOME_GLUCAGON_SIGNALING_IN_METABOLIC_REGULATION                                                                                 
REACTOME_ADORA2B_MEDIATED_ANTI_INFLAMMATORY_CYTOKINES_PRODUCTION                                                                   
REACTOME_HIGH_LAMINAR_FLOW_SHEAR_STRESS_ACTIVATES_SIGNALING_BY_PIEZENDOTHELIAL_CELLS 
                                                                                          NES
REACTOME_G_ALPHA_Z_SIGNALLING_EVENTS                                                 2.215561
REACTOME_ADENYLATE_CYCLASE_INHIBITORY_PATHWAY                                        2.127279
REACTOME_GLUCAGON_SIGNALING_IN_METABOLIC_REGULATION                                  2.124312
REACTOME_ADORA2B_MEDIATED_ANTI_INFLAMMATORY_CYTOKINES_PRODUCTION                     2.081539
REACTOME_HIGH_LAMINAR_FLOW_SHEAR_STRESS_ACTIVATES_SIGNALING_BY_PIEZENDOTHELIAL_CELLS 1.973407
                                                                                        p.adjust
REACTOME_G_ALPHA_Z_SIGNALLING_EVENTS                                                 0.006371968
REACTOME_ADENYLATE_CYCLASE_INHIBITORY_PATHWAY                                        0.006662178
REACTOME_GLUCAGON_SIGNALING_IN_METABOLIC_REGULATION                                  0.006662178
REACTOME_ADORA2B_MEDIATED_ANTI_INFLAMMATORY_CYTOKINES_PRODUCTION                     0.006662178
REACTOME_HIGH_LAMINAR_FLOW_SHEAR_STRESS_ACTIVATES_SIGNALING_BY_PIEZENDOTHELIAL_CELLS 0.013356248
                                                                                     setSize
REACTOME_G_ALPHA_Z_SIGNALLING_EVENTS                                                      27
REACTOME_ADENYLATE_CYCLASE_INHIBITORY_PATHWAY                                             12
REACTOME_GLUCAGON_SIGNALING_IN_METABOLIC_REGULATION                                       20
REACTOME_ADORA2B_MEDIATED_ANTI_INFLAMMATORY_CYTOKINES_PRODUCTION                          26
REACTOME_HIGH_LAMINAR_FLOW_SHEAR_STRESS_ACTIVATES_SIGNALING_BY_PIEZENDOTHELIAL_CELLS      32
```

## Explanation of the output
1. **Plot 1: MSigDB Reactome GSEA Enrichment (Dotplot)**
- `Description` (or `ID`): This maps to the labels on the y-axis (e.g., PKA MEDIATED PHOSPHORYLATION OF CREB, ADENYLATE CYCLASE INHIBITORY PATHWAY). The "`REACTOME_`" prefix was successfully stripped off by our `gsub` code block.
- `GeneRatio`: This maps to the x-axis positioning. In our CSV, this value represents the number of core enrichment genes found in our dataset divided by the total size of the pathway. For instance, *PKA MEDIATED PHOSPHORYLATION OF CREB* has a high GeneRatio of `0.7`.
- `p.adjust`: This maps directly to the dot color spectrum. Highly significant pathways like *GLUCAGON SIGNALING* appear bright red (closer to `0.007`), while pathways like *GABA B RECEPTOR ACTIVATION* appear dark blue (closer to `0.020` in our table).
- `Count` (derived from `setSize`): This dictates the size of the dot. Pathways with larger dots (like *GLUCAGON SIGNALING* or *G ALPHA Z SIGNALLING EVENTS*) contain a higher absolute count of genes shifting together in our ranked data.

2. **Plot 2: Running Score Plot (Gseaplot2)**
- `pvalue` and `p.adjust`: The exact numbers inside the inset table of the plot (pvalue 1.159e-05 and p.adjust 0.006372) match the floating-point values saved in the corresponding rows of your CSV file.
- `core_enrichment`: This is a critical column in your CSV containing a list of gene symbols separated by slashes (e.g., `GNAZ/GNAI1/...`). Visually, these genes are represented by the **dense cluster of vertical black barcode** ticks bunched up on the far-left side of the chart. They are the specific genes that successfully pushed the green enrichment line up to its maximum peak.












