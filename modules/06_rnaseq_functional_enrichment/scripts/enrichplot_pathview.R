suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(pathview)
    library(ggplot2)
    library(ggnewscale)
    library(dplyr)
})

set.seed(42)

# ==============================================================================
# 1. Simulate DE Data with Expression Fold Changes
# ==============================================================================

# Fetch Entrez IDs for human genes
all_entrez <- keys(org.Hs.eg.db, keytype = "ENTREZID")[1:3000]

# Simulate log2 fold changes (log2FC)
log2fc_vector <- rnorm(length(all_entrez), mean = 0, sd = 2.0)
names(log2fc_vector) <- all_entrez

# Extract significant DEGs for ORA and keep named log2FC vector for mapping
sig_degs <- names(log2fc_vector[abs(log2fc_vector) > 1.5])

cat("Total background genes:", length(all_entrez), "\n")
cat("Total significant DEGs (|log2FC| > 1.5):", length(sig_degs), "\n")


# ==============================================================================
# 2. Run KEGG & GO Functional Enrichment
# ==============================================================================

# 2A. KEGG Over-Representation Analysis
kegg_res <- enrichKEGG(
    gene          = sig_degs,
    universe      = all_entrez,
    organism      = "hsa",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    pAdjustMethod = "BH"
)

# 2. Convert Entrez IDs to Gene Symbols
kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")


# Gene-Concept Network Plot (cnetplot)
p_cnet <- cnetplot(
    kegg_res,
    showCategory   = 5,
    foldChange     = log2fc_vector,
    layout         = igraph::layout_in_circle,
    color_edge     = "category",
    node_label     = "all",
    size_category  = 1.2,
    size_item      = 0.8
) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(title = "Gene-Concept Network (cnetplot) with Fold Change")

print(p_cnet)


# Save the plot with specified dimensions (adjust width/height if your labels get cut off)
ggsave(
    filename = "kegg_cnetplot.pdf",
    plot     = p_cnet,
    device   = "pdf",
    width    = 10, # Width in inches
    height   = 10, # Height in inches
    units    = "in",
    dpi      = 300 # High resolution production quality
)


# 2B. GO Biological Process Enrichment
go_res <- enrichGO(
    gene          = sig_degs,
    universe      = all_entrez,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 1, # Relax to 1 to bypass filtering
    qvalueCutoff  = 1, # Prevent strict q-value dropouts
    pAdjustMethod = "BH",
    readable      = TRUE
)


# Enrichment Map Plot (emapplot)
# Clusters overlapping gene sets into a network based on pairwise term similarity
go_sim <- pairwise_termsim(go_res)

library(enrichplot)
library(ggplot2)

# Generate Enrichment Map Plot
library(enrichplot)
library(ggplot2)

p_emap <- emapplot(
    go_sim,
    showCategory    = 20,
    node_label_size = 3 # Standardized parameter in version 1.32.1+ (default is 5)
) +
    labs(title = "GO Enrichment Map (emapplot)")

print(p_emap)

ggsave("go_emapplot.pdf", plot = p_emap, width = 10, height = 10, units = "in")


# Functional Treeplot (treeplot)
# Hierarchical clustering of enriched terms based on semantic similarity
p_tree <- treeplot(
    go_sim,
    showCategory = 15,
    nCluster     = 4,
    label_format = 30
) + labs(title = "Hierarchical Functional Treeplot")

print(p_tree)

ggsave("go_treeplot.pdf", plot = p_tree, width = 12, height = 10, units = "in")

# Heatmap-like Plot (heatplot)
# Displays gene-pathway relationships in a compact grid matrix
p_heat <- heatplot(
    kegg_res,
    showCategory = 10,
    foldChange   = log2fc_vector
) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))

print(p_heat)

ggsave(
    filename = "kegg_heatplot.pdf",
    plot     = p_heat,
    device   = "pdf",
    width    = 16,
    height   = 8,
    units    = "in",
    dpi      = 300
)


# ==============================================================================
# 4. Native KEGG Pathway Mapping with pathview
# ==============================================================================

cat("\n--- Rendering Pathview Native KEGG Diagrams ---\n")

# Pathview requires named vectors with ENTREZ IDs as names
if (nrow(as.data.frame(kegg_res)) > 0) {
    top_kegg_id <- kegg_res@result$ID[1] # e.g., "hsa04110"

    cat("Generating Pathview diagram for top KEGG term:", top_kegg_id, "\n")

    # Renders and saves PNG + XML files locally in the working directory
    pathview_out <- pathview(
        gene.data    = log2fc_vector,
        pathway.id   = top_kegg_id,
        species      = "hsa",
        limit        = list(gene = max(abs(log2fc_vector)), cpd = 1),
        bins         = list(gene = 10, cpd = 10),
        low          = list(gene = "blue", cpd = "green"),
        mid          = list(gene = "gray", cpd = "gray"),
        high         = list(gene = "red", cpd = "yellow"),
        kegg.native  = TRUE, # Native KEGG PNG vs Graphviz rendering
        same.layer   = FALSE
    )

    cat("Pathview output generated successfully:", paste0(top_kegg_id, ".pathview.png"), "\n")
}
