# ==============================================================================
# Module 06: GSEA for GO and KEGG Pathways
# File: modules/06_enrichment_gsea/06_gsea_go_kegg.R
# Description: Threshold-free Gene Set Enrichment Analysis (GSEA) using
#              gseGO() and gseKEGG() from clusterProfiler.
# ==============================================================================

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(ggplot2)
    library(tidyverse)
})

set.seed(42)

# ==============================================================================
# 1. Synthetic Ranked Gene Vector Preparation
# ==============================================================================

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
print(head(ranked_genes, 3))


# ==============================================================================
# 2. GSEA with Gene Ontology (GO)
# ==============================================================================

cat("\n--- Running GSEA for GO (Biological Process) ---\n")

gse_go_res <- gseGO(
    geneList      = ranked_genes,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP", # Options: "BP" (Biological Process), "MF", "CC", "ALL"
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

saveRDS(gse_go_res, file = "./outputs/gse_go_kegg/gse_go_results.rds")


# ==============================================================================
# 3. GSEA with KEGG Pathways
# ==============================================================================

cat("\n--- Running GSEA for KEGG Pathways ---\n")

# Note: gseKEGG requires Entrez IDs by default for human ("hsa")
gse_kegg_res <- gseKEGG(
    geneList      = ranked_genes,
    organism      = "hsa", # "hsa" for Homo sapiens.
    keyType       = "ncbi-geneid", # Matches Entrez IDs
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

saveRDS(gse_kegg_res, file = "./outputs/gse_go_kegg/gse_kegg_results.rds")


# ==============================================================================
# 4. Visualizations
# ==============================================================================

# 4A. Dotplot comparing Normalized Enrichment Scores (NES) across top GO terms
if (nrow(gse_go_df) > 0) {
    p_go_dot <- dotplot(gse_go_res, showCategory = 10, split = ".sign") +
        facet_grid(. ~ .sign) +
        ggtitle("GO BP GSEA: Activated vs Suppressed Pathways") +
        theme_minimal()

    print(p_go_dot)

    ggsave("gsea_go_score.pdf", plot = p_go_dot, width = 8, height = 6)
}

# 4B. GSEA Running Score Plot for the top KEGG pathway
if (nrow(gse_kegg_df) > 0) {
    top_kegg_id <- gse_kegg_df$ID[1]

    p_kegg_run <- gseaplot2(
        gse_kegg_res,
        geneSetID = top_kegg_id,
        title     = paste("KEGG Pathway GSEA:", gse_kegg_df$Description[1])
    )

    print(p_kegg_run)

    ggsave("gsea_kegg_score.pdf", plot = p_kegg_run, width = 8, height = 6)
}
