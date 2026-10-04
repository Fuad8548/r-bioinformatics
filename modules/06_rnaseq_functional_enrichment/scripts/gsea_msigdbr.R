# ==============================================================================
# Module 07: Custom Gene Sets & MSigDB GSEA (Fixed Production Script)
# ==============================================================================

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(msigdbr)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(ggplot2)
    library(tidyverse)
})

set.seed(42)

# --------------------------------------------------------------------------
# 1. Construct Ranked Gene Vector (Symbol-Based)
# --------------------------------------------------------------------------
# Fetch universe of human gene symbols
all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")[1:2500]

# Simulate differential expression statistics (e.g., DESeq2 Wald stats)
stat_values <- rnorm(length(all_symbols), mean = 0, sd = 2.0)
names(stat_values) <- all_symbols

# Inject synthetic biological signal for Hallmark/Reactome enrichment
stat_values[1:180] <- stat_values[1:180] + 4.5 # Up-regulated module
stat_values[181:300] <- stat_values[181:300] - 3.5 # Down-regulated module

# Format as strictly decreasing ranked vector
ranked_genes <- sort(stat_values, decreasing = TRUE)
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]

cat("Total ranked genes for GSEA:", length(ranked_genes), "\n")


# --------------------------------------------------------------------------
# 2. Run GSEA: MSigDB Hallmark Collection (Category "H")
# --------------------------------------------------------------------------
cat("\n--- Running Hallmark GSEA ---\n")
msig_h <- msigdbr(species = "Homo sapiens", collection = "H") %>%
    dplyr::select(gs_name, gene_symbol)

gsea_hallmark <- GSEA(
    geneList      = ranked_genes,
    TERM2GENE     = msig_h,
    minGSSize     = 10,
    maxGSSize     = 500,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    eps           = 1e-10,
    verbose       = FALSE
)

# Clean up Hallmark labels DEFENSIVELY
if (!is.null(gsea_hallmark) && nrow(as.data.frame(gsea_hallmark)) > 0) {
    gsea_hallmark@result$Description <- gsub("^HALLMARK_", "", gsea_hallmark@result$Description)
    res_h_df <- as.data.frame(gsea_hallmark)
} else {
    res_h_df <- data.frame() # Keep as empty dataframe safely
}


# --------------------------------------------------------------------------
# 3. Run GSEA: MSigDB Reactome Collection (Category "C2", Subcategory "CP:REACTOME")
# --------------------------------------------------------------------------
cat("\n--- Running Reactome GSEA ---\n")
msig_reactome <- msigdbr(
    species     = "Homo sapiens",
    category    = "C2",
    subcategory = "CP:REACTOME"
) %>%
    dplyr::select(gs_name, gene_symbol)

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

# Clean up Reactome labels DEFENSIVELY
if (!is.null(gsea_reactome) && nrow(as.data.frame(gsea_reactome)) > 0) {
    gsea_reactome@result$Description <- gsub("^REACTOME_", "", gsea_reactome@result$Description)
    res_reactome_df <- as.data.frame(gsea_reactome)
} else {
    res_reactome_df <- data.frame()
}


# --------------------------------------------------------------------------
# 4. Extract Top Results Safely
# --------------------------------------------------------------------------
cat("\n=================== Hallmark GSEA Summary ===================\n")
cat("Total Significant Pathways:", nrow(res_h_df), "\n\n")
if (nrow(res_h_df) > 0) {
    print(head(res_h_df[, c("ID", "NES", "p.adjust", "setSize")], 5))
} else {
    cat("No significant Hallmark pathways detected.\n")
}

cat("\n=================== Reactome GSEA Summary ===================\n")
cat("Total Significant Pathways:", nrow(res_reactome_df), "\n\n")
if (nrow(res_reactome_df) > 0) {
    print(head(res_reactome_df[, c("ID", "NES", "p.adjust", "setSize")], 5))
} else {
    cat("No significant Reactome pathways detected.\n")
}

# Save the significant Reactome pathways to a CSV file
if (nrow(res_reactome_df) > 0) {
    write.csv(
        x         = res_reactome_df,
        file      = "Reactome_GSEA_Results.csv",
        row.names = FALSE
    )
    cat("Saved Tabular Data: gsea_outputs/Reactome_GSEA_Results.csv\n")
} else {
    cat("Skipping CSV Export: No significant Reactome pathways to save.\n")
}

# --------------------------------------------------------------------------
# 5. Visualization: Dotplot & Running Score Comparison
# --------------------------------------------------------------------------
cat("\n--- Generating Visualizations ---\n")

# Comparative Dotplot: Use Reactome since Hallmark is empty
if (nrow(res_reactome_df) > 0) {
    p_dot_r <- dotplot(
        gsea_reactome,
        showCategory = 10,
        font.size    = 8,
        title        = "MSigDB Reactome GSEA Enrichment"
    ) + theme_minimal()

    print(p_dot_r)
} else {
    cat("Skipping Dotplot: No significant data found.\n")
}

ggsave(
    filename = "Reactome_GSEA_Dotplot.pdf",
    plot     = p_dot_r,
    width    = 7,
    height   = 6
)

# Running score plot for top Reactome pathway
if (nrow(res_reactome_df) > 0) {
    p_gsea_reactome <- gseaplot2(
        gsea_reactome,
        geneSetID    = 1,
        title        = paste("Top Reactome Set:", res_reactome_df$Description[1]),
        pvalue_table = TRUE
    )
    print(p_gsea_reactome)
}

ggsave(
    filename = "Reactome_Top_Running_Score.pdf",
    plot     = p_gsea_reactome,
    width    = 8,
    height   = 7
)
