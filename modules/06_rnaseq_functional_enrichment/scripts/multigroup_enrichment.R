# ==============================================================================
# Module 08: Comparative & Multi-Group Enrichment using compareCluster
# File: modules/08_comparative_enrichment/08_compare_cluster.R
# ==============================================================================

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(msigdbr)
})

set.seed(42)

# Ensure an output directory exists to organize all exports
dir.create("gsea_outputs", showWarnings = FALSE)

# ==============================================================================
# 1. Simulate Multi-Group Gene Lists (e.g., Time-Series DEGs)
# ==============================================================================
# Fetch Entrez IDs for human genes
all_entrez <- keys(org.Hs.eg.db, keytype = "ENTREZID")[1:3000]

# Simulate 3 experimental timepoints with overlapping/distinct gene sets
deg_time_6h <- sample(all_entrez, 180)
deg_time_12h <- c(sample(deg_time_6h, 60), sample(all_entrez[301:1000], 140))
deg_time_24h <- c(sample(deg_time_12h, 40), sample(all_entrez[1001:2000], 160))

# Package into a named list of character vectors
multi_group_genes <- list(
    "Time_6h"  = deg_time_6h,
    "Time_12h" = deg_time_12h,
    "Time_24h" = deg_time_24h
)

cat("--- Gene List Counts per Cluster ---\n")
print(sapply(multi_group_genes, length))


# ==============================================================================
# 2. Comparative ORA across Groups using GO Biological Process
# ==============================================================================
cat("\n--- Running compareCluster (GO:BP) across 3 Timepoints ---\n")

ck_go <- compareCluster(
    geneClusters  = multi_group_genes,
    fun           = "enrichGO",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    readable      = TRUE
)

# Convert to data frame safely using a defensive row-count shield
if (!is.null(ck_go)) {
    ck_go_df <- as.data.frame(ck_go)
} else {
    ck_go_df <- data.frame()
}

cat("Total comparative GO terms identified:", nrow(ck_go_df), "\n")
if (nrow(ck_go_df) > 0) {
    print(head(ck_go_df[, c("Cluster", "ID", "Description", "p.adjust")], 6))

    # Export multi-group ORA table to CSV
    write.csv(ck_go_df, "gsea_outputs/Comparative_GO_ORA_Results.csv", row.names = FALSE)
    cat("Saved: gsea_outputs/Comparative_GO_ORA_Results.csv\n")
}


# ==============================================================================
# 3. Comparative ORA using Data Frame Formula Notation
# ==============================================================================
cat("\n--- Running compareCluster using Formula Syntax (Entrez ~ Cluster) ---\n")

# Reformat list into a single long data frame (standard tidy output)
df_long <- data.frame(
    Entrez = unlist(multi_group_genes, use.names = FALSE),
    Cluster = rep(names(multi_group_genes), lengths(multi_group_genes)),
    stringsAsFactors = FALSE
)

# Run compareCluster with formula interface
ck_formula <- compareCluster(
    Entrez ~ Cluster,
    data          = df_long,
    fun           = "enrichGO",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 0.05,
    readable      = TRUE
)


# ==============================================================================
# ==============================================================================
# 4. Multi-Group GSEA (Threshold-Free Comparison) - FIXED WITH SIGNAL INJECTION
# ==============================================================================
cat("\n--- Running Multi-Group GSEA across Conditions ---\n")

gene_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")[1:2000]

# Updated function to inject simulated biological signal so GSEA succeeds
set_ranked_vector_with_signal <- function(seed_val, signal_type) {
    set.seed(seed_val)
    stats <- rnorm(length(gene_symbols), mean = 0, sd = 2)
    names(stats) <- gene_symbols

    # Fetch real Hallmark genes to targetedly inject signal
    msig_sample <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")

    if (signal_type == "A") {
        # Inject signal into Hypoxia-associated hallmark genes
        hypoxia_genes <- unique(
            msig_sample$gene_symbol[msig_sample$gs_name == "HALLMARK_HYPOXIA"]
        )
        hypoxia_hits <- names(stats) %in% hypoxia_genes
        stats[hypoxia_hits] <- stats[hypoxia_hits] + 5.0
    } else if (signal_type == "B") {
        # Inject signal into Apoptosis-associated hallmark genes
        apoptosis_genes <- unique(
            msig_sample$gene_symbol[msig_sample$gs_name == "HALLMARK_APOPTOSIS"]
        )
        apoptosis_hits <- names(stats) %in% apoptosis_genes
        stats[apoptosis_hits] <- stats[apoptosis_hits] + 5.0
    } else if (signal_type == "C") {
        # Inject signal into Glycolysis-associated hallmark genes
        glyco_genes <- msig_sample %>%
            dplyr::filter(gs_name == "HALLMARK_GLYCOLYSIS") %>%
            dplyr::pull(gene_symbol) %>%
            unique()
        glyco_hits <- names(stats) %in% glyco_genes
        stats[glyco_hits] <- stats[glyco_hits] + 5.0
    }

    sort(stats, decreasing = TRUE)
}

# Construct treatments carrying independent biological footprints
multi_ranked_list <- list(
    "Treatment_A" = set_ranked_vector_with_signal(101, "A"),
    "Treatment_B" = set_ranked_vector_with_signal(202, "B"),
    "Treatment_C" = set_ranked_vector_with_signal(303, "C")
)

# Fetch database mapping
msig_h <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") %>%
    dplyr::select(gs_name, gene_symbol)

# Run comparative GSEA
ck_gsea <- compareCluster(
    geneClusters = multi_ranked_list,
    fun          = "GSEA",
    TERM2GENE    = msig_h,
    pvalueCutoff = 0.05, # Can keep strict now due to strong signal injection
    verbose      = FALSE
)

if (!is.null(ck_gsea)) {
    ck_gsea_df <- as.data.frame(ck_gsea)
} else {
    ck_gsea_df <- data.frame()
}

cat("Multi-Group GSEA Complete. Total enriched pathways:", nrow(ck_gsea_df), "\n")

# This will now execute successfully because nrow > 0
if (nrow(ck_gsea_df) > 0) {
    write.csv(ck_gsea_df, "gsea_outputs/Comparative_MSigDB_GSEA_Results.csv", row.names = FALSE)
    cat("Saved: gsea_outputs/Comparative_MSigDB_GSEA_Results.csv\n")
}


# ==============================================================================
# 5. Visualizations & Comparative Plotting
# ==============================================================================
cat("\n--- Generating and Saving Multi-Group Visualizations ---\n")

# 5A. Multi-Cluster Dotplot
if (!is.null(ck_go) && nrow(ck_go_df) > 0) {
    p_dot <- dotplot(
        ck_go,
        showCategory = 5,
        title        = "Comparative GO:BP Enrichment Across Timepoints",
        font.size    = 9
    ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    print(p_dot)

    # Save Dotplot as Vector PDF
    ggsave(
        filename = "Comparative_Timecourse_Dotplot.pdf",
        plot     = p_dot,
        width    = 8,
        height   = 7
    )
    cat("Saved Vector PDF: gsea_outputs/Comparative_Timecourse_Dotplot.pdf\n")
} else {
    cat("Skipping Dotplot: No significant data found.\n")
}

# --------------------------------------------------------------------------
# 5B. Functional Treeplot / Clustered Terms across Groups (Optimized)
# --------------------------------------------------------------------------
if (!is.null(ck_go) && nrow(ck_go_df) > 0) {
    cat("\n--- Generating Clean, High-Resolution Treeplot ---\n")

    # 1. Recalculate semantic similarity matrix
    ck_go_sim <- pairwise_termsim(ck_go)

    # 2. Build the optimized treeplot lawet
    p_tree <- treeplot(
        ck_go_sim,
        showCategory     = 3, # Number of top terms to show per cluster
        nCluster         = 5, # Targets the 5 major operational groups
        font.size        = 3, # Reduces font scale to prevent collisions
        label_words_as_y = TRUE, # Aligns category keywords cleanly on the margins
        label_format     = 30 # CRITICAL: Automatically wraps long labels at 30 characters
    ) +
        theme_minimal() +
        theme(
            text = element_text(family = "sans"),
            plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
        ) +
        labs(title = "Clustered Functional Pathway Map Across Timepoints")

    # Display in our active R session window
    print(p_tree)

    # 3. Save as a large canvas vector PDF to prevent any overlapping
    ggsave(
        filename = "Comparative_Functional_Treeplot_Clean.pdf",
        plot     = p_tree,
        width    = 14, # Expanded horizontally to accommodate long names
        height   = 11, # Expanded vertically to completely separate rows
        device   = "pdf"
    )
    cat("Saved Clean Vector PDF: Comparative_Functional_Treeplot_Clean.pdf\n")
} else {
    cat("Skipping Treeplot: No significant data found.\n")
}


# --------------------------------------------------------------------------
# 5C. Multi-Treatment GSEA Dotplot
# --------------------------------------------------------------------------
if (!is.null(ck_gsea) && nrow(ck_gsea_df) > 0) {
    p_dot_gsea <- dotplot(
        ck_gsea,
        showCategory = 5,
        title        = "Comparative MSigDB Hallmark GSEA Across Treatments",
        font.size    = 9
    ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    print(p_dot_gsea)

    # Save GSEA Dotplot as an independent Vector PDF
    ggsave(
        filename = "Comparative_Treatment_GSEA_Dotplot.pdf",
        plot     = p_dot_gsea,
        width    = 8,
        height   = 6
    )
    cat("Saved Vector PDF: Comparative_Treatment_GSEA_Dotplot.pdf\n")
} else {
    cat("Skipping GSEA Dotplot: No significant GSEA data found.\n")
}
