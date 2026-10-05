suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
})

set.seed(123)

# ------------------------------------------------------------------------------
# 1. Simulate Directional DEGs across 3 Treatment Conditions
# ------------------------------------------------------------------------------
all_entrez <- keys(org.Hs.eg.db, keytype = "ENTREZID")[1:4000]

# Generate 6 distinct clusters: (Drug_A, Drug_B, Drug_C) x (Up, Down)
deg_clusters <- list(
    "Drug_A_Up"   = sample(all_entrez[1:600], 160),
    "Drug_A_Down" = sample(all_entrez[601:1200], 130),
    "Drug_B_Up"   = sample(all_entrez[100:700], 180),
    "Drug_B_Down" = sample(all_entrez[1201:1800], 145),
    "Drug_C_Up"   = sample(all_entrez[300:900], 170),
    "Drug_C_Down" = sample(all_entrez[1500:2100], 125)
)

cat("Cluster gene counts:\n")
print(sapply(deg_clusters, length))


# ------------------------------------------------------------------------------
# 2. Execute compareCluster (GO Biological Process)
# ------------------------------------------------------------------------------
cat("\n--- Running compareCluster across 6 Directional Clusters ---\n")

ck_directional <- compareCluster(
    geneClusters  = deg_clusters,
    fun           = "enrichGO",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    readable      = TRUE
)

ck_df <- as.data.frame(ck_directional)

cat("Total comparative terms identified across all clusters:", nrow(ck_df), "\n\n")


# ------------------------------------------------------------------------------
# 3. Extract Top Enriched Pathways per Cluster Direction
# ------------------------------------------------------------------------------
top_per_cluster <- ck_df %>%
    group_by(Cluster) %>%
    slice_min(order_by = p.adjust, n = 2, with_ties = FALSE) %>%
    select(Cluster, ID, Description, GeneRatio, p.adjust, geneID)

print(top_per_cluster[, c("Cluster", "Description", "GeneRatio", "p.adjust")])


# ------------------------------------------------------------------------------
# 4. Customized Dotplot Visualization
# ------------------------------------------------------------------------------
p_dot <- dotplot(
    ck_directional,
    showCategory = 3,
    font.size    = 8,
    title        = "Directional Comparative Enrichment: Up vs Down across Treatments"
) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
        axis.text.y = element_text(size = 8),
        panel.grid.major = element_line(color = "grey92")
    )

print(p_dot)


# ------------------------------------------------------------------------------
# 5. Export Data Table and High-Resolution Plot
# ------------------------------------------------------------------------------

# Export the entire comparative data frame to a CSV file
write.csv(
    ck_df,
    file      = "directional_GO_comparative_results.csv",
    row.names = FALSE
)
cat("Saved enrichment table to 'directional_GO_comparative_results.csv'\n")


# Save the customized dotplot using ggplot2's ggsave
# Adjust width/height as needed depending on the length of our GO terms
ggsave(
    filename = "directional_comparative_dotplot.pdf",
    plot     = p_dot,
    device   = "pdf",
    width    = 9, # Width in inches
    height   = 7, # Height in inches
    bg       = "white" # Forces solid background
)
cat("Saved publication-ready plot to 'directional_comparative_dotplot.pdf'\n")
