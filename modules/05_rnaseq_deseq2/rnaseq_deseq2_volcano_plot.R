# ===============================================================================
# 5. Visualization of Results
# ==============================================================================

# Ensure ggrepel is loaded for smart text positioning
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)
library(ggplot2)
library(dplyr)

# 1. Classify genes and flag the top 5 for labeling
plot_data <- de_results_df %>%
    mutate(
        status = case_when(
            log2FoldChange >= 1 & padj < 0.05 ~ "Up-regulated",
            log2FoldChange <= -1 & padj < 0.05 ~ "Down-regulated",
            TRUE ~ "Not Significant"
        ),
        # Create a column that only contains names for the top 5 genes
        label_gene = if_else(gene_id %in% head(de_results_df$gene_id, 5), gene_id, "")
    )

# 2. Build the enhanced plot
p_volcano <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj))) +
    # Background points
    geom_point(aes(color = status), alpha = 0.6, size = 2) +
    scale_color_manual(values = c(
        "Up-regulated" = "#E41A1C",
        "Down-regulated" = "#377EB8",
        "Not Significant" = "#999999"
    )) +

    # CRITICAL FEATURE 1: Significance Cutoff Grid Lines
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgray", linewidth = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgray", linewidth = 0.6) +

    # CRITICAL FEATURE 2: Smart Gene Labels
    geom_text_repel(
        aes(label = label_gene),
        size = 4,
        fontface = "bold",
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "black",
        max.overlaps = Inf
    ) +

    # Styling and clean layout
    theme_classic() +
    labs(
        title = "Volcano Plot: Differential Gene Expression",
        subtitle = "Treated vs Control (Top 5 labeled; Thresholds: |LFC| > 1, padj < 0.05)",
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-value (padj)",
        color = "Expression Status"
    ) +
    theme(
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "right"
    )

# 3. Render and export the plot
print(p_volcano)
ggsave("deseq2_volcano_plot.pdf", plot = p_volcano, width = 8, height = 6)
