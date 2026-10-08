# ==============================================================================
# Module 01: Tidyverse Foundations for Bioinformatics
# File: modules/01_tidyverse_foundations/01_tidyverse_foundations.R
# Description: Data manipulation (dplyr), reshaping (tidyr), and plotting (ggplot2)
#              using synthetic biological expression and metadata tables.
# ==============================================================================

# 0. Load Core Libraries & Source Project Utilities
library(tidyverse)
source("R/plotting_theme.R") # Loads theme_bio()

# Set reproducible seed for synthetic dataset generation
set.seed(42)

# ==============================================================================
# 1. Synthetic Dataset Creation (RNA-Seq Metadata & Expression)
# ==============================================================================

# Sample metadata frame
metadata <- tibble(
    sample_id  = paste0("SMP_", 1:12),
    condition  = rep(c("Control", "Treated"), each = 6),
    cell_line  = rep(c("HeLa", "HEK293"), times = 6),
    batch      = sample(c("Batch_A", "Batch_B"), 12, replace = TRUE),
    rin_score  = round(runif(12, min = 6.8, max = 9.9), 1) # RNA Integrity Number
)

# Gene expression data in wide matrix format (log2 TPM)
expression_wide <- tibble(
    sample_id = paste0("SMP_", 1:12),
    TP53      = rnorm(12, mean = 8.5, sd = 0.6),
    BRCA1     = rnorm(12, mean = ifelse(metadata$condition == "Treated", 11.2, 7.8), sd = 0.8),
    EGFR      = rnorm(12, mean = ifelse(metadata$condition == "Treated", 5.4, 9.1), sd = 0.7),
    GAPDH     = rnorm(12, mean = 13.0, sd = 0.2) # Housekeeping control
)

# ==============================================================================
# 2. Reshaping Data with tidyr (Wide <-> Long)
# ==============================================================================

# Pivot wide matrix into tidy/long format (essential for ggplot2 & Bioconductor)
expression_long <- expression_wide %>%
    pivot_longer(
        cols      = -sample_id,
        names_to  = "gene",
        values_to = "log2_tpm"
    )

# Merge expression values with sample metadata
tidy_genomics_df <- expression_long %>%
    inner_join(metadata, by = "sample_id")

# ==============================================================================
# 3. Data Wrangling with dplyr
# ==============================================================================

# 3.1 Filter high-quality samples and target genes
filtered_df <- tidy_genomics_df %>%
    filter(rin_score >= 7.5, gene %in% c("BRCA1", "EGFR"))

# 3.2 Compute mean expression & fold-change across conditions
summary_stats <- tidy_genomics_df %>%
    group_by(gene, condition) %>%
    summarize(
        mean_exp = mean(log2_tpm),
        sd_exp = sd(log2_tpm),
        n_samples = n(),
        .groups = "drop"
    ) %>%
    pivot_wider(
        names_from  = condition,
        values_from = c(mean_exp, sd_exp)
    ) %>%
    mutate(
        log2_fold_change = mean_exp_Treated - mean_exp_Control
    )

print(summary_stats)

# ==============================================================================
# 4. Publication-Ready Visualizations with ggplot2
# ==============================================================================

# 4.1 Comparative Expression Boxplot with Jitter
p1 <- ggplot(tidy_genomics_df, aes(x = gene, y = log2_tpm, fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(0.8)) +
    geom_jitter(position = position_dodge(0.8), size = 2, alpha = 0.8) +
    scale_fill_manual(values = c("Control" = "#2b5c8f", "Treated" = "#d95f02")) +
    labs(
        title    = "Target Gene Expression Across Experimental Groups",
        subtitle = "Comparing Control vs. Treated samples across cell lines",
        x        = "Target Gene",
        y        = "Expression (log2 TPM)",
        fill     = "Condition"
    ) +
    theme_bio()

# Save output plot
ggsave("outputs/figures/m01_gene_expression_boxplot.png", plot = p1, width = 8, height = 5, dpi = 300)
