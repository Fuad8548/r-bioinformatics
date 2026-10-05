# ==============================================================================

# File: modules/09_microbiome/01_phyloseq_standard_workflow.R


suppressPackageStartupMessages({
    library(phyloseq)
    library(ggplot2)
    library(dplyr)
})

set.seed(42)

# ==============================================================================
# 1. Simulate Realistic ASV Count Table, Taxonomy, and Metadata
# ==============================================================================
cat("\n--- Step 1: Simulating Microbiome Multi-Table Data ---\n")

num_asvs <- 60
num_samples <- 12

asv_names <- paste0("ASV_", sprintf("%03d", 1:num_asvs))
sample_names <- paste0("Sample_", sprintf("%02d", 1:num_samples))

# Generate count matrix (ASVs x Samples)
asv_counts <- matrix(
    rpois(num_asvs * num_samples, lambda = 20),
    nrow = num_asvs,
    ncol = num_samples,
    dimnames = list(asv_names, sample_names)
)

# Introduce synthetic abundance shift in Treatment group (Samples 7-12)
asv_counts[1:10, 7:12] <- asv_counts[1:10, 7:12] + rpois(10 * 6, lambda = 100)

# Build Taxonomy Table
phyla_list <- c("Bacteroidetes", "Firmicutes", "Proteobacteria", "Actinobacteria")
tax_matrix <- matrix(
    c(
        rep("Bacteria", num_asvs),
        sample(phyla_list, num_asvs, replace = TRUE),
        paste0("Class_", sample(1:4, num_asvs, replace = TRUE)),
        paste0("Order_", sample(1:5, num_asvs, replace = TRUE)),
        paste0("Family_", sample(1:5, num_asvs, replace = TRUE)),
        paste0("Genus_", sample(1:10, num_asvs, replace = TRUE))
    ),
    nrow = num_asvs,
    ncol = 6,
    dimnames = list(asv_names, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
)

# Sample Metadata
sample_metadata <- data.frame(
    SampleID  = sample_names,
    Group     = rep(c("Control", "Treatment"), each = 6),
    Batch     = rep(c("Batch_1", "Batch_2"), times = 6),
    row.names = sample_names
)


# ==============================================================================
# 2. Assemble phyloseq S4 Container
# ==============================================================================
cat("\n--- Step 2: Assembling phyloseq S4 Object ---\n")

ps <- phyloseq(
    otu_table(asv_counts, taxa_are_rows = TRUE),
    tax_table(tax_matrix),
    sample_data(sample_metadata)
)

cat("Phyloseq Summary:\n")
print(ps)


# ==============================================================================
# 3. Alpha Diversity Metrics
# ==============================================================================
cat("\n--- Step 3: Calculating Within-Sample Alpha Diversity ---\n")

# Calculate Observed Richness, Shannon, and Simpson indices
alpha_df <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
cat("Sample Alpha Diversity Head:\n")
print(head(alpha_df, 4))

p_alpha <- plot_richness(ps, x = "Group", measures = c("Observed", "Shannon")) +
    geom_boxplot(aes(fill = Group), alpha = 0.4) +
    theme_bw() +
    labs(title = "Microbial Alpha Diversity across Conditions")

print(p_alpha)


# ==============================================================================
# 4. Relative Abundance Transformation & Beta Diversity Ordination
# ==============================================================================
cat("\n--- Step 4: Beta Diversity & PCoA Ordination ---\n")

# Transform counts to relative proportions (0 to 1)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Perform Principal Coordinate Analysis (PCoA) using Bray-Curtis distance
bray_pcoa <- ordinate(ps_rel, method = "PCoA", distance = "bray")

p_beta <- plot_ordination(ps_rel, bray_pcoa, color = "Group") +
    geom_point(size = 4, alpha = 0.8) +
    theme_bw() +
    labs(title = "PCoA Ordination (Bray-Curtis Distance)")

print(p_beta)


# ==============================================================================
# 5. Taxonomic Aggregation & Composition Barplots
# ==============================================================================
cat("\n--- Step 5: Taxonomic Composition at Phylum Rank ---\n")

# Collapse ASVs to Phylum level
ps_phylum <- tax_glom(ps_rel, taxrank = "Phylum")

p_bar <- plot_bar(ps_phylum, x = "SampleID", fill = "Phylum") +
    facet_wrap(~Group, scales = "free_x") +
    theme_bw() +
    labs(title = "Phylum-Level Relative Abundance Composition", y = "Relative Abundance")

print(p_bar)

cat("\nModule 09 Pipeline Completed Successfully!\n")


# 1. Save Alpha Diversity Boxplot
ggsave(
    filename = "microbiome_alpha_diversity.pdf",
    plot     = p_alpha,
    width    = 7,
    height   = 5
)

# 2. Save Beta Diversity PCoA Plot
ggsave(
    filename = "microbiome_beta_pcoa.pdf",
    plot     = p_beta,
    width    = 7,
    height   = 5
)

# 3. Save Taxonomic Phylum Composition Barplot
# Note: Barplots with many samples benefit from a wider canvas
ggsave(
    filename = "microbiome_phylum_composition.pdf",
    plot     = p_bar,
    width    = 10,
    height   = 6
)
