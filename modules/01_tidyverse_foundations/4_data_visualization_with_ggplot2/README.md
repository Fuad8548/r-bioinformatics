
# 4. Publication-Ready Visualizations with ggplot2

## 4.1 Comparative Expression Boxplot with Jitter

```r
library(ggplot2)

# Make sure we use 'tidy_genomics_df' for plotting individual points
ggplot(tidy_genomics_df, aes(x = condition, y = log2_tpm, fill = condition)) +
  # 1. Add cleaner boxplots with transparency
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5) +
  
  # 2. Add individual sample points to see actual distribution spread
  geom_jitter(width = 0.15, size = 2.5, shape = 21, color = "black", alpha = 0.8) +
  
  # 3. Split the plot into 4 sub-panels (one for each gene)
  facet_wrap(~ gene, scales = "free_y") +
  
  # 4. Apply clean journal aesthetics
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("Control" = "#4A90E2", "Treated" = "#E15759")) +
  
  # 5. Label axes cleanly
  labs(
    title = "RNA-Seq Synthetic Gene Expression Profile",
    subtitle = "Comparing Treatment vs. Control Across target and control genes",
    x = "Experimental Condition",
    y = "Expression Level (log2 TPM)",
    fill = "Condition"
  )

# Save output plot
ggsave(
  filename = "outputs/figures/m01_gene_expression_boxplot.pdf",
  width = 8, height = 5, units = "in"
)
```

**Explanation:**
1. **Target Genes (Clear Biological Effects)**
   - `BRCA1` **(Up-regulated by Treatment)**:
     - There is a major upward shift in expression from the Control group (centering around 7.2–8.0) to the Treated group (centering around 9.5–11.0).
     - The individual sample dots show zero overlap between the two conditions. This indicates a highly consistent, strong induction of `BRCA1` activity under treatment.
   - `EGFR` **(Down-regulated by Treatment)**:We see the exact opposite trend. Under control conditions, EGFR maintains a high baseline around 9.0–10.0.
   - Following treatment, its expression plummets sharply to a tight cluster around 5.0–6.0. The treatment strongly suppresses or down-regulates this gene.

2. **Control Genes (Experimental Benchmarks)**
   - `GAPDH` **(Housekeeping Control)**:
     - Notice the y-axis scale for this panel: it spans an incredibly narrow range (from 12.75 to 13.25).
     - Both the Control and Treated boxplots sit at almost identical heights (~13.0). This confirms our lab processing worked perfectly: `GAPDH` expression remains stable and is unaffected by the drug treatment.
   - `TP53` **(Unchanged Dynamic Baseline)**:
     - Unlike `GAPDH`, the data points for TP53 are much more vertically spread out (ranging from 7.5 to 9.5), meaning this gene naturally has higher sample-to-sample biological variance.
     - However, because the median lines (the dark horizontal bars inside the boxes) sit at roughly the same level (~8.2) for both groups, we can conclude that the treatment had no actual effect on `TP53`.