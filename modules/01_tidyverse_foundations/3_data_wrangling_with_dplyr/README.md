
# 3. Data Wrangling with dplyr

## 3.1 Filter high-quality samples and target genes

```r
filtered_df <- tidy_genomics_df |>
  filter(rin_score >= 7.5, gene %in% c("BRCA1", "EGFR"))

# 3.2 Compute mean expression & fold-change across conditions
summary_stats <- tidy_genomics_df |>
  group_by(gene, condition) |>
  summarize(
    mean_exp = mean(log2_tpm),
    sd_exp = sd(log2_tpm),
    n_samples = n(),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from  = condition,
    values_from = c(mean_exp, sd_exp)
  ) |>
  mutate(
    log2_fold_change = mean_exp_Treated - mean_exp_Control
  )

print(summary_stats)
```

**Explanation:**
It performs two distinct tasks: 
  - It throws away degraded biological samples and isolates specific genes of interest;
  - It calculates exactly how much gene expression changed between the **Control** and **Treated** groups.

1. **Filtering Quality and Selecting Targets (filtered_df)**
```r
filtered_df <- tidy_genomics_df |>
  filter(rin_score >= 7.5, gene %in% c("BRCA1", "EGFR"))
```
Think of `filter()` as a sieve. It drops any data rows that do not meet both of these conditions:
- `rin_score >= 7.5`: In RNA sequencing, a low RNA Integrity Number (RIN) means the sample degraded in the lab. This keeps only high-quality samples.
- `gene %in% c("BRCA1", "EGFR")`: It ignores `TP53` and `GAPDH`, keeping rows only for these two specific target genes.

2. **Calculating Differential Expression (`summary_stats`)**
- **Grouping and Aggregating (`group_by` + `summarize`)**
```r
  group_by(gene, condition) |>
  summarize(
    mean_exp = mean(log2_tpm),
    sd_exp = sd(log2_tpm),
    n_samples = n(),
    .groups = "drop"
  )
```
  - `group_by(gene, condition)`: R splits the dataset into subsets (e.g., all `BRCA1` rows that are `Control`, all `BRCA1` rows that are Treated, etc.).
  - `summarize(...)`: R calculates the **average** expression (`mean_exp`), the standard deviation/variance (`sd_exp`), and counts how many samples are in that group (`n_samples`).

- **Spreading the Data Side-by-Side (pivot_wider)**
```r
  pivot_wider(
    names_from  = condition,
    values_from = c(mean_exp, sd_exp)
  )
```
To calculate mathematical differences between groups, it is easier if the Control and Treated averages sit right next to each other on the same row. `pivot_wider` reshapes the data to look like this:
| Operation |   Description    | Header 4         |
| :-------: | :--------------: | ---------------- |
|   gene    | mean_exp_Control | mean_exp_Treated |
|   BRCA1   |       7.70       | 11.2             |

- **Finding the Log-Fold Change (mutate)**
```r
  mutate(
    log2_fold_change = mean_exp_Treated - mean_exp_Control
  )
```
  - `log2_fold_change`: Because the data is already on a $\[\log _{2}\]$ scale, subtracting the Control average from the Treated average yields the log-fold change.
    - For example: $\(11.2 - 7.8 = +3.4\)$. This positive number proves that BRCA1 was significantly up-regulated by the treatment.


**Summary Table**
| Operation | CIGAR Code |   Description    | Header 4         | Header 5       |
| :-------: | :--------: | :--------------: | ---------------- | -------------- |
|   gene    | n_samples  | mean_exp_Control | mean_exp_Treated | sd_exp_Control |
|   BRCA1   |     6      |       7.70       | 11.2             | 0.979          |
|   EGFR    |     6      |       8.61       | 5.37             | 0.660          |
|   GAPDH   |     6      |       13.1       | 12.9             | 0.169          |
|   TP53    |     6      |       8.09       | 8.28             | 0.697          |

**Explanation**:
This table summarizes our RNA-Seq experiment across 12 total samples (6 Control vs. 6 Treated) for four distinct genes. Since the data is in $\[\log _{2}\]$ TPM, a 1-unit difference in `mean_exp` equates to a **2-fold change** in biological expression. 
1. **BRCA1 (Strong Up-regulation)**:
   - **Control**: 7.70 $\[\rightarrow \]$ **Treated**: 11.23.
   - **Meaning**: The mean expression increased by roughly 3.5 units. Because $\(2^{3.5} \approx 11.3\)$, `BRCA1` is expressed over **11 times higher** in our treated samples compared to controls.
2. **EGFR (Strong Down-regulation)**: 
   - **Control**: 8.61 $\[\rightarrow \]$ **Treated**: 5.37.
   - Meaning: Expression dropped by more than 3 units. This indicates that treatment significantly **suppresses** `EGFR` expression.
3. **GAPDH (Stable Housekeeping Control)**:
   - **Control**: 13.06 $\[\rightarrow \]$ **Treated**: 12.92.
   - **Meaning**: The values are nearly identical, and the Standard Deviation (`sd_exp_Control` = 0.169) is exceptionally low. This confirms it is an excellent control gene because the treatment did not alter it. 
4. **TP53 (Unchanged Variable Biomarker):**
   - **Control**: 8.09 $\[\rightarrow \]$ Treated: 8.28.
   - **Meaning**: The minor increase is just random background statistical noise (within the range of normal variation), meaning `TP53` expression is unaffected by this treatment.