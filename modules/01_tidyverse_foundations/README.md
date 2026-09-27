# Data reshaping (tidyr), manipulation (dplyr) and plotting (ggplot2) using synthetic biological expression and metadata tables.


### Summary
While the data itself was randomly generated using R's `rnorm()` function to test our code, **the biological logic holds up perfectly**: our **Treatment** (let's say chemotherapy) acted exactly like an aggressive cancer intervention—shutting down a primary growth driver (`EGFR`) while triggering a massive DNA-damage distress response (`BRCA1`), all while leaving basic cellular survival genes (`GAPDH`, `TP53`) untouched.


## Exercises:

**1: Filtering & Quality Control**
Using `tidy_genomics_df`:

- Filter the dataset to include only samples where rin_score is strictly greater than `8.0` and `cell_line` is "`HeLa`".
- Select only the columns `sample_id`, `gene`, `log2_tpm`, and `condition`.
- Sort the resulting table in descending order of `log2_tpm`.


**Exercise 2: Fold Change & Normalization Calculations**
Housekeeping genes like `GAPDH` should remain constant across conditions.

- Group `tidy_genomics_df` by `sample_id` and compute a normalized expression column `relative_to_gapdh` by subtracting each sample's `GAPDH` expression value from all other genes in that sample ($log_2(A) - log_2(B) = log_2(A/B)$).
- Filter out `GAPDH` from the output.
- Calculate the average `relative_to_gapdh` for each remaining gene split by `condition`.


**Exercise 3: Faceted Visualization**
Create a publication-quality `ggplot2` visualization:

- Plot `log2_tpm` (y-axis) against `condition` (x-axis) for all genes.
- Add a `geom_violin(alpha = 0.5)` layer with individual sample points overlaid using `geom_point()`.
- Use `facet_wrap(~ gene, scales = "free_y")` to render each gene in its own panel.
- Apply the custom `theme_bio()` theme and export the figure to `outputs/figures/m01_exercise_violin_facets.png`.


**Exercise 4: Matrix Reconstruction**
- Take the filtered dataset from Exercise 2 (excluding `GAPDH`).
- Reshape it back into a wide expression matrix where rows are `gene` names and columns are `sample_ids`, with `relative_to_gapdh` as cell values.
- Convert the resulting tibble into a standard R matrix with row names matching the gene identifiers.

