# 2. Reshaping Data with tidyr (Wide <-> Long)

```r
# Pivot wide matrix into tidy/long format (essential for ggplot2 & Bioconductor)
expression_long <- expression_wide |>
  pivot_longer(
    cols      = -sample_id,
    names_to  = "gene",
    values_to = "log2_tpm"
  )

# Merge expression values with sample metadata
tidy_genomics_df <- expression_long |>
  inner_join(metadata, by = "sample_id")
```

**Explanation:**
This step is mandatory in modern bioinformatics workflows because data visualization tools like `ggplot2` and downstream statistical models require data to be grouped by features rather than spread across wide matrix columns.
1. Converting "Wide" to "Long" Format (`pivot_longer`)
`expression_wide` was structured like a spreadsheet matrix: each gene was its own column (`TP53`, `BRCA1`, etc.). While humans read this easily, plotting software struggles with it because the genes are separated.
The `pivot_longer()` function stacks those individual gene columns into a single, structured column.
How the code parameters work:
   - `cols = -sample_id`: This tells R to pivot every single column except `sample_id`. It leaves the sample names alone and collapses the rest.
   - `names_to = "gene"`: This creates a brand new column named "`gene`" to hold the names of the columns that are being collapsed (`TP53`, `BRCA1`, `EGFR`, `GAPDH`).
   - `values_to = "log2_tpm"`: This creates a new column to store the actual numerical expression numbers associated with those genes.

  The Transformation Visualized:
  Before `expression_wide`:
  | Operation | Description | Header 4 |
  | :-------: | :---------: | :------: |
  | sample_id |    TP53     |  BRCA1   |
  |   SMP_1   |     8.2     |   7.9    |
  |   SMP_2   |     8.6     |   11.4   |

  After `expression_long`: 
  | Operation | Description | Header 4  |
  | :-------: | :---------: | :-------: |
  | sample_id |    gene     | long2_tpm |
  |   SMP_1   |    TP53     |    8.2    |
  |   SMP_1   |    BRCA1    |    7.9    |
  |   SMP_2   |    TP53     |    8.6    |
  |   SMP_2   |    BRCA1    |   11.4    |

2. **Merging Expression with Experimental Metadata (`inner_join`)**
Right now, our gene data (`expression_long`) knows which sample it belongs to (e.g., `SMP_1`), but it doesn't know if that sample was a Control or a Treated cell line. That information is trapped inside our `metadata` data frame.
- `inner_join(metadata, by = "sample_id")`: R looks at both datasets, matches the rows based on the shared `sample_id` column, and merges all the metadata features (`condition`, `cell_line`, `batch`, `rin_score`) directly onto the corresponding expression rows.