# Module 07: Single-Cell RNA-Seq Processing Pipeline using Seurat
## Description: End-to-end scRNA-seq workflow covering QC, Log-Normalization, PCA, Louvain Clustering, UMAP, and Cell Marker Discovery.

In single-cell RNA sequencing (scRNA-seq), instead of measuring the average expression across thousands of cells as in bulk RNA-seq, we profile transcriptomes at individual cell resolution.
The fundamental object in Seurat is the `Seurat` S4 object, which holds the raw count matrix, normalized expression layers, cell metadata (e.g., QC metrics, cluster assignments), and dimensional reduction coordinates (PCA, UMAP, t-SNE).

1. Quality Control & Cell Filtering:
   - Removes dead cells, doublets, and low-quality barcodes.
   - Filter out low-quality cells based on three key thresholds:
     - Low `nFeature_RNA`: Empty droplets or poor sequencing depth.
     - High `nFeature_RNA` / `nCount_RNA`: Cell doublets or multiplets.
     - High percent.mt: Dying/apoptotic cells whose cytoplasmic mRNA leaked out, leaving elevated mitochondrial transcripts.

2. Data Normalization & Variable Feature Selection:
   - `NormalizeData()`: Applies global log-normalization ($log(1 + \frac{\text{gene count}}{\text{total counts}} \times 10000)$) to account for sequencing depth variations across cells.
   - `FindVariableFeatures()`: Identifies high-variance genes (typically top 2,000) that drive biological heterogeneity, ignoring stably expressed housekeeping genes.

3. Scaling & Linear Dimensionality Reduction (PCA):
   - `ScaleData()`: Shifts mean expression to $0$ and variance to $1$ for each gene, ensuring highly expressed genes do not dominate downstream PCA.
   - `RunPCA()`: Compresses high-dimensional space (2,000 genes) into orthogonal principal components (PCs) capturing maximal variance.

4. **Graph-Based Clustering & Non-Linear Visualization** (UMAP):
   - `FindNeighbors()` & `FindClusters()`: Constructs a $K$-nearest neighbor ($KNN$) graph based on top PCA dimensions and uses the Louvain algorithm to group cells into discrete clusters.
   - `RunUMAP()`: Projects the multi-dimensional PCA space into 2D coordinates for visual exploration while preserving cell-cell relationships.

5. **Differential Marker Analysis & Cell-Type Identification:**
- `FindAllMarkers()`: Performs differential expression tests (Wilcoxon rank-sum test by default) comparing each cluster against all other cells to extract cluster-defining marker genes.

### Essential Seurat v5 Function Reference
|             Function             |                      Primary Purpose                      |                 Key Arguments / Parameters                  |
| :------------------------------: | :-------------------------------------------------------: | :---------------------------------------------------------: |
|       CreateSeuratObject()       |    Encapsulates raw count matrix into Seurat container    |               counts, min.cells, min.features               |
|      PercentageFeatureSet()      | Calculates proportions of gene sets (e.g., mitochondrial) |            pattern = "^MT-" or pattern = "^HB-"             |
|         NormalizeData()          |    Accounts for sequencing depth via log-normalization    | normalization.method = "LogNormalize", scale.factor = 10000 |
|      FindVariableFeatures()      |     Extracts highly variable genes for downstream PCA     |         selection.method = "vst", nfeatures = 2000          |
|           ScaleData()            |       Z-score standardizes expression across cells        |  features, vars.to.regress = c("percent.mt", "nCount_RNA")  |
|             RunPCA()             |      Linear dimension reduction on variable features      |              npcs = 50, reduction.name = "pca"              |
| FindNeighbors() / FindClusters() |       Builds KNN graph and runs Louvain clustering        | dims = 1:15, resolution = 0.5 (higher res = more clusters)  |
|            RunUMAP()             |       2D manifold projection for cell visualization       |                dims = 1:15, n.neighbors = 30                |
|         FindAllMarkers()         |      Differential expression testing across clusters      |   only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25   |


This R script implements a standard workflow for single-cell RNA sequencing (scRNA-seq) data analysis using the Seurat ecosystem.
The ultimate goal of this pipeline is to take raw, messy gene expression data from thousands of individual cells, filter out dead or low-quality cells, cluster similar cells together, and identify the specific "marker genes" that define those cell types.


1. **Data Simulation & Sparse Matrix Conversion**

```r
raw_counts <- matrix(rpois(num_genes * num_cells, lambda = 0.3), ...)
# ... Intentionally spiking in cluster and mitochondrial signals ...
sparse_counts <- as(raw_counts, "dgCMatrix")
```
   - The Purpose: ScRNA-seq data is notoriously zero-inflated (sparse) because a sequencing machine cannot capture every single mRNA molecule inside a tiny cell.
   - Significance: Converting the raw matrix into a `dgCMatrix` (sparse format) via the `Matrix` package is a computational necessity. It stores only the non-zero values, reducing the computer memory footprint by up to 90%. The script also intentionally spikes in heavy mitochondrial counts into a subset of cells (cells 750 to 800) to mimic dying, low-quality cells for the next Quality Control stage.

2. **Quality Control (QC) & Filtering**

```r
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 150 & nFeature_RNA < 800 & percent.mt < 10)
```
- The Purpose: To clean the dataset so that downstream clusters represent genuine biological cell types rather than technical artifacts.
- Significance:
  - `nFeature_RNA` counts how many distinct genes are detected in a cell. Cells with very low counts (<150) are usually empty droplets or dead cell fragments. Cells with abnormally high counts (>800) are often "doublets" (two cells accidentally trapped in a single reaction droplet).
  - `percent.mt` calculates the ratio of reads mapping to mitochondrial genes (`MT-`). High mitochondrial content (>10%) indicates a dying cell whose outer membrane has ruptured, leaking its cytoplasmic RNA and leaving behind mostly mitochondrial transcripts.

3. **Normalization & Variance Stabilization**

```r
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 300)
```

- The Purpose: To make gene expression levels directly comparable across cells and isolate the biological signals from background noise.
- Significance:
	- `NormalizeData`: Some cells are sequenced deeper than others simply due to machine variation. This function scales the counts of each cell to a standard total of 10,000 molecules, followed by a $\(log_{2}\)$ transformation to minimize the mathematical dominance of ultra-highly expressed genes.
	- `FindVariableFeatures`: Out of thousands of genes, most are "housekeeping genes" that stay at the same level in every cell. This function flags a specific subset (300 in this script; usually 2,000 in real datasets) that show high variation across cells. These Highly Variable Genes (HVGs) drive the mathematical separation of different cell types.
   [**Housekeeping genes** are genes that turned on all the time in every single cell of an organism, regardless of what that cell's specialized job is. They provide essential baseline references for normalizing high-throughput gene expression data and correcting technical variations across samples]

4. **Linear & Non-Linear Dimensionality Reduction**
```r
pbmc <- ScaleData(pbmc, features = all_genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
```

- The Purpose: Single-cell datasets suffer from the "curse of dimensionality" (thousands of genes acting as independent mathematical dimensions). This step reduces the data to a clean, visualizable space.
- Significance:
	- `ScaleData`: Shifts the mean expression of each gene to 0 and variance to 1. This prevents highly abundant genes from entirely eclipsing low-abundance genes during modeling.
	- `RunPCA` (Linear): Compresses the 300 variable genes into a handful of "Principal Components" (PCs) that capture the absolute core variations of the dataset.
	- `RunUMAP` (Non-Linear): Takes those top 10 PCs and compresses them down further into a 2D scatter plot coordinate layout. UMAP preserves local relationships, meaning cells that are biologically similar end up packed close together on our screen.

5. **Louvain Graph-Based Clustering**
```r
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
```

- The Purpose: To partition individual cells into distinct, discrete groups (clusters) without assuming beforehand what cell types exist.
- Significance: `FindNeighbors` constructs a K-Nearest Neighbor (KNN) graph connecting cells with similar PCA expression profiles. `FindClusters` then applies the Louvain community detection algorithm to optimize the density of connections, grouping heavily interconnected cells into numbered clusters (0, 1, 2, ...). The resolution parameter controls granularity: higher values yield more clusters, while lower values yield fewer, larger clusters.

6. **Biomarker Identification & Diagnostics**
```r
cluster_markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
```

- The Purpose: To uncover the unique genetic signature of each newly generated cluster so that a biologist can identity what specific cell type each cluster represents (e.g., T-cells, B-cells, Monocytes).
- Significance: `FindAllMarkers` acts as an automated differential expression machine. It compares Cluster 0 against all other clusters combined, then repeats this for Cluster 1, and so on.
	- `only.pos` = TRUE ignores genes that are downregulated, isolating genes that are explicitly turned on uniquely in that group.
	- `min.pct = 0.25` speeds up calculations by only testing genes that are expressed in at least 25% of the cells within that target cluster.

### Output Explanations:
```bash
   p_val avg_log2FC pct.1 pct.2     p_val_adj    cluster gene     
1     0       3.86 0.968 0.224         0 0       GENE-0045
2     0       3.85 0.976 0.238         0 0       GENE-0006
3     0       3.84 0.948 0.224         0 0       GENE-0009
4     0       4.00 0.976 0.218         0 1       GENE-0084
5     0       4.00 0.98  0.226         0 1       GENE-0059
6     0       3.97 0.992 0.234         0 1       GENE-0055
7     0       3.97 0.98  0.21          0 2       GENE-0112
8     0       3.90 0.968 0.224         0 2       GENE-0141
9     0       3.88 0.984 0.244         0 2       GENE-0137
```

Breaking down what the numerical output table means using our top marker, `GENE-0045`:
- `cluster = 0`: This gene is flagged as a signature marker specifically for Cluster 0.
- `avg_log2FC = 3.86`: The expression of this gene is heavily upregulated (by a $\(log_{2}\)$ factor of 3.86) in Cluster 0 compared to all other cells combined.
- `pct.1 = 0.968`: 96.8% of the cells inside Cluster 0 express this gene.
- `pct.2 = 0.224`: Only 22.4% of the cells outside Cluster 0 express this gene.
- `p_val_adj = 0`: The difference is highly statistically significant.
   
**Plot 1: Single-Cell Clustering (`DimPlot`)**
- This plot maps the global landscape of our single-cell dataset. It has successfully grouped our 800 simulated cells into three distinct, non-overlapping clusters (Cluster 0 in red/salmon, Cluster 1 in green, and Cluster 2 in blue) based on their shared transcriptional profiles.
- Why we need it: It provides cellular identity and context. Without this plot, we have no way of knowing how many cell types or states exist in our sample, how crisp the boundaries between them are, or how many cells belong to each group. It acts as the "map" of the tissue.

**Plot 2: Expression Feature Plot (`FeaturePlot`)**
- This plot maps the expression of `GENE-0045` back onto the spatial UMAP coordinates. We can see a dense cluster of dark purple dots concentrated heavily on the left island (Cluster 0). The top-right (Cluster 1) and bottom-right (Cluster 2) islands remain mostly light gray, showing zero or trace expression.
- Why we need it: It provides spatial validation of markers. It bridges the raw statistics from our table to the UMAP coordinates, visually confirming that the gene's activity perfectly aligns with the physical boundaries of Cluster 0.


**Plot 3: Marker Violin Plot (`VlnPlot`)**
- This plot shows the exact expression distribution of `GENE-0045` across the three identities. Cluster 0 shows a high, bulbous density cloud hovering between expression levels 3 and 5. Meanwhile, Cluster 1 and Cluster 2 show flat lines tightly resting on the 0 baseline, with only a handful of minor outlier points.
- Why we need it: It provides expression resolution and continuous granularity. While a `FeaturePlot` can hide overlapping cells or background noise, the `VlnPlot` directly displays the exact range, shape, density, and baseline expression of every cell inside those clusters. It visually proves that the background expression (`pct.2 = 0.224`) is extremely weak, validating `GENE-0045` as a highly specific, clean biomarker.

**Summary: Why do we need all 3 together?**
|    Plot     |         Core Question it Answers         |                      Biological Insight                       |
| :---------: | :--------------------------------------: | :-----------------------------------------------------------: |
|   DimPlot   |  What groups of cells are in my sample?  |            Identifies the cell populations present            |
| FeaturePlot |        Where is my marker active?        | Confirms if a marker physically matches a specific cell group |
|   VlnPlot   | How strongly and cleanly is it expressed |    Shows the background noise versus true signal threshold    |






