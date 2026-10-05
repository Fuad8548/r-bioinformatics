suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(Matrix)
})

set.seed(42)

# ==============================================================================
# 1. Simulate Realistic Single-Cell Count Matrix with Mitochondrial Genes
# ==============================================================================
cat("\n--- Step 1: Simulating scRNA-seq Sparse Count Matrix ---\n")

num_genes <- 1000
num_cells <- 800

# Create realistic gene names including Mitochondrial (MT-) genes
regular_genes <- paste0("GENE-", sprintf("%04d", 1:(num_genes - 10)))
mt_genes <- paste0(
    "MT-",
    c(
        "ND1", "ND2", "CO1", "CO2", "ATP6", "ND4",
        "CYB", "ND5", "CO3", "ATP8"
    )
)
all_gene_names <- c(regular_genes, mt_genes)
cell_barcodes <- paste0("CELL_", sprintf("%04d", 1:num_cells))

# Generate sparse expression matrix with zero inflation
raw_counts <- matrix(rpois(num_genes * num_cells, lambda = 0.3),
    nrow = num_genes, ncol = num_cells
)
# Introduce synthetic cluster signal
raw_counts[1:50, 1:250] <- raw_counts[1:50, 1:250] +
    rpois(50 * 250, lambda = 3)
raw_counts[51:100, 251:500] <-
    raw_counts[51:100, 251:500] +
    rpois(50 * 250, lambda = 4)
raw_counts[101:150, 501:800] <- raw_counts[101:150, 501:800] +
    rpois(50 * 300, lambda = 3.5)

# High mitochondrial content in a subset of dying cells (for QC demo)
raw_counts[991:1000, 750:800] <- raw_counts[991:1000, 750:800] +
    rpois(10 * 51, lambda = 12)

rownames(raw_counts) <- all_gene_names
colnames(raw_counts) <- cell_barcodes

# Convert to Matrix dgCMatrix (sparse format)
sparse_counts <- as(raw_counts, "dgCMatrix")


# ==============================================================================
# 2. Create Seurat Object & Calculate QC Metrics
# ==============================================================================
cat("\n--- Step 2: Initializing Seurat Object & Quality Control ---\n")

pbmc <- CreateSeuratObject(
    counts       = sparse_counts,
    project      = "scRNA_Tutorial",
    min.cells    = 3,
    min.features = 100
)

# Calculate percentage of mitochondrial transcripts
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

cat("Initial cells:", ncol(pbmc), "| Initial genes:", nrow(pbmc), "\n")

# Filter out low-quality barcodes
pbmc <- subset(
    pbmc,
    subset = nFeature_RNA > 150 & nFeature_RNA < 800 & percent.mt < 10
)

cat("Cells remaining after QC filtering:", ncol(pbmc), "\n")


# ==============================================================================
# 3. Normalization, Feature Selection & Scaling
# ==============================================================================
cat("\n--- Step 3: Normalization & Variable Feature Selection ---\n")

# Log-normalization
pbmc <- NormalizeData(
    pbmc,
    normalization.method = "LogNormalize",
    scale.factor         = 10000
)

# Identify top highly variable genes (HVGs)
pbmc <- FindVariableFeatures(
    pbmc,
    selection.method = "vst",
    nfeatures        = 300
)

# Scale data across all genes
all_genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all_genes)


# ==============================================================================
# 4. Dimensionality Reduction (PCA & UMAP) & Clustering
# ==============================================================================
cat("\n--- Step 4: PCA, Louvain Clustering & UMAP ---\n")

# Linear dimension reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)

# Construct KNN graph & cluster cells (Louvain algorithm)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Non-linear dimension reduction for 2D visualization
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)


# ==============================================================================
# 5. Find Cluster Marker Genes & Visualization
# ==============================================================================
cat("\n--- Step 5: Identifying Cluster Markers ---\n")

# Find markers for every cluster compared to all remaining cells
cluster_markers <- FindAllMarkers(
    pbmc,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
)

top3_markers <- cluster_markers %>%
    group_by(cluster) %>%
    slice_max(n = 3, order_by = avg_log2FC)

print(head(top3_markers, 10))


# ==============================================================================
# 6. Generate Diagnostic & Publication Plots
# ==============================================================================
cat("\n--- Step 6: Rendering Visualizations ---\n")

# 6A. UMAP Plot colored by cluster
p_umap <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.8) +
    labs(title = "Single-Cell Clustering (UMAP)")

# 6B. Feature Plot of top marker gene on UMAP
top_gene <- top3_markers$gene[1]
p_feat <- FeaturePlot(pbmc, features = top_gene) +
    labs(title = paste("Expression of Marker:", top_gene))

# 6C. Violin Plot of marker expression across clusters
p_vln <- VlnPlot(pbmc, features = top_gene)

print(p_umap)
print(p_feat)
print(p_vln)

cat("\n--- Saving Plots to Local Directory ---\n")

# 1. Save the UMAP Cluster Plot (Square dimensions work best for UMAPs)
ggsave(
    filename = "scRNA_umap_clusters.pdf",
    plot     = p_umap,
    device   = "pdf",
    width    = 7,
    height   = 6,
    units    = "in",
    dpi      = 300
)

# 2. Save the Feature Plot (Tracks spatial gene expression)
ggsave(
    filename = "scRNA_marker_feature_plot.pdf",
    plot     = p_feat,
    device   = "pdf",
    width    = 7,
    height   = 6,
    units    = "in",
    dpi      = 300
)

# 3. Save the Violin Plot (Wider dimensions allow clusters to breath horizontally)
ggsave(
    filename = "scRNA_marker_violin_plot.pdf",
    plot     = p_vln,
    device   = "pdf",
    width    = 8,
    height   = 5,
    units    = "in",
    dpi      = 300
)

cat("All single-cell plots saved successfully as high-resolution PDFs!\n")
