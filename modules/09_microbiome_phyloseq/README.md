# Module 09: Microbiome & Metagenomics Pipeline using phyloseq
## Description: End-to-end 16S / Metagenomics workflow: object assembly, alpha diversity profiling, Bray-Curtis ordination, and taxonomy.

Microbiome profiling (16S rRNA gene sequencing, ITS, or shotgun metagenomics) evaluates microbial community structure: identifying taxonomic composition ("who is there") and relative abundance ("how many are there").
In R, the standard infrastructure is built on the `phyloseq` S4 object, which unifies four core biological data tables into a single synchronized container:
- `otu_table`: Matrix of count abundances (Taxa/ASVs $\times$ Samples).
- `tax_table`: Taxonomic assignment matrix (Kingdom $\rightarrow$ Species).
- `sample_data`: Experimental phenotype metadata (`data.frame`).
- `phy_tree`: Phylogenetic tree connecting ASV sequences (`phylo` object).


### Core phyloseq Function Reference
|         Function          |                     Primary Purpose                      |                         Key Parameters                          |
| :-----------------------: | :------------------------------------------------------: | :-------------------------------------------------------------: |
|        phyloseq()         | Constructor uniting all individual microbiome components |       otu_table(), tax_table(), sample_data(), phy_tree()       |
|        otu_table()        |           Instantiates abundance count matrix            |                  object, taxa_are_rows = TRUE                   |
|        tax_table()        |     Instantiates matrix of taxonomy classifications      |                 matrix with named column ranks                  |
|       sample_data()       |             Instantiates metadata container              |          data.frame with rownames matching sample IDs           |
|    estimate_richness()    |            Calculates alpha diversity metrics            | physeq, measures = c("Observed", "Shannon", "Simpson", "Chao1") |
| transform_sample_counts() |     Applies mathematical functions to sample columns     |              physeq, fun = function(x) x / sum(x)               |
|        distance()         |        Computes ecological dissimilarity matrices        |         physeq, method = "bray" / "unifrac" / "jaccard"         |
|        ordinate()         |  Performs dimension reduction on dissimilarity matrices  |       physeq, method = "PCoA" / "NMDS", distance = "bray"       |
|        tax_glom()         |   Aggregates taxa counts to a specified taxonomic rank   |                   physeq, taxrank = "Phylum"                    |


Instead of looking at a single genome (like human or mouse), metagenomics looks at the DNA of an entire community of microorganisms (bacteria, fungi, viruses) sampled from an environment (e.g., the human gut, soil, or ocean water).

The script outlines an end-to-end microbiome processing workflow:
- **Object Assembly**: High-throughput sequencing generates millions of reads. After a pipeline identifies which bacteria are present and counts them (creating an OTU/ASV table), we need to stitch this table together with sample metadata (e.g., Patient Age, Disease vs. Healthy) and a phylogenetic tree. This combined structure is called a **Phyloseq Object**.
- **Alpha Diversity Profiling**: This measures diversity within a single sample (How rich is this ecosystem?). For example, a healthy gut usually has high alpha diversity (many different bacterial species), whereas a gut treated with antibiotics might show low alpha diversity (dominated by just one or two resilient species).
- **Bray-Curtis Ordination (Beta Diversity)**: This measures diversity between different samples (How different is Sample A from Sample B?). Bray-Curtis dissimilarity is a mathematical formula that looks at both the presence and abundance of species. Ordination plots these differences on a 2D grid (like a PCA or PCoA plot). If our "Treated" samples cluster on the left and "Control" samples cluster on the right, it proves the treatment significantly altered the microbiome structure.
- **Taxonomy**: This tracks the family tree of the microbes. It lets you break down the samples by Kingdom, Phylum, Class, Order, Family, Genus, or Species to see exactly which bugs are shifting (e.g., an increase in inflammatory Proteobacteria).

**The Role of Each Package**
- `phyloseq`: The undisputed standard framework in R for handling microbiome data. It acts as a wrapper that binds our abundance tables, taxonomic classifications, metadata, and phylogenetic trees into one easily manipulable object.
- `dplyr`: A data manipulation library used to clean metadata, filter out low-abundance taxa, or group data by experimental cohorts.


1. **Simulating the Multi-Table Core**
```r
asv_counts <- matrix(rpois(num_asvs * num_samples, lambda = 20), ...)
asv_counts[1:10, 7:12] <- asv_counts[1:10, 7:12] + rpois(10 * 6, lambda = 100)
```
- It creates a count matrix of 60 ASVs (Amplicon Sequence Variants) across 12 samples. It deliberately seeds a biological signal by injecting a massive spike of extra bacterial counts (lambda = 100) into the Treatment group (Samples 7 to 12) for the first 10 ASVs.
- Purpose: In a real-world pipeline, we extract these from raw sequencing data using algorithms like DADA2 or Deblur. ASVs have completely replaced the older concept of OTUs (Operational Taxonomic Units) because ASVs track exact sequences down to single-nucleotide differences, providing much cleaner biological resolution.

2. **Assembling the `phyloseq` S4 Container**
```r
ps <- phyloseq(
    otu_table(asv_counts, taxa_are_rows = TRUE),
    tax_table(tax_matrix),
    sample_data(sample_metadata)
)
```

- It binds three completely different structural layouts—a raw count matrix (`otu_table`), a taxonomic dictionary tracking names (`tax_table`), and a clinical metadata spreadsheet (`sample_data`)—into a unified data class called an S4 object.
- Purpose: Microbiome data is notoriously hard to manage because if we remove a single sample or rename a bacterium in one spreadsheet, all our tables lose alignment. The `phyloseq` object works like an automated safe-deposit box: if we filter out low-quality samples, `phyloseq` automatically removes those exact same samples from the metadata and count matrices simultaneously, preventing indexing errors.


3. **Alpha Diversity Metrics (Within-Sample Diversity)**
```r
alpha_df <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
```
- It calculates how rich and balanced the microbial ecosystem is inside each distinct sample box.
	- Observed: Simply counts how many unique bacterial types are present.
	- Shannon / Simpson: Factoring in evenness. If a sample has 50 species but 99% of the population belongs to just one invasive species, these mathematical metrics drop significantly to reflect that the ecosystem is highly unstable.
- Bioinformatics Purpose: Because our script injected massive amounts of new reads into the Treatment group, the Treatment group will show highly inflated microbial counts, making it distinctly stand out during boxplot evaluations.

4. **Normalization and Beta Diversity (Between-Sample Comparison)**
```r
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
bray_pcoa <- ordinate(ps_rel, method = "PCoA", distance = "bray")
```
- What it does:
	1. `transform_sample_counts` scales the raw integer sequences into fractions/percentages between 0.0 and 1.0.
	2. `ordinate` measures the similarity landscape using Bray-Curtis distance and flattens that multi-dimensional variance landscape onto a 2D PCoA grid (Principal Coordinate Analysis).
- Purpose: Raw DNA sequencing output fluctuates wildy between samples simply due to machine loading errors. Transforming counts to Relative Abundance levels scales the profiles to a common baseline. Because our script created a drastic count difference between Control and Treatment groups, our PCoA plot will show a distinct separation, with the 6 Control dots clustering tightly together on one side and the 6 Treatment dots grouping on the opposite side.

5. **Taxonomic Aggregation** (`tax_glom`)
```r
ps_phylum <- tax_glom(ps_rel, taxrank = "Phylum")
```
- It collapses (agglomerates) the 60 raw individual ASV rows down by summing up their relative abundance metrics based entirely on their high-level classification category (`Phylum`).
- Purpose: Examining 60 to 10,000 individual ASVs on a single bar chart is impossible to visually interpret. Aggregating the data up to the Phylum level condenses the biological view so a scientist can quickly see broad shifts (e.g., whether a drug treatment caused an expansion of inflammatory Proteobacteria relative to healthy Firmicutes).

## Explanation of the output:
```bash
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 60 taxa and 12 samples ]
sample_data() Sample Data:       [ 12 samples by 3 sample variables ]
tax_table()   Taxonomy Table:    [ 60 taxa by 6 taxonomic ranks ]
```

**The Data Blueprint (`print(ps)`)**
The summary confirms our **S4 phyloseq object** assembled perfectly. It has synchronized 60 distinct bacterial types (taxa/ASVs) across 12 individual samples. The structural links are completely intact, ensuring that any downstream filtering will safely update all tables at once.

```bash
          Observed Shannon  Simpson
Sample_01       60 4.067311 0.9824729
Sample_02       60 4.074729 0.9826639
Sample_03       60 4.069010 0.9824828
Sample_04       60 4.061285 0.9821783
```

### Plot 1: Microbial Alpha Diversity Boxplot (p_alpha)
Our data snapshot shows that the Control samples (Sample_01 to Sample_04) have an Observed richness of exactly 60, a Shannon index around 4.06, and a Simpson index around 0.98.
- Observed Richness Metric: Because the simulation used a random distribution (rpois) where every single ASV received at least a few reads, all 60 bacteria are detected in every sample. As a result, the "Observed" boxplot will look completely flat at the 60 line for both the Control and Treatment groups.
- Shannon & Simpson Indices (Evenness): The Control samples have a very high Shannon (~4.06) and Simpson (~0.98) score because the baseline counts are highly uniform and evenly balanced. However, because the code injected a massive count spike (lambda = 100) into the first 10 ASVs for the Treatment group, the Treatment boxplot will drop significantly. Those 10 spiked bacteria now completely dominate the ecosystem, destroying the biological "evenness" and causing the diversity indices to plummet.

### Plot 2: PCoA Ordination Plot (`p_beta`)
This plot calculates the Bray-Curtis dissimilarity to map how similar the whole community structures are to one another.
- Visual Layout: Our plot will show a stark, undeniable visual separation. The 6 Control samples (Samples 1–6) will form a tight cluster on one side of the chart because their microbial communities are nearly identical. The 6 Treatment samples (Samples 7–12) will cluster tightly on the exact opposite side.
- Mathematical Cause: The primary horizontal axis (Axis 1) is capturing the massive, synthetic variation introduced by our treatment grouping, proving that the treatment completely remodeled the microbiome landscape.

### Plot 3: Phylum-Level Relative Abundance Composition (`p_bar`)
This stacked bar chart displays what percentage of the total community belongs to the four simulated phyla (Bacteroidetes, Firmicutes, Proteobacteria, and Actinobacteria).
- Control Facet: The 6 bars under the Control group will look incredibly uniform and balanced. Each phylum will occupy roughly equal vertical space because their initial generation counts were completely random and evenly spread.
- Treatment Facet: The 6 bars under the Treatment group will look drastically distorted. Whichever phyla were randomly assigned to the first 10 ASVs will visually balloon to take up the vast majority of the bar space, shifting the relative abundance totals and squishing the remaining non-spiked phyla down into thin slivers.










































