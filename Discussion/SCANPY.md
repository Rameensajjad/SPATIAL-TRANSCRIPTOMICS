
# Scanpy Analysis and Visualization of Spatial Transcriptomics Data

---

## Introduction

This analysis presents a comprehensive workflow for processing, analyzing, and visualizing spatial transcriptomics data using the Scanpy framework in Python. The primary dataset used is a 10x Genomics Visium spatial transcriptomics sample derived from a **human lymph node** tissue section, publicly available through the 10x Genomics data portal. A secondary analysis using **MERFISH** data is also presented to demonstrate that the same analytical framework generalizes across different spatial transcriptomics technologies.

Spatial transcriptomics extends conventional single-cell RNA sequencing by preserving the **physical location** of gene expression measurements within tissue. This enables not only the identification of transcriptionally distinct cell populations, but also the mapping of those populations to specific anatomical regions — a capability that is central to understanding tissue organization, cellular communication, and disease microenvironments.

The computational environment relied on four core libraries: **Scanpy** as the primary analysis engine, **AnnData** as the underlying data structure, **Matplotlib** and **Seaborn** for visualization, and **Pandas** for data manipulation. Global figure parameters were configured at the outset — white background, 8×8 inch default figure size — and Scanpy's verbosity was set to level 3 to ensure detailed progress logging throughout the pipeline.

---

## 1. Data Loading

The human lymph node Visium dataset was loaded using Scanpy's built-in dataset loader, which downloads the data directly from 10x Genomics and returns it in AnnData format. The AnnData object is the central data structure of the entire workflow, storing the gene expression count matrix, spatial barcoded spot coordinates, the high-resolution H&E tissue image, and all associated metadata in a single unified object.

Upon loading, gene names were made unique to prevent downstream identifier conflicts. Mitochondrial genes were then flagged by identifying all gene names beginning with the prefix "MT-". Quality control metrics were computed across all spots, yielding per-spot statistics including total counts, number of detected genes, log-transformed equivalents of these metrics, and the percentage of reads attributable to mitochondrial genes.

The resulting AnnData object comprised **4,035 spots and 36,601 genes**, with comprehensive observation-level and gene-level metadata populated across its annotation slots, and spatial image data stored in the unstructured metadata layer.

---

## 2. Quality Control and Preprocessing

### 2.1 Quality Control

Prior to any analytical steps, a thorough quality control assessment was performed to identify and remove unreliable spots and uninformative genes. Four histograms were generated in a single figure to inspect the distribution of total counts and detected genes per spot. Two of the four panels displayed full distributions, while the remaining two displayed zoomed-in lower-range views of the same metrics — a technique that reveals fine-grained structure that would otherwise be obscured at full scale. These visualizations informed the selection of appropriate filtering thresholds.

The following four filtering steps were applied sequentially:

- Spots with **fewer than 5,000 total counts** were removed, as these represent sparsely captured locations likely reflecting poor tissue coverage or technical failure
- Spots with **more than 35,000 total counts** were removed, as these may represent doublets or sequencing artifacts
- Spots in which **more than 20% of total counts** originated from mitochondrial genes were removed, as elevated mitochondrial read proportions are indicative of damaged or necrotic tissue
- Genes detected in **fewer than 10 spots** were removed, as these are too rare across the dataset to contribute reliable signal

After filtering, **3,861 spots** were retained and approximately 16,916 low-information genes were eliminated, resulting in a substantially leaner and higher-quality dataset.

### 2.2 Normalization and Feature Selection

Total-count normalization was applied across all spots so that each location was brought to a comparable expression scale, removing the confounding effects of variable sequencing depth. This was followed by a **log-plus-one transformation**, which compresses extreme expression values and stabilizes variance across the dynamic range of the data.

The **top 2,000 highly variable genes** were then selected using the Seurat method. These are genes that exhibit the greatest transcriptional variability across spots relative to their mean expression level, and they serve as the informative feature set for all subsequent dimensionality reduction, clustering, and visualization steps.

---

## 3. Dimensionality Reduction and Clustering

The standard Scanpy single-cell analysis pipeline was applied to embed and cluster the transcriptional manifold.

**Principal Component Analysis (PCA)** was first performed, compressing the high-dimensional gene expression matrix — spanning thousands of genes — down to 50 principal components. These components capture the dominant axes of transcriptional variation across the dataset while discarding noise.

A **neighbor graph** was then constructed using the PCA coordinates. In this graph, spots that are transcriptionally similar — close together in PCA space — are connected as neighbors. This graph structure underpins both the clustering algorithm and the UMAP embedding.

**UMAP (Uniform Manifold Approximation and Projection)** was computed to generate a 2D representation of the transcriptional landscape for visual exploration. UMAP attempts to preserve local neighborhood relationships from high-dimensional space in its 2D layout.

**Leiden clustering** was applied to the neighbor graph to detect communities of transcriptionally similar spots. The algorithm identified **10 distinct clusters** in this dataset, each assigned a numerical label stored as per-spot metadata.

To verify that the clustering structure reflects biological signal rather than technical artifacts, three UMAP plots were generated side by side — colored by total counts, by number of detected genes, and by cluster identity respectively. The absence of strong technical gradients aligning with cluster boundaries confirmed that the identified clusters are driven by genuine transcriptional differences.

---

## 4. Visualization in Spatial Coordinates

A defining feature of spatial transcriptomics analysis is the ability to project analytical results back onto the physical tissue. Using Scanpy's `spatial()` plotting function, colored circular spots were overlaid directly onto the **Hematoxylin and Eosin (H&E) stained tissue image** at high resolution.

The `spatial()` function offers several parameters that were actively used throughout this section:

- **`img_key`** — specifies which image resolution to render from the stored image data
- **`crop_coord`** — defines a rectangular crop region using pixel coordinates (left, right, top, bottom) for zooming into regions of interest
- **`alpha_img`** — controls the transparency of the tissue image beneath the spots
- **`bw`** — converts the tissue image to grayscale when needed
- **`size`** — acts as a scaling factor for spot radius rather than an absolute size

**Total counts and gene counts** were first mapped onto the tissue to assess whether sequencing depth varied spatially across the section — an important check for technical bias before interpreting biological patterns.

**Leiden cluster labels** were then overlaid on the tissue, revealing a critical result: spots belonging to the same transcriptional cluster tend to co-localize in space. This demonstrates that the gene expression-based clusters correspond to actual anatomical regions within the lymph node. A notable spatial relationship was observed between clusters 5 and 0 — spots from cluster 5 were frequently surrounded by cluster 0 spots — suggesting a spatial adjacency between two distinct cell populations in the tissue microenvironment.

A **zoomed-in crop** of clusters 5 and 9 was produced by specifying pixel crop coordinates and highlighting only those two cluster groups. Spot transparency was reduced to 50% so that the underlying H&E tissue morphology remained visible beneath the expression overlay, enabling qualitative interpretation of the anatomical context.

---

## 5. Cluster Marker Gene Analysis

To characterize each cluster biologically, a **t-test-based differential expression analysis** was performed, ranking genes by the specificity of their expression in each cluster relative to all other spots. This produced, for each cluster, a ranked list of marker genes along with associated log fold changes, p-values, and adjusted p-values.

The analysis focused on **cluster 9**, for which a heatmap was generated displaying the top 10 marker genes across all clusters. The clusters were arranged along the heatmap's axis according to a **dendrogram** computed from PCA coordinates, which hierarchically groups clusters by overall transcriptional similarity and provides additional context for interpreting relationships between populations.

The most notable finding from this analysis was the gene **CR2**. When its expression was mapped spatially on the tissue section — displayed alongside the cluster label overlay for direct comparison — the CR2 expression pattern closely recapitulated the spatial distribution of cluster 9. This gene-to-space correspondence validates that the transcriptional cluster identity has a direct anatomical footprint in the lymph node tissue, and demonstrates the power of combining differential expression analysis with spatial visualization.

Two additional genes — **COL1A2** and **SYPL1** — were mapped spatially at 70% spot transparency to illustrate how individual gene expression can be explored independently of cluster-level assignments, enabling finer-grained spatial interrogation of the tissue.

---

## 6. MERFISH Secondary Analysis

To demonstrate the generalizability of the Scanpy spatial workflow, the pipeline was extended to a **MERFISH dataset** from Xia et al. (2019), which profiled gene expression in cultured **U2-OS cells** — a human osteosarcoma cell line. MERFISH is a fluorescence in situ hybridization-based technology that measures spatial gene expression at single-cell resolution but without an accompanying tissue image.

Since MERFISH data has no dedicated Scanpy loader, the data was ingested manually. A coordinate table was read from an Excel file and a counts matrix was read from a CSV file. The counts matrix was transposed to orient cells as rows and genes as columns, and the two tables were aligned by matching cell indices. Spatial coordinates were then manually assigned to the AnnData object's spatial embedding slot.

A preprocessing pipeline was applied with parameters appropriate for this smaller and simpler dataset: per-cell normalization scaling to one million counts, log-plus-one transformation, PCA reduced to 15 components, neighbor graph construction, UMAP, and Leiden clustering at a resolution of 0.5. The algorithm identified **6 clusters** across **645 cells and 12,903 genes**.

The biological interpretation here differs from the Visium analysis. Since all cells in this experiment are of the same type (U2-OS), the 6 clusters do not represent distinct cell types. Rather, they reflect **cells at different stages of the cell cycle**, which is a primary driver of transcriptional variability in monotypic cell culture systems.

Two visualizations were produced: a UMAP colored by cluster to show the transcriptional structure, and a spatial embedding plot using the MERFISH physical coordinates as axes, also colored by cluster. As expected from a cell culture experimental setup — where cells are distributed randomly rather than arranged in tissue — no meaningful spatial organization was observed in the spatial plot. The clusters appeared scattered without anatomical structure, which is fully consistent with the experimental design.

---

## Summary and Conclusions

This analysis demonstrated a complete spatial transcriptomics workflow from raw data ingestion through quality control, normalization, dimensionality reduction, unsupervised clustering, spatial visualization, and marker gene identification.

The central finding from the Visium human lymph node analysis is that **transcriptional clusters identified purely from gene expression space have clear and reproducible spatial identities within the tissue**. Clusters do not merely reflect abstract molecular groupings — they correspond to physically distinct regions of the lymph node, as confirmed both by direct spatial overlay visualization and by the spatial mapping of individual marker genes. The gene CR2 in particular provided a clean, spatially coherent validation of cluster 9's anatomical footprint.

The MERFISH secondary analysis confirmed that the Scanpy spatial framework accommodates alternative spatial technologies with straightforward adaptations to the data loading and preprocessing steps. The pipeline is modular and transferable across platforms, making it a robust foundation for spatial transcriptomics research.

The MERFISH secondary analysis confirmed that the Scanpy spatial framework accommodates alternative spatial technologies with straightforward adaptations to the data loading and preprocessing steps. The pipeline is modular and transferable across platforms, making it a robust foundation for spatial transcriptomics research.
