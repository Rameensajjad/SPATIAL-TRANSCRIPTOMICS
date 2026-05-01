---

# Analysis of Visium H&E Data Using Squidpy


---

## Introduction

This analysis presents a comprehensive spatial transcriptomics workflow applied to a **Hematoxylin and Eosin (H&E) stained Visium dataset** of a coronal section of the mouse brain, using the Squidpy framework. Unlike the companion fluorescence Visium tutorial — which focuses on multi-channel image segmentation and intensity-based feature extraction — this analysis places emphasis on **spatial graph construction, neighborhood statistics, ligand-receptor interaction analysis, and spatially variable gene detection**. Together, these methods move beyond simple visualization to extract quantitative biological insights about tissue organization, intercellular communication, and spatially patterned gene expression.

The dataset is publicly available through the 10x Genomics data portal and is provided here in a pre-processed, pre-annotated form. The preprocessing pipeline follows the standard Scanpy spatial transcriptomics workflow. Cluster annotations were established using several neuroanatomical references including the **Allen Brain Atlas**, the **Mouse Brain gene expression atlas** from the Linnarsson lab, and related literature.

The computational environment uses **NumPy** and **Pandas** for numerical and tabular data handling, **Scanpy** for gene expression analysis, **Squidpy** for spatial graph analysis and image feature extraction, and **AnnData** as the unified data container.

---

## 1. Data Loading and Initial Visualization

Two objects are loaded at the start:

- **`adata`** — a pre-processed AnnData object containing normalized gene expression counts and pre-annotated cluster labels for **2,688 Visium spots**
- **`img`** — a Squidpy `ImageContainer` object holding the high-resolution H&E tissue image

The pre-annotated cluster labels — assigned using established neuroanatomical resources — are immediately visualized in spatial context using `squidpy.pl.spatial_scatter()`. This initial map provides anatomical orientation and confirms that the tissue section captures several distinct brain regions, including the Hippocampus (with its structurally distinct sub-layers), cortical areas, fiber tracts, and lateral ventricles.

---

## 2. Image Feature Extraction and Image-Based Clustering

### 2.1 Rationale

Even with H&E images — which encode general tissue morphology through differential staining of nuclei and cytoplasm rather than specific molecular markers — quantitative image features can yield biologically relevant information. Regions with different cell types often have visually distinct morphological appearances, and those differences can be captured numerically. The goal is to extract a compact, per-spot representation of tissue morphology that can be analyzed alongside gene expression data.

### 2.2 Summary Feature Extraction at Multiple Scales

Summary features — statistical summaries of pixel intensity values (e.g., mean, standard deviation, percentiles across color channels) — are extracted for each Visium spot using `squidpy.im.calculate_image_features()`. Critically, this is performed **at two different spatial scales**:

- **Scale 1.0:** Features are computed from the tissue area corresponding exactly to the size of each Visium spot, capturing local morphology directly beneath the spot
- **Scale 2.0:** Features are computed from a larger crop centered on each spot, incorporating more surrounding tissue context — equivalent to zooming out and seeing a broader neighborhood

A higher scale value means more spatial context is incorporated into each spot's feature vector. This multi-scale approach is designed to capture morphological information at different levels of spatial granularity simultaneously, since some tissue structures are only apparent at coarser resolutions.

The two sets of features — one per scale — are computed independently for all 2,688 spots and then **concatenated into a single combined feature matrix** stored in `adata.obsm["features"]`. Column names are made unique to resolve any naming conflicts introduced by merging.

### 2.3 Image-Based Leiden Clustering

A helper function is defined to perform Leiden clustering on any subset of extracted image features. The function constructs a temporary AnnData object from the feature matrix, **scales the features** (necessary because image feature values are not inherently on comparable numeric scales), runs PCA with up to 10 components, builds a neighbor graph, and applies Leiden clustering. The resulting cluster labels are stored in the main AnnData observation metadata as `features_cluster`.

### 2.4 Comparison of Image and Gene Clusters

Image-based feature clusters are visualized side by side with gene expression-based clusters on the tissue. The comparison reveals an important pattern:

- In some regions — particularly the **Fiber Tract** and areas surrounding the **Hippocampus** — the image feature clusters closely resemble the gene expression clusters, indicating that tissue morphology and transcriptional identity are well-aligned in these anatomically distinct regions
- In the **cortex**, the two clustering approaches diverge: gene expression clusters capture the **laminar layered structure** of the cortex, whereas image feature clusters tend to delineate different **spatial regions** of the cortex rather than its layers

This divergence is biologically informative. It suggests that while cortical layers are transcriptionally distinct, their H&E morphological appearance may be more similar across layers within a given cortical area than between cortical areas. The tutorial notes that a more sophisticated joint analysis — such as computing a shared neighbor graph from concatenated PCA representations of both the gene expression and image feature spaces — could produce an integrated clustering that leverages both modalities simultaneously.

---

## 3. Spatial Statistics and Graph Analysis

### 3.1 Spatial Neighbor Graph Construction

All graph-based spatial analyses in Squidpy operate on a **spatial connectivity matrix** — a graph where Visium spots are nodes and edges connect spatially adjacent spots. This graph is constructed using `squidpy.gr.spatial_neighbors()`, which identifies neighboring spots based on their physical coordinates in the tissue. The resulting connectivity and distance matrices are stored in the AnnData object and serve as the foundation for all subsequent spatial statistics.

---

### 3.2 Neighborhood Enrichment Analysis

**Concept:** Neighborhood enrichment analysis quantifies whether spots belonging to two different clusters are spatially adjacent to each other more often — or less often — than would be expected by chance. If two clusters are frequently neighbors in tissue space, they receive a high enrichment score (they are "enriched" neighbors). If they are rarely adjacent, they receive a low or negative score (they are "depleted" neighbors).

**Method:** The enrichment score is computed using a **permutation-based test** via `squidpy.gr.nhood_enrichment()`. The cluster labels are randomly shuffled across spots a specified number of times (1,000 permutations by default), and the observed co-adjacency frequency between each pair of clusters is compared against this null distribution. This generates a score that reflects genuine spatial proximity rather than random chance.

**Visualization:** Results are displayed as a symmetric matrix heatmap using `squidpy.pl.nhood_enrichment()`, where each cell represents the enrichment score between a pair of clusters.

**Findings:** The analysis reveals strong neighborhood enrichment within the Hippocampus region. Specifically, the *Pyramidal_layer_dentate_gyrus* and *Pyramidal_layer* clusters are frequently spatially adjacent to the broader *Hippocampus* cluster. This is anatomically expected — the pyramidal cell layers are sub-structures embedded within the larger hippocampal formation — and serves as a validation that the enrichment score correctly captures known neuroanatomical spatial relationships.

---

### 3.3 Co-occurrence Analysis Across Spatial Dimensions

**Concept:** While neighborhood enrichment operates on the discrete connectivity graph (binary adjacency between spots), co-occurrence analysis works directly from the **continuous spatial coordinates** and evaluates how the probability of observing one cluster changes as a function of distance from another cluster. This provides a spatially resolved, scale-dependent view of cluster association.

**Score definition:** The co-occurrence score is defined as the conditional probability of observing cluster *exp* given the presence of cluster *cond* within a given radius, divided by the unconditional probability of observing cluster *exp* at that radius. A score greater than 1 indicates that cluster *exp* occurs more frequently near cluster *cond* than expected by chance; a score near 1 indicates no association; a score below 1 indicates spatial exclusion.

**Method:** The score is computed across a range of increasing radii around each spot using `squidpy.gr.co_occurrence()`, producing a curve that shows how cluster associations strengthen or weaken with spatial scale.

**Visualization and findings:** The analysis is focused on the **Hippocampus** cluster as the conditioning cluster. The resulting co-occurrence curves confirm that the *Pyramidal_layer* cluster shows elevated co-occurrence scores at **short distances** from the Hippocampus cluster, consistent with the neighborhood enrichment findings and with the known anatomy of the hippocampal region. The tutorial notes that distance units correspond to pixels of the source tissue image, directly linked to the spatial coordinate system stored in the AnnData object.

---

### 3.4 Ligand-Receptor Interaction Analysis

**Concept:** Having established that certain clusters are spatially proximal to one another, the logical next question is whether there is molecular evidence for **intercellular communication** between these neighboring populations. Ligand-receptor interaction analysis attempts to identify gene pairs — one encoding a secreted ligand and the other encoding its cognate receptor — that are co-expressed in adjacent cell populations, suggesting potential signaling events.

**Method:** Squidpy implements a fast reimplementation of **CellPhoneDB**, a widely used method for cell-cell communication inference, extended with the **Omnipath** database of annotated ligand-receptor pairs. The analysis is run using `squidpy.gr.ligrec()` with 100 permutation iterations. For each pair of clusters and each annotated ligand-receptor pair in the database, a permutation test evaluates whether the mean expression of the ligand in one cluster and the receptor in the other cluster is higher than expected under random cluster label shuffling. The output includes mean expression values and adjusted p-values for every cluster pair and interaction pair tested.

**Visualization:** The results are visualized as a **dot plot** using `squidpy.pl.ligrec()`, filtered to show only interactions with high mean expression (means greater than 3) and high statistical significance (adjusted p-value below 0.0001). The visualization is further focused on the **Hippocampus** cluster as the source, with *Pyramidal_layer* and *Pyramidal_layer_dentate_gyrus* as the target clusters — the same pairs identified as spatially proximal in the preceding analyses.

**Findings:** The dot plot reveals a set of candidate ligand-receptor pairs that may drive cellular communication between the Hippocampus and its pyramidal sub-layers. Each dot encodes both the mean co-expression level (dot size) and statistical significance (color or shading). The tutorial notes that these results represent candidate interactions rather than confirmed communication events, and that further refinement — such as integration with cell-type deconvolution results to determine the precise cell-type composition of each spot — would strengthen the interpretation.

---

### 3.5 Spatially Variable Gene Detection with Moran's I

**Concept:** A spatially variable gene is one whose expression is not randomly distributed across the tissue but instead follows a spatial pattern — with high expression concentrated in certain regions and low expression in others. Identifying such genes is important for understanding which molecular processes are spatially organized in the tissue, which often aligns with functional compartmentalization.

**Method:** Moran's I is a classical spatial autocorrelation statistic borrowed from spatial statistics and geography. It measures whether nearby observations (here, nearby spots) tend to have similar values (positive spatial autocorrelation, I close to +1), dissimilar values (negative autocorrelation, I close to -1), or no spatial structure (I close to 0). A gene with a Moran's I close to +1 has expression that is strongly spatially clustered — spots with high expression tend to be surrounded by other spots with high expression.

The analysis is performed using `squidpy.gr.spatial_autocorr()` in Moran's I mode (`mode="moran"`), applied to a subset of the 1,000 most highly variable genes to keep computation time manageable. Statistical significance is assessed using 100 permutations, and both nominal and permutation-based p-values with Benjamini-Hochberg FDR correction are reported. Results are stored in `adata.uns["moranI"]`, sorted by Moran's I statistic in descending order.

**Top spatially variable genes:** The top 10 genes ranked by Moran's I are:

| Gene | Moran's I |
|---|---|
| Olfm1 | 0.763 |
| Plp1 | 0.748 |
| Itpka | 0.727 |
| Snap25 | 0.721 |
| Nnat | 0.709 |
| Ppp3ca | 0.693 |
| Chn1 | 0.685 |
| Mal | 0.680 |
| Tmsb4x | 0.677 |
| Cldn11 | 0.674 |

All ten genes show very high Moran's I values (all above 0.67) with adjusted p-values of zero, indicating extremely strong and statistically robust spatial patterning.

**Spatial visualization:** The top three genes — *Olfm1*, *Plp1*, and *Itpka* — are mapped spatially alongside the cluster labels. Their expression patterns reveal associations with specific anatomical structures: notably, these genes appear enriched in the **pyramidal layers** of the Hippocampus and in **fiber tract** regions, consistent with their known biological roles. *Plp1* (Proteolipid Protein 1) is a well-established myelin gene enriched in oligodendrocytes, which are abundant in white matter fiber tracts. *Snap25* is a synaptic protein gene strongly expressed in neurons. The spatial patterns of these genes therefore reflect the underlying cellular composition and anatomical architecture of the mouse brain section.

---

## Summary and Conclusions

This analysis demonstrated a complete spatial transcriptomics workflow for H&E Visium data using Squidpy, spanning image feature extraction, spatial graph construction, neighborhood enrichment, co-occurrence scoring, ligand-receptor interaction analysis, and spatially variable gene detection.

The key findings can be summarized across four analytical dimensions:

**Image features** showed that H&E morphology captures meaningful biological structure. Image-based clustering broadly recapitulated gene-space clusters in anatomically distinct regions such as fiber tracts and the Hippocampus, while diverging in the cortex — where gene expression resolves laminar organization that image morphology alone cannot distinguish at equivalent cluster resolution.

**Neighborhood enrichment and co-occurrence** consistently identified strong spatial associations between the Hippocampus cluster and its pyramidal sub-layer clusters, validating that the spatial graph analysis correctly recovers known neuroanatomical adjacency relationships.

**Ligand-receptor analysis** extended the spatial proximity findings into the molecular domain, identifying candidate intercellular signaling interactions between the Hippocampus and the Pyramidal layer clusters, providing a basis for hypothesis generation about communication mechanisms in this brain region.

**Spatially variable genes** identified by Moran's I revealed strong spatial patterning across a subset of highly variable genes, with top-ranked genes such as *Plp1*, *Olfm1*, and *Snap25* showing expression patterns consistent with the known cellular composition of hippocampal and fiber tract regions.

Together, these results illustrate how Squidpy's spatial analysis toolkit — combining image processing, graph-based spatial statistics, and molecular interaction inference — enables a multi-layered interpretation of spatial transcriptomics data that substantially exceeds what gene expression analysis alone can provide.
