# Spatial Transcriptomics 

> A comprehensive guide covering four spatial transcriptomics tutorials spanning the **Scanpy** and **Squidpy** ecosystems from foundational Visium analysis to single-cell resolution Xenium data.

---

## Table of Contents

- [Environment & Installation](#environment--installation)
- [Overview Comparison](#overview-comparison)
- [Tutorial 1 — Scanpy: Basic Visium Analysis](#tutorial-1--scanpy-basic-visium-analysis)
- [Tutorial 2 — Squidpy: Visium Fluorescence Analysis](#tutorial-2--squidpy-visium-fluorescence-analysis)
- [Tutorial 3 — Squidpy: Visium H&E Analysis](#tutorial-3--squidpy-visium-he-analysis)
- [Tutorial 4 — Squidpy: Xenium Analysis](#tutorial-4--squidpy-xenium-analysis)
- [Method Coverage at a Glance](#method-coverage-at-a-glance)
- [When to Use Which Tutorial](#when-to-use-which-tutorial)
- [Glossary](#glossary)

---

## Environment & Installation

All four tutorials share a common Python environment rooted in the **scverse** ecosystem. Setting up a clean conda environment is strongly recommended to avoid dependency conflicts, particularly between `squidpy`, `spatialdata`, and legacy `scanpy` versions.

### Core Dependencies

| Package | Role |
|---|---|
| `scanpy` | Single-cell / spatial expression analysis, preprocessing, clustering, visualization |
| `squidpy` | Spatial graph statistics, image analysis, neighborhood enrichment, ligand-receptor |
| `anndata` | Core data structure (expression matrix + metadata + spatial coordinates) |
| `spatialdata` | Unified multi-modal spatial data format (required for Xenium, Tutorial 4) |
| `spatialdata-io` | I/O readers for Xenium, Visium HD, MERFISH into SpatialData format |
| `leidenalg` | Leiden clustering algorithm |
| `python-igraph` | Graph backend for Leiden and neighbor graph computations |

### Additional Dependencies by Tutorial

| Tutorial | Extra Packages |
|---|---|
| Tutorial 1 | `SpatialDE` (external; for spatially variable genes — now superseded by Squidpy's Moran's I) |
| Tutorial 2 | `Dask`, `xarray` (for ImageContainer lazy loading); optionally `Cellpose` or `StarDist` for custom segmentation |
| Tutorial 3 | No additional packages beyond core stack |
| Tutorial 4 | `spatialdata`, `spatialdata-io` (mandatory); `Dask` (for large Xenium datasets) |

###Setup Summary

```
conda create -n spatial_env python=3.10
conda activate spatial_env
pip install scanpy squidpy spatialdata spatialdata-io
pip install leidenalg python-igraph dask xarray
# For Tutorial 1 only (optional):
pip install SpatialDE
# For Tutorial 2 advanced segmentation (optional):
pip install cellpose stardist
```

> **Note:** `spatialdata` and `spatialdata-io` are actively developed; always install the latest version to ensure compatibility with current Xenium output formats from 10x Genomics.

---

## Overview Comparison

| Property | Tutorial 1 (Scanpy Basic) | Tutorial 2 (Squidpy Fluo) | Tutorial 3 (Squidpy H&E) | Tutorial 4 (Squidpy Xenium) |
|---|---|---|---|---|
| **Library** | Scanpy | Squidpy + Scanpy | Squidpy + Scanpy | Squidpy + SpatialData |
| **Technology** | 10x Visium | 10x Visium | 10x Visium | 10x Xenium In Situ |
| **Image type** | H&E (histology) | Fluorescence (DAPI, NEUN, GFAP) | H&E (histology) | None / subcellular |
| **Resolution** | Spot-level (~55 µm spots) | Spot-level (~55 µm spots) | Spot-level (~55 µm spots) | Single-cell resolution |
| **Tissue** | Human lymph node | Mouse brain (cropped) | Mouse brain (coronal) | Mouse brain (Xenium panel) |
| **Primary focus** | End-to-end preprocessing + clustering | Image feature extraction + segmentation | Spatial graph statistics + ligand-receptor | QC + graph analysis + spatial autocorrelation |
| **Key output** | Cluster visualization on tissue | Segmentation features mapped to gene space | Neighborhood enrichment + spatially variable genes | Leiden clusters + Moran's I spatial gene patterns |
| **Difficulty** | Beginner | Intermediate | Intermediate | Intermediate–Advanced |

---

## Tutorial 1 — Scanpy: Basic Visium Analysis

**Source:** https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html

### Purpose

This is the **foundational tutorial** for spatial transcriptomics in the Scanpy ecosystem. It walks through the entire pipeline from raw data loading to clustering and spatial visualization, using a human lymph node Visium dataset. A brief MERFISH example is included at the end.

### Dataset

- **Tissue:** Human lymph node (10x Visium)
- **Loaded via:** `sc.datasets.visium_sge('V1_Human_Lymph_Node')` — auto-downloads from 10x Genomics
- **Also covered:** MERFISH (U2-OS cells, demonstrating cell cycle clustering in spatial coordinates)

### Analysis Steps

**1. Environment Setup**
Configure Scanpy with appropriate verbosity and white-background figure settings for spatial overlays.
```python
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3
```

**2. Data Loading**
The AnnData object stores the gene expression count matrix (spots × genes), 2D spot coordinates in pixel space, the tissue image and scale factors, and per-spot metadata.
```python
adata = sc.datasets.visium_sge('V1_Human_Lymph_Node')
adata.var_names_make_unique()
```

**3. Quality Control**
Key QC metrics include: `total_counts` (total UMI per spot), `n_genes_by_counts` (number of expressed genes), and `mt_frac` (mitochondrial read fraction as a proxy for spot quality). These are visualized in spatial coordinates to identify low-quality tissue regions.
```python
sc.pp.calculate_qc_metrics(adata, inplace=True)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])
```

**4. Filtering**
Spots and genes with insufficient data are removed before normalization.
```python
sc.pp.filter_genes(adata, min_cells=10)
```

**5. Normalization & Highly Variable Gene Selection**
`normalize_total` scales each spot to the same total count, followed by log-transformation. 2000 highly variable genes are selected for dimensionality reduction using the Seurat flavor.
```python
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
```

**6. Dimensionality Reduction & Clustering**
Standard single-cell workflow: PCA → KNN graph → UMAP embedding → Leiden clustering. Clustering is performed in **gene expression space**, not spatial space.
```python
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters")
```

**7. Spatial Visualization of Clusters**
Clusters computed in expression space are projected back onto spatial coordinates overlaid on the H&E tissue image.
```python
sc.pl.spatial(adata, img_key="hires", color="clusters")
```

**8. Marker Gene Visualization in Space**
Individual marker genes (e.g., `CR2` for follicular B-cells) are overlaid on tissue to confirm that spatial patterns match known biology.
```python
sc.pl.spatial(adata, img_key="hires", color="CR2", size=1.5)
```

**9. Spatially Variable Genes with SpatialDE**
SpatialDE tests each gene for spatial autocorrelation using a Gaussian Process model, ranking genes by spatial variability. This is the original approach; Squidpy's Moran's I is now the preferred built-in alternative.
```python
results = SpatialDE.run(coord, counts)
results.sort_values("qvalue").head(10)
```

**10. MERFISH Bonus Example**
Brief demonstration that the same Scanpy clustering pipeline applies to MERFISH subcellular data.
```python
sc.pl.embedding(adata_merfish, basis="spatial", color="clusters")
```

### Key Outputs

- **UMAP plot** colored by Leiden cluster — shows transcriptionally distinct populations
- **Spatial scatter plots** of `total_counts` and `n_genes_by_counts` overlaid on H&E — reveals tissue regions with low capture or poor quality
- **Spatial cluster map** — cluster IDs projected onto the H&E image, revealing spatially coherent tissue zones (germinal centers, mantle zones, T-cell zones in lymph node)
- **Single gene spatial maps** — e.g., `CR2` expression patterns confirming known follicular B-cell localization
- **SpatialDE result table** — ranked list of spatially variable genes with p-values and spatial pattern types (linear, periodic, etc.)

### Biological & Downstream Significance

Projecting transcriptomic clusters onto tissue reveals that gene expression clusters correspond to known anatomical compartments — germinal centers, mantle zones, and paracortical T-cell regions in the lymph node. This validates the clustering approach and demonstrates that **transcriptional identity and spatial organization are coupled**. Spatially variable genes identified by SpatialDE are more informative than standard HVGs because they capture genes whose expression is structured by tissue architecture rather than random cell-to-cell variability. Downstream, these clusters and spatial gene patterns serve as the input for deconvolution methods (e.g., Tangram, RCTD) to estimate cell-type proportions per spot, and for trajectory or interaction analyses in Squidpy.

---

## Tutorial 2 — Squidpy: Visium Fluorescence Analysis

**Source:** https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_fluo.html

### Purpose

This tutorial focuses on **image analysis** capabilities of Squidpy's `ImageContainer` object using a fluorescence (not H&E) Visium dataset. The goal is to extract morphological features from the tissue image and integrate them with gene expression data to create a richer multimodal representation of each spot.

### Dataset

- **Tissue:** Mouse brain, coronal section — cropped subset for computational speed
- **Image:** Multi-channel fluorescence with 3 stains: DAPI (nuclei), anti-NEUN (neurons), anti-GFAP (glial cells)
- **Pre-processing:** Already clustered; gene expression preprocessing mirrors Tutorial 1

### Analysis Steps

**1. Initial Visualization**
`sq.pl.spatial_scatter()` — Squidpy's modern spatial plotting function — displays the pre-clustered spots on tissue.

**2. Image Channel Inspection**
The fluorescence image has three biologically distinct channels: DAPI (Channel 0) for nuclear staining and cell segmentation, NEUN (Channel 1) as a neuron-specific marker, and GFAP (Channel 2) as a glial cell marker. Each is visualized independently using `ImageContainer.show()`.

**3. Image Preprocessing (Smoothing)**
Gaussian smoothing reduces noise before segmentation. Results are saved as a new named layer within the `ImageContainer`, preserving the original image.

**4. Cell/Nucleus Segmentation**
Watershed segmentation applied to the DAPI channel detects individual nuclei. The output is a label image where each nucleus is assigned a unique integer ID. For more complex cases, custom segmentation models such as Cellpose or StarDist can be plugged in.

**5. Segmentation Feature Calculation**
`sq.im.calculate_image_features()` is the central function. For each Visium spot, it crops the image to that spot's area and computes: number of nuclei (cells) per spot, average DAPI signal intensity, and average NEUN signal (neuronal density). Features are stored in `adata.obsm["features_segmentation"]` as a spots × features matrix.

**6. Other Feature Types**
Summary features (mean, std, quantiles of pixel intensities), texture features (Haralick descriptors from co-occurrence matrices), and custom user-defined features are also supported. Features can be computed at multiple image scales to capture information from different spatial contexts.

**7. Integrated Visualization**
A multi-panel spatial plot compares cell counts and signal intensities against gene expression clusters — the key integration step between image and transcriptomic data.

### Key Outputs

- **Fluorescence channel maps** — separate spatial maps for DAPI, NEUN, and GFAP channels per spot area
- **Segmentation overlay image** — side-by-side crop of raw DAPI channel vs. watershed label image, allowing visual QC of segmentation quality
- **Feature matrix** (`adata.obsm["features_segmentation"]`) — a spots × features table including nucleus count per spot, mean DAPI intensity, and mean NEUN/GFAP intensity per spot
- **Multi-panel integrated plot** — spatial scatter showing cell count, NEUN signal, GFAP signal, and transcriptomic cluster side by side for direct comparison

### Biological & Downstream Significance

Fluorescence channels provide orthogonal biological information — NEUN staining marks neurons and GFAP marks astrocytes — that can validate or augment clustering derived purely from gene expression. For example, a Visium spot cluster identified as neuronal by transcriptomics should show high NEUN signal; discordance may indicate poor segmentation or doublet spots. Cell density per spot (nucleus count from segmentation) is critical input for **deconvolution methods** like Tangram or RCTD, which require estimated cell counts to correctly proportionate cell types within multi-cell spots. The texture features (Haralick descriptors) can capture histological texture differences between tissue regions, enabling morphology-informed clustering independent of gene expression. Downstream, this combined image + transcriptome feature space opens pathways to multimodal dimensionality reduction (e.g., MOFA+) and morphological biomarker discovery.

---

## Tutorial 3 — Squidpy: Visium H&E Analysis

**Source:** https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html

### Purpose

This is the **most analytically complete** Visium tutorial. Using a pre-clustered mouse brain coronal section with H&E staining, it demonstrates Squidpy's full suite of spatial graph statistics — neighborhood enrichment, co-occurrence, ligand-receptor interactions, Ripley's statistics, and spatially variable gene detection.

### Dataset

- **Tissue:** Mouse brain, coronal section (full slide, not cropped)
- **Image:** H&E staining
- **Cluster annotation:** Performed using Allen Brain Atlas and Mouse Brain gene expression atlas (Linnarsson lab)
- **Pre-processing:** Identical pipeline to Tutorial 1

### Analysis Steps

**1. Spatial Scatter Visualization**
Reveals the anatomical structure of the mouse brain: distinct clusters map to hippocampus, cortex, cerebellum, etc.

**2. Building the Spatial Neighbor Graph**
A spatial graph is constructed where each Visium spot is connected to its physical neighbors. For grid-based Visium data, the hexagonal spot layout is used, connecting each spot to its 6 immediate neighbors. The graph is stored as a sparse connectivity matrix and distance matrix in the AnnData object. **This graph is a prerequisite for all subsequent spatial statistics.**

**3. Neighborhood Enrichment Score**
Tests whether pairs of cluster types are spatially co-localized more than expected by chance via a permutation test (default 1000 permutations). A high positive score for clusters A and B means spots of those types are frequently adjacent in tissue. Results are visualized as an annotated heatmap with hierarchical clustering of clusters.

**4. Co-occurrence Score**
Extends neighborhood enrichment by computing the conditional probability that a spot of cluster B is found within increasing distance thresholds from a spot of cluster A, plotted as a curve of probability vs. distance. This quantifies **how close** (not just **whether**) two clusters are associated.

**5. Ligand-Receptor Interaction Analysis**
A re-implementation of CellPhoneDB's permutation test, integrated with the Omnipath ligand-receptor database. For each cluster pair and each annotated ligand-receptor pair, the method tests whether co-expression of the ligand (in source cluster) and receptor (in target cluster) is significantly higher than a permuted null. Results are displayed as a dotplot where dot size encodes mean expression and color encodes p-value significance after FDR correction.

**6. Ripley's Statistics**
Ripley's L (a normalized form of Ripley's K) quantifies the **spatial distribution pattern** of each cluster independently: clustered, dispersed, or random. High L values indicate tight spatial clustering of a cell type across increasing radius.

**7. Spatially Variable Genes with Moran's I**
Moran's I is a global spatial autocorrelation statistic computed per gene. Values near +1 indicate strong positive spatial autocorrelation (spatially clustered expression). Example high-scoring genes in mouse brain include `Mbp` (myelin basic protein, oligodendrocytes) and `Nrgn` (neurogranin, neurons).

### Key Outputs

- **Anatomical spatial scatter** — Leiden or annotated clusters overlaid on the H&E brain section, matching Atlas-defined regions
- **Neighborhood enrichment heatmap** — square matrix of cluster pairs with color-coded enrichment z-scores; clusters with high mutual enrichment are adjacent in tissue
- **Co-occurrence curves** — line plots showing conditional proximity probability vs. distance for selected cluster pairs (e.g., Pyramidal_layer vs. Hippocampus)
- **Ligand-receptor dotplot** — interaction pairs between source and target cluster groups, filtered by mean expression and FDR threshold
- **Ripley's L curves** — per-cluster spatial aggregation curves across distance radii
- **Moran's I ranked gene table** (`adata.uns["moranI"]`) — genes sorted by spatial autocorrelation score with associated p-values

### Biological & Downstream Significance

The spatial neighbor graph-based statistics in this tutorial move beyond simple visualization to **quantitative spatial biology**. Neighborhood enrichment reveals anatomically expected co-localizations (e.g., Pyramidal_layer enriched with Hippocampus) as well as unexpected spatial associations that may point to functional interactions. The co-occurrence score adds a distance dimension, distinguishing direct contact from longer-range proximity. Ligand-receptor analysis translates these spatial co-localizations into molecular hypotheses: which signaling pathways might mediate communication between adjacent cell types? Moran's I spatially variable genes are particularly valuable because, unlike standard HVGs, they are selected based on **tissue-structure relevance** — genes like `Mbp` and `Nrgn` reflect white matter vs. grey matter organization. Downstream, these SVGs can be used for spatial trajectory analysis, spatial domain calling, or as interpretable features for spatial deep learning models.

---

## Tutorial 4 — Squidpy: Xenium Analysis

**Source:** https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_xenium.html

### Purpose

This tutorial demonstrates working with **single-cell resolution** spatial transcriptomics data from the 10x Xenium In Situ platform. Unlike Visium (which captures ~10–100 cells per spot), Xenium assigns transcripts to individual segmented cells, providing a fundamentally different data structure and set of analytical considerations.

### Dataset

- **Platform:** 10x Xenium In Situ (mouse brain)
- **Resolution:** Single-cell — each observation is one cell, not a spot
- **Data format:** SpatialData object (loaded via `spatialdata-io`)
- **Panel:** Targeted gene panel (~hundreds of genes, vs. ~33,000 for Visium)
- **Coordinates:** Expressed in microns (µm), not pixel space

### Analysis Steps

**1. Data Loading via SpatialData**
Xenium output is loaded as a `SpatialData` object containing cell boundaries, transcript coordinates, and cell-by-gene count tables. The AnnData table is extracted for downstream analysis.

**2. Quality Control**
QC at single-cell resolution uses `total_counts` (transcript count per cell) and `n_genes_by_counts` (number of detected genes per cell). Cells with very low counts may represent empty segmentation masks or debris. Because Xenium uses a targeted panel, QC thresholds differ significantly from Visium whole-transcriptome data.

**3. Preprocessing & Clustering**
Standard Scanpy single-cell workflow: normalization → log-transform → PCA → KNN graph → Leiden clustering. Because Xenium captures single cells, the clustering output **directly represents cell types** — no deconvolution step is needed.

**4. Spatial Scatter Visualization**
Because Xenium has millions of cells at high density, direct visualization requires subsampling. Cells are displayed as points rather than circles, appropriate for single-cell data.

**5. Building the Spatial Neighbor Graph**
Xenium cells are irregularly distributed in tissue, so **Delaunay triangulation** is used instead of the hexagonal grid approach. Delaunay triangulation connects each cell to its geometrically nearest neighbors without a fixed radius, naturally adapting to local cell density variations. This is a critical methodological difference from the Visium tutorials.

**6. Centrality Scores**
Graph centrality metrics computed per cluster: closeness centrality (how close each node is to all others), betweenness centrality (how often a node lies on shortest paths between others), and degree centrality (number of neighbors relative to maximum possible). These describe the **topological role** of each cell type in the tissue graph.

**7. Neighborhood Enrichment**
Same approach as Tutorial 3, now at single-cell resolution. Because there are far more observations, the permutation test has higher statistical power. Clusters with high mutual enrichment represent cell types that physically interact in tissue.

**8. Moran's I Spatial Autocorrelation**
Moran's I is computed on the Delaunay-based spatial graph. At single-cell resolution, a high Moran's I gene means that neighboring cells tend to have similar expression levels — capturing **cell-type-specific** spatial patterns at finer scale than Visium.

### Key Outputs

- **Single-cell spatial scatter** — subsampled point cloud of cells colored by Leiden cluster, overlaid on tissue coordinates in µm
- **QC metric distributions** — histograms of transcript count per cell and detected genes per cell, with filtering thresholds
- **Centrality score bar plots** — per-cluster closeness, betweenness, and degree centrality, revealing topological hubs in the tissue cellular graph
- **Neighborhood enrichment heatmap** — at single-cell resolution, showing which cell types are spatially co-localized
- **Moran's I SVG table** — top spatially variable genes with spatial scatter maps at single-cell resolution showing fine-grained expression domains

### Biological & Downstream Significance

Xenium's single-cell resolution unlocks spatial analyses that are impossible with spot-based data. Centrality scores identify cell types that act as **spatial hubs** — for example, astrocytes may have high betweenness centrality because they bridge neuronal and vascular compartments, consistent with their architectural role. Neighborhood enrichment at single-cell resolution reflects **direct cell-cell contact** rather than co-occupancy of a 55 µm spot, making the biological interpretation more precise. Moran's I at single-cell resolution can capture fine-grained expression domains — e.g., layer-specific gene expression in cortex — that would be blurred at spot resolution. Because the Xenium panel is targeted, the SVGs reflect biology captured by the selected markers, making panel choice a critical experimental design decision. The SpatialData format enables future integration of transcripts, cell boundaries, and tissue images in a common coordinate system, setting the stage for truly multimodal spatial analyses.

---

## Method Coverage at a Glance

| Method | Tutorial 1 | Tutorial 2 | Tutorial 3 | Tutorial 4 |
|---|---|---|---|---|
| Spatial scatter plot | `sc.pl.spatial` | `sq.pl.spatial_scatter` | `sq.pl.spatial_scatter` | `sq.pl.spatial_scatter` |
| Spatial neighbor graph | ❌ | ❌ | ✅ (grid) | ✅ (Delaunay) |
| Neighborhood enrichment | ❌ | ❌ | ✅ | ✅ |
| Co-occurrence score | ❌ | ❌ | ✅ | ❌ |
| Ripley's statistics | ❌ | ❌ | ✅ | ❌ |
| Ligand-receptor (`ligrec`) | ❌ | ❌ | ✅ | ❌ |
| Moran's I (SVGs) | ❌ (SpatialDE) | ❌ | ✅ | ✅ |
| Centrality scores | ❌ | ❌ | ❌ | ✅ |
| Image segmentation | ❌ | ✅ | ❌ | ❌ |
| Image feature extraction | ❌ | ✅ | ❌ | ❌ |

---

## When to Use Which Tutorial

- **Tutorial 1** — Start here if you are new to spatial transcriptomics. Covers the full Visium preprocessing pipeline (QC → normalization → clustering → spatial visualization). Everything else builds on this foundation.
- **Tutorial 2** — Use if your Visium experiment includes a fluorescence image (DAPI, immunofluorescence). The image segmentation and feature extraction pipeline is unique to this tutorial and is the entry point for spot deconvolution with Tangram.
- **Tutorial 3** — Use for quantitative spatial analysis of Visium data: neighborhood enrichment, co-occurrence, Ripley's L, ligand-receptor interactions, and Moran's I. The primary reference for spatial graph statistics on spot-based data.
- **Tutorial 4** — Use for Xenium or other single-cell resolution platforms. Covers SpatialData format, Delaunay triangulation for neighbor graphs, and centrality scores — none of which appear in the Visium tutorials.

---

## Glossary

**SpatialData** — Unified multi-modal spatial omics format (scverse). Stores images, transcript coordinates, cell boundaries, and tabular data in a common coordinate system. Required for Xenium import in Tutorial 4.

**Deconvolution** — Computational methods (e.g., Tangram, RCTD) that estimate cell-type proportions within each multi-cell Visium spot. Requires reference single-cell data and is informed by nucleus counts from image segmentation (Tutorial 2).

**Co-occurrence Score** — Conditional probability of finding cluster B within increasing distance radii from cluster A. Unlike neighborhood enrichment (which flags *whether* clusters are adjacent), this quantifies *how close* they are as a function of distance.

**Ripley's L** — Normalized spatial point process statistic measuring whether a cluster's cells are more spatially clustered, dispersed, or random than expected by chance, computed across a range of distance radii.

**Delaunay Triangulation** — Graph construction method connecting each cell to its geometrically nearest neighbors without a fixed radius, adapting naturally to local cell density. Used in Tutorial 4 for irregular single-cell layouts; contrasts with the fixed hexagonal grid used for Visium spots.

---

