# Spatial Transcriptomics

> A comprehensive guide covering four spatial transcriptomics tutorials spanning the **Scanpy** and **Squidpy** ecosystems — from foundational Visium analysis to single-cell resolution Xenium data.

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
- [Biological Interpretation Guide](#biological-interpretation-guide)
- [Downstream Analysis Pointers](#downstream-analysis-pointers)
- [Glossary](#glossary)

---

## Environment & Installation

All four tutorials share a common Python environment rooted in the **scverse** ecosystem. A clean conda environment is strongly recommended to avoid dependency conflicts, particularly between `squidpy`, `spatialdata`, and legacy `scanpy` versions.

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

| Tutorial | Extra Packages | Notes |
|---|---|---|
| Tutorial 1 | `SpatialDE` | For spatially variable genes — now superseded by Squidpy's Moran's I |
| Tutorial 2 | `Dask`, `xarray` | For ImageContainer lazy loading; optionally `Cellpose` or `StarDist` for custom segmentation |
| Tutorial 3 | None | No additional packages beyond the core stack |
| Tutorial 4 | `spatialdata`, `spatialdata-io` | Mandatory; `Dask` recommended for large Xenium datasets |

### Setup

```bash
conda create -n spatial_env python=3.10
conda activate spatial_env
pip install scanpy squidpy spatialdata spatialdata-io
pip install leidenalg python-igraph dask xarray

# Tutorial 1 only (optional — Moran's I in Squidpy is now preferred):
pip install SpatialDE

# Tutorial 2 advanced segmentation (optional):
pip install cellpose stardist
```

### Version Reference

| Package | Minimum Version | Notes |
|---|---|---|
| `python` | 3.9 (3.10 recommended) | Runtime |
| `scanpy` | ≥ 1.9 | `sc.pl.spatial()` deprecated in 1.11+; replaced by `sq.pl.spatial_scatter()` |
| `squidpy` | ≥ 1.2 | ImageContainer API stable from 1.2 |
| `anndata` | ≥ 0.8 | Required by both scanpy and squidpy |
| `spatialdata-io` | ≥ 0.1 | Required only for Tutorial 4 |
| `SpatialDE` | 1.1.3 | Required only for Tutorial 1 |

> **Note:** `spatialdata` and `spatialdata-io` are actively developed. Always install the latest version to ensure compatibility with current Xenium output formats from 10x Genomics.

---

## Overview Comparison

| Property | Tutorial 1 — Scanpy Basic | Tutorial 2 — Squidpy Fluo | Tutorial 3 — Squidpy H&E | Tutorial 4 — Squidpy Xenium |
|---|---|---|---|---|
| **Library** | Scanpy | Squidpy + Scanpy | Squidpy + Scanpy | Squidpy + SpatialData |
| **Technology** | 10x Visium | 10x Visium | 10x Visium | 10x Xenium In Situ |
| **Image type** | H&E histology | Fluorescence (DAPI, NEUN, GFAP) | H&E histology | None / subcellular |
| **Resolution** | Spot-level (~55 µm) | Spot-level (~55 µm) | Spot-level (~55 µm) | Single-cell |
| **Tissue** | Human lymph node | Mouse brain (cropped) | Mouse brain (coronal) | Mouse brain (Xenium panel) |
| **Primary focus** | End-to-end preprocessing + clustering | Image feature extraction + segmentation | Spatial graph statistics + ligand-receptor | QC + graph analysis + spatial autocorrelation |
| **Key output** | Cluster visualization on tissue | Segmentation features mapped to gene space | Neighborhood enrichment + spatially variable genes | Leiden clusters + Moran's I spatial gene patterns |
| **Difficulty** | Beginner | Intermediate | Intermediate | Intermediate–Advanced |

---

## Tutorial 1 — Scanpy: Basic Visium Analysis

**Source:** https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html

### Purpose

This is the **foundational tutorial** for spatial transcriptomics in the Scanpy ecosystem. It walks through the entire pipeline from raw data loading to clustering and spatial visualization, using a human lymph node Visium dataset. A brief MERFISH example is included at the end to demonstrate multi-format compatibility.

### Dataset

- **Tissue:** Human lymph node (10x Visium)
- **Loaded via:** `sc.datasets.visium_sge('V1_Human_Lymph_Node')` — auto-downloads from 10x Genomics
- **Also covered:** MERFISH (U2-OS cells, demonstrating cell cycle clustering in spatial coordinates)

### AnnData Structure for Visium

When a Visium dataset is loaded, the AnnData object is organized as follows:

| Slot | Contents |
|---|---|
| `adata.X` | Gene expression count matrix (spots × genes) |
| `adata.obsm['spatial']` | 2D spot coordinates in pixel space |
| `adata.uns['spatial']` | Tissue image and scale factors |
| `adata.obs` | Per-spot metadata (QC metrics, cluster labels) |

### Analysis Steps

---

**Step 1 — Environment Setup**

Configure Scanpy with appropriate verbosity and white-background figure settings for spatial overlays. Verbosity level 3 prints informative status messages during computation, supporting reproducibility.

```python
import scanpy as sc
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3
sc.logging.print_versions()
```

---

**Step 2 — Data Loading**

The dataset is fetched and parsed into an AnnData object in a single call. Gene names are made unique before any processing.

```python
adata = sc.datasets.visium_sge('V1_Human_Lymph_Node')
adata.var_names_make_unique()
print(adata)  # Inspect shape and structure
```

---

**Step 3 — Quality Control**

Three standard QC metrics are computed per spot: `total_counts` (total UMI), `n_genes_by_counts` (detected gene count), and `mt_frac` (mitochondrial read fraction — a proxy for spot quality and cellular integrity). These are visualized in spatial coordinates to detect tissue regions with poor capture or cell damage.

```python
sc.pp.calculate_qc_metrics(adata, inplace=True)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.obs['mt_frac'] = (
    adata[:, adata.var['mt']].X.sum(1).A.squeeze() / adata.obs['total_counts']
)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])
```

> `sc.pl.spatial()` parameters: `img_key` selects `"hires"` or `"lowres"` tissue image; `alpha_img` sets transparency; `size` scales spot rendering; `crop_coord` restricts the view to a sub-region.

---

**Step 4 — Filtering**

Genes detected in fewer than 10 spots are removed to eliminate low-confidence transcripts. Spot-level filtering by minimum total count thresholds removes damaged or empty spots.

```python
sc.pp.filter_genes(adata, min_cells=10)
# Optionally filter spots: adata = adata[adata.obs['total_counts'] > 1000]
```

---

**Step 5 — Normalization & Highly Variable Gene Selection**

Library-size normalization scales each spot to the same total count. Log-transformation compresses the dynamic range. The top 2,000 highly variable genes (HVGs) are selected using the Seurat flavor for downstream dimensionality reduction.

```python
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
```

---

**Step 6 — Dimensionality Reduction & Clustering**

Standard single-cell workflow: PCA reduces to major axes of variance → KNN graph built in PCA space → UMAP projects to 2D → Leiden partitions the graph into expression-based communities. Clustering happens in **gene expression space**, not spatial space.

```python
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters")
sc.pl.umap(adata, color="clusters")
```

---

**Step 7 — Spatial Visualization of Clusters**

Expression-derived cluster labels are projected back onto spot coordinates overlaid on the H&E tissue image, revealing spatially coherent tissue zones.

```python
sc.pl.spatial(adata, img_key="hires", color="clusters")
```

---

**Step 8 — Marker Gene Visualization in Space**

Individual genes can be plotted the same way as cluster labels — spots are colored by expression level, overlaid on the tissue image. Confirms that cluster identity matches known spatial biology.

```python
sc.pl.spatial(adata, img_key="hires", color="CR2", size=1.5)
# CR2 (CD21) marks follicular B-cells and should localize to germinal centers
```

---

**Step 9 — Spatially Variable Genes with SpatialDE**

SpatialDE uses a Gaussian Process model to score each gene for spatial autocorrelation — whether its expression is more spatially structured than expected by chance, given spot coordinates. Returns a ranked table of spatially variable genes (SVGs).

```python
import SpatialDE
import pandas as pd
counts = pd.DataFrame(adata.X.todense(), columns=adata.var_names, index=adata.obs_names)
coord = pd.DataFrame(adata.obsm['spatial'], columns=['x_coord', 'y_coord'], index=adata.obs_names)
results = SpatialDE.run(coord, counts)
results.sort_values("qvalue").head(10)
```

> Squidpy's Moran's I (Tutorial 3) is now the preferred built-in alternative — faster, no separate install required.

---

**Step 10 — MERFISH Bonus Example**

Demonstrates that the same Scanpy clustering workflow applies to MERFISH subcellular data by using `sc.pl.embedding()` with `basis="spatial"` instead of `sc.pl.spatial()`.

```python
sc.pl.embedding(adata_merfish, basis="spatial", color="clusters")
```

---

### Key Outputs

- **UMAP plot** colored by Leiden cluster — distinguishes transcriptionally distinct populations
- **Spatial QC maps** of `total_counts` and `n_genes_by_counts` on the H&E image — identifies low-quality tissue regions
- **Spatial cluster map** — cluster IDs projected onto tissue, revealing germinal centers, mantle zones, and T-cell zones in the lymph node
- **Single gene spatial maps** — e.g., `CR2` confirming follicular B-cell localization
- **SpatialDE result table** — genes ranked by spatial variability with q-values and pattern type annotations

### Biological & Downstream Significance

Projecting transcriptomic clusters onto tissue reveals that gene expression communities correspond directly to histological compartments — germinal centers, mantle zones, and paracortical T-cell regions in the lymph node. This validates clustering and demonstrates that **transcriptional identity and spatial organization are tightly coupled**. SVGs from SpatialDE capture genes whose expression is structured by tissue architecture rather than random cell-to-cell noise — these are more biologically meaningful for spatial interpretation than standard HVGs. Downstream, the cluster assignments and SVGs serve as inputs for deconvolution (Tangram, RCTD), spatial interaction analysis (Tutorial 3), and trajectory methods.

---

## Tutorial 2 — Squidpy: Visium Fluorescence Analysis

**Source:** https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_fluo.html

### Purpose

This tutorial focuses on **image analysis** as a co-equal data modality alongside gene expression. Using a fluorescence-stained Visium dataset, it demonstrates Squidpy's `ImageContainer` object and `squidpy.im` image processing module — from nucleus segmentation through per-spot feature extraction and integration with transcriptomic data.

### Dataset

- **Tissue:** Mouse brain, coronal section — cropped subset for computational speed
- **Image:** Multi-channel fluorescence with 3 stains: DAPI (nuclei), anti-NEUN (neurons), anti-GFAP (glial cells)
- **Pre-processing:** Already clustered; gene expression preprocessing mirrors Tutorial 1

```python
import squidpy as sq
img = sq.datasets.visium_fluo_image_crop()    # ImageContainer object
adata = sq.datasets.visium_fluo_adata_crop()  # Pre-processed AnnData
```

### Analysis Steps

---

**Step 1 — Initial Visualization**

Squidpy's modern spatial plotting function displays pre-clustered spots on tissue coordinates.

```python
sq.pl.spatial_scatter(adata, color="cluster")
```

---

**Step 2 — Image Channel Inspection**

The fluorescence image has three biologically distinct channels: **DAPI** (Channel 0) for nuclear staining, **NEUN** (Channel 1) as a neuron-specific marker, and **GFAP** (Channel 2) as a glial cell marker. Visualized independently to assess image quality and guide channel selection for segmentation.

```python
img.show(layer="image", channel=0)  # DAPI
img.show(layer="image", channel=1)  # NEUN
img.show(layer="image", channel=2)  # GFAP
```

---

**Step 3 — Image Preprocessing (Smoothing)**

Gaussian smoothing reduces pixel-level noise before segmentation. The smoothed image is saved as a new named layer within the `ImageContainer` — the original is preserved. This non-destructive multi-layer design means multiple processed versions coexist without overwriting.

```python
sq.im.process(img=img, layer="image", method="smooth")
# Result stored as img["image_smooth"]
```

---

**Step 4 — Cell/Nucleus Segmentation**

Watershed segmentation applied to the smoothed DAPI channel detects individual nuclei. The output is a **label image** — each nucleus is assigned a unique integer ID; background pixels are zero. A cropped region is inspected side by side to assess segmentation quality.

```python
sq.im.segment(img=img, layer="image_smooth", method="watershed", channel=0, chunks=1000)

# Inspect segmentation quality on a crop
crop = img.crop_corner(2000, 2000, size=500)
crop.show(layer="image", channel=0)            # Raw DAPI
crop.show(layer="segmented_watershed", channel=0)  # Label image
```

> For complex tissue where nuclei overlap, replace `method="watershed"` with a custom Cellpose or StarDist backend via `sq.im.SegmentationCustom`.

---

**Step 5 — Segmentation Feature Calculation**

`sq.im.calculate_image_features()` is the central function of this tutorial. For each Visium spot, it crops the image to that spot's area and computes summary features from the segmentation label image. Features are stored in `adata.obsm["features_segmentation"]`.

```python
sq.im.calculate_image_features(
    adata, img,
    features="segmentation",
    layer="image",
    key_added="features_segmentation",
    features_kwargs={"segmentation": {"label_layer": "segmented_watershed"}},
)
# adata.obsm["features_segmentation"] now contains per-spot features
```

Computed features include:

| Feature name | Meaning |
|---|---|
| `segmentation_label` | Number of detected nuclei (cells) per spot |
| `segmentation_ch-0_mean_intensity_mean` | Average DAPI signal intensity per cell in spot |
| `segmentation_ch-1_mean_intensity_mean` | Average NEUN signal (proxy for neuronal density) |

---

**Step 6 — Other Feature Types**

Three additional feature extractors are available and computed the same way:

```python
# Summary features: pixel intensity statistics per spot
sq.im.calculate_image_features(adata, img, features="summary", key_added="features_summary")

# Texture features: Haralick descriptors from co-occurrence matrices
sq.im.calculate_image_features(adata, img, features="texture", key_added="features_texture")
```

Features can be computed at multiple spatial scales with the `scale` parameter — smaller scale = fine-grained, larger scale = broad-area context per spot.

---

**Step 7 — Integrated Visualization**

`sq.pl.extract()` temporarily promotes features from `adata.obsm` into `adata.obs`, enabling direct color mapping in `spatial_scatter`. The multi-panel plot compares nucleus count, signal intensities, and transcriptomic clusters on the same spatial coordinates.

```python
sq.pl.spatial_scatter(
    sq.pl.extract(adata, "features_segmentation"),
    color=["segmentation_label", "cluster",
           "segmentation_ch-0_mean_intensity_mean",
           "segmentation_ch-1_mean_intensity_mean"],
    frameon=False, ncols=2,
)
```

---

### Key Outputs

- **Fluorescence channel maps** — separate spatial maps for DAPI, NEUN, and GFAP channels per spot area
- **Segmentation overlay** — side-by-side crop of raw DAPI vs. watershed label image for visual QC of nucleus detection quality
- **Feature matrix** (`adata.obsm["features_segmentation"]`) — a spots × features table including nucleus count, mean DAPI intensity, and mean NEUN/GFAP intensity per spot
- **Multi-panel integrated plot** — spatial scatter showing cell density, NEUN signal, GFAP signal, and gene expression cluster simultaneously for direct comparison

### Biological & Downstream Significance

Fluorescence channels provide orthogonal biological information — NEUN marks neurons, GFAP marks astrocytes — that can validate or challenge gene expression cluster assignments. A cluster labeled neuronal by transcriptomics should show high NEUN signal; discordance suggests poor segmentation, doublet spots, or unexpected biology. Cell density per spot (nucleus count from segmentation) is **critical input for deconvolution methods** like Tangram and RCTD, which require estimated cell counts to correctly proportion cell types within multi-cell Visium spots. Texture features (Haralick descriptors) can capture histological texture differences between tissue regions, enabling morphology-informed clustering independent of gene expression — a powerful orthogonal validation. Downstream, the combined image + transcriptome feature space enables multimodal dimensionality reduction (MOFA+) and morphological biomarker discovery.

---

## Tutorial 3 — Squidpy: Visium H&E Analysis

**Source:** https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html

### Purpose

This is the **most analytically complete** Visium tutorial and the primary reference for quantitative spatial analysis of spot-based data. Using a pre-clustered mouse brain coronal section with H&E staining, it demonstrates Squidpy's full suite of spatial graph statistics: neighborhood enrichment, co-occurrence, ligand-receptor interactions, Ripley's statistics, and Moran's I spatially variable gene detection.

### Dataset

- **Tissue:** Mouse brain, coronal section (full slide, not cropped)
- **Image:** H&E staining
- **Cluster annotation:** Performed using Allen Brain Atlas and Mouse Brain gene expression atlas (Linnarsson lab)
- **Pre-processing:** Identical pipeline to Tutorial 1

```python
import squidpy as sq
adata = sq.datasets.visium_hne_adata()
img   = sq.datasets.visium_hne_image()
```

### Analysis Steps

---

**Step 1 — Spatial Scatter Visualization**

Cluster annotations are plotted in spatial coordinates, revealing the anatomical structure of the mouse brain — hippocampus, cortex, cerebellum, and olfactory bulb each appear as spatially contiguous, transcriptionally distinct domains.

```python
sq.pl.spatial_scatter(adata, color="cluster", figsize=(10, 8))
```

---

**Step 2 — Building the Spatial Neighbor Graph**

A spatial connectivity graph is constructed using the hexagonal Visium spot layout. With `n_rings=1`, each spot connects to its 6 immediate hexagonal neighbors. The graph is stored as sparse connectivity and distance matrices in the AnnData object. **This graph is a prerequisite for all subsequent spatial statistics.**

```python
sq.gr.spatial_neighbors(adata, coord_type="grid", n_rings=1)
# Stored in: adata.obsp['spatial_connectivities']
#            adata.obsp['spatial_distances']
```

---

**Step 3 — Neighborhood Enrichment Score**

A permutation test (default 1000 permutations) asks: are spots of cluster A and cluster B adjacent in tissue more often than expected if cluster labels were randomly shuffled? Results are visualized as a symmetric heatmap with hierarchical ordering.

```python
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(adata, cluster_key="cluster", figsize=(8, 8))
```

---

**Step 4 — Co-occurrence Score**

Extends neighborhood enrichment by computing the conditional probability that a spot of cluster B is found within increasing distance thresholds from cluster A — plotted as probability vs. distance (pixels).

```python
sq.gr.co_occurrence(adata, cluster_key="cluster")
sq.pl.co_occurrence(adata, cluster_key="cluster",
                    clusters="Hippocampus", figsize=(8, 4))
```

---

**Step 5 — Ligand-Receptor Interaction Analysis**

A CellPhoneDB-style permutation test integrated with the Omnipath database. Tests whether mean co-expression of a ligand (in source cluster) and its cognate receptor (in target cluster) exceeds a permuted null. Results are displayed as a dotplot filtered by expression level and FDR threshold.

```python
sq.gr.ligrec(adata, n_perms=1000, cluster_key="cluster", use_raw=False)
sq.pl.ligrec(
    adata, cluster_key="cluster",
    source_groups="Hippocampus",
    target_groups=["Pyramidal_layer", "Pyramidal_layer dentate gyrus"],
    means_range=(3, float("inf")), alpha=1e-4, swap_axes=True,
)
```

---

**Step 6 — Ripley's Statistics**

Ripley's L (normalized Ripley's K) measures the spatial aggregation pattern of each cluster independently — clustered, dispersed, or random — across a range of increasing radii.

```python
sq.gr.ripley(adata, cluster_key="cluster", mode="L")
sq.pl.ripley(adata, cluster_key="cluster", mode="L")
```

---

**Step 7 — Spatially Variable Genes with Moran's I**

Moran's I is computed for every gene. Values near +1 = strong positive spatial autocorrelation (expression follows tissue structure). Results are stored in `adata.uns["moranI"]` as a ranked table.

```python
sq.gr.spatial_autocorr(adata, mode="moran")
top_svgs = adata.uns["moranI"].sort_values("I", ascending=False).head(10)
print(top_svgs[["I", "pval_norm", "pval_z_sim"]])

# Visualize top SVGs on tissue
sq.pl.spatial_scatter(adata, color=top_svgs.index[:4].tolist(), ncols=2)
```

---

### Key Outputs

- **Anatomical spatial scatter** — cluster IDs overlaid on the H&E brain section, matching Allen Brain Atlas-defined regions (hippocampus, cortex, cerebellum)
- **Neighborhood enrichment heatmap** — square cluster-by-cluster matrix with color-coded enrichment z-scores; anatomically adjacent cluster pairs show high enrichment
- **Co-occurrence curves** — line plots of conditional proximity probability vs. distance for selected cluster pairs (e.g., `Pyramidal_layer` vs. `Hippocampus`)
- **Ligand-receptor dotplot** — significant interaction pairs between source and target cluster groups, filtered by mean expression and FDR; dot size = mean expression, color = significance
- **Ripley's L curves** — per-cluster spatial aggregation curves showing compactness across distance radii
- **Moran's I gene table** (`adata.uns["moranI"]`) — genes sorted by spatial autocorrelation score with p-values; high scorers include `Mbp` (myelin, oligodendrocytes) and `Nrgn` (neurons)

### Biological & Downstream Significance

The spatial neighbor graph-based statistics move beyond visualization into **quantitative spatial biology**. Neighborhood enrichment reveals anatomically expected co-localizations (e.g., `Pyramidal_layer` enriched with `Hippocampus`) and unexpected associations that may indicate novel functional interactions. The co-occurrence score adds a distance dimension — distinguishing direct cellular contact from longer-range regional co-distribution. Ligand-receptor analysis translates spatial co-localizations into molecular hypotheses: which signaling pathways might mediate communication between adjacent cell types? These are **correlative, hypothesis-generating** results that point toward experimental validation targets. Moran's I SVGs are particularly valuable because they capture genes whose expression is structured by **tissue architecture** (`Mbp` reflecting white vs. grey matter organization) rather than global expression variance. Downstream, SVGs serve as interpretable features for spatial domain calling, spatial trajectory analysis, and spatial deep learning models.

---

## Tutorial 4 — Squidpy: Xenium Analysis

**Source:** https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_xenium.html

### Purpose

This tutorial demonstrates working with **single-cell resolution** spatial transcriptomics data from the 10x Xenium In Situ platform. Unlike Visium (which captures ~10–100 cells per 55 µm spot), Xenium assigns transcripts to individual segmented cells — requiring a different data format, a different neighbor graph construction method, and a different analytical interpretation at every step.

### Dataset

- **Platform:** 10x Xenium In Situ (mouse brain)
- **Resolution:** Single-cell — each observation (`adata.obs` row) is one segmented cell
- **Data format:** SpatialData object (loaded via `spatialdata-io`)
- **Panel:** Targeted gene panel (~hundreds of genes vs. ~33,000 for Visium)
- **Coordinates:** Expressed in microns (µm), not pixel space

```python
import spatialdata_io
import squidpy as sq
sdata = spatialdata_io.xenium("path/to/xenium_output/")
adata = sdata["table"]  # Extract AnnData from SpatialData container
print(adata)            # Each row is one cell
```

### Analysis Steps

---

**Step 1 — Quality Control**

QC metrics are computed per cell rather than per spot. Because Xenium uses a targeted panel, absolute count thresholds are much lower than Visium — a typical cell may have 50–500 total transcripts. Very low count cells likely represent empty segmentation objects or debris.

```python
import scanpy as sc
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].hist(adata.obs['total_counts'], bins=50)
axes[0].set_xlabel('Total counts per cell')
axes[1].hist(adata.obs['n_genes_by_counts'], bins=50)
axes[1].set_xlabel('Genes detected per cell')
plt.tight_layout()
```

---

**Step 2 — Preprocessing & Clustering**

The same Scanpy normalization pipeline from Tutorial 1 is applied. Because Xenium captures single cells, the output **directly represents cell types** — no deconvolution step is needed.

```python
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata)
sc.pl.umap(adata, color="leiden")
```

---

**Step 3 — Spatial Scatter Visualization**

Because Xenium datasets contain millions of cells at high density, direct plotting requires subsampling. Cells are rendered as points (`shape=None`) rather than circles — appropriate for single-cell data with no predefined spot diameter.

```python
adata_sub = sc.pp.subsample(adata, fraction=0.1, copy=True)
sq.pl.spatial_scatter(adata_sub, color="leiden", shape=None, size=2)
```

---

**Step 4 — Building the Spatial Neighbor Graph**

Xenium cells are irregularly distributed in tissue — a fixed hexagonal grid is inappropriate. **Delaunay triangulation** connects each cell to its geometrically nearest neighbors without a fixed radius, naturally adapting to local cell density variations. This is the critical methodological difference from the Visium tutorials.

```python
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
# Stored in: adata.obsp['spatial_connectivities']
#            adata.obsp['spatial_distances']
```

---

**Step 5 — Centrality Scores**

Three graph centrality metrics are computed per cluster within the spatial cell graph. These describe the **topological role** of each cell type in the tissue architecture.

```python
sq.gr.centrality_scores(adata, cluster_key="leiden")
sq.pl.centrality_scores(adata, cluster_key="leiden", figsize=(12, 4))
```

| Metric | Biological meaning |
|---|---|
| Closeness centrality | How quickly a cluster can reach all others — spatially central cell types |
| Betweenness centrality | How often a cluster lies on shortest paths — spatial gatekeepers or bridges |
| Degree centrality | Number of direct cell neighbors — locally connected hub populations |

---

**Step 6 — Neighborhood Enrichment**

Same permutation-based approach as Tutorial 3, but computed at single-cell resolution. With far more observations (individual cells rather than spots), the permutation test has substantially higher statistical power. High enrichment = cell types in direct physical contact.

```python
sq.gr.nhood_enrichment(adata, cluster_key="leiden")
sq.pl.nhood_enrichment(adata, cluster_key="leiden", figsize=(8, 8),
                        title="Neighborhood enrichment — single-cell")
```

---

**Step 7 — Moran's I Spatial Autocorrelation**

Moran's I is computed on the Delaunay-based cell graph. At single-cell resolution, a high I-value gene means **neighboring cells** (not spots) share similar expression — capturing fine-grained, cell-type-specific spatial patterns impossible to resolve at spot resolution.

```python
sq.gr.spatial_autocorr(adata, mode="moran")
top_svgs = adata.uns["moranI"].sort_values("I", ascending=False).head(10)
print(top_svgs)

# Map top SVGs onto single-cell tissue coordinates
sq.pl.spatial_scatter(adata_sub,
                       color=top_svgs.index[:6].tolist(),
                       shape=None, size=2, ncols=3)
```

---

### Key Outputs

- **Single-cell spatial scatter** — subsampled point cloud of cells colored by Leiden cluster, overlaid on tissue coordinates in µm; reveals fine-grained tissue organization at cellular resolution
- **QC metric histograms** — transcript count and detected gene distributions per cell, with filtering thresholds annotated
- **Centrality score bar plots** — per-cluster closeness, betweenness, and degree centrality, identifying spatial hub cell types in the tissue cellular graph
- **Neighborhood enrichment heatmap** — at single-cell resolution, showing which Leiden cell types are spatially co-localized in direct contact
- **Moran's I SVG table + spatial maps** — top spatially variable genes with single-cell expression maps showing fine-grained expression domains (e.g., layer-specific cortical genes)

### Biological & Downstream Significance

Xenium's single-cell resolution unlocks analyses impossible with spot-based data. **Centrality scores** identify cell types that act as spatial hubs — for example, astrocytes may show high betweenness centrality because they bridge neuronal and vascular compartments, consistent with their architectural scaffolding role in brain tissue. **Neighborhood enrichment at single-cell resolution** reflects direct cell-cell contact rather than co-occupancy of a 55 µm spot, making biological interpretations far more precise — two cell types with high enrichment are likely physically interacting or in signaling proximity. **Moran's I at single-cell resolution** can resolve layer-specific gene expression patterns in cortex or subfield-specific patterns in hippocampus that would be blurred or averaged at spot resolution. Because the Xenium panel is targeted, SVG results are constrained by panel design — making panel choice a critical experimental design decision upstream of data collection. The SpatialData format enables future integration of transcriptomics with cell boundaries, subcellular transcript point clouds, and tissue images in a common coordinate system, setting the stage for truly multimodal spatial analyses.

---

## Method Coverage at a Glance

| Method | Tutorial 1 | Tutorial 2 | Tutorial 3 | Tutorial 4 |
|---|---|---|---|---|
| Spatial scatter / tissue overlay | `sc.pl.spatial` | `sq.pl.spatial_scatter` | `sq.pl.spatial_scatter` | `sq.pl.spatial_scatter` |
| Image smoothing | ❌ | ✅ | ❌ | ❌ |
| Nucleus / cell segmentation | ❌ | ✅ (watershed) | ❌ | Pre-done by instrument |
| Image feature extraction | ❌ | ✅ (segmentation, summary, texture) | ❌ | ❌ |
| Spatial neighbor graph | ❌ | ❌ | ✅ (hexagonal grid) | ✅ (Delaunay triangulation) |
| Centrality scores | ❌ | ❌ | ❌ | ✅ |
| Neighborhood enrichment | ❌ | ❌ | ✅ | ✅ |
| Co-occurrence score | ❌ | ❌ | ✅ | ❌ |
| Ripley's L statistics | ❌ | ❌ | ✅ | ❌ |
| Ligand-receptor (`ligrec`) | ❌ | ❌ | ✅ | ❌ |
| Spatially variable genes | ✅ SpatialDE (external) | ❌ | ✅ Moran's I (built-in) | ✅ Moran's I (built-in) |

---

## When to Use Which Tutorial

- **Tutorial 1** — Start here if new to spatial transcriptomics. Covers the complete Visium preprocessing pipeline (QC → normalization → clustering → spatial visualization). Every subsequent tutorial builds on this foundation and skips these steps.

- **Tutorial 2** — Use when your Visium experiment includes a fluorescence image (DAPI, immunofluorescence markers). The image segmentation and feature extraction pipeline is exclusive to this tutorial and is the entry point for spot deconvolution with Tangram or RCTD.

- **Tutorial 3** — Use for quantitative spatial analysis of Visium data: neighborhood enrichment, co-occurrence, Ripley's L, ligand-receptor interactions, and Moran's I. The primary reference for spatial graph statistics on spot-based data.

- **Tutorial 4** — Use for Xenium or other single-cell resolution platforms. Covers SpatialData format, Delaunay triangulation for neighbor graphs, and centrality scores — none of which appear in the Visium tutorials.

---

## Biological Interpretation Guide

### Clusters on Tissue (Tutorial 1)

When expression-derived clusters form spatially contiguous regions, they reflect underlying tissue architecture — anatomical layers, cell-type niches, or functional compartments. In the human lymph node, clusters map to germinal centers, mantle zones, and T-cell paracortex. Spatially scattered clusters suggest transcriptional states distributed throughout tissue — typical of infiltrating immune cells or stress-response populations. Always overlay known marker genes to confirm cluster identities before assigning biological labels.

### Image Features vs. Gene Clusters (Tutorial 2)

When cell density (nucleus count from segmentation) correlates with high total-count expression clusters, the correlation may be technical — more cells simply produce more transcripts. Conversely, when NEUN or GFAP mean intensity per spot aligns with gene expression clusters, this is strong evidence that the cluster captures a real cell-type difference. Discordance between image features and gene clusters may reveal transcriptional heterogeneity within a morphologically uniform cell population — itself an important biological finding.

### Neighborhood Enrichment (Tutorial 3)

A high positive score between two clusters means those cell types are consistently found in direct physical adjacency across tissue — a necessary condition for short-range signaling. Negative scores indicate spatial segregation, which may reflect developmental boundary maintenance or antagonistic cell fates. Use strongly enriched pairs as candidates for ligand-receptor analysis.

### Co-occurrence Score (Tutorial 3)

A curve rising sharply at short distances then decaying to 1.0 indicates tight local pairing — the types are immediate neighbors but not broadly co-distributed. A broadly elevated curve indicates shared regional occupancy. Values below 1.0 at any distance indicate active spatial exclusion at that scale.

### Ligand-Receptor Results (Tutorial 3)

A statistically significant pair (after FDR correction) means the data are consistent with a source cluster expressing a signaling molecule that could bind to a receptor expressed by an adjacent target cluster. Combined with high neighborhood enrichment, this constitutes a strong mechanistic hypothesis. These results are always **correlative and hypothesis-generating** — experimental validation is required to confirm active signaling.

### Moran's I — Spatially Variable Genes (Tutorials 3 & 4)

High Moran's I genes have expression that mirrors tissue structure. Unlike HVGs (selected for global expression variance), SVGs are selected for **spatial relevance** to tissue organization. A Moran's I > 0.2 with FDR-corrected p < 0.05 is a conventional significance threshold. Always visualize top candidates on tissue before interpretation.

### Centrality Scores at Single-Cell Resolution (Tutorial 4)

Cell types with high **betweenness centrality** spatially bridge different populations — potential signaling relays or structural organizers. High **closeness centrality** indicates a centrally located, hub-like position in the tissue microenvironment. Low centrality values identify peripheral, edge-associated, or isolated populations.

---

## Downstream Analysis Pointers

**Cell-Type Deconvolution** (follows Tutorial 2)
Visium spots contain mixtures of cells. **Tangram** maps single-cell RNA-seq reference atlas cell types to each spot, constrained by nucleus count from image segmentation. **Cell2Location** offers a hierarchical Bayesian alternative that does not require segmentation-derived cell counts. Both are integrated with Squidpy and have dedicated tutorials in the scverse documentation.

**Multi-Sample Integration** (follows Tutorial 1 or 3)
For experiments with multiple tissue sections or conditions, **Harmony** or **scVI** correct batch effects before cross-sample clustering. After integration, neighborhood enrichment and ligand-receptor analyses can be run per condition to identify context-dependent spatial interactions.

**Spatially Informed Differential Expression** (follows Tutorial 3)
**SPARK-X** and **nnSVG** identify spatially variable genes while explicitly modeling spatial autocorrelation — more sensitive than Moran's I for comparative SVG testing between conditions or tissue regions.

**Trajectory Analysis in Space** (follows Tutorial 1 or 4)
**PAGA** traces developmental or differentiation trajectories through the KNN graph in expression space. Overlaying PAGA trajectories on tissue coordinates reveals where cell state transitions occur physically — powerful for developmental tissues, tumor margins, or wound-healing contexts.

**Interactive Visualization**
All tutorial outputs are static matplotlib figures. For interactive exploration of large spatial datasets, **napari-spatialdata** enables real-time browsing of SpatialData objects with linked expression and spatial views directly from a Jupyter environment.

**Advanced Image Segmentation** (extends Tutorial 2)
Replace the default watershed backend with **Cellpose** or **StarDist** for improved segmentation in tissue with overlapping nuclei or low fluorescence signal. These are integrated via `sq.im.SegmentationCustom`, requiring no changes to the downstream feature extraction pipeline.

---

## Glossary

| Term | Definition |
|---|---|
| **AnnData** | Core data structure. Holds expression matrix (`X`), spot/cell metadata (`obs`), gene metadata (`var`), embeddings (`obsm`), spatial graph (`obsp`), and images (`uns`). |
| **SpatialData** | scverse unified format for multi-modal spatial omics — images, transcript coordinates, cell boundaries, and AnnData in one coordinate system. Required for Tutorial 4. |
| **ImageContainer** | Squidpy's xarray/Dask-backed image object. Stores multi-channel tissue images with lazy loading and non-destructive multi-layer editing. Used in Tutorials 2 & 3. |
| **Spatial Neighbor Graph** | Sparse connectivity matrix of physically adjacent spots or cells. Built by `sq.gr.spatial_neighbors()`. Prerequisite for all graph-based spatial statistics. |
| **Neighborhood Enrichment** | Permutation test for cluster co-localization. High score = types frequently adjacent; negative score = spatial segregation. |
| **Co-occurrence Score** | Conditional probability of finding cluster B within increasing radii from cluster A — adds a distance dimension to neighborhood enrichment. |
| **Moran's I** | Per-gene spatial autocorrelation statistic (−1 dispersed → +1 clustered). High-I genes are spatially variable genes (SVGs) reflecting tissue structure. |
| **Ripley's L** | Normalized point process statistic measuring whether a cluster is more spatially aggregated, dispersed, or random than expected across a range of radii. |
| **`ligrec`** | CellPhoneDB-style permutation test (Omnipath database). Tests significant co-expression of a ligand in a source cluster and its receptor in a target cluster. |
| **SVGs vs. HVGs** | SVGs are selected for spatial patterning across tissue; HVGs are selected for global expression variance. A gene can be one without being the other. |
| **Delaunay Triangulation** | Connects each cell to its geometrically nearest neighbors without a fixed radius — used for Xenium's irregular cell layouts (Tutorial 4) vs. the hexagonal grid used for Visium. |
| **Deconvolution** | Estimation of cell-type proportions within multi-cell Visium spots using a single-cell reference (e.g., Tangram, Cell2Location). Not needed for Xenium single-cell data. |
| **Leiden Clustering** | Modularity-optimizing community detection on the KNN expression graph. Standard in both single-cell and spatial transcriptomics workflows. |

---

## References

### Primary Software Publications

1. **Scanpy** — Wolf, F.A., Angerer, P. & Theis, F.J. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19, 15. https://doi.org/10.1186/s13059-017-1382-0

2. **Squidpy** — Palla, G., Spitzer, H., Klein, M. et al. (2022). Squidpy: a scalable framework for spatial omics analysis. *Nature Methods*, 19, 171–178. https://doi.org/10.1038/s41592-021-01358-2

3. **AnnData** — Virshup, I., Rybakov, S., Theis, F.J. et al. (2021). anndata: Annotated data. *bioRxiv*. https://doi.org/10.1101/2021.12.16.473007

4. **SpatialData** — Marconato, L., Palla, G., Yamauchi, K.A. et al. (2024). SpatialData: an open and universal data framework for spatial omics. *Nature Methods*, 21, 1352–1360. https://doi.org/10.1038/s41592-024-02212-x

### Spatial Transcriptomics Technology

5. **10x Visium** — Ståhl, P.L., Salmén, F., Vickovic, S. et al. (2016). Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. *Science*, 353(6294), 78–82. https://doi.org/10.1126/science.aaf2403

6. **10x Xenium In Situ** — 10x Genomics Xenium Platform white paper. https://www.10xgenomics.com/platforms/xenium

7. **MERFISH** — Chen, K.H., Boettiger, A.N., Moffitt, J.R. et al. (2015). Spatially resolved, highly multiplexed RNA profiling in single cells. *Science*, 348(6233), aaa6090. https://doi.org/10.1126/science.aaa6090

### Statistical Methods

8. **SpatialDE** — Svensson, V., Teichmann, S.A. & Stegle, O. (2018). SpatialDE: identification of spatially variable genes. *Nature Methods*, 15, 343–346. https://doi.org/10.1038/nmeth.4636

9. **Moran's I (original)** — Moran, P.A.P. (1950). Notes on continuous stochastic phenomena. *Biometrika*, 37(1/2), 17–23. https://doi.org/10.2307/2332142

10. **Leiden Algorithm** — Traag, V.A., Waltman, L. & van Eck, N.J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*, 9, 5233. https://doi.org/10.1038/s41598-019-41695-z

11. **CellPhoneDB / ligrec** — Efremova, M., Vento-Tormo, M., Teichmann, S.A. & Vento-Tormo, R. (2020). CellPhoneDB: inferring cell–cell communication from combined expression of multi-subunit ligand–receptor complexes. *Nature Protocols*, 15, 1484–1506. https://doi.org/10.1038/s41596-020-0292-x

12. **Omnipath (ligand-receptor database)** — Türei, D., Korcsmáros, T. & Saez-Rodriguez, J. (2016). OmniPath: guidelines and gateway for literature-curated signaling pathway resources. *Nature Methods*, 13, 966–967. https://doi.org/10.1038/nmeth.4077

### Downstream Tools Referenced

13. **Tangram (deconvolution)** — Biancalani, T., Scalia, G., Buffoni, L. et al. (2021). Deep learning and alignment of spatially resolved single-cell transcriptomes with Tangram. *Nature Methods*, 18, 1352–1362. https://doi.org/10.1038/s41592-021-01264-7

14. **Cell2Location (deconvolution)** — Kleshchevnikov, V., Shmatko, A., Dann, E. et al. (2022). Cell2location maps fine-grained cell types in spatial transcriptomics. *Nature Biotechnology*, 40, 661–671. https://doi.org/10.1038/s41587-021-01139-4

15. **Harmony (batch integration)** — Korsunsky, I., Millard, N., Fan, J. et al. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods*, 16, 1289–1296. https://doi.org/10.1038/s41592-019-0619-0

16. **SPARK-X (spatially variable genes)** — Zhu, J., Sun, S. & Zhou, X. (2021). SPARK-X: non-parametric modeling enables scalable and robust detection of spatial expression patterns for large spatial transcriptomic studies. *Genome Biology*, 22, 184. https://doi.org/10.1186/s13059-021-02404-0

17. **PAGA (trajectory analysis)** — Wolf, F.A., Hamey, F.K., Plass, M. et al. (2019). PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells. *Genome Biology*, 20, 59. https://doi.org/10.1186/s13059-019-1663-x


