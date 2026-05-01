---

# Analysis of Visium Fluorescence Data Using Squidpy


---

## Introduction

This analysis demonstrates how Squidpy's image analysis capabilities can be applied to fluorescence-based Visium spatial transcriptomics data. Unlike standard Visium workflows that rely on Hematoxylin and Eosin (H&E) brightfield images, this dataset is accompanied by a **multi-channel fluorescence image**, which enables a richer layer of morphological and cell-type-specific information to be extracted directly from the tissue image and integrated with gene expression data.

The dataset is a coronal section of the **mouse brain**, originally available through the 10x Genomics data portal. A cropped, pre-processed version is used here to keep computation times manageable. The pre-processing pipeline follows the standard Scanpy spatial workflow, and cluster annotations were performed using several established neuroanatomical references including the **Allen Brain Atlas**, the **Mouse Brain gene expression atlas** from the Linnarsson lab, and relevant literature.

The primary aim of this analysis is to extract quantitative image features from the fluorescence tissue image — including segmentation-based cell counts, channel intensities, texture, and summary statistics — and use these features to generate a complementary clustering of the tissue that can be compared against and validated by gene expression-based cluster assignments.

The computational environment uses **Scanpy** for gene expression analysis, **Squidpy** for spatial image processing, **AnnData** as the data container, **Pandas** for tabular data manipulation, and **Matplotlib** for visualization.

---

## 1. Data Loading and Initial Visualization

Two objects are loaded at the start of the analysis:

- **`adata`** — a pre-processed AnnData object containing gene expression counts and pre-annotated cluster labels for each Visium spot
- **`img`** — a Squidpy `ImageContainer` object holding the high-resolution fluorescence tissue image

The AnnData object contains **704 spots**, representing a crop of the full brain section. Pre-annotated cluster labels — assigned using neuroanatomical references — are already stored in the observation metadata.

The cluster annotations are immediately visualized in spatial context using Squidpy's `spatial_scatter()` function, which overlays the spot-level cluster labels directly onto the tissue image. This provides an initial anatomical orientation before any image processing begins.

The fluorescence image is then displayed channel by channel using the `ImageContainer.show()` method with the `channelwise=True` argument. The image contains **three fluorescence channels**:

- **Channel 0 — DAPI:** A DNA-binding stain that labels all cell nuclei, useful for identifying cell locations and counting cells
- **Channel 1 — anti-NEUN:** A marker specific to neurons, indicating the presence and density of neuronal cells
- **Channel 2 — anti-GFAP:** A marker specific to glial cells, indicating the presence of astrocytes and other glial populations

This multi-channel structure is what fundamentally distinguishes fluorescence Visium from H&E Visium — each channel encodes specific biological information about cell identity, not just overall morphology.

---

## 2. Image Segmentation

### 2.1 Purpose

Before image features can be extracted, the tissue image must be **segmented** — a process that delineates individual objects (in this case, cell nuclei) from the background. Segmentation enables the quantification of discrete cellular structures rather than just continuous pixel intensities.

### 2.2 Image Smoothing (Pre-processing)

As a best practice prior to segmentation, the image is first smoothed using Squidpy's `process()` function with the smoothing method. Smoothing reduces high-frequency noise in the image, which would otherwise cause the segmentation algorithm to over-fragment nuclei or produce spurious detections. The smoothed image is saved as a new layer within the ImageContainer, leaving the original image intact.

### 2.3 Watershed Segmentation

Segmentation is performed on the **DAPI channel** (channel 0) of the smoothed image using Squidpy's `segment()` function with the **watershed algorithm**. Watershed is a classical image segmentation method that treats pixel intensity values as a topographical surface and identifies individual objects by finding the "basins" — regions that would fill with water from local intensity minima — separated by "ridges." It is well-suited to nucleus segmentation from DAPI images because nuclei appear as bright, compact, well-separated objects against a dark background.

The `chunks` parameter controls how the image is divided into blocks for memory-efficient processing. The resulting segmentation is stored as a new layer in the ImageContainer named `segmented_watershed`. This is a **label image** — an array of the same dimensions as the original image, where each pixel belonging to a distinct segmented nucleus is assigned a unique integer identifier, and background pixels receive a value of zero.

### 2.4 Segmentation Visualization

A 500×500 pixel crop is taken from a corner of the image to visually validate the segmentation quality. Two panels are shown side by side — the raw DAPI channel and the corresponding watershed label image — allowing direct comparison of the detected nuclei against the original signal.

---

## 3. Segmentation Features

### 3.1 Concept

With the nucleus segmentation in place, **segmentation features** are extracted for each Visium spot using the `calculate_image_features()` function with the `segmentation` feature type specified. These features summarize what is happening within the segmentation mask underneath each spot and include:

- **Segmentation label count:** The number of distinct segmented objects (nuclei) detected within each spot — serving as a proxy for **cell count per spot**
- **Mean channel intensity within segmented objects:** For each fluorescence channel, the average pixel intensity within the segmented nuclei regions — providing an estimate of **cell-type composition** per spot based on which markers are expressed

### 3.2 Execution

The feature calculation processes all 704 spots and stores the results in a dedicated slot in the AnnData observation embedding (`adata.obsm["features_segmentation"]`). To make these features accessible for spatial visualization, the `extract()` helper function temporarily promotes the stored features into the main observation metadata (`adata.obs`), enabling them to be passed directly to the plotting function.

### 3.3 Results and Interpretation

Four spatial maps are plotted in a 2×2 grid:

**Cell count per spot (segmentation label):** The number of nuclei per Visium spot is mapped spatially. A notable observation is that the **pyramidal layer of the Hippocampus** (visible in the upper-left region of the crop) contains significantly more cells than the surrounding areas. This fine-grained cellular density pattern is not visible in the gene expression clusters, where the entire Hippocampus is grouped into a single cluster — demonstrating that image-derived features provide information that gene expression alone cannot resolve at standard clustering resolution.

**Gene-space cluster labels:** The pre-annotated cluster assignments are shown alongside for direct comparison.

**Channel 0 (DAPI) mean intensity within segmented objects:** Reflects overall nuclear content and cell density.

**Channel 1 (anti-NEUN) mean intensity within segmented objects:** Mapped spatially, this reveals that the regions labeled *Cortex_1* and *Cortex_3* exhibit higher NEUN signal, indicating a **higher density of neurons** in these cortical areas relative to the rest of the crop.

**Channel 2 (anti-GFAP) mean intensity within segmented objects:** Regions labeled *Fiber_tracts* and *Lateral Ventricles* show elevated GFAP signal, indicating enrichment of **glial cells** in these anatomical structures — a result consistent with known neurobiology, as white matter fiber tracts and ventricular linings are indeed rich in astrocytes and other glial populations.

These results demonstrate that fluorescence channel intensities within segmented nuclei provide direct, biologically interpretable readouts of local cell-type composition that complement and extend the transcriptional clustering.

---

## 4. Multi-Scale Feature Extraction and Image-Based Clustering

### 4.1 Feature Types

Beyond segmentation, three additional categories of image features are extracted to provide a comprehensive characterization of each spot's image content:

- **Summary features:** Statistical summaries of pixel intensity values within each spot (e.g., mean, standard deviation, percentiles). These capture the overall brightness and intensity distribution of the image patch.
- **Histogram features:** Pixel intensity histograms within each spot. These encode how pixel values are distributed across the intensity range, capturing more distributional detail than simple summary statistics.
- **Texture features:** Gray-level co-occurrence matrix (GLCM) based texture descriptors. These capture the spatial arrangement and patterns of pixel intensities — quantifying properties like homogeneity, contrast, and entropy — which reflect the microscopic structural organization of the tissue.

### 4.2 Multi-Scale Extraction Strategy

To provide richer contextual information, features are extracted under **three different parameter configurations**, each reflecting a different spatial scale:

- **`features_orig`:** Summary, texture, and histogram features extracted only from the tissue directly beneath each spot (a circular mask is applied to exclude tissue outside the Visium spot boundary), at full resolution. This gives the most precise, spot-specific readout.
- **`features_context`:** Summary and histogram features extracted from a slightly larger area around each spot at full resolution, incorporating some spatial context from neighboring tissue.
- **`features_lowres`:** Summary and histogram features extracted at a quarter of the original resolution (scale = 0.25), providing a broader, lower-detail view of the surrounding tissue. This captures large-scale tissue organization patterns that higher-resolution features may miss.

Each of these three feature sets is computed independently for all 704 spots and stored in separate slots within the AnnData embedding. They are then **concatenated into a single combined feature matrix** stored in `adata.obsm["features"]`. Duplicate column names that may arise from combining multiple feature sets are resolved by making the combined index unique.

### 4.3 Image-Based Clustering

A helper function is defined to perform Leiden clustering on any subset of the extracted image features. The function operates as follows: it optionally filters the combined feature matrix by a keyword to select a specific feature type, creates a temporary AnnData object from those features, **scales the features** (a critical step since image feature values are not inherently on comparable scales), computes PCA with up to 10 components, builds a neighbor graph, and runs Leiden clustering. The resulting cluster labels are returned and stored back into the main AnnData observation metadata.

This clustering is performed three times — once for summary features, once for histogram features, and once for texture features — yielding three independent image-based cluster assignments per spot.

### 4.4 Comparison of Image-Based and Gene-Based Clusters

All four cluster assignments (summary, histogram, texture, and gene-space) are visualized together in a spatial scatter plot with three columns, enabling direct spatial comparison.

The key findings from this comparison are:

- All three image-based feature cluster maps are **spatially coherent** — nearby spots tend to receive the same cluster label — demonstrating that the image features capture genuine tissue structure rather than noise.
- All three image-based cluster maps **resolve the Hippocampus into multiple sub-clusters**, reflecting the distinct structural layers (e.g., the pyramidal cell layer, the dentate gyrus, the stratum radiatum) that are well-known in hippocampal neuroanatomy. The gene-space clustering assigns the entire Hippocampus to a single cluster, representing a lower anatomical resolution.
- Similarly, the image feature clusters **subdivide the cortex into more fine-grained regions** than the gene expression clustering, potentially capturing laminar cortical organization.
- The three image-based methods differ from one another in their specific cluster boundaries, reflecting that each feature type encodes a different aspect of tissue morphology — intensity distribution, spatial texture, and histogram shape — and therefore segments the tissue differently.

The overall conclusion is that image-derived features provide **complementary and additive information** relative to gene expression-based clustering, and that integrating both modalities yields a more complete and anatomically detailed picture of tissue organization than either approach achieves independently.

---

## Summary and Conclusions

This analysis established a complete pipeline for extracting, analyzing, and interpreting image features from fluorescence Visium spatial transcriptomics data using Squidpy. The workflow proceeded through four major stages: data loading and initial cluster visualization, nucleus segmentation via the watershed algorithm applied to the DAPI channel, extraction of segmentation-based features including per-spot cell counts and fluorescence channel intensities, and multi-scale extraction of summary, histogram, and texture features followed by image-based Leiden clustering.

The central conclusion is that fluorescence tissue images in Visium datasets contain substantial biological information that is distinct from and complementary to gene expression measurements. Segmentation-based cell counts revealed cellular density patterns — particularly in the Hippocampus — invisible to gene-space clustering. Channel-specific intensity features identified cortical regions enriched in neurons and subcortical regions enriched in glial cells, consistent with established neuroanatomy. Image-based clustering produced spatially coherent and anatomically fine-grained cluster maps that subdivided brain structures at higher resolution than the corresponding transcriptional clusters.

Together, these results argue for the routine integration of image feature analysis alongside gene expression analysis in spatial transcriptomics workflows, particularly when high-quality tissue images with biologically informative staining are available.
