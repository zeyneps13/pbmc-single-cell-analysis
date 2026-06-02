# Single-Cell RNA-seq Analysis Pipeline — PBMC3K

![Python](https://img.shields.io/badge/python-3.11-blue)
![Scanpy](https://img.shields.io/badge/scanpy-1.12.1-green)
![License](https://img.shields.io/badge/license-MIT-lightgrey)
[![nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.org/github/YOUR_USERNAME/YOUR_REPO/blob/main/notebooks/Single-Cell_RNA-seq_Analysis_of_PBMC_Dataset.ipynb)

## Overview

This project presents a comprehensive end-to-end single-cell RNA sequencing (scRNA-seq) analysis pipeline applied to the PBMC3K benchmark dataset. The workflow was developed as a reproducible and biologically interpretable reference pipeline for single-cell transcriptomics, and serves as a bioinformatics portfolio project.

The pipeline covers the full analytical stack: from raw count matrix quality control to cell type annotation, trajectory inference, pathway enrichment, and cell-cell communication analysis.

---

## Biological Background

Peripheral Blood Mononuclear Cells (PBMCs) are a heterogeneous mixture of immune cells isolated from peripheral blood:

| Cell Type | Key Markers |
|-----------|-------------|
| CD4+ T cells | CD3D, LDHB |
| CD8+ Cytotoxic T cells | CD3D, CCL5, GZMK |
| B cells | CD79A, MS4A1 (CD20) |
| NK cells | NKG7, PRF1 |
| CD14+ Classical Monocytes | LYZ, S100A9 |
| FCGR3A+ Non-classical Monocytes | FCGR3A, LST1 |

PBMC datasets are widely used in single-cell transcriptomics benchmarking because they contain well-characterised, transcriptomically distinct immune populations that can be robustly separated by computational analysis.

---

## Dataset

| Property | Value |
|----------|-------|
| Source | `sc.datasets.pbmc3k()` — built-in Scanpy dataset |
| Organism | *Homo sapiens* |
| Tissue | Peripheral blood |
| Technology | 10x Genomics Chromium v1 |
| Cells | ~2,700 (after QC) |
| GEO accession | GSE96772 |
| Format | AnnData (.h5ad) |

---

## Pipeline Overview

```
Raw counts (2,700 cells × 32,738 genes)
    │
    ├── 1. Doublet detection (Scrublet)
    ├── 2. Gene filtering (min_cells = 3)
    ├── 3. QC metrics (mt%, ribo%)
    ├── 4. Adaptive MAD-based filtering (±3 MADs)
    │
    ├── 5. Normalization (CPM) + log1p
    ├── 6. Cell cycle scoring (Tirosh et al. 2015)
    ├── 7. HVG selection (top 2,000, Seurat flavor)
    ├── 8. Regression (mt%, S_score, G2M_score) + scaling
    │
    ├── 9.  PCA (50 components, elbow at PC15)
    ├── 10. kNN graph (n_neighbors=10, n_pcs=15)
    ├── 11. UMAP + t-SNE embeddings
    ├── 12. Leiden clustering (res=0.5 → 7 clusters)
    │
    ├── 13. Wilcoxon marker genes (BH correction)
    ├── 14. Manual annotation (marker-based)
    ├── 15. CellTypist automated annotation (Immune_All_Low)
    │
    ├── 16. PAGA connectivity graph
    ├── 17. Diffusion pseudotime (DPT) — monocyte axis
    ├── 18. Pseudotime gene dynamics heatmap
    │
    ├── 19. gProfiler GO:BP enrichment (monocyte markers)
    └── 20. LIANA cell-cell communication (monocyte → T/NK)
```

---

## Repository Structure

```
.
├── notebooks/
│   └── Single-Cell_RNA-seq_Analysis_of_PBMC_Dataset.ipynb
├── figures/          # all generated figures (.pdf)
├── results/
│   └── pbmc3k_analyzed.h5ad  # final annotated AnnData object
├── environment.yml   # conda environment (recommended)
├── requirements.txt  # pip fallback
└── README.md
```

> **Note:** `results/pbmc3k_analyzed.h5ad` is excluded from version control (see `.gitignore`).
> The dataset is loaded programmatically via `sc.datasets.pbmc3k()` — no manual download required.

---

## Installation

### Option 1 — Conda (recommended)

```bash
git clone https://github.com/zeyneps13/pbmc-single-cell-analysis-.git
cd pbmc_single_cell_analysis-
conda env create -f environment.yml
conda activate pbmc3k-scrna
jupyter lab
```

### Option 2 — pip

```bash
git clone https://github.com/zeyneps13/pbmc-single-cell-analysis-.git
cd pbmc_single_cell_analysis-
pip install -r requirements.txt
jupyter lab
```

Then open `notebooks/Single-Cell_RNA-seq_Analysis_of_PBMC_Dataset.ipynb`.

---

## Key Libraries

| Library | Version | Purpose |
|---------|---------|---------|
| [Scanpy](https://scanpy.readthedocs.io) | 1.12.1 | Core scRNA-seq analysis |
| [CellTypist](https://celltypist.readthedocs.io) | 1.7.1 | Automated cell type annotation |
| [LIANA](https://liana-py.readthedocs.io) | 1.7.1 | Cell-cell communication inference |
| [gProfiler](https://biit.cs.ut.ee/gprofiler) | 1.0.0 | Pathway enrichment (GO:BP) |
| [Scrublet](https://github.com/swolock/scrublet) | 0.2.3 | Doublet detection |
| NumPy / Pandas / SciPy | latest | Numerical and statistical computing |
| Matplotlib / Seaborn | latest | Visualization |

---

## Analytical Steps

### Quality Control

- Doublets removed using Scrublet (expected rate = 6%)
- Ribosomal (RPS*, RPL*) and mitochondrial (MT-*) gene fractions computed
- Adaptive **MAD-based thresholds** (±3 MADs) applied per metric to avoid hardcoded cutoffs
- Pre- and post-filtration QC distributions visualised

### Preprocessing

- Library-size normalization to 10,000 counts + log1p transformation
- Cell cycle scoring using Tirosh et al. (2015) S/G2M gene lists
- Top 2,000 highly variable genes selected (Seurat flavor)
- Confounders regressed out: `pct_counts_mt`, `S_score`, `G2M_score`
- Scaled to unit variance (max_value = 10)

### Dimensionality Reduction & Clustering

- PCA with 50 components; n_pcs = 15 selected from elbow plot
- kNN graph: n_neighbors = 10
- UMAP and t-SNE computed and compared
- Leiden clustering at resolution 0.5 → 7 clusters

### Annotation

- **Manual:** Canonical marker gene dotplot, matrixplot, stacked violin, tracksplot
- **Automated:** CellTypist `Immune_All_Low.pkl` with majority voting
- Cross-validation via heatmap and hierarchical clustermap

### Trajectory Inference

- PAGA connectivity graph on annotated cell types
- PAGA-initialized UMAP (original coordinates preserved)
- Diffusion Pseudotime (DPT) rooted at CD14+ Classical Monocytes
- Pseudotime gene dynamics heatmap (LYZ → FCGR3A maturation axis)

### Downstream Analysis

- GO:BP enrichment for CD14+ monocyte markers (gProfiler, FDR correction)
- Ligand-receptor communication inference (LIANA, consensus resource)
- Monocyte → T/NK cell signalling axes visualised

---

## Example Results

Cell types identified at Leiden resolution 0.5:

| Cluster | Cell Type | n cells |
|---------|-----------|---------|
| 0 | CD4+ T-Cells | ~1,100 |
| 1 | CD8+ Cytotoxic T-Cells | ~380 |
| 2 | CD14+ Classic Monocytes | ~480 |
| 3 | B-Cells | ~340 |
| 4 | NK Cells | ~155 |
| 5 | Cytotoxic NK Cells | ~80 |
| 6 | FCGR3A+ Non-classic Monocytes | ~65 |

---

## Reproducibility

- Random seed fixed at 42 throughout all stochastic steps
- All package versions pinned in `environment.yml` and `requirements.txt`
- Dataset loaded programmatically (no local data files required)
- All figures saved to `figures/` as PDF (300 dpi)
- Final AnnData object saved to `results/pbmc3k_analyzed.h5ad`

---

## Future Directions

- RNA velocity analysis (scVelo)
- Multi-sample batch correction (Harmony, scVI)
- Differential abundance testing across conditions
- Deep learning-based annotation (scANVI)
- Spatial transcriptomics integration

---

## References

1. Wolf FA, Angerer P, Theis FJ. **SCANPY: large-scale single-cell gene expression data analysis.** *Genome Biology* 19:15 (2018). https://doi.org/10.1186/s13059-017-1382-0

2. Domínguez Conde C et al. **Cross-tissue immune cell analysis reveals tissue-specific features in humans.** *Science* 376(6594) (2022). https://doi.org/10.1126/science.abl5197

3. Dimitrov D et al. **LIANA provides an all-in-one framework for cell-cell communication inference.** *Nature Communications* 13:5755 (2022). https://doi.org/10.1038/s41467-022-30755-0

4. Wolf FA, Hamey FK, Plass M et al. **PAGA: graph abstraction reconciles clustering with trajectory inference.** *Genome Biology* 20:59 (2019). https://doi.org/10.1186/s13059-019-1663-x

5. Tirosh I et al. **Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq.** *Science* 352(6282) (2015). https://doi.org/10.1126/science.aad0501

---

## Author

**Zeynep Sude Kırlı**
Master's Student — Computational Biology & Bioinformatics

This project was developed as a bioinformatics portfolio focused on single-cell transcriptomics analysis using the Python scientific computing ecosystem.
