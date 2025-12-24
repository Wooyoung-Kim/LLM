# Seurat - Single-Cell RNA-Seq Analysis

Comprehensive R skill for single-cell RNA-seq analysis using Seurat.

## Overview

Seurat is the most widely-used R package for scRNA-seq analysis, providing tools for quality control, normalization, clustering, differential expression, visualization, and integration.

## Quick Start

```r
library(Seurat)
source('scripts/seurat_helpers.R')

# Load data
counts <- Read10X(data.dir = "filtered_feature_bc_matrix/")
seurat_obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)

# QC
seurat_obj <- add_qc_metrics(seurat_obj, species = "human")
plot_qc_metrics(seurat_obj)
suggest_qc_thresholds(seurat_obj)

# Filter
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & 
                     nFeature_RNA < 2500 & percent.mt < 5)

# Standard workflow
seurat_obj <- run_standard_workflow(seurat_obj, n_pcs = 30, resolution = 0.5)

# Visualize
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

# Find markers
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25)
top_markers <- export_top_markers(markers, top_n = 10)

# Save
saveRDS(seurat_obj, "seurat_analysis.rds")
```

## Key Features

### Quality Control
- Automatic QC metrics calculation
- QC visualization and threshold suggestions
- Cell and gene filtering

### Normalization
- Log-normalization (standard)
- SCTransform (recommended for variance stabilization)
- Batch effect regression

### Clustering
- Graph-based clustering
- Multiple resolution testing
- Cluster validation

### Differential Expression
- FindMarkers for pairwise comparisons
- FindAllMarkers for cluster markers
- Multiple statistical tests (Wilcoxon, MAST, DESeq2)

### Visualization
- DimPlot, FeaturePlot, VlnPlot
- DotPlot, DoHeatmap, RidgePlot
- Custom plotting functions

### Integration
- Standard CCA integration
- Harmony integration
- SCTransform integration

### Data Export
- Export to AnnData/H5AD format
- Save/load between sessions
- Metadata export

## Directory Structure

```
seurat/
├── SKILL.md                    # Main documentation
├── README.md                   # This file
├── scripts/
│   └── seurat_helpers.R       # QC and workflow helpers
└── references/
    ├── seurat_examples.md     # Complete workflow examples
    └── markers_database.md    # Cell type marker genes
```

## Complete Workflows

See `references/seurat_examples.md` for 8 complete examples:
1. Standard PBMC analysis
2. Multi-sample integration with Harmony
3. Differential expression between conditions
4. Marker discovery and visualization
5. Subset and re-clustering
6. SCTransform workflow
7. Export to AnnData
8. Cell cycle scoring and regression

## Cell Type Markers

See `references/markers_database.md` for comprehensive lists of:
- Human PBMC markers (T, B, NK, monocytes, DCs)
- Mouse immune cell markers
- Tissue-specific markers (lung, gut, etc.)
- Cancer-related markers
- QC markers

## Common Use Cases

### Basic scRNA-seq Analysis
```r
# Complete workflow from counts to clusters
seurat_obj <- run_standard_workflow(seurat_obj)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
```

### Find Cell Type Markers
```r
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
DotPlot(seurat_obj, features = unique(top10$gene)) + RotatedAxis()
```

### Compare Conditions
```r
Idents(seurat_obj) <- "condition"
de_genes <- FindMarkers(seurat_obj, ident.1 = "treated", ident.2 = "control")
```

### Integrate Datasets
```r
library(harmony)
seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "orig.ident")
```

## Requirements

```r
# Core
install.packages("Seurat")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("patchwork")

# Optional but recommended
install.packages("harmony")         # Integration
install.packages("clustree")        # Clustering vis
install.packages("MAST")             # DE testing
BiocManager::install("SingleR")     # Auto-annotation
install.packages("MuDataSeurat")    # H5AD export
```

## Seurat Version

This skill is designed for **Seurat v5**.

Key differences from v4:
- Layer-based data structure
- Use `JoinLayers()` before export
- Improved SCTransform (vst.flavor = "v2")
- Sketch-based analysis for large datasets

## Related Skills

- **nichenet**: Ligand-receptor analysis and cell communication prediction
- **cellchat**: Cell-cell communication network analysis
- **scientific-visualization-r**: Publication-quality figure creation

## Citation

When using Seurat, cite:
```
Hao et al. (2021). Integrated analysis of multimodal single-cell data. 
Cell. doi:10.1016/j.cell.2021.04.048
```

## Documentation

- Main skill: [SKILL.md](SKILL.md)
- Examples: [references/seurat_examples.md](references/seurat_examples.md)
- Markers: [references/markers_database.md](references/markers_database.md)
- Official docs: https://satijalab.org/seurat/
