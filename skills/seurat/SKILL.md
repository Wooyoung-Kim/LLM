---
name: seurat
description: "Single-cell RNA-seq analysis with Seurat: QC, normalization, clustering, differential expression, marker discovery, integration, visualization (UMAP/tSNE), cell type annotation, and data export."
---

# Seurat - Single-Cell RNA-Seq Analysis

## Overview

Seurat is the most widely used R package for single-cell RNA-seq analysis. This skill covers the complete workflow: quality control, normalization, dimensionality reduction, clustering, differential expression analysis, visualization, data integration, and cell type annotation.

## When to Use This Skill

This skill should be used when:
- Analyzing single-cell or single-nucleus RNA-seq data
- Performing quality control on scRNA-seq datasets
- Identifying cell populations and clusters
- Finding marker genes and differential expression
- Integrating multiple datasets or batches
- Visualizing single-cell data (UMAP, tSNE, feature plots)
- Annotating cell types
- Converting between Seurat and AnnData formats
- Preparing data for downstream analysis (trajectory, velocity, etc.)

## Quick Start Guide

### Basic Seurat Workflow

```r
library(Seurat)
library(dplyr)
library(ggplot2)

# 1. Load data
counts <- Read10X(data.dir = "path/to/filtered_feature_bc_matrix/")
seurat_obj <- CreateSeuratObject(counts = counts, project = "my_project",
                                  min.cells = 3, min.features = 200)

# 2. Quality control
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

# Filter cells
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & 
                     nFeature_RNA < 2500 & percent.mt < 5)

# 3. Normalize and find variable features
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", 
                                    nfeatures = 2000)

# 4. Scale and PCA
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# 5. Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# 6. UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# 7. Visualization
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

# 8. Find markers
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, 
                          min.pct = 0.25, logfc.threshold = 0.25)

# 9. Save object
saveRDS(seurat_obj, "my_seurat_object.rds")
```

## Core Workflow Steps

### 1. Data Loading and Object Creation

#### Load 10X Genomics Data

```r
library(Seurat)

# Single sample
counts <- Read10X(data.dir = "filtered_feature_bc_matrix/")
seurat_obj <- CreateSeuratObject(counts = counts, 
                                  project = "sample1",
                                  min.cells = 3,      # Include features detected in >= 3 cells
                                  min.features = 200)  # Include cells with >= 200 features

# Multiple samples
sample_dirs <- list(
  ctrl1 = "data/ctrl1/filtered_feature_bc_matrix/",
  ctrl2 = "data/ctrl2/filtered_feature_bc_matrix/",
  treat1 = "data/treat1/filtered_feature_bc_matrix/"
)

seurat_list <- lapply(names(sample_dirs), function(sample_name) {
  counts <- Read10X(data.dir = sample_dirs[[sample_name]])
  CreateSeuratObject(counts = counts, project = sample_name, 
                     min.cells = 3, min.features = 200)
})

# Merge multiple samples
seurat_obj <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)],
                    add.cell.ids = names(sample_dirs))
```

#### Load from Count Matrix

```r
# From matrix
counts <- as.matrix(read.csv("counts.csv", row.names = 1))
seurat_obj <- CreateSeuratObject(counts = counts)

# From sparse matrix
library(Matrix)
counts <- readMM("matrix.mtx")
genes <- read.csv("genes.csv", header = FALSE)$V1
barcodes <- read.csv("barcodes.csv", header = FALSE)$V1
rownames(counts) <- genes
colnames(counts) <- barcodes

seurat_obj <- CreateSeuratObject(counts = counts)
```

### 2. Quality Control (QC)

#### Calculate QC Metrics

```r
# Mitochondrial percentage (human: ^MT-, mouse: ^mt-)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Ribosomal percentage
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

# Hemoglobin percentage (for blood samples)
seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^HB[^(P)]")

# View metadata
head(seurat_obj@meta.data)
```

#### Visualize QC Metrics

```r
# Violin plots
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.1)

# Scatter plots
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")

# Summary statistics
summary(seurat_obj$nFeature_RNA)
summary(seurat_obj$nCount_RNA)
summary(seurat_obj$percent.mt)
```

#### Filter Cells

```r
# Standard filtering
seurat_obj <- subset(seurat_obj, subset = 
                       nFeature_RNA > 200 & 
                       nFeature_RNA < 2500 & 
                       percent.mt < 5)

# More stringent filtering
seurat_obj <- subset(seurat_obj, subset = 
                       nFeature_RNA > 500 & 
                       nFeature_RNA < 6000 & 
                       nCount_RNA > 1000 &
                       nCount_RNA < 30000 &
                       percent.mt < 10)

# Remove doublets (after running DoubletFinder)
seurat_obj <- subset(seurat_obj, subset = DF.classifications == "Singlet")
```

### 3. Normalization

#### Log-Normalization (Default)

```r
# Standard log-normalization
seurat_obj <- NormalizeData(seurat_obj, 
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)
```

#### SCTransform (Recommended for Better Variance Stabilization)

```r
# SCTransform (replaces NormalizeData, ScaleData, FindVariableFeatures)
seurat_obj <- SCTransform(seurat_obj, 
                          vars.to.regress = "percent.mt",
                          verbose = FALSE)

# For large datasets
seurat_obj <- SCTransform(seurat_obj, 
                          vars.to.regress = c("percent.mt", "nCount_RNA"),
                          vst.flavor = "v2",
                          verbose = FALSE)
```

### 4. Feature Selection and Scaling

#### Find Variable Features

```r
# Find highly variable features
seurat_obj <- FindVariableFeatures(seurat_obj, 
                                    selection.method = "vst",
                                    nfeatures = 2000)

# Identify top 10 most variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)
print(top10)

# Plot variable features
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
```

#### Scale Data

```r
# Scale all genes
all_genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all_genes)

# Scale with regression (removes unwanted variation)
seurat_obj <- ScaleData(seurat_obj, 
                        vars.to.regress = c("percent.mt", "nCount_RNA"))

# Scale only variable features (faster)
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
```

### 5. Dimensionality Reduction

#### PCA

```r
# Run PCA
seurat_obj <- RunPCA(seurat_obj, 
                     features = VariableFeatures(seurat_obj),
                     npcs = 50)

# Visualize PCA
DimPlot(seurat_obj, reduction = "pca")
DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)

# Elbow plot to determine dimensionality
ElbowPlot(seurat_obj, ndims = 50)

# Print top genes for each PC
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
```

#### UMAP

```r
# Run UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# Visualize
DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 5)

# Colored by metadata
DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident")
DimPlot(seurat_obj, reduction = "umap", split.by = "condition")
```

#### tSNE

```r
# Run tSNE
seurat_obj <- RunTSNE(seurat_obj, dims = 1:30)

# Visualize
DimPlot(seurat_obj, reduction = "tsne", label = TRUE)
```

### 6. Clustering

```r
# Find neighbors
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)

# Find clusters
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Try multiple resolutions
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.1, 0.3, 0.5, 0.8, 1.0))

# View clustering results
head(Idents(seurat_obj))
table(Idents(seurat_obj))

# Visualize clusters
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

# Set active identity to specific resolution
Idents(seurat_obj) <- "RNA_snn_res.0.5"
```

### 7. Differential Expression Analysis

#### Find Markers for All Clusters

```r
# Find markers for all clusters
all_markers <- FindAllMarkers(seurat_obj, 
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25,
                              test.use = "wilcox")  # wilcox, MAST, DESeq2

# View top markers
top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
print(top_markers)

# Save markers
write.csv(all_markers, "all_cluster_markers.csv")
```

#### Find Markers for Specific Cluster

```r
# Markers for cluster 0
cluster0_markers <- FindMarkers(seurat_obj, ident.1 = 0, min.pct = 0.25)

# Markers for cluster 0 vs cluster 1
cluster0_vs_1 <- FindMarkers(seurat_obj, ident.1 = 0, ident.2 = 1, 
                              min.pct = 0.25)

# Conserved markers (across conditions)
conserved_markers <- FindConservedMarkers(seurat_obj, 
                                           ident.1 = 0,
                                           grouping.var = "condition")
```

#### Compare Groups

```r
# Set identity to condition
Idents(seurat_obj) <- "condition"

# Find DEGs between conditions
de_genes <- FindMarkers(seurat_obj, 
                        ident.1 = "treated", 
                        ident.2 = "control",
                        test.use = "MAST",
                        min.pct = 0.10,
                        logfc.threshold = 0.25)

# Within specific cell type
tcells <- subset(seurat_obj, subset = celltype == "T cells")
tcell_de <- FindMarkers(tcells, 
                        ident.1 = "treated",
                        ident.2 = "control",
                        group.by = "condition")
```

### 8. Visualization

#### DimPlot - Dimensionality Reduction Plots

```r
# Basic UMAP
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

# Color by metadata
DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type")

# Split by condition
DimPlot(seurat_obj, reduction = "umap", split.by = "condition", ncol = 2)

# Highlight specific clusters
DimPlot(seurat_obj, reduction = "umap", cells.highlight = WhichCells(seurat_obj, idents = c(0, 1, 2)))
```

#### FeaturePlot - Gene Expression

```r
# Single gene
FeaturePlot(seurat_obj, features = "CD3D")

# Multiple genes
FeaturePlot(seurat_obj, features = c("CD3D", "CD8A", "CD4", "IL7R"))

# Split by condition
FeaturePlot(seurat_obj, features = "CD3D", split.by = "condition")

# Blended feature plot
FeaturePlot(seurat_obj, features = c("CD3D", "CD8A"), blend = TRUE)
```

#### VlnPlot - Violin Plots

```r
# Single feature
VlnPlot(seurat_obj, features = "CD3D")

# Multiple features
VlnPlot(seurat_obj, features = c("CD3D", "CD8A", "CD4"), ncol = 3)

# Split by group
VlnPlot(seurat_obj, features = "CD3D", split.by = "condition")

# No points (cleaner)
VlnPlot(seurat_obj, features = "CD3D", pt.size = 0)
```

#### DotPlot - Marker Expression Heatmap

```r
# Top markers per cluster
markers_to_plot <- c("CD3D", "CD8A", "CD4", "IL7R", "CCR7", 
                     "CD14", "FCGR3A", "MS4A1", "CD79A")

DotPlot(seurat_obj, features = markers_to_plot) + 
  RotatedAxis()

# With customization
DotPlot(seurat_obj, features = markers_to_plot, 
        cols = c("lightgrey", "red"),
        dot.scale = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

#### DoHeatmap - Expression Heatmap

```r
# Top markers
top10_markers <- all_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

DoHeatmap(seurat_obj, features = top10_markers$gene) + 
  NoLegend()

# Subset of cells
DoHeatmap(subset(seurat_obj, downsample = 100), 
          features = top10_markers$gene, 
          size = 3)
```

#### RidgePlot - Distribution of Expression

```r
RidgePlot(seurat_obj, features = c("CD3D", "CD8A", "CD4"), ncol = 1)
```

### 9. Cell Type Annotation

#### Manual Annotation

```r
# View cluster markers
DotPlot(seurat_obj, features = canonical_markers) + RotatedAxis()

# Assign cell types
new_cluster_ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", 
                     "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")

names(new_cluster_ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new_cluster_ids)

# Save as metadata
seurat_obj$celltype <- Idents(seurat_obj)

# Visualize
DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 4.5)
```

#### Automated Annotation with SingleR

```r
library(SingleR)
library(celldex)

# Get reference
ref <- celldex::HumanPrimaryCellAtlasData()

# Extract counts
counts <- GetAssayData(seurat_obj, slot = "data")

# Run SingleR
pred <- SingleR(test = counts, ref = ref, 
                labels = ref$label.main)

# Add to Seurat object
seurat_obj$singler_labels <- pred$labels

# Visualize
DimPlot(seurat_obj, reduction = "umap", group.by = "singler_labels", label = TRUE)
```

#### Automated Annotation with Azimuth

```r
library(Azimuth)

# Run Azimuth (requires internet connection)
seurat_obj <- RunAzimuth(seurat_obj, reference = "pbmcref")

# View predicted cell types
DimPlot(seurat_obj, group.by = "predicted.celltype.l2", label = TRUE)
```

### 10. Data Integration

#### Standard Integration

```r
# List of Seurat objects
seurat_list <- list(sample1 = obj1, sample2 = obj2, sample3 = obj3)

# Normalize and find variable features
seurat_list <- lapply(seurat_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30)

# Integrate
seurat_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

# Switch to integrated assay
DefaultAssay(seurat_integrated) <- "integrated"

# Run standard workflow
seurat_integrated <- ScaleData(seurat_integrated)
seurat_integrated <- RunPCA(seurat_integrated, npcs = 30)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:30)
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:30)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.5)

# Visualize
DimPlot(seurat_integrated, reduction = "umap", group.by = "orig.ident")
DimPlot(seurat_integrated, reduction = "umap", label = TRUE)
```

#### Harmony Integration

```r
library(harmony)

# Merge datasets first
seurat_obj <- merge(obj1, y = c(obj2, obj3), 
                    add.cell.ids = c("sample1", "sample2", "sample3"))

# Standard workflow
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

# Run Harmony
seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "orig.ident")

# UMAP on Harmony embedding
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:30)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Visualize
DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident")
```

### 11. Data Export and Conversion

#### Export to AnnData (for Python/Scanpy)

```r
library(MuDataSeurat)

# Prepare Seurat object for export
seurat_obj <- DietSeurat(
  seurat_obj,
  counts = TRUE,          # Raw counts -> adata.layers['counts']
  data = TRUE,            # Log-normalized -> adata.X
  scale.data = FALSE,     # Don't export scaled data
  features = rownames(seurat_obj),  # All genes
  assays = "RNA",
  dimreducs = c("pca", "umap"),
  graphs = c("RNA_nn", "RNA_snn")
)

# Join layers (required for Seurat v5)
seurat_obj <- JoinLayers(seurat_obj)

# Export
MuDataSeurat::WriteH5AD(seurat_obj, "seurat_object.h5ad", assay = "RNA")
```

#### Save and Load Seurat Objects

```r
# Save
saveRDS(seurat_obj, "my_seurat_object.rds")

# Load
seurat_obj <- readRDS("my_seurat_object.rds")

# Update old Seurat objects
seurat_obj <- UpdateSeuratObject(seurat_obj)
```

## Advanced Workflows

See `references/advanced_workflows.md` for:
- Trajectory analysis preparation
- RNA velocity with velocyto
- Cell-cell communication (CellChat, NicheNet)
- Spatial transcriptomics
- Multi-modal analysis (CITE-seq, ATAC-seq)
- Custom scoring and signatures

## Common Tasks

### Task 1: Standard scRNA-seq Analysis

Complete workflow from raw counts to annotated clusters.

See `references/seurat_examples.md` Example 1.

### Task 2: Batch Integration

Integrate multiple samples while preserving biological variation.

See `references/seurat_examples.md` Example 2.

### Task 3: Differential Expression Analysis

Find DEGs between conditions within specific cell types.

See `references/seurat_examples.md` Example 3.

### Task 4: Marker Discovery and Visualization

Identify and visualize cell type-specific markers.

See `references/seurat_examples.md` Example 4.

### Task 5: Subset and Re-cluster

Extract specific cell populations for detailed analysis.

See `references/seurat_examples.md` Example 5.

## Best Practices

### Quality Control
- Always check QC metrics before filtering
- Adjust thresholds based on your data type
- Consider tissue-specific QC (e.g., higher mt% for heart/muscle)
- Remove doublets using DoubletFinder or Scrublet

### Normalization
- Use SCTransform for better variance stabilization
- Consider using sctransform_v2 for very large datasets
- Regress out technical variation (mt%, cell cycle) when needed

### Clustering
- Try multiple resolutions (0.1-1.0)
- Use Clustree to visualize clustering stability
- Validate clusters with known markers
- Consider overclustering then merging similar clusters

### Differential Expression
- Use appropriate statistical tests (MAST for complex designs)
- Set reasonable thresholds (min.pct, logfc.threshold)
- Multiple testing correction is applied automatically
- Validate key findings with additional methods

### Visualization
- Use consistent color schemes across figures
- Include scale bars and labels
- Show both clustering and biological groupings
- Export high-resolution figures for publication

## Common Pitfalls to Avoid

1. **Over-filtering**: Too stringent QC removes real cell types
2. **Wrong mitochondrial pattern**: Use ^MT- for human, ^mt- for mouse
3. **Not scaling data**: Required before PCA
4. **Too few PCs**: Check elbow plot, typically use 20-50
5. **Single resolution**: Always try multiple clustering resolutions
6. **Wrong assay**: Switch between "RNA" and "integrated" as needed
7. **Ignoring batch effects**: Use integration methods when necessary
8. **Not saving intermediate results**: Save objects at key steps
9. **Memory issues**: Use sketch-based methods for >1M cells
10. **Version incompatibility**: Update old objects with UpdateSeuratObject

## Resources

### Scripts Directory
- `seurat_qc.R` - Quality control helper functions
- `seurat_integration.R` - Integration workflows
- `seurat_visualization.R` - Publication-quality plots

### References Directory
- `seurat_examples.md` - Complete workflow examples
- `markers_database.md` - Common cell type markers
- `troubleshooting.md` - Common issues and solutions
- `advanced_workflows.md` - Specialized analyses

### Assets Directory
- `marker_genes.csv` - Curated marker gene lists
- `color_schemes.R` - Consistent color palettes

## Seurat Version Notes

This skill is designed for **Seurat v5**. Key differences from v4:
- Layer-based data structure (use `JoinLayers()` before certain operations)
- Improved SCTransform (`vst.flavor = "v2"`)
- Better integration methods
- Sketch-based analysis for large datasets

For legacy code, use `UpdateSeuratObject()` to convert v3/v4 objects.

## Citation

If you use Seurat in your research, please cite:

```
Hao et al. (2021). Integrated analysis of multimodal single-cell data. 
Cell. doi:10.1016/j.cell.2021.04.048

Stuart*, Butler*, et al. (2019). Comprehensive Integration of Single-Cell Data. 
Cell. doi:10.1016/j.cell.2019.05.031
```

Use this skill to perform comprehensive, reproducible single-cell RNA-seq analysis with Seurat.
