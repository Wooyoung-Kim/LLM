---
name: celltype-annotation
description: "Automated cell type annotation: ScType, Scibet, scMayoMap for marker-based, reference-based, and database-driven cell type identification in scRNA-seq."
---

# Cell Type Annotation

## Overview

Automated cell type annotation assigns cell identities to clusters in single-cell RNA-seq data based on gene expression patterns. This skill covers three complementary approaches: ScType (marker-based), Scibet (reference-based machine learning), and scMayoMap (Mayo Clinic tissue atlas-based).

## When to Use This Skill

This skill should be used when:
- You have clustered scRNA-seq data that needs cell type labels
- Manual annotation is time-consuming or inconsistent
- You need reference-based validation of cell types
- You want to compare multiple annotation methods
- You're working with well-characterized tissues (blood, immune, brain, etc.)
- You need tissue-specific cell type databases

## Quick Start Guide

### Basic ScType Annotation

```r
library(Seurat)
library(dplyr)

# Load annotation functions
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/I anevskiAleksandr/sc-type/master/R/sctype_score_.R")

# Load Seurat object
seurat_obj <- readRDS("seurat_clustered.rds")

# Load ScType database
db <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
tissue <- "Immune system"  # Options: Immune system, Brain, Lung, Liver, etc.

# Prepare gene sets
gs_list <- gene_sets_prepare(db, tissue)

# Extract expression matrix
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(seurat_obj[["RNA"]])))
scRNAseqData_scaled <- if (seurat_package_v5) {
  as.matrix(seurat_obj[["RNA"]]$scale.data)
} else {
  as.matrix(seurat_obj[["RNA"]]@scale.data)
}

# Run ScType scoring
es.max <- sctype_score(
  scRNAseqData = scRNAseqData_scaled,
  scaled = TRUE,
  gs = gs_list$gs_positive,
  gs2 = gs_list$gs_negative
)

# Assign cell types per cluster
cL_results <- do.call("rbind", lapply(unique(seurat_obj$seurat_clusters), function(cl) {
  es.max.cl <- sort(rowSums(es.max[, rownames(seurat_obj@meta.data[seurat_obj$seurat_clusters == cl, ])]), decreasing = TRUE)
  head(data.frame(
    cluster = cl,
    type = names(es.max.cl),
    scores = es.max.cl,
    ncells = sum(seurat_obj$seurat_clusters == cl)
  ), 10)
}))

sctype_scores <- cL_results %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = scores)

# Filter low-confidence predictions
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"

# Add to Seurat object
seurat_obj$sctype <- ""
for(j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[sctype_scores$cluster == j, ]
  seurat_obj$sctype[seurat_obj$seurat_clusters == j] <- as.character(cl_type$type[1])
}

# Visualize
DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = "sctype")
```

## Core Workflow Steps

### 1. ScType: Marker Gene-Based Annotation

#### Load ScType Functions

```r
# Load required libraries
library(Seurat)
library(dplyr)
library(HGNChelper)

# Load ScType functions from GitHub
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
```

#### Prepare Gene Sets

```r
# ScType database URL
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"

# Available tissues:
# "Immune system", "Pancreas", "Liver", "Eye", "Kidney", "Brain", 
# "Lung", "Adrenal", "Heart", "Intestine", "Muscle", "Placenta", 
# "Spleen", "Stomach", "Thymus"

tissue <- "Immune system"

# Prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)

# Check gene sets
names(gs_list)  # $gs_positive, $gs_negative
```

#### Run ScType Annotation

```r
# Check Seurat version (v4 vs v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(seurat_obj[["RNA"]])))
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# Extract scaled RNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) {
  as.matrix(seurat_obj[["RNA"]]$scale.data)
} else {
  as.matrix(seurat_obj[["RNA"]]@scale.data)
}

# Run ScType scoring
# For SCTransform: use seurat_obj[["SCT"]]$scale.data
es.max <- sctype_score(
  scRNAseqData = scRNAseqData_scaled,
  scaled = TRUE,
  gs = gs_list$gs_positive,
  gs2 = gs_list$gs_negative
)

# Aggregate scores by cluster
cL_results <- do.call("rbind", lapply(unique(seurat_obj$seurat_clusters), function(cl) {
  es.max.cl <- sort(
    rowSums(es.max[, rownames(seurat_obj@meta.data[seurat_obj$seurat_clusters == cl, ])]),
    decreasing = TRUE
  )
  head(data.frame(
    cluster = cl,
    type = names(es.max.cl),
    scores = es.max.cl,
    ncells = sum(seurat_obj$seurat_clusters == cl)
  ), 10)
}))

# Get top cell type per cluster
sctype_scores <- cL_results %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = scores)

# Set low-confidence clusters to "Unknown"
# Threshold: score must be >= ncells/4
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"

print(sctype_scores[, 1:3])
```

#### Add Annotations to Seurat

```r
# Create annotation column
seurat_obj$sctype <- ""

# Assign cell types
for(j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[sctype_scores$cluster == j, ]
  seurat_obj$sctype[seurat_obj$seurat_clusters == j] <- as.character(cl_type$type[1])
}

# Visualize
DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype')
```

### 2. Scibet: Machine Learning-Based Annotation

#### Load Reference Model

```r
library(scibet)
library(dplyr)

# Download reference model from Scibet website
# http://scibet.cancer-pku.cn/download_references.html

# Available reference models:
# - major_human_cell_types.csv (broad cell types)
# - human_cell_atlas.csv (detailed annotation)
# - mouse_brain.csv (brain-specific)
# - etc.

model <- readr::read_csv("major_human_cell_types.csv")
model <- pro.core(model)

head(model)
```

#### Prepare Query Data

```r
# Join layers (Seurat v5)
seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")

# Extract normalized counts
query <- seurat_obj@assays$RNA$data %>%
  t() %>%
  as.data.frame()

# Add cluster labels for gene selection
query_meta <- seurat_obj@meta.data %>%
  select(seurat_clusters)
colnames(query_meta) <- "label"
query <- cbind(query, query_meta)
```

#### Feature Selection and Prediction

```r
# E(ntropy)-test for supervised gene selection
etest_gene <- SelectGene(query, k = 50)
print(etest_gene)

# Visualize marker genes
Marker_heatmap(query, etest_gene)

# Remove label for prediction
ori_label <- query$label
query <- query[, -ncol(query)]

# Load model and predict
prd <- LoadModel(model)
predicted_labels <- prd(query)

# Add to Seurat object
seurat_obj$scibet <- predicted_labels

# Visualize
DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = "scibet")
```

#### Consensus Annotation per Cluster

```r
# Get most common cell type per cluster
cluster_celltype_map <- seurat_obj@meta.data %>%
  group_by(seurat_clusters, scibet) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(seurat_clusters, desc(count)) %>%
  group_by(seurat_clusters) %>%
  slice(1) %>%
  select(seurat_clusters, scibet)

colnames(cluster_celltype_map) <- c("seurat_clusters", "scibet_consensus")

# Add consensus to metadata
seurat_meta <- seurat_obj@meta.data %>%
  rownames_to_column("cells") %>%
  full_join(cluster_celltype_map, by = "seurat_clusters") %>%
  column_to_rownames("cells")

seurat_obj@meta.data <- seurat_meta

# Visualize consensus
DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = "scibet_consensus")
```

### 3. scMayoMap: Mayo Clinic Atlas-Based Annotation

#### Load scMayoMap

```r
library(scMayoMap)
library(dplyr)
library(tibble)

# Check available tissues
scMayoMapDatabase$tissue %>% unique()

# Options include:
# "blood", "bone marrow", "spleen", "brain", "liver", "lung", "pancreas", etc.
```

#### Run scMayoMap on Marker Genes

```r
# First, find all markers
Idents(seurat_obj) <- "seurat_clusters"
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25)

# Save markers
write.csv(markers, "cluster_markers.csv", row.names = FALSE)

# Run scMayoMap
obj <- scMayoMap(data = markers, tissue = 'blood')

# Get annotations
mayo_markers <- obj$markers

# Get top annotation per cluster
top_markers <- mayo_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = score, n = 1) %>%
  ungroup() %>%
  select(cluster, celltype, genes)

colnames(top_markers) <- c("seurat_clusters", "mayo", "mayo_gene")

# Add to Seurat metadata
seurat_meta <- seurat_obj@meta.data %>%
  rownames_to_column("cells") %>%
  full_join(top_markers, by = "seurat_clusters") %>%
  column_to_rownames("cells")

seurat_obj@meta.data <- seurat_meta

# Visualize
DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = "mayo")
```

### 4. Celltypist: Python-Based Deep Learning Annotation

#### Overview

Celltypist is a Python tool that uses machine learning for automated cell type annotation. It's particularly effective for immune cells and provides hierarchical predictions.

**Note**: Celltypist runs in Python, not R. You'll need to:
1. Export Seurat object to h5ad format
2. Run celltypist in Python
3. Import results back to R

#### Export Seurat to h5ad

```r
library(MuDataSeurat)
library(Seurat)

# Prepare Seurat object
seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")

# Must use LogNormalize with scale.factor = 10000
# (Celltypist assumes this normalization)
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# Export to h5ad
MuDataSeurat::WriteH5AD(seurat_obj, "seurat_data.h5ad", assay = "RNA")
```

#### Run Celltypist in Python

```python
import celltypist
import scanpy as sc

# Load data
adata = sc.read_h5ad("seurat_data.h5ad")

# Download model (first time only)
# Available models: 'Immune_All_Low.pkl', 'Immune_All_High.pkl', 
# 'AIFI_L1.pkl', 'AIFI_L2.pkl', 'AIFI_L3.pkl', etc.
celltypist.models.download_models(model = 'Immune_All_High.pkl')

# Load model
model = celltypist.models.Model.load(model = 'Immune_All_High.pkl')

# Run prediction
predictions = celltypist.annotate(
    adata,
    model = 'Immune_All_High.pkl',
    majority_voting = True  # Use cluster-level consensus
)

# Add predictions to adata
adata.obs['celltypist_prediction'] = predictions.predicted_labels
adata.obs['celltypist_conf_score'] = predictions.conf_scores

# Save results
adata.write_h5ad("seurat_celltypist.h5ad")
```

#### Import Results to R

```r
library(MuDataSeurat)
library(Seurat)

# Read h5ad with celltypist annotations
seurat_annotated <- ReadH5AD("seurat_celltypist.h5ad")

# Or add predictions to existing Seurat object
predictions <- read.csv("celltypist_predictions.csv", row.names = 1)
seurat_obj$celltypist <- predictions$celltypist_prediction
seurat_obj$celltypist_score <- predictions$celltypist_conf_score

# Create cluster-level consensus
cluster_consensus <- seurat_obj@meta.data %>%
  group_by(seurat_clusters, celltypist) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  slice_max(count, n = 1) %>%
  select(seurat_clusters, celltypist)

colnames(cluster_consensus)[2] <- "celltypist_consensus"

# Add to metadata
seurat_meta <- seurat_obj@meta.data %>%
  rownames_to_column("cells") %>%
  left_join(cluster_consensus, by = "seurat_clusters") %>%
  column_to_rownames("cells")

seurat_obj@meta.data <- seurat_meta

# Visualize
DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = "celltypist_consensus")
```

#### Available Celltypist Models

```python
# List all available models
celltypist.models.models_description()

# Common models:
# - 'Immune_All_Low.pkl': Low-resolution immune cells
# - 'Immune_All_High.pkl': High-resolution immune cells
# - 'AIFI_L1.pkl': AIFI Level 1 (broad categories)
# - 'AIFI_L2.pkl': AIFI Level 2 (medium resolution)
# - 'AIFI_L3.pkl': AIFI Level 3 (fine resolution)
# - 'Healthy_Adult_Heart.pkl': Heart-specific
# - 'Cells_Intestinal_Tract.pkl': Intestine-specific
# - 'Healthy_COVID19_PBMC.pkl': PBMC with COVID samples
```

### 5. Compare Multiple Annotation Methods

#### Side-by-Side Comparison

```r
library(patchwork)

# Create comparison plots
p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = "sctype") +
  ggtitle("ScType") + NoLegend()

p2 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = "scibet_consensus") +
  ggtitle("Scibet") + NoLegend()

p3 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = "mayo") +
  ggtitle("scMayoMap") + NoLegend()

# Combined plot
(p1 | p2 | p3)
```

#### Create Consensus Annotation

```r
# Create summary table
annotation_summary <- seurat_obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    sctype_annot = names(which.max(table(sctype))),
    scibet_annot = names(which.max(table(scibet_consensus))),
    mayo_annot = names(which.max(table(mayo))),
    n_cells = n()
  )

print(annotation_summary)

# Export for manual curation
write.csv(annotation_summary, "annotation_comparison.csv", row.names = FALSE)

# Manual final annotation
# Read curated annotations
final_annot <- read.csv("annotation_final.csv")

# Add to Seurat
seurat_meta <- seurat_obj@meta.data %>%
  rownames_to_column("cells") %>%
  full_join(final_annot, by = "seurat_clusters") %>%
  column_to_rownames("cells")

seurat_obj@meta.data <- seurat_meta

# Final visualization
DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = "final_celltype")
```

## Advanced Features

### Custom Marker Database for ScType

```r
# Create custom marker database in Excel format
# Columns: tissueType, cellName, geneSymbolmore1, geneSymbolmore2

# Or load from local file
db_custom <- "/path/to/custom_markers.xlsx"
gs_list_custom <- gene_sets_prepare(db_custom, "Custom Tissue")

# Run ScType with custom database
es.max_custom <- sctype_score(
  scRNAseqData = scRNAseqData_scaled,
  scaled = TRUE,
  gs = gs_list_custom$gs_positive,
  gs2 = gs_list_custom$gs_negative
)
```

### Subset-Specific Annotation

```r
# Annotate T cells only
tcells <- subset(seurat_obj, subset = broad_celltype == "T cells")

# Use Immune system database with T cell focus
tissue <- "Immune system"
gs_list <- gene_sets_prepare(db_, tissue)

# Run annotation on subset
# ... (same workflow as above)

# Add back to main object
seurat_obj$tcell_subtype <- NA
seurat_obj$tcell_subtype[colnames(tcells)] <- tcells$sctype
```

### Confidence Scoring

```r
# Calculate annotation confidence
# Higher score = more confident

calculate_confidence <- function(sctype_scores, threshold_multiplier = 4) {
  sctype_scores %>%
    mutate(
      confidence = case_when(
        scores >= ncells * threshold_multiplier ~ "High",
        scores >= ncells * threshold_multiplier / 2 ~ "Medium",
        TRUE ~ "Low"
      )
    )
}

sctype_scores_conf <- calculate_confidence(sctype_scores)

# Visualize confidence
seurat_obj$annotation_confidence <- "Unknown"
for(j in unique(sctype_scores_conf$cluster)) {
  conf <- sctype_scores_conf$confidence[sctype_scores_conf$cluster == j]
  seurat_obj$annotation_confidence[seurat_obj$seurat_clusters == j] <- conf
}

DimPlot(seurat_obj, reduction = "umap", group.by = "annotation_confidence")
```

## Best Practices

### Choosing Annotation Method

**Use ScType when**:
- You have limited reference data
- Working with standard tissues (immune, brain, lung, etc.)
- Need quick, automated annotation
- Want marker gene interpretability

**Use Scibet when**:
- You have high-quality reference datasets
- Need species-specific or custom references
- Want machine learning-based confidence
- Have similar experimental conditions to reference

**Use scMayoMap when**:
- Working with tissues in Mayo Clinic atlas
- Need tissue-specific annotations
- Want clinically-validated cell types

**Use Multiple Methods when**:
- You want consensus validation
- Annotation is challenging (mixed tissues, novel cell types)
- Need high confidence for critical conclusions

### Quality Control

```r
# 1. Check marker expression
# Validate that predicted cell types express expected markers

VlnPlot(seurat_obj, features = c("CD3D", "CD4", "CD8A"), group.by = "sctype")

# 2. Check cluster purity
# Clusters should be mostly one cell type

purity <- seurat_obj@meta.data %>%
  group_by(seurat_clusters, sctype) %>%
  summarise(count = n()) %>%
  group_by(seurat_clusters) %>%
  mutate(purity = max(count) / sum(count))

# 3. Compare methods
# High agreement = high confidence

agreement <- seurat_obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    sctype_scibet_agree = sum(sctype == scibet_consensus) / n(),
    sctype_mayo_agree = sum(sctype == mayo) / n()
  )
```

### Common Issues

1. **Low scores/Unknown predictions**:
   - Try different tissue databases
   - Check clustering resolution (too high = mixed clusters)
   - Validate with known markers

2. **Inconsistent annotations**:
   - Use consensus from multiple methods
   - Manually curate based on marker genes
   - Check for batch effects

3. **Novel cell types**:
   - ScType may mark as "Unknown"
   - Use de novo marker discovery
   - Create custom marker database

## Resources

### Scripts Directory
- `annotation_helpers.R` - Helper functions for all three methods

### References Directory
- `annotation_examples.md` - Complete workflow examples
- `marker_databases.md` - Available marker databases and tissues

## Citation

```bibtex
@article{ianevski2022sctype,
  title={Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data},
  author={Ianevski, Aleksandr and others},
  journal={Nature Communications},
  year={2022}
}

@article{li2020scibet,
  title={Reference-free cell type deconvolution of multi-cellular pixel-resolution spatially resolved transcriptomics data},
  author={Li, Huipeng and others},
  journal={Nature Communications},
  year={2020}
}

@article{mayo2024,
  title={scMayoMap: atlas-scale mapping for single-cell RNA-seq data},
  author={Mayo Clinic},
  year={2024}
}
```

Use this skill for rapid, accurate cell type annotation in scRNA-seq datasets.
