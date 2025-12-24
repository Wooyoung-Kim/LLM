# Seurat Workflow Examples

## Example 1: Complete Standard scRNA-seq Analysis

### PBMC 3k Dataset Analysis

```r
library(Seurat)
library(dplyr)
library(ggplot2)

# 1. Load data
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", 
                           min.cells = 3, min.features = 200)

# 2. QC and filtering
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
               percent.mt < 5)

# 3. Normalization
pbmc <- NormalizeData(pbmc)

# 4. Feature selection
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 5. Scaling
pbmc <- ScaleData(pbmc)

# 6. PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualize PCA
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Elbow plot
ElbowPlot(pbmc, ndims = 50)

# 7. Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# 8. UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize
DimPlot(pbmc, reduction = "umap", label = TRUE)

# 9. Find markers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold = 0.25)

# Top markers per cluster
top10 <- pbmc.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

# 10. Visualization of markers
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# 11. Cell type annotation
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", 
                     "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$celltype <- Idents(pbmc)

# Final visualization
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# 12. Save
saveRDS(pbmc, file = "pbmc_final.rds")
```

## Example 2: Integration of Multiple Samples with Harmony

### Integrating Treatment and Control Samples

```r
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)

# Load multiple samples
ctrl <- Read10X(data.dir = "data/ctrl/filtered_feature_bc_matrix/")
stim <- Read10X(data.dir = "data/stim/filtered_feature_bc_matrix/")

# Create Seurat objects
ctrl <- CreateSeuratObject(counts = ctrl, project = "CTRL", min.cells = 5)
stim <- CreateSeuratObject(counts = stim, project = "STIM", min.cells = 5)

# Add metadata
ctrl$condition <- "ctrl"
stim$condition <- "stim"

# Merge
immune.combined <- merge(ctrl, y = stim, add.cell.ids = c("CTRL", "STIM"))

# QC
immune.combined[["percent.mt"]] <- PercentageFeatureSet(immune.combined, pattern = "^MT-")
immune.combined <- subset(immune.combined, subset = nFeature_RNA > 200 & 
                          nFeature_RNA < 2500 & percent.mt < 5)

# Standard workflow
immune.combined <- NormalizeData(immune.combined)
immune.combined <- FindVariableFeatures(immune.combined, selection.method = "vst", 
                                         nfeatures = 2000)
immune.combined <- ScaleData(immune.combined)
immune.combined <- RunPCA(immune.combined, npcs = 30)

# Before integration
p1 <- DimPlot(object = immune.combined, reduction = "pca", group.by = "condition")
p2 <- VlnPlot(object = immune.combined, features = "PC_1", group.by = "condition")

# Run Harmony
immune.combined <- RunHarmony(immune.combined, 
                               group.by.vars = "condition",
                               plot_convergence = TRUE)

# After integration
p3 <- DimPlot(object = immune.combined, reduction = "harmony", group.by = "condition")
p4 <- VlnPlot(object = immune.combined, features = "harmony_1", group.by = "condition")

# UMAP and clustering on integrated data
immune.combined <- RunUMAP(immune.combined, reduction = "harmony", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "harmony", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Visualize
DimPlot(immune.combined, reduction = "umap", group.by = "condition")
DimPlot(immune.combined, reduction = "umap", label = TRUE)
DimPlot(immune.combined, reduction = "umap", split.by = "condition")

# Find conserved markers
markers <- FindConservedMarkers(immune.combined, 
                                ident.1 = 0,
                                grouping.var = "condition")

saveRDS(immune.combined, "immune_combined_harmony.rds")
```

## Example 3: Differential Expression Between Conditions

### Within-Cell-Type DE Analysis

```r
library(Seurat)
library(MAST)
library(dplyr)

# Load integrated object
seurat_obj <- readRDS("integrated_seurat.rds")

# Focus on specific cell type
tcells <- subset(seurat_obj, subset = celltype == "CD8 T")

# Set identity to condition
Idents(tcells) <- "condition"

# Find DEGs between treated and control
de_genes <- FindMarkers(tcells, 
                        ident.1 = "treated",
                        ident.2 = "control",
                        test.use = "MAST",
                        logfc.threshold = 0.25,
                        min.pct = 0.1)

# Add gene names
de_genes$gene <- rownames(de_genes)

# Filter significant genes
sig_genes <- de_genes %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5)

# Volcano plot
library(EnhancedVolcano)
EnhancedVolcano(de_genes,
                lab = rownames(de_genes),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Treated vs Control in CD8 T cells',
                pCutoff = 0.05,
                FCcutoff = 0.5)

# Visualize top DEGs
top_genes <- sig_genes %>%
  arrange(p_val_adj) %>%
  head(20) %>%
  pull(gene)

# Heatmap
DoHeatmap(tcells, features = top_genes, group.by = "condition")

# Violin plots
VlnPlot(tcells, features = top_genes[1:6], split.by = "condition", ncol = 3)

# Feature plots
FeaturePlot(tcells, features = top_genes[1:4], split.by = "condition")

# Save results
write.csv(sig_genes, "CD8_T_cells_DEGs.csv", row.names = FALSE)
```

## Example 4: Marker Discovery and Visualization

### Finding and Visualizing Cell Type Markers

```r
library(Seurat)
library(dplyr)
library(ggplot2)

seurat_obj <- readRDS("clustered_seurat.rds")

# Find all markers
all.markers <- FindAllMarkers(seurat_obj, 
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25,
                              test.use = "wilcox")

# Top 5 markers per cluster
top5 <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Dot plot
DotPlot(seurat_obj, features = unique(top5$gene)) + 
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Heatmap of top markers
top10 <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DoHeatmap(seurat_obj, features = top10$gene, size = 3) + NoLegend()

# Feature plots for key markers
canonical_markers <- c(
  "CD3D", "CD8A",     # T cells
  "CD14", "FCGR3A",   # Monocytes
  "MS4A1", "CD79A",   # B cells
  "GNLY", "NKG7",     # NK cells
  "FCER1A", "CST3"    # DCs
)

FeaturePlot(seurat_obj, features = canonical_markers, ncol = 4)

# Violin plots
VlnPlot(seurat_obj, features = canonical_markers, ncol = 5, pt.size = 0)

# Ridge plots for distribution
library(ggridges)
RidgePlot(seurat_obj, features = canonical_markers[1:6], ncol = 2)

# Percentage of cells expressing markers
marker_pct <- FoldChange(seurat_obj, ident.1 = levels(seurat_obj))

# Save markers
write.csv(all.markers, "all_cluster_markers.csv", row.names = FALSE)
write.csv(top10, "top10_markers_per_cluster.csv", row.names = FALSE)
```

## Example 5: Subset and Re-cluster Specific Cell Types

### Deep Dive into Myeloid Cells

```r
library(Seurat)
library(dplyr)

# Load annotated object
seurat_obj <- readRDS("annotated_seurat.rds")

# Subset myeloid cells
myeloid <- subset(seurat_obj, subset = celltype %in% c("CD14+ Mono", "FCGR3A+ Mono", "DC"))

# Re-normalize and re-cluster
myeloid <- NormalizeData(myeloid)
myeloid <- FindVariableFeatures(myeloid, selection.method = "vst", nfeatures = 2000)
myeloid <- ScaleData(myeloid)
myeloid <- RunPCA(myeloid, npcs = 30)

# Elbow plot
ElbowPlot(myeloid, ndims = 30)

# Clustering with higher resolution
myeloid <- FindNeighbors(myeloid, dims = 1:20)
myeloid <- FindClusters(myeloid, resolution = c(0.3, 0.5, 0.8, 1.0))

# UMAP
myeloid <- RunUMAP(myeloid, dims = 1:20)

# Compare resolutions
library(clustree)
clustree(myeloid, prefix = "RNA_snn_res.")

# Select appropriate resolution
myeloid <- FindClusters(myeloid, resolution = 0.5)

# Visualize
DimPlot(myeloid, reduction = "umap", label = TRUE)
DimPlot(myeloid, reduction = "umap", group.by = "celltype")
DimPlot(myeloid, reduction = "umap", split.by = "condition")

# Find markers for new clusters
myeloid.markers <- FindAllMarkers(myeloid, only.pos = TRUE, 
                                  min.pct = 0.25, logfc.threshold = 0.25)

# Myeloid-specific markers
myeloid_genes <- c("CD14", "FCGR3A", "CD16", "S100A8", "S100A9", 
                   "LYZ", "CD68", "CD163", "CLEC9A", "FCER1A")

DotPlot(myeloid, features = myeloid_genes) + RotatedAxis()

# Annotate refined clusters
new.ids <- c("CD14+ Mono 1", "CD14+ Mono 2", "CD16+ Mono", 
             "cDC1", "cDC2", "pDC")
names(new.ids) <- levels(myeloid)
myeloid <- RenameIdents(myeloid, new.ids)
myeloid$refined_celltype <- Idents(myeloid)

# Final visualization
DimPlot(myeloid, reduction = "umap", label = TRUE, repel = TRUE)

# Save refined myeloid object
saveRDS(myeloid, "myeloid_refined.rds")
```

## Example 6: SCTransform Workflow

### Using SCTransform for Normalization

```r
library(Seurat)
library(ggplot2)

# Load data
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc", 
                           min.cells = 3, min.features = 200)

# QC
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
               percent.mt < 5)

# SCTransform (replaces NormalizeData, FindVariableFeatures, ScaleData)
pbmc <- SCTransform(pbmc, 
                    vars.to.regress = "percent.mt",
                    verbose = FALSE)

# Note: After SCTransform, default assay is "SCT"
DefaultAssay(pbmc)  # "SCT"

# PCA
pbmc <- RunPCA(pbmc, npcs = 30, verbose = FALSE)

# Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# UMAP
pbmc <- RunUMAP(pbmc, dims = 1:30)

# Visualization
DimPlot(pbmc, reduction = "umap", label = TRUE)

# For visualization, switch to RNA assay for raw counts
DefaultAssay(pbmc) <- "RNA"
FeaturePlot(pbmc, features = c("CD3D", "CD14", "MS4A1"))

# For DE analysis, use SCT assay
DefaultAssay(pbmc) <- "SCT"
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25)

saveRDS(pbmc, "pbmc_sctransform.rds")
```

## Example 7: Exporting to AnnData

### Preparing Data for Python/Scanpy

```r
library(Seurat)
library(MuDataSeurat)

# Load Seurat object
seurat_obj <- readRDS("final_seurat_object.rds")

# Prepare for export - slim down object
seurat_obj <- DietSeurat(
  seurat_obj,
  counts = TRUE,        # Keep raw counts
  data = TRUE,          # Keep normalized data
  scale.data = FALSE,   # Don't export scaled data
  features = rownames(seurat_obj),  # All genes
  assays = "RNA",       # Which assay to export
  dimreducs = c("pca", "umap"),
  graphs = c("RNA_nn", "RNA_snn"),
  misc = TRUE
)

# For Seurat v5: Join layers first
seurat_obj <- JoinLayers(seurat_obj)

# Check what will be exported
DefaultAssay(seurat_obj)
Assays(seurat_obj)
names(seurat_obj@reductions)

# Export to h5ad
MuDataSeurat::WriteH5AD(seurat_obj, "seurat_export.h5ad", assay = "RNA")

# In Python, you can now read this:
# import scanpy as sc
# adata = sc.read_h5ad('seurat_export.h5ad')
# adata.X # normalized data
# adata.layers['counts'] # raw counts
# adata.obsm['X_pca'] # PCA
# adata.obsm['X_umap'] # UMAP
```

## Example 8: Cell Cycle Scoring and Regression

### Removing Cell Cycle Effects

```r
library(Seurat)

seurat_obj <- readRDS("seurat_obj.rds")

# Load cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Score cells for cell cycle
seurat_obj <- CellCycleScoring(seurat_obj, 
                                s.features = s.genes,
                                g2m.features = g2m.genes,
                                set.ident = TRUE)

# View cell cycle scores
head(seurat_obj@meta.data)

# Visualize cell cycle
DimPlot(seurat_obj, reduction = "umap", group.by = "Phase")

# Option 1: Regress out cell cycle during SCTransform
seurat_obj <- SCTransform(seurat_obj, 
                          vars.to.regress = c("S.Score", "G2M.Score"),
                          verbose = FALSE)

# Option 2: Regress out during ScaleData (if using standard workflow)
seurat_obj <- ScaleData(seurat_obj, 
                        vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))

# Continue with standard workflow
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Check if cell cycle effect is removed
DimPlot(seurat_obj, reduction = "umap", group.by = "Phase")

saveRDS(seurat_obj, "seurat_cc_regressed.rds")
```
