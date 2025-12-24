---
name: pseudobulk-deseq2
description: "Pseudobulk analysis: aggregate scRNA-seq to pseudobulk, run DESeq2 for differential expression, handle batch effects, compare conditions across cell types."
---

# Pseudobulk DESeq2 Analysis

## Overview

Pseudobulk analysis aggregates single-cell RNA-seq data by summing counts within biological replicates and cell types, then applies bulk RNA-seq differential expression methods like DESeq2. This approach provides better statistical power and controls false discovery rates compared to single-cell-level tests.

## When to Use This Skill

This skill should be used when:
- Comparing gene expression between conditions in scRNA-seq
- You have biological replicates for each condition
- You want robust statistical testing with better FDR control
- Analyzing treatment effects within specific cell types
- Dealing with complex experimental designs (batch, time course)
- Need to account for sample-level variation
- Publishing differential expression results (recommended by recent benchmarks)

## Quick Start Guide

### Basic Pseudobulk DESeq2 Workflow

```r
library(Seurat)
library(DESeq2)
library(dplyr)
library(tibble)

# 1. Load Seurat object
seurat_obj <- readRDS("seurat_annotated.rds")

# 2. Aggregate to pseudobulk
# Group by: condition + replicate + celltype
pseudobulk <- AggregateExpression(
  seurat_obj,
  assays = "RNA",
  slot = "counts",  # IMPORTANT: Use raw counts
  return.seurat = TRUE,
  group.by = c("condition", "replicate", "celltype")
)

# 3. Extract counts and metadata
counts <- GetAssayData(pseudobulk, layer = "counts")
samples <- colnames(counts)

# Parse sample names (e.g., "control_1_Tcells")
metadata <- tibble(sample = samples) %>%
  separate(sample, into = c("condition", "replicate", "celltype"), 
           sep = "_", remove = FALSE)

# 4. Run DESeq2 for specific cell type
celltype_of_interest <- "Tcells"

# Filter to cell type
keep_samples <- metadata$celltype == celltype_of_interest
counts_subset <- counts[, keep_samples]
meta_subset <- metadata[keep_samples, ]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts_subset,
  colData = meta_subset,
  design = ~ condition
)

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("condition", "treated", "control"))

# Convert to data frame
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

# View significant genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

print(head(sig_genes, 20))
```

## Core Workflow Steps

### 1. Prepare Seurat Object

```r
library(Seurat)

# Load data
seurat_obj <- readRDS("seurat_object.rds")

# Ensure you have necessary metadata
# Required columns: condition, replicate, celltype
head(seurat_obj@meta.data)

# Add replicate column if needed
seurat_obj$replicate <- seurat_obj$orig.ident

# Add condition column if needed
seurat_obj$condition <- ifelse(
  seurat_obj$treatment == "drug",
  "treated",
  "control"
)

# Check cell type annotations
table(seurat_obj$celltype)
table(seurat_obj$condition, seurat_obj$replicate)
```

### 2. Aggregate to Pseudobulk

```r
# Method 1: Using Seurat AggregateExpression
pseudobulk <- AggregateExpression(
  seurat_obj,
  assays = "RNA",
  slot = "counts",       # CRITICAL: Use raw counts, not normalized
  return.seurat = TRUE,
  group.by = c("condition", "replicate", "celltype"),
  verbose = TRUE
)

# Verify counts are integers
counts_matrix <- GetAssayData(pseudobulk, layer = "counts")
all(counts_matrix@x %% 1 == 0)  # Should be TRUE

# Method 2: Manual aggregation (more control)
library(Matrix)

# Create grouping variable
seurat_obj$sample_id <- paste(
  seurat_obj$condition,
  seurat_obj$replicate,
  seurat_obj$celltype,
  sep = "_"
)

# Get raw counts
raw_counts <- GetAssayData(seurat_obj, layer = "counts")

# Aggregate by sample_id
unique_samples <- unique(seurat_obj$sample_id)
pseudobulk_counts <- sapply(unique_samples, function(sample) {
  cells <- WhichCells(seurat_obj, expression = sample_id == sample)
  rowSums(raw_counts[, cells])
})

colnames(pseudobulk_counts) <- unique_samples
```

### 3. Create Metadata

```r
library(tidyr)
library(dplyr)

# Extract sample names
samples <- colnames(pseudobulk_counts)

# Parse sample names
metadata <- tibble(sample = samples) %>%
  separate(sample, 
           into = c("condition", "replicate", "celltype"),
           sep = "_",
           remove = FALSE) %>%
  mutate(
    replicate = as.integer(replicate),
    condition = factor(condition),
    celltype = factor(celltype)
  )

# Verify alignment
stopifnot(identical(metadata$sample, colnames(pseudobulk_counts)))

# Add additional metadata if needed
metadata$age_group <- ifelse(grepl("young", metadata$replicate), "young", "old")
```

### 4. Run DESeq2 for Each Cell Type

```r
library(DESeq2)

# Function to run DESeq2 for one cell type
run_deseq2_celltype <- function(counts, metadata, celltype, 
                                 contrast_var, ident1, ident2) {
  
  # Filter to cell type
  keep <- metadata$celltype == celltype
  if (sum(keep) < 3) return(NULL)  # Need at least 3 samples
  
  counts_sub <- counts[, keep, drop = FALSE]
  meta_sub <- metadata[keep, , drop = FALSE]
  
  # Remove zero-count genes
  keep_genes <- rowSums(counts_sub >= 1) >= 2
  counts_sub <- counts_sub[keep_genes, , drop = FALSE]
  
  if (nrow(counts_sub) < 10) return(NULL)
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData = meta_sub,
    design = as.formula(paste("~", contrast_var))
  )
  
  # Run DESeq2
  dds <- DESeq(dds, fitType = "local")
  
  # Get results
  res <- results(dds, contrast = c(contrast_var, ident1, ident2))
  
  # Shrink log2 fold changes
  res_shrink <- lfcShrink(dds, 
                          contrast = c(contrast_var, ident1, ident2),
                          res = res,
                          type = "ashr")
  
  # Convert to data frame
  res_df <- as.data.frame(res_shrink) %>%
    rownames_to_column("gene") %>%
    mutate(
      celltype = celltype,
      contrast = paste(ident1, "vs", ident2)
    )
  
  return(res_df)
}

# Run for all cell types
celltypes <- unique(metadata$celltype)

all_results <- lapply(celltypes, function(ct) {
  cat("Processing", ct, "...\n")
  run_deseq2_celltype(
    counts = pseudobulk_counts,
    metadata = metadata,
    celltype = ct,
    contrast_var = "condition",
    ident1 = "treated",
    ident2 = "control"
  )
})

# Combine results
combined_results <- bind_rows(all_results)

# Save
write.csv(combined_results, "pseudobulk_deseq2_results.csv", row.names = FALSE)
```

### 5. Filter and Visualize Results

```r
library(ggplot2)
library(EnhancedVolcano)

# Filter significant genes
sig_results <- combined_results %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

# Count DEGs per cell type
deg_counts <- sig_results %>%
  group_by(celltype) %>%
  summarise(
    n_up = sum(log2FoldChange > 1),
    n_down = sum(log2FoldChange < -1),
    total = n()
  )

print(deg_counts)

# Volcano plot for specific cell type
celltype_results <- combined_results %>%
  filter(celltype == "Tcells")

EnhancedVolcano(
  celltype_results,
  lab = celltype_results$gene,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'T cells: Treated vs Control',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2,
  labSize = 4
)

# Heatmap of top DEGs
library(ComplexHeatmap)
library(circlize)

top_genes <- sig_results %>%
  group_by(celltype) %>%
  slice_max(order_by = abs(log2FoldChange), n = 10) %>%
  pull(gene) %>%
  unique()

# Get normalized counts
dds_normalized <- vst(dds, blind = FALSE)
mat <- assay(dds_normalized)[top_genes, ]

Heatmap(
  mat,
  name = "Expression",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE
)
```

## Advanced Features

### Time Course Analysis

```r
# Design with time as continuous variable
metadata$time <- as.numeric(gsub("h", "", metadata$timepoint))

dds <- DESeqDataSetFromMatrix(
  countData = counts_sub,
  colData = meta_sub,
  design = ~ condition + time + condition:time
)

dds <- DESeq(dds, test = "LRT", reduced = ~ condition + time)

# Genes with significant interaction
res_time <- results(dds)
sig_interaction <- res_time %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(padj < 0.05)
```

### Batch Correction

```r
# Include batch in design
dds <- DESeqDataSetFromMatrix(
  countData = counts_sub,
  colData = meta_sub,
  design = ~ batch + condition
)

# Or use limma's removeBatchEffect for visualization
library(limma)
vsd <- vst(dds, blind = FALSE)
mat_corrected <- removeBatchEffect(
  assay(vsd),
  batch = vsd$batch,
  design = model.matrix(~ condition, colData(vsd))
)
```

### Multiple Comparisons

```r
# All pairwise comparisons
conditions <- unique(metadata$condition)
comparisons <- combn(conditions, 2, simplify = FALSE)

all_comparisons <- lapply(comparisons, function(pair) {
  run_deseq2_celltype(
    counts = pseudobulk_counts,
    metadata = metadata,
    celltype = "Tcells",
    contrast_var = "condition",
    ident1 = pair[1],
    ident2 = pair[2]
  )
})

combined_comparisons <- bind_rows(all_comparisons)
```

## Best Practices

### Sample Size
- Need at least 3 biological replicates per condition
- More replicates = better power and FDR control
- Technical replicates should be treated as biological

### Quality Control
- Remove low-quality cells before aggregation
- Check for outlier samples after aggregation
- Ensure balanced sample sizes across conditions

### Gene Filtering
- Remove genes with very low counts (< 10 total)
- Keep genes detected in at least 2-3 samples
- Consider minimum cell count per sample

### Statistical Considerations
- Use LFC shrinkage (ashr, apeglm) for more accurate effect sizes
- Set appropriate alpha (typically 0.05 or 0.01)
- Consider independent filtering for power
- Report FC threshold used (typically |log2FC| > 0.5 or 1)

## Common Pitfalls to Avoid

1. **Using normalized counts**: Always use RAW counts for DESeq2
2. **No replicates**: Need biological replicates for DESeq2
3. **Wrong slot**: Use `slot = "counts"`, not "data" or "scale.data"
4. **Ignoring batch effects**: Include batch in design if present
5. **Too few cells**: Need sufficient cells per sample (>10-20)
6. **No gene filtering**: Remove very low count genes
7. **Comparing wrong groups**: Verify factor levels and contrasts
8. **Not using shrinkage**: Always use LFC shrinkage for fold changes

## Resources

### Scripts Directory
- `pseudobulk_helpers.R` - Helper functions for aggregation and analysis

### References Directory
- `pseudobulk_examples.md` - Complete workflow examples
- `troubleshooting.md` - Common issues and solutions

## Citation

```
Love et al. (2014). Moderated estimation of fold change and dispersion 
for RNA-seq data with DESeq2. Genome Biology.

Squair et al. (2021). Confronting false discoveries in single-cell 
differential expression. Nature Communications.
```

Use this skill for robust differential expression analysis in scRNA-seq experiments.
