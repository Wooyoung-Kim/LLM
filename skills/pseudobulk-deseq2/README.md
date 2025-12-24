# Pseudobulk DESeq2 Analysis

Aggregate scRNA-seq to pseudobulk and run DESeq2 for differential expression.

## Quick Start

```r
library(Seurat)
library(DESeq2)

# Aggregate to pseudobulk
pseudobulk <- AggregateExpression(
  seurat_obj,
  assays = "RNA",
  slot = "counts",  # Important: raw counts
  return.seurat = TRUE,
  group.by = c("condition", "replicate", "celltype")
)

# Extract counts
counts <- GetAssayData(pseudobulk, layer = "counts")

# Run DESeq2 for specific cell type
# (See SKILL.md for complete function)
```

## Key Features

- **Pseudobulk aggregation**: Sum counts by replicate and cell type
- **DESeq2 analysis**: Robust statistical testing
- **Multiple cell types**: Analyze all cell types
- **Batch correction**: Account for batch effects
- **Time course**: Analyze temporal dynamics

## Why Pseudobulk?

- Better FDR control than single-cell tests
- Accounts for sample-level variation
- More statistical power with replicates
- Recommended by recent benchmarks

## Requirements

- At least 3 biological replicates per condition
- Raw counts (not normalized)
- Cell type annotations

## Citation

```
Love et al. (2014). DESeq2. Genome Biology.
Squair et al. (2021). Confronting false discoveries in single-cell DE. 
Nature Communications.
```

## Documentation

- Main skill: [SKILL.md](SKILL.md)
- Examples: [references/pseudobulk_examples.md](references/pseudobulk_examples.md)
