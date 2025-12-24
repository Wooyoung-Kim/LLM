# Cell Type Annotation

Automated cell type identification using ScType, Scibet, and scMayoMap.

## Quick Start

```r
library(Seurat)
library(dplyr)

# Load ScType
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# Prepare database
db <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
gs_list <- gene_sets_prepare(db, "Immune system")

# Get expression matrix
scData <- as.matrix(seurat_obj[["RNA"]]$scale.data)

# Score and assign
es.max <- sctype_score(scData, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# Add to Seurat
# (See SKILL.md for complete workflow)
```

## Key Features

- **ScType**: Marker gene database-based annotation
- **Scibet**: Machine learning reference-based prediction  
- **scMayoMap**: Mayo Clinic tissue atlas annotation
- **Celltypist**: Python-based deep learning (immune cells, hierarchical)
- **Consensus**: Compare multiple methods for validation

## Four Complementary Approaches

### 1. ScType (Marker-Based)
- Uses curated marker gene databases
- 15+ tissue types supported
- Fast, interpretable results

### 2. Scibet (ML Reference-Based)
- Trains on reference datasets
- Entropy-based feature selection
- Species-specific models available

### 3. scMayoMap (Atlas-Based)
- Mayo Clinic tissue atlases
- Clinically validated cell types
- Tissue-specific annotations

### 4. Celltypist (Deep Learning)
- Python-based machine learning
- Hierarchical predictions (L1/L2/L3)
- Excellent for immune cells
- 40+ pre-trained models

## Common Use Cases

### Immune System Annotation
```r
tissue <- "Immune system"
# Identifies: T cells, B cells, NK cells, Monocytes, etc.
```

### Brain Tissue Annotation
```r
tissue <- "Brain"
# Identifies: Neurons, Astrocytes, Oligodendrocytes, etc.
```

### Multi-Method Validation
```r
# Run all three methods
# Compare results
# Create consensus annotation
```

## Available Tissues

**ScType**: Immune system, Brain, Lung, Liver, Kidney, Eye, Pancreas, Heart, Intestine, Muscle, Placenta, Spleen, Stomach, Thymus

**Scibet**: Major cell types, Human Cell Atlas, Mouse brain, custom references

**scMayoMap**: Blood, Bone marrow, Spleen, Brain, Liver, Lung, Pancreas

## Citation

```
Ianevski et al. (2022). Fully-automated cell-type identification. Nature Communications.
Li et al. (2020). Reference-free cell type deconvolution. Nature Communications.
```

## Documentation

- Main skill: [SKILL.md](SKILL.md)
- Examples: [references/annotation_examples.md](references/annotation_examples.md)
