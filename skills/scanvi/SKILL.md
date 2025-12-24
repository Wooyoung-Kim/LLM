---
name: scanvi
description: "Semi-supervised cell type annotation using SCANVI: scVI → SCANVI workflow for reference-based annotation and label transfer in single-cell RNA-seq"
---

# SCANVI (Single-Cell ANnotation using Variational Inference)

## Overview

SCANVI is a semi-supervised deep generative model for single-cell RNA-seq data that performs automated cell type annotation by leveraging reference datasets with known labels. It builds on scVI (Variational Inference for scRNA-seq) and is particularly powerful for:

- **Reference-based annotation**: Transfer cell type labels from reference to query data
- **Semi-supervised learning**: Learn from both labeled (reference) and unlabeled (query) data
- **Batch correction**: Integrate data across batches while preserving biological variation
- **Uncertainty quantification**: Provide confidence scores for predictions

## When to Use This Skill

Use SCANVI when you need to:
- Annotate query cells using a well-annotated reference dataset
- Transfer labels across batches, technologies, or studies
- Integrate reference and query data while correcting batch effects
- Get probabilistic predictions with confidence scores
- Handle partially labeled datasets (semi-supervised scenario)

**Key advantage**: Unlike purely unsupervised methods, SCANVI leverages known cell type information to improve both integration and annotation quality.

## Quick Start Guide

### Basic Workflow: Train scVI → Convert to SCANVI

```python
import scanpy as sc
import scvi

# Load data (reference + query concatenated)
adata = sc.read_h5ad("combined_data.h5ad")

# Setup AnnData for scVI
scvi.model.SCVI.setup_anndata(
    adata,
    batch_key="batch",
    labels_key="cell_type"  # Has labels for reference, "Unknown" for query
)

# Step 1: Train scVI (unsupervised)
scvi_model = scvi.model.SCVI(adata)
scvi_model.train()

# Step 2: Initialize SCANVI from scVI model
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    labels_key="cell_type",
    unlabeled_category="Unknown"
)

# Step 3: Train SCANVI (semi-supervised)
scanvi_model.train(max_epochs=20)

# Get predictions for query cells
adata.obs["scanvi_prediction"] = scanvi_model.predict()
adata.obs["scanvi_confidence"] = scanvi_model.predict(soft=True).max(axis=1)

# Get latent representation
adata.obsm["X_scanvi"] = scanvi_model.get_latent_representation()

# Compute UMAP on SCANVI latent space
sc.pp.neighbors(adata, use_rep="X_scanvi")
sc.tl.umap(adata)

# Visualize
sc.pl.umap(adata, color=["scanvi_prediction", "scanvi_confidence"])
```

## Core Workflow Steps

### 1. Data Preparation

#### Combine Reference and Query Data

```python
import scanpy as sc
import pandas as pd

# Load reference (with cell type labels)
adata_ref = sc.read_h5ad("reference.h5ad")
adata_ref.obs["dataset"] = "reference"
adata_ref.obs["batch"] = "ref_batch1"

# Load query (no labels)
adata_query = sc.read_h5ad("query.h5ad")
adata_query.obs["dataset"] = "query"
adata_query.obs["cell_type"] = "Unknown"  # Mark as unlabeled
adata_query.obs["batch"] = "query_batch1"

# Concatenate
adata = sc.concat([adata_ref, adata_query], join="inner")

# Basic preprocessing
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

#### Critical: Store Raw Counts

SCANVI requires **raw counts**, not normalized data:

```python
# Method 1: Keep raw counts in a layer
adata.layers["counts"] = adata.raw.X.copy()

# Method 2: Use .raw attribute
adata.raw = adata.copy()
```

### 2. Setup AnnData for scVI/SCANVI

```python
import scvi

# Setup specifies which columns to use
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",  # Use raw counts
    batch_key="batch",  # Batch correction
    labels_key="cell_type",  # For SCANVI
    categorical_covariate_keys=["dataset"],  # Additional covariates
)

# Verify setup
adata.uns["_scvi"]["summary_stats"]
```

### 3. Train scVI Model (Unsupervised Phase)

```python
# Initialize scVI
scvi_model = scvi.model.SCVI(
    adata,
    n_latent=30,  # Latent dimension
    n_layers=2,  # Neural network layers
    gene_likelihood="nb"  # Negative binomial (default for UMI data)
)

# Train
scvi_model.train(
    max_epochs=400,
    early_stopping=True,
    early_stopping_patience=20,
    plan_kwargs={'lr': 1e-3}
)

# Check training history
scvi_model.history["elbo_train"].plot()

# Save model
scvi_model.save("scvi_model/", overwrite=True)
```

### 4. Convert to SCANVI (Semi-Supervised Phase)

#### Initialize from scVI Model

```python
# Convert scVI → SCANVI
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    labels_key="cell_type",
    unlabeled_category="Unknown"
)

# Train SCANVI (much faster since initialized from scVI)
scanvi_model.train(
    max_epochs=20,
    early_stopping=True,
    plan_kwargs={'lr': 5e-4}  # Lower learning rate
)

# Save SCANVI model
scanvi_model.save("scanvi_model/", overwrite=True)
```

#### Alternative: Train SCANVI from Scratch (Not Recommended)

```python
# Directly train SCANVI (slower and less stable)
scanvi_model = scvi.model.SCANVI(
    adata,
    n_latent=30,
    unlabeled_category="Unknown"
)
scanvi_model.train(max_epochs=400)
```

### 5. Predict Cell Types

#### Get Predictions

```python
# Hard predictions (most likely cell type)
predictions = scanvi_model.predict()
adata.obs["scanvi_prediction"] = predictions

# Soft predictions (probability matrix)
soft_predictions = scanvi_model.predict(soft=True)  # Shape: (n_cells, n_cell_types)

# Confidence score (max probability)
adata.obs["scanvi_confidence"] = soft_predictions.max(axis=1)

# Uncertainty (entropy)
from scipy.stats import entropy
adata.obs["scanvi_entropy"] = entropy(soft_predictions.T)
```

#### Filter by Confidence

```python
# Only keep high-confidence predictions
threshold = 0.7
adata.obs["scanvi_filtered"] = adata.obs["scanvi_prediction"].copy()
adata.obs.loc[adata.obs["scanvi_confidence"] < threshold, "scanvi_filtered"] = "Low_confidence"
```

### 6. Get Latent Representation

```python
# Extract SCANVI latent space
latent = scanvi_model.get_latent_representation()
adata.obsm["X_scanvi"] = latent

# Compute neighbors and UMAP on SCANVI space
sc.pp.neighbors(adata, use_rep="X_scanvi", n_neighbors=15)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="scanvi_leiden")
```

### 7. Visualization

```python
import matplotlib.pyplot as plt

# UMAP plots
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

sc.pl.umap(adata, color="scanvi_prediction", ax=axes[0], show=False)
axes[0].set_title("SCANVI Predictions")

sc.pl.umap(adata, color="scanvi_confidence", ax=axes[1], show=False, cmap="viridis")
axes[1].set_title("Prediction Confidence")

sc.pl.umap(adata, color="dataset", ax=axes[2], show=False)
axes[2].set_title("Dataset (Reference vs Query)")

plt.tight_layout()
plt.show()

# Compare with true labels (for reference cells)
ref_cells = adata.obs["dataset"] == "reference"
sc.pl.umap(adata[ref_cells], color=["cell_type", "scanvi_prediction"])
```

## Advanced Features

### 1. Differential Expression with SCANVI

```python
# DE between cell types in the latent space
de_df = scanvi_model.differential_expression(
    groupby="scanvi_prediction",
    group1="T_cells",
    group2="B_cells"
)

# Get top markers
top_markers = de_df.sort_values("lfc_mean", ascending=False).head(20)
```

### 2. Query-to-Reference Surgery

For very large reference datasets, use "surgery" to avoid retraining on reference:

```python
# Train scVI + SCANVI on reference only
scvi_ref = scvi.model.SCVI(adata_ref)
scvi_ref.train()

scanvi_ref = scvi.model.SCANVI.from_scvi_model(scvi_ref, unlabeled_category="Unknown")
scanvi_ref.train()

# Perform surgery on query
scanvi_query = scvi.model.SCANVI.load_query_data(
    adata_query,
    scanvi_ref
)
scanvi_query.train(max_epochs=100, plan_kwargs={"weight_decay": 0.0})

# Predict
adata_query.obs["prediction"] = scanvi_query.predict()
```

### 3. Multiple References

```python
# Combine multiple reference datasets
adata_ref1.obs["ref_source"] = "ref1"
adata_ref2.obs["ref_source"] = "ref2"

adata_combined = sc.concat([adata_ref1, adata_ref2, adata_query])

# Use ref_source as additional covariate
scvi.model.SCVI.setup_anndata(
    adata_combined,
    batch_key="batch",
    labels_key="cell_type",
    categorical_covariate_keys=["ref_source"]
)
```

### 4. Hierarchical Annotation

For multi-level annotations (broad → fine):

```python
# Level 1: Broad categories
scanvi_l1 = scvi.model.SCANVI(adata, labels_key="cell_type_l1", unlabeled_category="Unknown")
scanvi_l1.train()
adata.obs["prediction_l1"] = scanvi_l1.predict()

# Level 2: Fine-grained (subset by L1 prediction)
for broad_type in adata.obs["prediction_l1"].unique():
    adata_subset = adata[adata.obs["prediction_l1"] == broad_type].copy()
    
    scanvi_l2 = scvi.model.SCANVI(adata_subset, labels_key="cell_type_l2", unlabeled_category="Unknown")
    scanvi_l2.train()
    
    adata.obs.loc[adata.obs["prediction_l1"] == broad_type, "prediction_l2"] = scanvi_l2.predict()
```

## Best Practices

### Data Preparation

**DO**:
- Use raw counts (not normalized)
- Normalize reference and query the same way
- Include sufficient number of cells per cell type in reference (>50 recommended)
- Use high-quality, well-curated reference datasets

**DON'T**:
- Don't use heavily filtered genes (keep >2000)
- Don't use TPM or CPM (use counts)
- Don't mix UMI and non-UMI data without careful normalization

### Model Training

**DO**:
- Train scVI first, then initialize SCANVI (recommended)
- Use early stopping
- Check convergence (ELBO plot)
- Save models for reproducibility

**DON'T**:
- Don't skip scVI pre-training (much slower without it)
- Don't overtrain (watch for overfitting)
- Don't use too few epochs (<50 for scVI, <10 for SCANVI)

### Prediction Validation

**DO**:
- Check confidence scores
- Validate with known marker genes
- Compare predictions across multiple references
- Manually inspect low-confidence cells

**DON'T**:
- Don't blindly trust all predictions
- Don't ignore cells with low confidence
- Don't skip biological validation

### Batch Correction

**DO**:
- Always specify batch_key if batches exist
- Check batch mixing metrics (kBET, LISI)
- Visualize batches in latent space

**DON'T**:
- Don't over-correct (may lose biological variation)
- Don't ignore batch effects in reference itself

## Common Issues and Solutions

### Issue 1: Poor Predictions

**Symptoms**: Many "Unknown" or low-confidence predictions

**Solutions**:
```python
# Increase training epochs
scanvi_model.train(max_epochs=50)

# Lower learning rate
scanvi_model.train(plan_kwargs={'lr': 1e-4})

# Check if reference covers query cell types
ref_types = adata[adata.obs["dataset"] == "reference"]["cell_type"].unique()
print(f"Reference has {len(ref_types)} cell types")
```

### Issue 2: Batch Effects Not Corrected

**Symptoms**: Clear separation by batch in UMAP

**Solutions**:
```python
# Ensure batch_key is set
scvi.model.SCVI.setup_anndata(adata, batch_key="batch", ...)

# Increase latent dimension
scvi_model = scvi.model.SCVI(adata, n_latent=50)

# Check batch mixing
from scib import metrics
metrics.silhouette_batch(adata, batch_key="batch", label_key="cell_type", embed="X_scanvi")
```

### Issue 3: Memory Issues

**Solutions**:
```python
# Use smaller batch size
scvi_model.train(batch_size=64)

# Reduce latent dimension
scvi_model = scvi.model.SCVI(adata, n_latent=20)

# Use GPU if available
scvi_model.train(use_gpu=True)
```

### Issue 4: Convergence Issues

**Solutions**:
```python
# Reduce learning rate
scvi_model.train(plan_kwargs={'lr': 5e-4})

# Increase patience
scvi_model.train(early_stopping_patience=30)

# Check for NaNs in data
assert not adata.X.isnan().any()
```

## Comparison with Other Methods

| Method | Type | Strength | Limitation |
|--------|------|----------|------------|
| **SCANVI** | Probabilistic | Uncertainty, batch correction | Needs reference |
| **Celltypist** | ML classifier | Fast, pre-trained models | No batch correction |
| **scANVI surgery** | Transfer learning | Scales to huge references | More complex |
| **Azimuth** | Reference mapping | User-friendly, web-based | Limited references |
| **scArches** | Architecture surgery | Extremely scalable | Requires scVI/SCANVI knowledge |

## Integration Metrics

Evaluate integration quality:

```python
from scib.metrics import metrics

# Batch correction
batch_score = metrics.silhouette_batch(adata, batch_key="batch", label_key="cell_type", embed="X_scanvi")

# Bio conservation
bio_score = metrics.nmi(adata, cluster_key="scanvi_leiden", label_key="cell_type")

# Overall
metrics.metrics(
    adata,
    adata_int=adata,
    batch_key="batch",
    label_key="cell_type",
    embed="X_scanvi",
    cluster_key="scanvi_leiden",
    organism="human"
)
```

## Resources

### Documentation
- [SCANVI Paper](https://www.embopress.org/doi/full/10.15252/msb.20209620)
- [scvi-tools Documentation](https://docs.scvi-tools.org/)
- [scvi-tools Tutorials](https://docs.scvi-tools.org/en/stable/tutorials/)

### Example References
- [Tabula Sapiens](https://tabula-sapiens-portal.ds.czbiohub.org/)
- [Human Cell Atlas](https://data.humancellatlas.org/)
- [CELLxGENE Census](https://chanzuckerberg.github.io/cellxgene-census/)

## Citation

```bibtex
@article{xu2021scanvi,
  title={Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models},
  author={Xu, Chenling and others},
  journal={Molecular systems biology},
  year={2021}
}

@article{lopez2018scvi,
  title={Deep generative modeling for single-cell transcriptomics},
  author={Lopez, Romain and others},
  journal={Nature methods},
  year={2018}
}
```

Use SCANVI for principled, probabilistic cell type annotation with uncertainty quantification.
