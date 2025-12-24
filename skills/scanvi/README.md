# SCANVI - Semi-Supervised Cell Type Annotation

Automated cell type annotation using SCANVI (Single-Cell ANnotation using Variational Inference) for reference-based label transfer.

## Quick Start

```python
import scanpy as sc
import scvi

# Load data (reference + query)
adata = sc.read_h5ad("combined.h5ad")

# Setup
scvi.model.SCVI.setup_anndata(adata, batch_key="batch", labels_key="cell_type")

# Train scVI (unsupervised)
scvi_model = scvi.model.SCVI(adata)
scvi_model.train()

# Convert to SCANVI (semi-supervised)
scanvi_model = scvi.model.SCANVI.from_scvi_model(scvi_model, unlabeled_category="Unknown")
scanvi_model.train(max_epochs=20)

# Predict
adata.obs["prediction"] = scanvi_model.predict()
adata.obs["confidence"] = scanvi_model.predict(soft=True).max(axis=1)
```

## Key Features

- **scVI → SCANVI Workflow**: Recommended two-stage training
- **Semi-Supervised Learning**: Leverage both labeled and unlabeled data
- **Batch Correction**: Integrate across batches while preserving biology
- **Uncertainty Quantification**: Confidence scores for predictions
- **Reference-Based**: Transfer labels from curated references
- **Scalable**: Surgery mode for huge reference datasets

## When to Use

**Use SCANVI when**:
- You have a well-annotated reference dataset
- You need to annotate query data with uncertain identity
- Batch correction is important
- You want probabilistic predictions with confidence

**Use Celltypist when**:
- You want pre-trained models (no reference needed)
- Speed is critical
- You don't need batch correction

**Use scANVI surgery when**:
- Reference dataset is very large (>1M cells)
- You want to avoid retraining on reference

## Workflow

```
1. Combine reference (labeled) + query (unlabeled)
   ↓
2. Train scVI (unsupervised phase)
   - Learn manifold structure
   - Correct batch effects
   ↓
3. Initialize SCANVI from scVI
   - scvi.model.SCANVI.from_scvi_model()
   ↓
4. Train SCANVI (semi-supervised phase)
   - Refine latent space with labels
   - Much faster than from scratch
   ↓
5. Predict query cell types
   - Hard predictions
   - Soft predictions (probabilities)
   - Confidence scores
```

## Common Use Cases

### Case 1: Annotate New Dataset
```python
# Reference: PBMC with known cell types
# Query: New PBMC sample (unknown)

# Combine and mark query as "Unknown"
adata_query.obs["cell_type"] = "Unknown"
adata = sc.concat([adata_ref, adata_query])

# Run SCANVI workflow...
```

### Case 2: Cross-Study Integration
```python
# Multiple batches, some with labels
adata.obs["cell_type"].fillna("Unknown", inplace=True)

scvi.model.SCVI.setup_anndata(adata, batch_key="study", labels_key="cell_type")
# Continue with scVI → SCANVI...
```

### Case 3: Hierarchical Annotation
```python
# Level 1: Broad types (Immune, Epithelial, ...)
scanvi_l1.predict()

# Level 2: Fine types within each broad category
for cell_type in ["T_cells", "B_cells"]:
    subset = adata[adata.obs["pred_l1"] == cell_type]
    scanvi_l2 = train_scanvi(subset)
    # Predict subtypes...
```

## Outputs

- `adata.obs["scanvi_prediction"]`: Predicted cell type
- `adata.obs["scanvi_confidence"]`: Max probability (0-1)
- `adata.obsm["X_scanvi"]`: Latent representation
- Trained models saved to disk

## Requirements

```bash
pip install scvi-tools scanpy
```

## Documentation

- Main skill: [SKILL.md](SKILL.md)
- References: [references/](references/)
- Helper scripts: [scripts/](scripts/)

## Citation

```
Xu et al. (2021). Probabilistic harmonization and annotation of single-cell 
transcriptomics data with deep generative models. 
Molecular Systems Biology.
```

Advanced semi-supervised annotation with uncertainty quantification.
