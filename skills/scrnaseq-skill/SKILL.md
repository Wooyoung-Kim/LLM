---
name: scrnaseq-skill
description: scRNA-seq analysis with Scanpy/AnnData, including preprocessing, QC, normalization, dimensionality reduction, clustering, marker discovery, and visualization. Use when Codex needs to write Scanpy pipelines, explain Scanpy workflows, or troubleshoot Scanpy/AnnData tasks.
---

# Scanpy Skill

## Quick workflow

1. Load data into `AnnData` (e.g., `sc.read_10x_mtx`, `sc.read_h5ad`, or `anndata.AnnData`).
2. Run QC and filter: compute `sc.pp.calculate_qc_metrics`, filter cells/genes, and annotate QC thresholds.
3. Normalize and log-transform: `sc.pp.normalize_total`, `sc.pp.log1p`, optionally `sc.pp.highly_variable_genes`.
4. Scale, reduce, and neighbors: `sc.pp.scale`, `sc.tl.pca`, `sc.pp.neighbors`.
5. Cluster and embed: `sc.tl.leiden` or `sc.tl.louvain`, then `sc.tl.umap` or `sc.tl.tsne`.
6. Find markers and visualize: `sc.tl.rank_genes_groups` with plots like `sc.pl.umap`, `sc.pl.rank_genes_groups`.

## Guidance

- Prefer reproducible settings: set `random_state` where relevant and store key parameters in `adata.uns`.
- Keep raw counts in `adata.raw` before scaling, and use `adata.layers` for normalized or corrected matrices.
- Use clear metadata columns in `adata.obs` for batch, condition, and sample IDs.
- When integrating batches, use Scanpy tools or known workflows (e.g., BBKNN, Harmony via scanpy.external).
