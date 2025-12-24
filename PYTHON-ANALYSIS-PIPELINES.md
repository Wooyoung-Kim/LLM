# Python Skills Analysis Pipelines

Complete end-to-end analysis pipelines using Python skills for single-cell RNA-seq, genomics, and systems biology research.

---

## 📋 Table of Contents

- [Quick Start](#quick-start)
- [Complete Pipelines](#complete-pipelines)
  - [Single-Cell RNA-seq: From Raw Data to Publication](#1-single-cell-rna-seq-from-raw-data-to-publication)
  - [Cell Type Annotation with Deep Learning](#2-cell-type-annotation-with-deep-learning)
  - [Differential Expression Analysis](#3-differential-expression-analysis)
  - [Integration with Public Datasets](#4-integration-with-public-datasets)
  - [Trajectory Analysis](#5-trajectory-analysis)
  - [Multi-Modal Integration](#6-multi-modal-integration)
- [Skill Combinations](#skill-combinations)
- [Tips & Best Practices](#tips--best-practices)

---

## 🚀 Quick Start

각 Python 스킬을 사용하여 Claude에게 요청할 수 있는 분석입니다:

```python
# 예시: scRNA-seq 전체 분석 요청
"Use the Scanpy skill to analyze my 10X Genomics h5 file. 
Perform QC, normalization, clustering using Leiden algorithm. 
Then use Celltypist for automated cell type annotation,
run differential expression with PyDESeq2 in pseudobulk mode,
and create publication-quality figures with matplotlib and seaborn."
```

---

## 🔬 Complete Pipelines

### 1. Robust Single-Cell Analysis (Sample-by-Sample)

**목표**: 개별 샘플 QC 및 Celltypist Pre-annotation 후 SCANVI로 통합

**사용 스킬**: Scanpy → Celltypist → SCANVI → PyDESeq2 → Scientific-Visualization

**Claude 프롬프트**:
```
Use available Python skills you have access to whenever possible.

Step 1: Individual Sample Processing (Loop)
- For each sample file:
  1. Load data (sc.read_10x_h5)
  2. Create AnnData, add `sample_id`
  3. **Quality Control**:
     - Calculate metrics (n_genes, percent_mito, etc.)
     - Visualize with scientific-visualization (Violin)
     - Filter (standard cutoffs)
  4. **Pre-annotation (Celltypist)**:
     - Normalize to 10k, log1p
     - Run `celltypist.annotate()` (per-cell prediction)
     - Store predictions in `adata.obs['celltypist_pre']`
  5. Save individual processed files: `adata.write('results/01_individual/sample_id.h5ad')`

Step 2: Integration with SCANVI
- Concatenate all individual AnnDays
- Setup scVI:
  * `batch_key="sample_id"`
  * `layer="counts"` (use raw counts)
- **Train scVI** (Unsupervised integration)
- **Initialize SCANVI from scVI**:
  * Use `celltypist_pre` as labels_key (Semi-supervised!)
  * Treat low-confidence Celltypist predictions as "Unknown"
- Train SCANVI to refine integration and annotation
- Get latent representation (`X_scanvi`)
- Save: `adata.write('results/02_integrated.h5ad')`

Step 3: Clustering and Final Annotation
- Compute Neighbors and UMAP on `X_scanvi`
- Run Leiden clustering
- **Consensus Annotation**:
  * Compare SCANVI predictions vs Leiden clusters
  * Assign final cell types
- Visualize:
  * UMAP colored by Sample (Batch mixing)
  * UMAP colored by Final Cell Type
- Save: `adata.write('results/03_annotated.h5ad')`

Step 4: Marker Discovery and Downstream
- Rank genes groups (Wilcoxon)
- Dot plots and Heatmaps (Scientific-Visualization)
- **Save**: `marker_df.to_csv('results/04_markers.csv')`

Step 5: Pseudobulk DE Analysis (PyDESeq2)
- Aggregate by Sample + Cell Type
- Run PyDESeq2
- Filters and Volcano plots
- Save: `results.to_csv('results/05_degs.csv')`

Step 6: Publication Figures
- Use scientific-visualization skill for all plots
- High-res export (PDF/PNG)
```

**예상 결과**:
- 개별 샘플의 정밀한 품질 관리
- Celltypist와 SCANVI의 시너지 (Semi-supervised integration)
- 배치 효과가 완벽히 제거된 통합 데이터
- 검증된 최종 세포 타입 주석

Step 8: Pathway enrichment (optional)
- Use decoupler for PROgeny pathway activities
- Or export gene lists for enrichR/GSEA

Step 9: Publication-quality figures
- Use scientific-visualization skill to create:
  * Configure publication-ready matplotlib settings
  * High-resolution UMAP plots (300 DPI, PDF/PNG)
  * Styled dot plots with journal themes
  * Volcano plots for DE results (if applicable)
  * Heatmaps with proper color schemes and clustering
  * Multi-panel figures with consistent formatting
- Apply journal-specific requirements:
  * Nature: 89mm (single column) or 183mm (double column)
  * Science: 56mm or 120mm width
  * Cell: Custom dimensions with specific fonts
- Export formats: PDF (vector), PNG (raster), TIFF (high-res)
- Save with proper DPI (300-600) and embedded fonts

Step 10: Save results
- Save annotated AnnData as h5ad
- Export figures to figures/ directory
- Save markers and DEGs as CSV
- Create summary statistics report
```

**예상 결과**:
- 고품질 클러스터링 및 annotation
- 세포 타입별 마커 유전자
- Robust한 차별 발현 분석
- 논문 제출용 그림들

---

### 2. Cell Type Annotation with Deep Learning

**목표**: 다양한 방법으로 세포 타입 자동 주석 및 검증

**사용 스킬**: Scanpy → Celltypist → SCANVI → scvi-tools

**Claude 프롬프트**:
```
Use available Python skills you have access to whenever possible.

Step 1: Prepare data
- Load preprocessed AnnData object
- Ensure normalized and log-transformed
- Have clustering results available

Step 2: Celltypist annotation (Method 1)
- Download multiple models:
  * 'Immune_All_High.pkl' for high-resolution
  * 'Immune_All_Low.pkl' for broad categories
  * 'AIFI_L1.pkl', 'AIFI_L2.pkl', 'AIFI_L3.pkl' for hierarchical
- Run predictions with each model (per-cell)
- Compare results across models
- (Optional) Use majority_voting for cluster-level consensus
- **Save**: `adata.write('results/annotation/02_celltypist_raw.h5ad')`

Step 3: Semi-supervised annotation with scVI + SCANVI (Method 2)
- Concatenate Reference and Query data
- Train scVI model first (Unsupervised):
  * Learns the data manifold and corrects batch effects
  * `scvi.model.SCVI(adata)`
- Initialize SCANVI from scVI model:
  * `scvi.model.SCANVI.from_scvi_model(scvi_model, labels_key="cell_type")`
  * Treats query cells as "unlabeled"
- Train SCANVI (Semi-supervised)
- Predict labels for query cells
- Add to adata.obs['scanvi_prediction']
- **Save**:
  * `scanvi_model.save('models/scanvi_model/')`
  * `adata.write('results/annotation/03_scanvi_pred.h5ad')`

Step 4: Scanpy marker-based (Method 3)
- Use known marker genes for validation
- Create custom annotation based on:
  * Automated predictions
  * Marker gene expression
  * Literature knowledge
- Score cell types using sc.tl.score_genes()
- **Save**: `adata.write('results/annotation/04_marker_scored.h5ad')`

Step 5: Compare and validate
- Create comparison UMAP plots
- Calculate agreement between methods
- Generate confusion matrices
- Identify high-confidence vs ambiguous cells

Step 6: Create consensus annotation
- Combine predictions from multiple methods
- Assign final cell types based on:
  * Agreement between methods (>66%)
  * Confidence scores
  * Manual curation for ambiguous cases
- Add confidence levels to metadata

Step 7: Visualization
- Use scientific-visualization skill:
  * Multi-panel UMAPs showing different annotations
  * Publication-quality dot plots of marker genes
  * Sankey diagrams showing label transfer
  * Confidence score distributions with proper styling
  * Apply consistent color schemes (colorblind-safe)
  * Set journal-specific figure dimensions
- Export annotated AnnData

Step 8: Quality control report
- Generate annotation statistics
- Cell counts per type
- Confidence metrics
- Marker gene enrichment validation
```

**예상 결과**:
- 높은 신뢰도의 세포 타입 annotation
- 여러 방법 간 검증
- Confidence scores
- 상세한 품질 관리 보고서

---

### 3. Differential Expression Analysis

**목표**: Pseudobulk 방식의 robust한 차별 발현 분석

**사용 스킬**: Scanpy → PyDESeq2 → Matplotlib

**Claude 프롬프트**:
```
Use available Python skills you have access to whenever possible.

Step 1: Prepare annotated data
- Load AnnData with cell type annotations
- Verify metadata: condition, replicate, cell_type
- Check for at least 3 biological replicates per condition

Step 2: Pseudobulk aggregation
- Group cells by: condition + replicate + cell_type
- Sum raw counts for each pseudobulk sample
- Create pseudobulk count matrix
- Generate corresponding metadata DataFrame
- **Save**: `pseudobulk_adata.write('results/de/02_pseudobulk.h5ad')`

Step 3: Run PyDESeq2 per cell type
- For each cell type:
  * Filter pseudobulk matrix to cell type
  * Create DeseqDataSet
  * Set design formula: ~ condition
  * Run DESeq2 pipeline
  * Apply LFC shrinkage (lfcShrink)
  * Get results for contrast: treated vs control
- **Save**: `pickle.dump(dds_results, open('results/de/03_dds_objects.pkl', 'wb'))`

Step 4: Filter and annotate results
- For each cell type:
  * Filter significant genes (padj < 0.05, |log2FC| > 1)
  * Annotate with gene symbols
  * Add cell type information
  * Combine all results into single DataFrame
- **Save**: `all_degs.to_csv('results/de/04_all_significant_genes.csv')`

Step 5: Cross-cell-type analysis
- Identify shared DEGs across cell types
- Find cell-type-specific responses
- Create Venn diagrams or UpSet plots
- Calculate enrichment of shared genes

Step 6: Visualization
- Use scientific-visualization skill:
  * Create volcano plots for each cell type with proper styling
  * Generate publication-ready MA plots
  * Heatmap of top DEGs with journal themes
  * Bar plots of DEG counts with consistent colors
  * Multi-panel composite figures
  * Export high-resolution figures (PDF 300 DPI, PNG 600 DPI)

Step 7: Functional interpretation
- Prepare gene lists for pathway analysis
- Use decoupler for pathway activities (PROgeny, DoRothEA)
- Or export for external tools (GSEA, Enrichr)

Step 8: Export results
- Save results as CSV/TSV
- Create Excel file with multiple sheets (one per cell type)
- Generate summary statistics
- Export session info and parameters
```

**예상 결과**:
- 세포 타입별 robust한 DEG 목록
- 통계적으로 검증된 결과
- Publication-ready 그림들
- 추가 분석용 gene lists

---

### 4. Integration with Public Datasets

**목표**: CELLxGENE Census와 공개 데이터 통합 분석

**사용 스킬**: Scanpy → Cellxgene-Census → SCANVI → Celltypist

**Claude 프롬프트**:
```
Use available Python skills you have access to whenever possible.

Step 1: Load personal data
- Load your AnnData object
- Perform standard preprocessing
- Basic clustering and QC

Step 2: Query Cellxgene Census
- Use cellxgene_census to access data
- Query relevant datasets:
  * Filter by tissue, disease, organism
  * Select matching cell types
  * Download reference data
- Create reference AnnData object

Step 3: Data harmonization
- Align gene names between datasets
- Identify common highly variable genes
- Match metadata schemas
- Standardize normalization

Step 4: Integration with scvi-tools
- Concatenate datasets (adata_concat)
- Train scVI model on combined data:
  * Account for batch effects
  * Use categorical covariates (dataset, batch)
  * Train for sufficient epochs
- Get latent representation
- Compute UMAP on integrated space
- **Save**:
  * `scvi_model.save('models/integration_model/')`
  * `adata_concat.write('results/integration/04_integrated.h5ad')`

Step 5: Label transfer with SCANVI
- Initialize SCANVI from the trained scVI model:
  * `scvi.model.SCANVI.from_scvi_model(scvi_model)`
  * Use "Unknown" label for query cells
- Train SCANVI to refine latent space with labels
- Predict cell types for query data
- Get prediction probabilities (confidence)
- Validate transferred labels against clusters
- **Save**: `adata_concat.write('results/integration/05_label_transferred.h5ad')`

Step 6: Additional annotation with Celltypist
- Run Celltypist on query data
- Compare with transferred labels
- Create consensus annotation

Step 7: Batch effect assessment
- Calculate integration metrics:
  * Mixing scores
  * kBET (k-nearest neighbor batch effect test)
  * Local Inverse Simpson's Index (LISI)
- Visualize batch mixing in UMAP

Step 8: Comparative analysis
- Compare cell type compositions
- Identify dataset-specific populations
- Find conserved marker genes
- Detect differential states

Step 9: Visualization and export
- Create integrated UMAPs
- Split views by dataset/condition
- Dot plots of conserved markers
- Export integrated AnnData
```

**예상 결과**:
- 공개 데이터와 통합된 분석
- 검증된 cell type annotations
- 배치 효과 제거
- 더 큰 맥락에서의 해석

---

### 5. Trajectory Analysis

**목표**: 세포 분화 궤적 및 pseudotime 분석

**사용 스킬**: Scanpy → scvi-tools → Matplotlib

**Claude 프롬프트**:
```
Use available Python skills you have access to whenever possible.

Step 1: Prepare data for trajectory
- Load preprocessed AnnData
- Subset to cell types of interest (e.g., differentiation lineage)
- Ensure smooth gene expression manifold

Step 2: Trajectory inference with Scanpy
- Use sc.tl.paga() for trajectory structure:
  * Build PAGA graph
  * Visualize connectivity
  * Identify root cell type
- Compute diffusion pseudotime (DPT):
  * Set root cell
  * Calculate dpt values
  * Add to adata.obs
- **Save**: `adata.write('results/trajectory/02_paga_dpt.h5ad')`

Step 3: Alternative: scVelo for RNA velocity (if available)
- Prepare spliced/unspliced counts
- Estimate RNA velocity
- Project trajectories
- Infer directionality
- **Save**: `adata.write('results/trajectory/03_velocity.h5ad')`

Step 4: Identify trajectory-associated genes
- Correlate gene expression with pseudotime
- Use Spearman correlation
- Identify early, middle, late genes
- Find switch-like vs gradual changes
- **Save**: `traj_genes.to_csv('results/trajectory/04_trajectory_genes.csv')`

Step 5: Gene expression dynamics
- Plot gene expression along pseudotime
- Create heatmap of dynamic genes
- Cluster genes by temporal patterns
- Identify key regulators at transitions

Step 6: Functional characterization
- Run pathway analysis at different pseudotime windows
- Use decoupler for TF activities
- Identify stage-specific programs

Step 7: Visualization
- Plot trajectories on UMAP/Force-directed graph
- Color by pseudotime
- Expression dynamics plots
- Heatmaps of trajectory genes
- Export publication figures

Step 8: Export results
- Save trajectory-annotated AnnData
- Export pseudotime values
- Save dynamic gene lists
- Create summary report
```

**예상 결과**:
- 세포 분화 경로 맵
- Pseudotime 값
- 궤적 관련 유전자
- 발달 단계별 특성

---

### 6. Multi-Modal Integration

**목표**: scRNA-seq + ATAC-seq 또는 protein 데이터 통합

**사용 스킬**: Scanpy → AnnData → scvi-tools

**Claude 프롬프트**:
```
Use available Python skills you have access to whenever possible.

Step 1: Load multi-modal data
- Load RNA AnnData object
- Load ATAC or protein (ADT) data
- Verify matching cells between modalities

Step 2: Create MuData object (if applicable)
- Use mudata package for multi-modal
- Or store as separate layers in AnnData
- Align cell barcodes

Step 3: Individual modality processing
- RNA:
  * Standard Scanpy workflow
  * Normalization, HVG selection
- ATAC (if available):
  * Peak calling and filtering
  * TF-IDF transformation
  * LSI dimensionality reduction
- Protein (if ADT):
  * CLR normalization
  * Feature selection
- **Save**: `mdata.write('results/multimodal/03_preprocessed.h5mu')`

Step 4: Integration with scvi-tools
- Use totalVI for RNA+protein:
  * Joint latent representation
  * Model batch effects
  * Impute protein values
- Or use MultiVI for RNA+ATAC:
  * Learn shared latent space
  * Account for modality-specific variation
- **Save**:
  * `model.save('models/multimodal_model/')`
  * `mdata.write('results/multimodal/04_integrated.h5mu')`

Step 5: Joint clustering and visualization
- Compute UMAP on integrated latent space
- Leiden clustering
- Annotate cell types using:
  * RNA markers
  * Protein markers
  * Accessible chromatin patterns

Step 6: Multi-modal analysis
- Correlate RNA and protein expression
- Link ATAC peaks to genes
- Identify modality-specific patterns
- Regulatory element analysis

Step 7: Visualization
- Create multi-panel UMAPs (one per modality)
- Feature plots showing RNA + protein
- Correlation plots
- Integration quality metrics

Step 8: Export
- Save integrated MuData/AnnData
- Export modality-specific results
- Create comprehensive report
```

**예상 결과**:
- 통합된 multi-modal 분석
- RNA-protein/ATAC 관계
- 더 정확한 cell typing
- 조절 메커니즘 insight

---

## 🎯 Skill Combinations

### Core Workflows

#### Workflow 1: Basic scRNA-seq
```
Scanpy → Scientific-Visualization
```
- QC, normalization, clustering
- Publication-quality 시각화

#### Workflow 2: Annotated Analysis
```
Scanpy → Celltypist/SCANVI → PyDESeq2 → Scientific-Visualization
```
- 자동 annotation (Celltypist 또는 SCANVI)
- Robust DE analysis
- Publication figures

#### Workflow 3: Deep Integration
```
Scanpy → Cellxgene-Census → SCANVI → Celltypist → Scientific-Visualization
```
- 공개 데이터 통합
- Reference-based annotation (SCANVI)
- Validation with Celltypist
- 배치 효과 보정

#### Workflow 4: Complete Pipeline
```
Scanpy → scvi-tools → Celltypist → PyDESeq2 → decoupler → Scientific-Visualization
```
- 데이터 통합 및 배치 보정
- 자동 세포 타입 주석
- 전체 통합 분석 및 시각화

---

## 💡 Tips & Best Practices

### 1. 스킬 선택 가이드

**Scanpy를 사용하세요**:
- scRNA-seq 기본 분석
- QC, 정규화, 클러스터링
- 빠르고 메모리 효율적
- AnnData 생태계의 중심

**Celltypist를 사용하세요**:
- 자동 세포 타입 annotation
- 40+ pre-trained models
- 특히 immune cells에 강력
- 계층적 annotation (L1/L2/L3)

**SCANVI를 사용하세요**:
- Reference-based annotation
- Batch correction + label transfer
- 확률적 예측 및 신뢰도 점수
- scVI → SCANVI 2단계 학습
- Semi-supervised learning

**PyDESeq2를 사용하세요**:
- Pseudobulk 차별 발현 분석
- R의 DESeq2와 동일한 결과
- 생물학적 반복 필수 (최소 3개)
- FDR 제어 우수

**scvi-tools를 사용하세요**:
- 배치 효과 보정
- 데이터 통합
- Label transfer
- Probabilistic modeling

**Cellxgene-Census를 사용하세요**:
- 공개 데이터 접근
- Reference mapping
- Meta-analysis
- 대규모 atlas 활용

**Scientific-Visualization을 사용하세요**:
- Publication-quality 그림 생성
- 저널별 요구사항 자동 적용 (Nature, Science, Cell)
- Colorblind-safe palettes
- 일관된 스타일 및 테마
- 고해상도 export (PDF, PNG, TIFF)
- 멀티패널 레이아웃
- 폰트 임베딩 및 벡터 포맷

### 2. 데이터 구조

**AnnData 구조 이해**:
```python
adata.X          # Normalized expression matrix
adata.raw        # Raw counts (backup)
adata.obs        # Cell metadata
adata.var        # Gene metadata
adata.obsm       # Multi-dimensional annotations (PCA, UMAP)
adata.uns        # Unstructured annotations
adata.layers     # Alternative matrices (raw, scaled, etc.)
```

**필수 메타데이터**:
```python
adata.obs['condition']  # experimental condition
adata.obs['replicate']  # biological replicate
adata.obs['cell_type']  # cell type annotation
adata.obs['batch']      # batch information
```

### 3. 메모리 관리

**대용량 데이터 처리**:
```python
# 메모리 효율적 읽기
adata = sc.read_h5ad(filename, backed='r')

# Subset 후 메모리에 로드
adata_subset = adata[cells, :].to_memory()

# Raw 저장으로 메모리 절약
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
```

### 4. 재현성

**항상 포함하세요**:
```python
# Random seed 설정
import random
import numpy as np
random.seed(42)
np.random.seed(42)

# Session info
sc.logging.print_versions()

# 파라미터 저장
adata.uns['analysis_params'] = {
    'n_top_genes': 2000,
    'n_pcs': 40,
    'resolution': 0.5
}
```

### 5. 품질 관리

**QC 메트릭**:
```python
# 표준 QC
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# 추가 QC
adata.obs['percent_ribo'] = (
    adata[:, adata.var_names.str.startswith('RPS')].X.sum(1).A1 / 
    adata.obs['total_counts']
) * 100

# 시각화
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mt'], 
             multi_panel=True)
```

---

## 📚 Example Data Sources

### Public Datasets
- **CELLxGENE**: https://cellxgene.cziscience.com/
- **10X Genomics**: https://www.10xgenomics.com/resources/datasets
- **Human Cell Atlas**: https://data.humancellatlas.org/
- **GEO**: https://www.ncbi.nlm.nih.gov/geo/

### Tutorial Datasets
```python
# Scanpy built-in
import scanpy as sc
adata = sc.datasets.pbmc3k()
adata = sc.datasets.pbmc68k_reduced()

# 10X datasets
sc.datasets.paul15()  # hematopoiesis
sc.datasets.moignard15()  # blood development
```

---

## 🔗 Related Resources

### Python Package Documentation
- [Scanpy](https://scanpy.readthedocs.io/)
- [Celltypist](https://www.celltypist.org/)
- [PyDESeq2](https://pydeseq2.readthedocs.io/)
- [scvi-tools](https://scvi-tools.org/)
- [AnnData](https://anndata.readthedocs.io/)
- [decoupler-py](https://decoupler-py.readthedocs.io/)

### Tutorials
- [Scanpy Tutorials](https://scanpy-tutorials.readthedocs.io/)
- [Best Practices in scRNA-seq](https://www.sc-best-practices.org/)
- [scvi-tools Tutorials](https://docs.scvi-tools.org/en/stable/tutorials/)

---

## 📝 Citation

```bibtex
@article{wolf2018scanpy,
  title={SCANPY: large-scale single-cell gene expression data analysis},
  author={Wolf, F Alexander and others},
  journal={Genome biology},
  year={2018}
}

@article{dominguez2022celltypist,
  title={Cross-tissue immune cell analysis reveals tissue-specific features in humans},
  author={Domínguez Conde, C and others},
  journal={Science},
  year={2022}
}

@article{Lopez2018,
  title={Deep generative modeling for single-cell transcriptomics},
  author={Lopez, Romain and others},
  journal={Nature methods},
  year={2018}
}

@article{schaarschmidt2022pydeseq2,
  title={A Python implementation of DESeq2},
  author={Schaarschmidt, S and others},
  year={2022}
}
```

---

**Created**: 2025-12-24  
**Version**: 1.0  
**License**: MIT  

모든 스킬과 파이프라인은 `/home/kwy7605/LLM/skills/` 디렉토리에서 확인할 수 있습니다.
