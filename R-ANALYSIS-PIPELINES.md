# R Skills Analysis Pipelines

Complete end-to-end analysis pipelines using R skills for single-cell RNA-seq, genomics, and systems biology research.

---

## 📋 Table of Contents

- [Quick Start](#quick-start)
- [Complete Pipelines](#complete-pipelines)
  - [Single-Cell RNA-seq: From Raw Data to Publication](#1-single-cell-rna-seq-from-raw-data-to-publication)
  - [Cell-Cell Communication Analysis](#2-cell-cell-communication-analysis)
  - [Treatment Response Analysis](#3-treatment-response-analysis)
  - [Disease vs Healthy Comparison](#4-disease-vs-healthy-comparison)
  - [Time Course Analysis](#5-time-course-analysis)
  - [Multi-Sample Integration](#6-multi-sample-integration)
- [Skill Combinations](#skill-combinations)
- [Tips & Best Practices](#tips--best-practices)

---

## 🚀 Quick Start

각 R 스킬을 사용하여 Claude에게 요청할 수 있는 분석입니다:

```r
# 예시: scRNA-seq 전체 분석 요청
"Use the Seurat skill to analyze my 10X Genomics data. 
Perform QC, normalization, clustering, and find marker genes. 
Then use pathway-enrichment to run GSEA on the DEGs, 
and create publication-quality figures with scientific-visualization-r."
```

---

## 🔬 Complete Pipelines

### 1. Robust Single-Cell Analysis (Sample-by-Sample)

**목표**: 개별 샘플 QC 및 Pre-annotation 후 통합하여 정밀 분석 수행

**사용 스킬**: Seurat → Celltype-Annotation → Scientific-Visualization-R → Pathway-Enrichment

**Claude 프롬프트**:
```
Use available R skills you have access to whenever possible.

Step 1: Individual Sample Processing (Loop)
- For each sample in data directory:
  1. Load data (Read10X)
  2. Create Seurat object
  3. **Quality Control**:
     - Calculate percent.mt
     - Visualize with scientific-visualization-r (Violin plots)
     - Filter (nFeature 200-2500, percent.mt < 5%)
  4. **Pre-annotation (Automated)**:
     - Normalize & Scale
     - Run PCA & Clustering (low res 0.1)
     - Use **celltype-annotation** skill (ScType) for initial labels
     - Add 'pre_annotation' to metadata
  5. Save individual processed objects: `saveRDS(obj, "results/01_individual/sample_ID.rds")`

Step 2: Merge and Integration
- Merge all individual Seurat objects
- Perform standard processing:
  * NormalizeData
  * FindVariableFeatures
  * ScaleData
  * RunPCA
- Run **Harmony** integration (group.by.vars = "sample_id")
- Run UMAP and Clustering on integrated reduction
- Visualize:
  * UMAP colored by Sample (check mixing)
  * UMAP colored by Pre-annotation (check biological consistency)
- Save: `saveRDS(seurat_integrated, "results/02_integrated.rds")`

Step 3: Final Cell Type Annotation
- Use **celltype-annotation** skill to refine labels:
  * Compare Pre-annotation vs Cluster-based annotation
  * Run Scibet/scMayoMap on integrated clusters if needed
  * Create Consensus Annotation
- Manually curate final labels
- Visualize with scientific-visualization-r
- Save: `saveRDS(seurat_final, "results/03_annotated.rds")`

Step 4: Marker Discovery (Downstream)
- FindAllMarkers on final cell types
- Create DotPlots and Heatmaps (Scientific-Visualization-R)
- Save: `write.csv(markers, "results/04_markers.csv")`

Step 5: Pathway enrichment analysis
- For each cell type, get top DEGs
- Use pathway-enrichment skill (fgsea, clusterProfiler)
- Save: `saveRDS(enrichment, "results/05_pathway.rds")`

Step 6: Publication-ready figures
- Use scientific-visualization-r for all final plots
- Export as PDF (300 DPI)
```

**예상 결과**:
- 개별 샘플의 철저한 QC
- Integration 전 사전 주석 정보 확보
- 배치 효과가 제거된 고품질 통합 데이터
- 신뢰도 높은 최종 세포 타입 주석- Save final Seurat object as RDS
- Export all figures to figures/ directory
- Create summary CSV files for markers and pathways
- 고품질 클러스터링 결과
- 세포 타입별 마커 유전자
- 농축된 경로 분석
- 논문 제출용 그림들

---

### 2. Cell-Cell Communication Analysis

**목표**: 세포 간 통신 네트워크 분석 및 시그널링 경로 발견

**사용 스킬**: Seurat → NicheNet → CellChat → Scientific-Visualization-R

**Claude 프롬프트**:
```
Use available R skills you have access to whenever possible.

Step 1: Prepare Seurat object
- Load annotated Seurat object with cell types
- Ensure metadata has: celltype, condition, replicate

Step 2: NicheNet ligand-receptor analysis
- Load NicheNet networks (human or mouse)
- Define sender cells: "Fibroblast"
- Define receiver cells: "Epithelial"
- Get DEGs in receiver cells between conditions
- Run NicheNet to predict active ligands
- Identify top 20 ligands affecting target genes
- Find receptors for top ligands
- Create ligand activity plots and ligand-target heatmaps
- **Save**: `saveRDS(nichenet_output, "results/ccc/02_nichenet.rds")`

Step 3: CellChat network analysis
- Create CellChat object from Seurat
- Load CellChatDB (Secreted Signaling)
- Identify over-expressed genes and interactions
- Compute communication probability between all cell types
- Filter communications (min.cells=10)
- Compute pathway-level communication
- Aggregate network
- **Save**: `saveRDS(cellchat, "results/ccc/03_cellchat.rds")`

Step 4: Pattern discovery
- Identify outgoing communication patterns (k=3)
- Identify incoming communication patterns (k=3)
- Create river plots showing communication patterns
- Generate dot plots for signaling roles

Step 5: Specific pathway analysis
- Analyze TGFβ, WNT, NOTCH signaling pathways
- Create chord diagrams for each pathway
- Generate bubble plots showing ligand-receptor pairs
- Visualize gene expression for pathway components

Step 6: Publication figures
- Use scientific-visualization-r to create:
  * Circle plots showing overall communication
  * Heatmaps of pathway activities
  * Combined figures for multi-panel layouts
  * Export in journal-specific formats (Nature, Cell)

Step 7: Integration and summary
- Combine NicheNet and CellChat results
- Identify top ligand-receptor pairs from both methods
- Create comprehensive summary table
- Export all results as CSV files
```

**예상 결과**:
- 세포 간 통신 네트워크 맵
- 활성화된 ligand-receptor 쌍
- 주요 시그널링 경로
- 통신 패턴 분석

---

### 3. Treatment Response Analysis

**목표**: 치료 전후 또는 약물 처리 효과 분석

**사용 스킬**: Seurat → Pseudobulk-DESeq2 → Pathway-Enrichment → Scientific-Visualization-R

**Claude 프롬프트**:
```
Use available R skills you have access to whenever possible.

Step 1: Load and prepare data
- Load Seurat object with treatment and control samples
- Ensure metadata has: condition (treated/control), replicate, celltype
- Check sample distribution across conditions
- Verify you have at least 3 biological replicates per condition

Step 2: Seurat preprocessing
- For each condition separately:
  * Perform QC and filtering
  * Normalize and find variable features
- Integrate samples using Harmony
- Run clustering and UMAP
- Annotate cell types
- **Save**: `saveRDS(seurat_integrated, "results/treatment/02_integrated.rds")`

Step 3: Pseudobulk aggregation and DESeq2
- Use pseudobulk-deseq2 skill to:
  * Aggregate counts by condition + replicate + celltype
  * Run DESeq2 for each cell type
  * Compare treated vs control
  * Use LFC shrinkage (ashr method)
  * Get significant DEGs (padj < 0.05, |log2FC| > 1)
- **Save**:
  * `saveRDS(dds_list, "results/treatment/03_dds_objects.rds")`
  * `write.csv(all_degs, "results/treatment/03_all_degs.csv")`

Step 4: Cell type-specific responses
- For each cell type:
  * Count number of up/down-regulated genes
  * Create volcano plots
  * Identify cell types most responsive to treatment
- **Save**: `write.csv(response_summary, "results/treatment/04_summary.csv")`

Step 5: Pathway enrichment
- Use pathway-enrichment skill for each cell type:
  * Run fgsea with Hallmark gene sets
  * Run clusterProfiler for GO Biological Process
  * Find enriched pathways (padj < 0.05)
  * Identify common vs cell type-specific pathway responses

Step 6: Cross-cell-type analysis
- Find genes commonly affected across multiple cell types
- Create heatmap showing log2FC across all cell types
- Identify treatment response signatures

Step 7: Publication figures
- Use scientific-visualization-r to create:
  * Multi-panel volcano plots (one per cell type)
  * Heatmap of top DEGs across cell types
  * Pathway enrichment dot plots
  * UMAP showing treatment effects
  * Summary bar plots of DEG counts
- Export all figures in high resolution

Step 8: Report generation
- Create summary tables:
  * DEG counts per cell type
  * Top pathways per cell type
  * Common response genes
- Export as CSV and create markdown report
```

**예상 결과**:
- 세포 타입별 치료 반응
- 차별 발현 유전자 목록
- 농축된 경로 분석
- 치료 효과 시각화

---

### 4. Disease vs Healthy Comparison

**목표**: 질병과 정상 조직의 세포 수준 비교

**사용 스킬**: Seurat → Pseudobulk-DESeq2 → CellChat → Pathway-Enrichment → Scientific-Visualization-R

**Claude 프롬프트**:
```
Use available R skills you have access to whenever possible.

Step 1: Data integration
- Load disease and healthy Seurat objects
- Merge or integrate with Harmony
- Perform standard Seurat workflow
- Identify shared cell types
- **Save**: `saveRDS(seurat_merged, "results/disease/01_merged.rds")`

Step 2: Cell composition analysis
- Compare cell type proportions between conditions
- Create bar plots showing cell type frequencies
- Statistical test for composition changes
- **Save**: `write.csv(proportions, "results/disease/02_composition.csv")`

Step 3: Differential expression per cell type
- Use pseudobulk-deseq2 for robust DE testing:
  * For each cell type, compare disease vs healthy
  * Account for batch effects if present
  * Get significant DEGs for each cell type
- **Save**: `write.csv(disease_degs, "results/disease/03_disease_degs.csv")`

Step 4: Pathway dysregulation
- Use pathway-enrichment skill:
  * Run GSEA for each cell type
  * Identify disease-specific pathway alterations
  * Find common dysregulated pathways
  * Compare to Reactome and KEGG databases

Step 5: Communication network changes
- Use CellChat to compare networks:
  * Create CellChat objects for disease and healthy
  * Merge CellChat objects
  * Compare number and strength of interactions
  * Identify disease-specific communication changes
  * Find altered signaling pathways

Step 6: Disease signature discovery
- Identify genes consistently altered in disease
- Create disease signature gene sets
- Module score analysis in Seurat

Step 7: Visualization and reporting
- Use scientific-visualization-r to create:
  * Side-by-side UMAP comparisons
  * Heatmaps showing pathway activities
  * Communication network comparisons
  * Dot plots of disease markers
  * Multi-panel summary figures

Step 8: Clinical relevance
- Use pathway-enrichment with disease-related gene sets
- Identify druggable targets from altered pathways
- Create summary report with key findings
```

**예상 결과**:
- 질병 특이적 변화 규명
- 세포 타입별 병리학적 변화
- 변화된 통신 네트워크
- 잠재적 치료 타겟

---

### 5. Time Course Analysis

**목표**: 시간 경과에 따른 세포 반응 추적

**사용 스킬**: Seurat → Pseudobulk-DESeq2 → Pathway-Enrichment → Scientific-Visualization-R

**Claude 프롬프트**:
```
Use available R skills you have access to whenever possible.

Step 1: Load time series data
- Load Seurat objects for each time point (0h, 6h, 12h, 24h, 48h)
- Ensure consistent cell type annotations
- Add timepoint metadata

Step 2: Integration and clustering
- Integrate all time points with Harmony
- Maintain batch information
- Cluster and annotate cell types consistently
- **Save**: `saveRDS(seurat_time, "results/timecourse/02_integrated.rds")`

Step 3: Time course differential expression
- Use pseudobulk-deseq2 with time as continuous variable:
  * Design: ~ time + cell_type + time:cell_type
  * Identify genes with significant time effects
  * Find cell type-specific temporal responses
- **Save**: `saveRDS(lrt_results, "results/timecourse/03_lrt_results.rds")`

Step 4: Pairwise temporal comparisons
- For each consecutive time point pair:
  * Run DESeq2: time_n vs time_n-1
  * Identify early, middle, late response genes
  * Track gene expression trajectories

Step 5: Temporal pathway dynamics
- Use pathway-enrichment skill:
  * Run GSEA for each time point comparison
  * Track pathway activation over time
  * Identify temporal pathway patterns
  * Create time series of pathway scores

Step 6: Pattern identification
- Cluster genes by temporal expression patterns
- Identify:
  * Transiently activated genes
  * Sustained response genes
  * Late activation genes
- Create temporal heatmaps

Step 7: Cell type temporal dynamics
- Compare temporal responses across cell types
- Identify cell types with:
  * Early response
  * Delayed response
  * No response

Step 8: Visualization
- Use scientific-visualization-r to create:
  * Line plots showing gene expression over time
  * Heatmaps with temporal patterns
  * Pathway activity time series
  * Multi-panel figures for different cell types
  * Animated plots (optional)

Step 9: Summary and reporting
- Identify key temporal switches
- Find master regulators at each time point
- Create comprehensive temporal response summary
```

**예상 결과**:
- 시간에 따른 유전자 발현 패턴
- 경로 활성화 역학
- 세포 타입별 시간 반응
- 시간적 조절 인자

---

### 6. Multi-Sample Integration

**목표**: 여러 실험, 배치, 조건의 데이터 통합 분석

**사용 스킬**: Seurat → Pseudobulk-DESeq2 → NicheNet → Scientific-Visualization-R

**Claude 프롬프트**:
```
Use available R skills you have access to whenever possible.

Step 1: Load multiple datasets
- Load all Seurat objects (different experiments, batches, conditions)
- Check metadata consistency
- Ensure sample_id, batch, condition are properly labeled

Step 2: Quality control per sample
- Run QC for each sample individually
- Create QC summary table
- Identify and remove low-quality samples

Step 3: Integration strategy
- Use Seurat integration approaches:
  * Option 1: Standard CCA integration
  * Option 2: Harmony integration (faster, good for many samples)
  * Option 3: SCTransform integration
- Evaluate integration quality with batch mixing metrics
- **Save**:
  * `saveRDS(seurat_integrated, "results/integration/03_integrated.rds")`
  * `write.csv(metrics, "results/integration/03_metrics.csv")`

Step 4: Unified clustering
- Cluster integrated data
- Annotate cell types using consistent markers
- Verify cell types exist across batches
- **Save**: `saveRDS(seurat_annotated, "results/integration/04_annotated.rds")`

Step 5: Batch-aware differential expression
- Use pseudobulk-deseq2 with batch correction:
  * Design: ~ batch + condition
  * Account for batch effects
  * Run for each cell type
  * Get batch-corrected DEGs

Step 6: Cross-dataset validation
- Identify conserved markers across datasets
- Find dataset-specific effects
- Validate biological signals vs technical artifacts

Step 7: Communication analysis across conditions
- Use NicheNet for each condition separately
- Compare ligand activities between conditions
- Identify condition-specific ligand-receptor pairs

Step 8: Comprehensive visualization
- Use scientific-visualization-r to create:
  * UMAPs colored by batch, condition, cell type
  * Before/after integration comparisons
  * Batch effect assessment plots
  * Cross-dataset marker validation
  * Condition comparison figures

Step 9: Integration quality report
- Calculate integration metrics
- Create summary statistics
- Export integrated data
- Generate comprehensive analysis report
```

**예상 결과**:
- 통합된 다중 샘플 데이터
- 배치 효과 제거된 DE 분석
- 조건 간 검증된 생물학적 신호
- 품질 관리 보고서

---

## 🎯 Skill Combinations

### Core Workflows

#### Workflow 1: Basic scRNA-seq
```
Seurat → Scientific-Visualization-R
```
- QC, normalization, clustering
- 기본 시각화

#### Workflow 2: Advanced scRNA-seq
```
Seurat → Pseudobulk-DESeq2 → Pathway-Enrichment → Scientific-Visualization-R
```
- 완전한 발현 분석
- 통계적으로 robust한 DE
- 경로 분석

#### Workflow 3: Communication Analysis
```
Seurat → NicheNet + CellChat → Scientific-Visualization-R
```
- 세포 간 통신 규명
- 시그널링 경로 분석

#### Workflow 4: Complete Pipeline
```
Seurat → Celltype-Annotation → Pseudobulk-DESeq2 → Pathway-Enrichment → NicheNet → CellChat → Scientific-Visualization-R
```
- 자동 세포 타입 주석
- 전체 통합 분석
- 모든 측면의 세포 생물학

---

## 💡 Tips & Best Practices

### 1. 스킬 선택 가이드

**Seurat를 사용하세요**:
- scRNA-seq 기본 분석
- QC, 정규화, 클러스터링
- 세포 타입 주석
- 탐색적 데이터 분석

**Pseudobulk-DESeq2를 사용하세요**:
- 조건 간 비교 (disease vs healthy, treated vs control)
- 생물학적 반복이 있을 때 (최소 3개)
- 통계적으로 robust한 결과 필요시
- Batch effect 고려 필요시

**Pathway-Enrichment를 사용하세요**:
- DEG 리스트의 생물학적 의미 파악
- GSEA 분석
- GO/KEGG/Reactome 경로 농축
- 빠른 탐색 (enrichR)

**NicheNet을 사용하세요**:
- Ligand-receptor 예측
- 특정 유전자 변화의 upstream 조절자 찾기
- 기계론적 가설 생성

**CellChat을 사용하세요**:
- 전체 통신 네트워크 분석
- 패턴 발견
- 조건 간 통신 비교
- 시그널링 경로 농축

**Scientific-Visualization-R을 사용하세요**:
- 논문 제출용 그림
- 저널 특정 형식 (Nature, Science, Cell)
- 일관된 색상과 테마
- 고해상도 내보내기

**Celltype-Annotation을 사용하세요**:
- 자동 세포 타입 주석
- 여러 방법 비교 (ScType, Scibet, scMayoMap, Celltypist)
- 마커 기반, 레퍼런스 기반, Atlas 기반, 딥러닝 기반
- Consensus annotation으로 신뢰도 향상

### 2. 프롬프트 작성 팁

**구체적으로 요청하세요**:
```
❌ "Analyze my scRNA-seq data"
✅ "Use Seurat to load 10X data, perform QC (filter nFeature 200-2500, 
    percent.mt < 5%), normalize with LogNormalize, find 2000 variable 
    features, run PCA and UMAP with 30 dims, cluster with resolution 0.5"
```

**단계별로 나누세요**:
```
Step 1: Data loading and QC with Seurat
Step 2: Clustering and annotation
Step 3: Differential expression with Pseudobulk-DESeq2
Step 4: Pathway enrichment
Step 5: Visualization
```

**필요한 출력을 명시하세요**:
```
"Save results as:
- seurat_final.rds
- figures/umap.pdf (300 DPI)
- results/markers.csv
- results/pathways.csv"
```

### 3. 데이터 준비

**메타데이터 필수 항목**:
- `condition`: 조건 (control, treated, disease, healthy 등)
- `replicate`: 생물학적 반복 번호 (1, 2, 3, ...)
- `celltype`: 세포 타입 주석
- `batch`: 배치 정보 (있다면)

**파일 구조**:
```
project/
├── data/
│   ├── sample1/filtered_feature_bc_matrix/
│   ├── sample2/filtered_feature_bc_matrix/
│   └── ...
├── results/
├── figures/
└── analysis.R
```

### 4. 리소스 관리

**메모리 고려사항**:
- 큰 데이터셋: DietSeurat() 사용
- Integration: Harmony 추천 (메모리 효율적)
- 저장 전: 불필요한 assay 제거

**병렬 처리**:
- Pseudobulk-DESeq2: 여러 cell type 동시 분석 가능
- future 패키지 활용

### 5. 재현성

**항상 포함하세요**:
```r
# Session info
sessionInfo()

# Random seed
set.seed(42)

# Software versions
packageVersion("Seurat")
packageVersion("DESeq2")
```

**Save intermediate results**:
```r
saveRDS(seurat_qc, "intermediate/01_qc.rds")
saveRDS(seurat_clustered, "intermediate/02_clustered.rds")
saveRDS(deg_results, "intermediate/03_deg.rds")
```

---

## 📚 Example Data Sources

### Public Datasets
- **10X Genomics**: https://www.10xgenomics.com/resources/datasets
- **Human Cell Atlas**: https://data.humancellatlas.org/
- **GEO**: https://www.ncbi.nlm.nih.gov/geo/
- **Single Cell Portal**: https://singlecell.broadinstitute.org/

### Tutorial Datasets
- PBMC 3k (Seurat tutorial)
- Mouse brain (10X)
- COVID-19 BALF
- Cancer atlases (TCGA, CPTAC)

---

## 🔗 Related Resources

### R Package Documentation
- [Seurat](https://satijalab.org/seurat/)
- [DESeq2](https://bioconductor.org/packages/DESeq2/)
- [clusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/)
- [fgsea](https://bioconductor.org/packages/fgsea/)
- [NicheNet](https://github.com/saeyslab/nichenetr)
- [CellChat](https://github.com/sqjin/CellChat)

### Tutorials
- [Seurat Tutorials](https://satijalab.org/seurat/articles/)
- [Orchestrating Single-Cell Analysis with Bioconductor](http://bioconductor.org/books/OSCA/)
- [Current best practices in single‐cell RNA‐seq analysis](https://www.embopress.org/doi/full/10.15252/msb.20188746)

---

## 📝 Citation

사용한 스킬들을 논문에 인용하세요:

```bibtex
@article{hao2021seurat,
  title={Integrated analysis of multimodal single-cell data},
  author={Hao, Yuhan and others},
  journal={Cell},
  year={2021}
}

@article{love2014deseq2,
  title={Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2},
  author={Love, Michael I and others},
  journal={Genome biology},
  year={2014}
}

@article{browaeys2020nichenet,
  title={NicheNet: modeling intercellular communication by linking ligands to target genes},
  author={Browaeys, Robin and others},
  journal={Nature methods},
  year={2020}
}

@article{jin2021cellchat,
  title={Inference and analysis of cell-cell communication using CellChat},
  author={Jin, Suoqin and others},
  journal={Nature communications},
  year={2021}
}

@article{korotkevich2021fgsea,
  title={Fast gene set enrichment analysis},
  author={Korotkevich, Gennady and others},
  journal={bioRxiv},
  year={2021}
}

@article{yu2012clusterprofiler,
  title={clusterProfiler: an R package for comparing biological themes among gene clusters},
  author={Yu, Guangchuang and others},
  journal={Omics},
  year={2012}
}

@article{ianevski2022sctype,
  title={Fully-automated and ultra-fast cell-type identification using specific marker combinations},
  author={Ianevski, Aleksandr and others},
  journal={Nature Communications},
  year={2022}
}

@article{dominguez2022celltypist,
  title={Cross-tissue immune cell analysis reveals tissue-specific features in humans},
  author={Domínguez Conde, C and others},
  journal={Science},
  year={2022}
}
```

---

**Created**: 2025-12-24  
**Version**: 1.0  
**License**: MIT  

모든 스킬과 파이프라인은 `/home/kwy7605/LLM/skills/` 디렉토리에서 확인할 수 있습니다.
