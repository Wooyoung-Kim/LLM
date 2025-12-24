# Pathway Enrichment Analysis

Gene set enrichment and pathway analysis using fgsea, clusterProfiler, and enrichR.

## Quick Start

```r
library(fgsea)
library(msigdbr)
library(dplyr)

# Prepare ranked gene list
gene_ranks <- deg_results$avg_log2FC
names(gene_ranks) <- deg_results$gene
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

# Load Hallmark gene sets
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

# Run GSEA
fgsea_results <- fgsea(pathways = pathways, stats = gene_ranks)

# View significant
fgsea_results %>%
  filter(padj < 0.05) %>%
  arrange(desc(NES))
```

## Key Features

- **GSEA**: Fast gene set enrichment with fgsea
- **ORA**: Over-representation analysis with clusterProfiler
- **Multiple databases**: GO, KEGG, Reactome, WikiPathways, MSigDB
- **Visualization**: Dot plots, bar plots, enrichment maps
- **Quick access**: enrichR for 100+ databases

## Common Use Cases

### GO Enrichment
```r
library(clusterProfiler)
library(org.Hs.eg.db)

go_results <- enrichGO(
  gene = sig_genes,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05
)

dotplot(go_results)
```

### Quick enrichR
```r
library(enrichR)
dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human")
enriched <- enrichr(sig_genes, dbs)
```

## Citation

```
Korotkevich et al. (2021). Fast gene set enrichment analysis.
Yu et al. (2012). clusterProfiler. OMICS.
```

## Documentation

- Main skill: [SKILL.md](SKILL.md)
- Examples: [references/enrichment_examples.md](references/enrichment_examples.md)
