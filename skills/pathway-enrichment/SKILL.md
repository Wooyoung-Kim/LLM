---
name: pathway-enrichment
description: "Pathway enrichment analysis: GSEA, fgsea, clusterProfiler, enrichR for gene set enrichment, GO/KEGG/Reactome analysis, visualization of enrichment results."
---

# Pathway Enrichment Analysis

## Overview

Pathway enrichment analysis identifies biological pathways and processes that are significantly enriched in your gene lists. This skill covers multiple R packages for comprehensive pathway analysis: fgsea for fast GSEA, clusterProfiler for GO/KEGG/Reactome, and enrichR for access to hundreds of gene set libraries.

## When to Use This Skill

This skill should be used when:
- Analyzing differential expression gene (DEG) lists
- Identifying enriched biological pathways
- Performing Gene Set Enrichment Analysis (GSEA)
- Comparing pathway enrichment across conditions
- Visualizing enrichment results for publication
- Accessing multiple gene set databases (MSigDB, GO, KEGG, Reactome, WikiPathways)
- Understanding biological functions of gene clusters

## Quick Start Guide

### Basic GSEA with fgsea

```r
library(fgsea)
library(dplyr)
library(ggplot2)

# 1. Prepare ranked gene list
deg_results <- read.csv("deg_results.csv")
deg_results$gene <- rownames(deg_results)

# Create ranked list (by log2FC or stat)
gene_ranks <- deg_results$avg_log2FC
names(gene_ranks) <- deg_results$gene
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

# 2. Load pathways (Hallmark)
library(msigdbr)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

# 3. Run fgsea
fgsea_results <- fgsea(
  pathways = pathways,
  stats = gene_ranks,
  minSize = 15,
  maxSize = 500,
  nperm = 10000
)

# 4. View results
fgsea_results %>%
  filter(padj < 0.05) %>%
  arrange(desc(NES))

# 5. Visualize top pathways
top_pathways <- fgsea_results %>%
  filter(padj < 0.05) %>%
  top_n(20, wt = abs(NES))

ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES)) +
  geom_col(aes(fill = padj < 0.01)) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score") +
  theme_minimal()
```

### Over-Representation Analysis with clusterProfiler

```r
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# 1. Get significant DEGs
deg_up <- deg_results %>%
  filter(avg_log2FC > 0.5 & p_val_adj < 0.05) %>%
  pull(gene)

# 2. GO enrichment
go_results <- enrichGO(
  gene = deg_up,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# 3. Visualize
dotplot(go_results, showCategory = 20)
barplot(go_results, showCategory = 20)
enrichplot::cnetplot(go_results, categorySize = "pvalue")

# 4. KEGG enrichment
# Convert to Entrez IDs
gene_entrez <- bitr(deg_up, 
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

kegg_results <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

dotplot(kegg_results, showCategory = 20)
```

### Quick Enrichment with enrichR

```r
library(enrichR)

# 1. List available databases
dbs <- listEnrichrDbs()
head(dbs)

# 2. Select databases
selected_dbs <- c(
  "GO_Biological_Process_2023",
  "KEGG_2021_Human",
  "WikiPathways_2023_Human",
  "Reactome_2022"
)

# 3. Run enrichment
enriched <- enrichr(deg_up, selected_dbs)

# 4. View results
enriched$GO_Biological_Process_2023 %>%
  filter(Adjusted.P.value < 0.05) %>%
  head(20)

# 5. Plot
plotEnrich(enriched[[1]], showTerms = 20, 
           title = "GO Biological Process")
```

## Core Workflow Steps

### 1. Prepare Gene Lists

#### From DEG Results

```r
library(dplyr)

# Load DEG results
deg_results <- read.csv("deg_results.csv", row.names = 1)
deg_results$gene <- rownames(deg_results)

# Option 1: Significant genes only (for ORA)
sig_genes <- deg_results %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>%
  pull(gene)

# Option 2: Upregulated only
up_genes <- deg_results %>%
  filter(avg_log2FC > 0.5 & p_val_adj < 0.05) %>%
  pull(gene)

# Option 3: Downregulated only
down_genes <- deg_results %>%
  filter(avg_log2FC < -0.5 & p_val_adj < 0.05) %>%
  pull(gene)

# Option 4: Ranked list (for GSEA)
gene_ranks <- deg_results$avg_log2FC
names(gene_ranks) <- deg_results$gene
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

# Remove NA values
gene_ranks <- gene_ranks[!is.na(gene_ranks)]
```

### 2. fgsea - Fast Gene Set Enrichment Analysis

#### Load Gene Sets

```r
library(fgsea)
library(msigdbr)

# Hallmark gene sets
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways_hallmark <- split(hallmark$gene_symbol, hallmark$gs_name)

# C2: Curated gene sets (KEGG, Reactome, BioCarta)
c2 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP")
pathways_c2 <- split(c2$gene_symbol, c2$gs_name)

# C5: GO gene sets
c5_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
pathways_go_bp <- split(c5_bp$gene_symbol, c5_bp$gs_name)

# Or load from GMT file
pathways_custom <- gmtPathways("custom_pathways.gmt")
```

#### Run fgsea

```r
# Run GSEA
fgsea_results <- fgsea(
  pathways = pathways_hallmark,
  stats = gene_ranks,
  minSize = 15,      # Minimum pathway size
  maxSize = 500,     # Maximum pathway size
  nperm = 10000      # Number of permutations
)

# Filter significant pathways
sig_pathways <- fgsea_results %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(NES)))

print(sig_pathways, n = 20)
```

#### Visualize fgsea Results

```r
library(ggplot2)

# Bar plot of top pathways
top_pathways <- fgsea_results %>%
  filter(padj < 0.05) %>%
  top_n(20, wt = abs(NES)) %>%
  mutate(pathway = gsub("^HALLMARK_", "", pathway),
         pathway = gsub("_", " ", pathway))

ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES)) +
  geom_col(aes(fill = -log10(padj))) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       fill = "-log10(padj)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Enrichment plot for specific pathway
plotEnrichment(pathways_hallmark[["HALLMARK_INFLAMMATORY_RESPONSE"]], 
               gene_ranks) +
  labs(title = "Inflammatory Response")

# Table plot
plotGseaTable(
  pathways_hallmark[top_pathways$pathway[1:10]], 
  gene_ranks,
  fgsea_results,
  gseaParam = 0.5
)
```

### 3. clusterProfiler - Comprehensive Pathway Analysis

#### GO Enrichment

```r
library(clusterProfiler)
library(org.Hs.eg.db)

# Biological Process
go_bp <- enrichGO(
  gene = sig_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# Molecular Function
go_mf <- enrichGO(
  gene = sig_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "MF",
  pvalueCutoff = 0.05
)

# Cellular Component
go_cc <- enrichGO(
  gene = sig_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "CC",
  pvalueCutoff = 0.05
)

# View results
head(go_bp, 20)
```

#### KEGG Pathway Enrichment

```r
# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(sig_genes,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

# KEGG enrichment
kegg <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# Convert IDs back to symbols
kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
```

#### Reactome Pathway Enrichment

```r
library(ReactomePA)

# Reactome enrichment
reactome <- enrichPathway(
  gene = gene_entrez$ENTREZID,
  organism = "human",
  pvalueCutoff = 0.05,
  readable = TRUE
)
```

#### GSEA with clusterProfiler

```r
# Prepare ranked list with Entrez IDs
gene_list <- deg_results$avg_log2FC
names(gene_list) <- deg_results$gene

# Convert to Entrez
gene_df <- bitr(names(gene_list),
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

gene_list_entrez <- gene_list[gene_df$SYMBOL]
names(gene_list_entrez) <- gene_df$ENTREZID
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

# GO GSEA
gsea_go <- gseGO(
  geneList = gene_list_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# KEGG GSEA
gsea_kegg <- gseKEGG(
  geneList = gene_list_entrez,
  organism = "hsa",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)
```

#### Visualization

```r
library(enrichplot)

# Dot plot
dotplot(go_bp, showCategory = 20) +
  ggtitle("GO Biological Process")

# Bar plot
barplot(go_bp, showCategory = 20)

# Gene-concept network
cnetplot(go_bp, categorySize = "pvalue", foldChange = gene_list)

# Heatmap plot
heatplot(go_bp, foldChange = gene_list)

# Enrichment map
emapplot(pairwise_termsim(go_bp), showCategory = 30)

# GSEA plot
gseaplot2(gsea_go, geneSetID = 1:3)

# Ridge plot
ridgeplot(gsea_go) + labs(x = "Enrichment Distribution")

# Upset plot
upsetplot(go_bp)
```

### 4. enrichR - Access Multiple Databases

#### Quick Enrichment

```r
library(enrichR)

# List all available databases
dbs <- listEnrichrDbs()
View(dbs)

# Popular databases
databases <- c(
  "GO_Biological_Process_2023",
  "GO_Molecular_Function_2023",
  "GO_Cellular_Component_2023",
  "KEGG_2021_Human",
  "WikiPathways_2023_Human",
  "Reactome_2022",
  "MSigDB_Hallmark_2020",
  "BioPlanet_2019",
  "Panther_2016"
)

# Run enrichment
enriched <- enrichr(sig_genes, databases)

# View results for each database
enriched$GO_Biological_Process_2023 %>%
  filter(Adjusted.P.value < 0.05) %>%
  arrange(Adjusted.P.value) %>%
  head(20)

# Plot specific database
plotEnrich(enriched$GO_Biological_Process_2023, 
           showTerms = 20,
           numChar = 50,
           y = "Count", 
           orderBy = "P.value",
           title = "GO BP Enrichment")
```

## Visualization

See `references/enrichment_examples.md` for detailed visualization examples.

## Best Practices

### Gene List Preparation
- Remove duplicate gene names
- Use official gene symbols or Entrez IDs
- Filter out low-quality genes (low expression, high dropout)
- Consider background gene set (all detected genes)

### Statistical Considerations
- Use FDR/BH correction for multiple testing
- Set appropriate p-value cutoffs (typically 0.05)
- Consider pathway size (min 10-15, max 500 genes)
- Use ranked lists for GSEA when possible

### Choosing Methods
- **GSEA/fgsea**: When you have full ranked gene lists
- **ORA (enrichGO, enrichKEGG)**: When you have discrete gene sets
- **enrichR**: For quick exploration of many databases
- **clusterProfiler**: For comprehensive, reproducible analysis

### Interpretation
- Focus on biological coherence, not just statistics
- Consider pathway overlap and redundancy
- Validate key findings with literature
- Report exact p-values and gene sets used

## Common Pitfalls to Avoid

1. **Small gene lists**: Need at least 10-20 genes for meaningful enrichment
2. **Wrong organism**: Ensure gene names match the species database
3. **ID mismatch**: Convert between gene symbols and Entrez IDs correctly
4. **No background**: Use appropriate background for ORA
5. **Over-interpretation**: Don't read too much into marginally significant pathways
6. **Ignoring pathway size**: Very small/large pathways can be unreliable
7. **Multiple testing**: Always use FDR correction
8. **Outdated databases**: Use recent versions of gene set databases

## Resources

### Scripts Directory
- `enrichment_helpers.R` - Helper functions for common workflows

### References Directory
- `enrichment_examples.md` - Complete workflow examples
- `databases_guide.md` - Guide to gene set databases

## Citation

When using these tools, cite:

```
# fgsea
Korotkevich et al. (2021). Fast gene set enrichment analysis. bioRxiv.

# clusterProfiler
Yu et al. (2012). clusterProfiler: an R package for comparing biological themes 
among gene clusters. OMICS.

# enrichR
Kuleshov et al. (2016). Enrichr: a comprehensive gene set enrichment analysis 
web server 2016 update. Nucleic Acids Research.

# MSigDB
Liberzon et al. (2015). The Molecular Signatures Database Hallmark Gene Set Collection. 
Cell Systems.
```

Use this skill for comprehensive pathway enrichment analysis of your genomics data.
