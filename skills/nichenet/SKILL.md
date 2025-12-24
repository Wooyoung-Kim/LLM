---
name: nichenet
description: "NicheNet for ligand-receptor analysis: predict ligand-target links, identify active ligand-receptor pairs, prioritize ligands, cell-cell communication inference from scRNA-seq data."
---

# NicheNet - Ligand-Receptor and Cell Communication Analysis

## Overview

NicheNet predicts which ligands from sender cells are most likely to affect gene expression in receiver cells, based on prior knowledge of signaling and gene regulatory networks. It's particularly useful for understanding cell-cell communication mechanisms in scRNA-seq data.

## When to Use This Skill

This skill should be used when:
- Predicting which ligands affect gene expression in target cells
- Identifying active ligand-receptor pairs between cell types
- Prioritizing ligands based on target gene expression changes
- Understanding intercellular communication mechanisms
- Analyzing condition-specific communication (e.g., disease vs. healthy)
- Validating cell-cell communication hypotheses
- Finding upstream regulators of DEGs

## Quick Start Guide

### Basic NicheNet Workflow

```r
library(nichenetr)
library(Seurat)
library(tidyverse)

# 1. Load NicheNet networks (human)
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

# 2. Define cell populations
sender_cells <- "Fibroblast"
receiver_cells <- "Epithelial"

# 3. Define gene set of interest (e.g., DEGs in receiver cells)
geneset_oi <- c("GENE1", "GENE2", "GENE3")  # Your genes of interest

# 4. Define background genes (all expressed genes)
background_expressed_genes <- rownames(seurat_obj)[rowSums(seurat_obj@assays$RNA@counts) > 0]

# 5. Get expressed genes in sender and receiver
sender_expressed <- get_expressed_genes(sender_cells, seurat_obj, pct = 0.10)
receiver_expressed <- get_expressed_genes(receiver_cells, seurat_obj, pct = 0.10)

# 6. Get potential ligands
ligands <- lr_network %>% pull(from) %>% unique()
expressed_ligands <- intersect(ligands, sender_expressed)

# 7. Get receptors
receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors <- intersect(receptors, receiver_expressed)

# 8. Get potential ligands targeting genes of interest
potential_ligands <- expressed_ligands %>%
  intersect(colnames(ligand_target_matrix))

# 9. Perform ligand activity analysis
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

# 10. Get top ligands
best_upstream_ligands <- ligand_activities %>%
  top_n(20, aupr_corrected) %>%
  arrange(-aupr_corrected) %>%
  pull(test_ligand)

# 11. Infer target genes of top ligands
active_ligand_target_links <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links, 
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 250) %>%
  bind_rows()

# 12. Get receptors for top ligands
lr_network_top <- lr_network %>%
  filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>%
  distinct(from, to) %>%
  rename(ligand = from, receptor = to)

# Visualize
p <- make_mushroom_plot(ligand_activities, top_n = 20)
print(p)
```

## Core Workflow Steps

### 1. Load NicheNet Networks

#### Human Networks

```r
library(nichenetr)

# Option 1: Load from Zenodo (requires internet)
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

# Option 2: Load from local files (faster)
lr_network <- readRDS("data/nichenet/lr_network_human_21122021.rds")
ligand_target_matrix <- readRDS("data/nichenet/ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS("data/nichenet/weighted_networks_nsga2r_final.rds")

# Inspect networks
head(lr_network)
dim(ligand_target_matrix)
```

#### Mouse Networks

```r
# Mouse networks
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
```

### 2. Define Sender and Receiver Cell Populations

```r
library(Seurat)
library(dplyr)

# Load Seurat object
seurat_obj <- readRDS("seurat_object.rds")

# Option 1: Define by cell type annotation
sender_cells <- "Fibroblast"
receiver_cells <- "Epithelial"

# Option 2: Define by cluster
sender_cells <- WhichCells(seurat_obj, idents = c(0, 1, 2))
receiver_cells <- WhichCells(seurat_obj, idents = c(3, 4))

# Option 3: Define by condition and cell type
sender_cells <- seurat_obj@meta.data %>%
  filter(celltype == "Fibroblast" & condition == "diseased") %>%
  rownames()

receiver_cells <- seurat_obj@meta.data %>%
  filter(celltype == "Epithelial" & condition == "diseased") %>%
  rownames()
```

### 3. Define Gene Set of Interest

```r
# Option 1: From differential expression analysis
Idents(seurat_obj) <- "condition"
receiver_subset <- subset(seurat_obj, subset = celltype == "Epithelial")

de_results <- FindMarkers(receiver_subset,
                          ident.1 = "diseased",
                          ident.2 = "healthy",
                          min.pct = 0.10)

geneset_oi <- rownames(de_results %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.25))

# Option 2: Manually defined gene set
geneset_oi <- c("COL1A1", "FN1", "ACTA2", "TGFB1", "VIM")

# Option 3: From gene set enrichment
# Genes enriched in specific pathway
geneset_oi <- c("genes", "from", "pathway", "analysis")
```

### 4. Get Expressed Genes

```r
#' Get expressed genes in specific cell population
#'
#' @param cell_ids Cell identities (cell type name or cell barcodes)
#' @param seurat_obj Seurat object
#' @param pct Minimum percentage of cells expressing gene (default: 0.10)
#' @return Vector of expressed gene names
get_expressed_genes <- function(cell_ids, seurat_obj, pct = 0.10) {
  if (length(cell_ids) == 1 && is.character(cell_ids)) {
    # Cell type name provided
    cells <- WhichCells(seurat_obj, idents = cell_ids)
  } else {
    # Cell barcodes provided
    cells <- cell_ids
  }
  
  # Get expression data
  expression_data <- GetAssayData(seurat_obj, slot = "data")[, cells]
  
  # Find genes expressed in > pct of cells
  expressed_genes <- rownames(expression_data)[
    rowSums(expression_data > 0) >= pct * ncol(expression_data)
  ]
  
  return(expressed_genes)
}

# Use function
sender_expressed <- get_expressed_genes(sender_cells, seurat_obj, pct = 0.10)
receiver_expressed <- get_expressed_genes(receiver_cells, seurat_obj, pct = 0.10)

# Background: all genes expressed in dataset
background_expressed_genes <- rownames(seurat_obj)[
  rowSums(GetAssayData(seurat_obj, slot = "counts")) > 0
]
```

### 5. Define Potential Ligands and Receptors

```r
# Get all ligands and receptors from network
ligands <- lr_network %>% pull(from) %>% unique()
receptors <- lr_network %>% pull(to) %>% unique()

# Filter by expression in sender and receiver
expressed_ligands <- intersect(ligands, sender_expressed)
expressed_receptors <- intersect(receptors, receiver_expressed)

# Filter ligands that are in the ligand-target matrix
potential_ligands <- expressed_ligands %>%
  intersect(colnames(ligand_target_matrix))

cat(sprintf("Expressed ligands: %d\n", length(expressed_ligands)))
cat(sprintf("Expressed receptors: %d\n", length(expressed_receptors)))
cat(sprintf("Potential ligands in matrix: %d\n", length(potential_ligands)))
```

### 6. Ligand Activity Analysis

```r
# Predict which ligands can regulate the gene set of interest
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

# View results
ligand_activities <- ligand_activities %>%
  arrange(-aupr_corrected)

head(ligand_activities, 20)

# Select top ligands
best_upstream_ligands <- ligand_activities %>%
  top_n(20, aupr_corrected) %>%
  arrange(-aupr_corrected) %>%
  pull(test_ligand)

cat("\nTop predicted upstream ligands:\n")
print(best_upstream_ligands)
```

### 7. Infer Target Genes of Top Ligands

```r
# Get weighted ligand-target links
active_ligand_target_links <- best_upstream_ligands %>%
  lapply(function(ligand) {
    get_weighted_ligand_target_links(
      best_upstream_ligands = ligand,
      expressed_genes = receiver_expressed,
      ligand_target_matrix = ligand_target_matrix,
      n = 250
    )
  }) %>%
  bind_rows()

# View top targets
active_ligand_target_links %>%
  arrange(-weight) %>%
  head(20)

# Prepare for visualization
active_ligand_target_links_df <- active_ligand_target_links %>%
  filter(target %in% geneset_oi) %>%
  prepare_ligand_target_visualization(
    best_upstream_ligands = best_upstream_ligands,
    order_ligands = TRUE
  )
```

### 8. Identify Active Ligand-Receptor Pairs

```r
# Get ligand-receptor links for top ligands
lr_network_top <- lr_network %>%
  filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>%
  distinct(from, to) %>%
  rename(ligand = from, receptor = to)

# Add expression data
lr_network_top_df <- lr_network_top %>%
  inner_join(ligand_activities %>% select(test_ligand, aupr_corrected), 
             by = c("ligand" = "test_ligand"))

# View results
lr_network_top_df %>%
  arrange(-aupr_corrected)
```

## Visualization

### 1. Mushroom Plot - Ligand Activity

```r
library(ggplot2)

# Create mushroom plot
p <- make_mushroom_plot(ligand_activities, top_n = 20)
print(p)

# Save
ggsave("ligand_activities.pdf", p, width = 6, height = 8)
```

### 2. Heatmap - Ligand-Target Links

```r
library(ComplexHeatmap)
library(circlize)

# Prepare matrix for heatmap
ligand_target_matrix_plot <- active_ligand_target_links %>%
  filter(target %in% geneset_oi) %>%
  spread(target, weight, fill = 0) %>%
  column_to_rownames("ligand") %>%
  as.matrix()

# Create heatmap
ht <- Heatmap(
  ligand_target_matrix_plot,
  name = "Regulatory\nPotential",
  col = colorRamp2(c(0, max(ligand_target_matrix_plot)), c("white", "red")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8)
)

pdf("ligand_target_heatmap.pdf", width = 10, height = 8)
draw(ht)
dev.off()
```

### 3. Dot Plot - Ligand-Receptor Expression

```r
library(Seurat)

# Prepare ligand-receptor pairs for visualization
lr_pairs <- lr_network_top_df %>%
  mutate(pair = paste(ligand, receptor, sep = " - ")) %>%
  pull(pair)

# Create expression data frame
sender_subset <- subset(seurat_obj, cells = sender_cells)
receiver_subset <- subset(seurat_obj, cells = receiver_cells)

# Dot plot of top pairs
ligands_to_plot <- lr_network_top_df$ligand %>% unique() %>% head(10)
receptors_to_plot <- lr_network_top_df %>%
  filter(ligand %in% ligands_to_plot) %>%
  pull(receptor) %>%
  unique()

# Sender ligands
p1 <- DotPlot(sender_subset, features = ligands_to_plot) +
  coord_flip() +
  ggtitle("Ligands (Sender cells)")

# Receiver receptors
p2 <- DotPlot(receiver_subset, features = receptors_to_plot) +
  coord_flip() +
  ggtitle("Receptors (Receiver cells)")

library(patchwork)
combined <- p1 | p2
ggsave("ligand_receptor_expression.pdf", combined, width = 12, height = 6)
```

### 4. Circos Plot - Ligand-Receptor Interactions

```r
# Prepare for circos plot
lr_network_top_for_circos <- lr_network_top_df %>%
  head(20)

# Create circos plot
p <- make_circos_lr(lr_network_top_for_circos, 
                    ligand_colors = "red",
                    receptor_colors = "blue")

pdf("ligand_receptor_circos.pdf", width = 8, height = 8)
print(p)
dev.off()
```

## Advanced Analyses

See `references/nichenet_examples.md` for:
- Condition-specific analysis (disease vs. healthy)
- Multi-sender analysis
- Custom prioritization weights
- Integration with other tools
- Validation strategies

## Best Practices

### Data Preparation
- Use high-quality, well-annotated Seurat objects
- Ensure cell type annotations are accurate
- Filter low-quality cells before analysis
- Consider using SCTransform-normalized data

### Gene Set Definition
- Use stringent thresholds for DEGs (p < 0.05, |log2FC| > 0.5)
- Focus on biologically relevant genes
- Consider both upregulated and downregulated genes separately
- Validate gene sets with pathway analysis

### Ligand Selection
- Set appropriate expression thresholds (typically 10%)
- Consider tissue-specific expression patterns
- Validate top ligands with literature
- Check for biological plausibility

### Interpretation
- Combine with experimental validation
- Consider timing of signaling events
- Account for spatial context
- Validate with orthogonal methods

## Common Pitfalls to Avoid

1. **Too many genes of interest**: Keep geneset focused (<500 genes)
2. **Low expression threshold**: May include spurious ligands
3. **Wrong species**: Ensure networks match your data
4. **Ignoring expression**: Top ligands must be expressed in sender cells
5. **Over-interpretation**: NicheNet predictions are hypotheses, not facts
6. **Batch effects**: Correct batch effects before analysis
7. **Cell contamination**: Remove doublets and low-quality cells
8. **Wrong cell populations**: Ensure sender/receiver are correctly defined

## Resources

### Scripts Directory
- `nichenet_helpers.R` - Helper functions for common tasks

### References Directory
- `nichenet_examples.md` - Complete workflow examples
- `ligand_databases.md` - Common ligand families

## NicheNet Networks

Pre-computed networks available for:
- **Human**: GRCh38
- **Mouse**: GRCm38

Networks include:
- **lr_network**: Ligand-receptor interactions
- **ligand_target_matrix**: Prior knowledge of ligand-target regulatory potential
- **weighted_networks**: Weighted signaling and gene regulatory networks

## Citation

If you use NicheNet in your research, cite:

```
Browaeys et al. (2020). NicheNet: modeling intercellular communication 
by linking ligands to target genes. Nature Methods. 
doi:10.1038/s41592-019-0667-5
```

Use this skill to predict and analyze ligand-mediated cell-cell communication in your scRNA-seq data.
