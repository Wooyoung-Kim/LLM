---
name: cellchat
description: "CellChat for cell-cell communication: infer communication networks, identify signaling pathways, compare groups, visualize interactions with circle/chord/bubble plots, pattern analysis."
---

# CellChat - Cell-Cell Communication Network Analysis

## Overview

CellChat quantitatively infers and analyzes intercellular communication networks from scRNA-seq data. It identifies biologically significant cell-cell interactions, infers communication patterns, and compares communication networks across different datasets or conditions.

## When to Use This Skill

This skill should be used when:
- Inferring cell-cell communication networks from scRNA-seq
- Quantifying signaling between specific cell types
- Identifying dominant senders and receivers
- Comparing communication networks across conditions (disease vs. healthy)
- Finding enriched signaling pathways
- Discovering communication patterns
- Visualizing intercellular interactions
- Analyzing signaling roles of cell populations

## Quick Start Guide

### Basic CellChat Workflow

```r
library(CellChat)
library(Seurat)
library(patchwork)
library(dplyr)

# 1. Create CellChat object from Seurat
cellchat <- createCellChat(object = seurat_obj, 
                           meta = seurat_obj@meta.data,
                           group.by = "celltype",
                           assay = "RNA")

# 2. Load database
CellChatDB <- CellChatDB.human  # or CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# Use subset (optional)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

# 3. Preprocessing
cellchat <- subsetData(cellchat)  # Subset to signaling genes
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 4. Infer communication
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# 5. Infer pathway-level communication
cellchat <- computeCommunProbPathway(cellchat)

# 6. Aggregate network
cellchat <- aggregateNet(cellchat)

# 7. Visualize
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = TRUE, label.edge = FALSE)

# 8. Identify patterns
cellchat <- identifyCommunicationPatterns(cellchat, 
                                          pattern = "outgoing", k = 3)
cellchat <- identifyCommunicationPatterns(cellchat, 
                                          pattern = "incoming", k = 3)

# 9. Save
saveRDS(cellchat, "cellchat_analysis.rds")
```

## Core Workflow Steps

### 1. Create CellChat Object

#### From Seurat Object

```r
library(CellChat)
library(Seurat)

# Load Seurat object
seurat_obj <- readRDS("seurat_object.rds")

# Create CellChat object
cellchat <- createCellChat(object = seurat_obj,
                           meta = seurat_obj@meta.data,
                           group.by = "celltype",  # Column for cell grouping
                           assay = "RNA")

# Check object
cellchat
slotNames(cellchat)
```

#### From Data Matrix

```r
# From count matrix
data.input <- GetAssayData(seurat_obj, slot = "data")  # Normalized data
meta <- seurat_obj@meta.data

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
```

### 2. Load CellChat Database

```r
# Human database
CellChatDB <- CellChatDB.human

# Mouse database
# CellChatDB <- CellChatDB.mouse

# View database structure
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

# View available categories
unique(CellChatDB$interaction$annotation)

# Subset database (optional)
# Option 1: Secreted signaling only
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

# Option 2: ECM-Receptor interactions only
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor")

# Option 3: Cell-Cell Contact
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact")

# Option 4: Use full database
CellChatDB.use <- CellChatDB

# Assign to CellChat object
cellchat@DB <- CellChatDB.use
```

### 3. Preprocessing

```r
# Set the used database
cellchat@DB <- CellChatDB.use

# Subset data to signaling genes (reduces computational cost)
cellchat <- subsetData(cellchat)

# Identify overexpressed genes and interactions
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Optional: Project gene expression onto protein-protein interaction network
# cellchat <- projectData(cellchat, PPI.human)
```

### 4. Infer Cell-Cell Communication

```r
# Compute communication probability
cellchat <- computeCommunProb(cellchat)

# Optional: Compute using a different method
# cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)

# Filter out communications with few cells
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract communication probability
df.net <- subsetCommunication(cellchat)
head(df.net)

# Extract specific interactions
df.net.specific <- subsetCommunication(cellchat, 
                                        sources.use = c("Fibroblast"),
                                        targets.use = c("Epithelial"))
```

### 5. Infer Pathway-Level Communication

```r
# Compute pathway-level communication
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate network
cellchat <- aggregateNet(cellchat)

# View aggregated network
cellchat@net$count  # Number of interactions
cellchat@net$weight  # Interaction strength
```

### 6. Calculate Communication Network Centrality

```r
# Compute centrality measures
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Visualize signaling roles
netAnalysis_signalingRole_network(cellchat, 
                                   signaling = "TGFb",
                                   width = 8, height = 2.5)
```

## Visualization

### 1. Circle Plot - Overall Communication

```r
# Prepare data
groupSize <- as.numeric(table(cellchat@idents))

# Number of interactions
par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(cellchat@net$count, 
                 vertex.weight = groupSize,
                 weight.scale = TRUE, 
                 label.edge = FALSE,
                 title.name = "Number of interactions")

# Interaction strength
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Interaction weights/strength")
```

### 2. Circle Plot - Specific Sources/Targets

```r
# Visualize from specific senders to specific receivers
par(mfrow = c(1,1))

# Example: From Fibroblasts to Epithelial cells
sender_ids <- c(1, 2)    # Cluster IDs
receiver_ids <- c(3, 4)  # Cluster IDs

netVisual_circle(cellchat@net$weight,
                 sources.use = sender_ids,
                 targets.use = receiver_ids,
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Fibroblast to Epithelial")
```

### 3. Chord Diagram - Specific Pathway

```r
# Visualize specific signaling pathway
pathways.show <- "TGFb"

# Chord diagram for ligand-receptor pairs
netVisual_chord_gene(cellchat, 
                     signaling = pathways.show,
                     slot.name = "netP",
                     legend.pos.x = 10)

# From specific sources
netVisual_chord_gene(cellchat,
                     sources.use = c(1, 2),
                     targets.use = c(3, 4),
                     signaling = pathways.show,
                     legend.pos.x = 10)
```

### 4. Bubble Plot - L-R Pairs

```r
# Show all significant interactions
netVisual_bubble(cellchat, 
                 sources.use = c(1, 2),
                 targets.use = c(3, 4),
                 remove.isolate = FALSE,
                 angle.x = 45)

# Specific pathway
netVisual_bubble(cellchat,
                 sources.use = c(1),
                 targets.use = c(3),
                 signaling = c("TGFb", "EGF"),
                 remove.isolate = FALSE)
```

### 5. Heatmap - Communication Probability

```r
# Heatmap showing all interactions
netVisual_heatmap(cellchat, 
                  signaling = "TGFb",
                  color.heatmap = "Reds")

# Heatmap for multiple pathways
pathways.show <- c("TGFb", "EGF", "WNT", "NOTCH")
par(mfrow = c(2,2))
for (pathway in pathways.show) {
  netVisual_heatmap(cellchat, signaling = pathway, 
                    color.heatmap = "Reds")
}
```

### 6. Dot Plot - Signaling Patterns

```r
# Outgoing signaling patterns
netAnalysis_dot(cellchat, pattern = "outgoing")

#Incoming signaling patterns
netAnalysis_dot(cellchat, pattern = "incoming")
```

### 7. River Plot - Communication Patterns

```r
# Requires pattern identification first
cellchat <- identifyCommunicationPatterns(cellchat, 
                                          pattern = "outgoing", k = 3)

# Visualize
netAnalysis_river(cellchat, pattern = "outgoing")

# Incoming patterns
cellchat <- identifyCommunicationPatterns(cellchat, 
                                          pattern = "incoming", k = 3)
netAnalysis_river(cellchat, pattern = "incoming")
```

### 8. Gene Expression - Signaling Genes

```r
# Plot gene expression for specific pathway
plotGeneExpression(cellchat, signaling = "TGFb")

# Multiple pathways
plotGeneExpression(cellchat, 
                   signaling = c("TGFb", "EGF", "WNT"),
                   enriched.only = FALSE)
```

## Comparing Multiple Datasets

### Load and Prepare Multiple Datasets

```r
# Create CellChat objects for each condition
cellchat_control <- createCellChat(...)  # Control
cellchat_treated <- createCellChat(...)  # Treated

# Process both
process_cellchat <- function(cc) {
  cc@DB <- CellChatDB.use
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc)
  return(cc)
}

cellchat_control <- process_cellchat(cellchat_control)
cellchat_treated <- process_cellchat(cellchat_treated)
```

### Merge and Compare

```r
# Create list
cellchat_list <- list(Control = cellchat_control, 
                      Treated = cellchat_treated)

# Merge
cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

# Compare number of interactions
compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1,2))

# Compare interaction strength
compareInteractions(cellchat_merged, show.legend = FALSE, 
                    group = c(1,2), measure = "weight")

# Differential number and strength
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE)
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight")

# Heatmap comparison
par(mfrow = c(1,2))
netVisual_heatmap(cellchat_control, color.heatmap = "Reds")
netVisual_heatmap(cellchat_treated, color.heatmap = "Reds")

# Identify changed signaling
pos.dataset <- "Treated"
features.name <- netAnalysis_signalingChanges_scatter(cellchat_merged)
```

## Pattern Analysis

### Identify Communication Patterns

```r
# Determine number of patterns
selectK(cellchat, pattern = "outgoing")  # Suggests K value
selectK(cellchat, pattern = "incoming")

# Identify outgoing patterns
npatterns <- 3
cellchat <- identifyCommunicationPatterns(cellchat, 
                                          pattern = "outgoing",
                                          k = npatterns,
                                          height = 9)

# Identify incoming patterns
cellchat <- identifyCommunicationPatterns(cellchat,
                                          pattern = "incoming", 
                                          k = npatterns,
                                          height = 9)

# Visualize
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "incoming")
```

## Best Practices

### Data Preparation
- Use well-annotated Seurat objects
- Remove low-quality cells before analysis
- Ensure cell type annotations are accurate
- Consider using SCTransform-normalized data

### Parameter Selection
- **min.cells**: Set to 10 or higher to filter noise
- **trim**: Adjust if using truncatedMean (typically 0.1)
- **Database subset**: Choose appropriate interaction types
- **Pattern number (k)**: Use selectK() to determine optimal value

### Interpretation
- Focus on biologically relevant pathways
- Validate with literature and experiments
- Consider spatial context in tissues
- Account for temporal dynamics
- Cross-reference with NicheNet for mechanism

### Visualization
- Use appropriate plot types for your question
- Combine multiple visualization approaches
- Export high-resolution figures for publication
- Use consistent color schemes

## Common Pitfalls to Avoid

1. **Too few cells**: Need sufficient cells per population (>20)
2. **Wrong database**: Ensure human/mouse matches your data
3. **Over-filtering**: Too stringent thresholds remove real interactions
4. **Ignoring expression**: Communication requires receptor expression
5. **Batch effects**: Correct before CellChat analysis
6. **Wrong normalization**: Use appropriate assay (RNA, not SCT)
7. **Missing metadata**: Ensure cell type annotations are present
8. **Single condition**: Compare multiple conditions when possible

## Resources

### Scripts Directory
- `cellchat_helpers.R` - Helper functions for common tasks

### References Directory
- `cellchat_examples.md` - Complete workflow examples
- `signaling_pathways.md` - Common signaling pathways

## CellChat Database

Pre-compiled databases available for:
- **Human**: CellChatDB.human
- **Mouse**: CellChatDB.mouse

Database categories:
- **Secreted Signaling**: Cytokines, growth factors, etc.
- **ECM-Receptor**: Extracellular matrix interactions
- **Cell-Cell Contact**: Membrane-bound interactions
- **Non-protein Signaling**: Metabolic signaling

## Citation

If you use CellChat in your research, cite:

```
Jin et al. (2021). Inference and analysis of cell-cell communication 
using CellChat. Nature Communications. 
doi:10.1038/s41467-021-21246-9
```

Use this skill to comprehensively analyze intercellular communication networks in your scRNA-seq data.
