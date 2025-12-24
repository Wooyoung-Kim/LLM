# CellChat - Cell-Cell Communication Network Analysis

Infer and analyze intercellular communication networks from scRNA-seq data.

## Overview

CellChat quantifies cell-cell communication, identifies signaling pathways, and compares networks across conditions using a curated database of ligand-receptor interactions.

## Quick Start

```r
library(CellChat)
library(Seurat)
library(patchwork)

# Create object
cellchat <- createCellChat(object = seurat_obj,
                           group.by = "celltype",
                           assay = "RNA")

# Load database and preprocess
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Infer communication
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Visualize
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize)

# Identify patterns
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 3)
netAnalysis_river(cellchat, pattern = "outgoing")
```

## Key Features

- **Communication Inference**: Quantify signaling between cell types
- **Pathway Analysis**: Identify active signaling pathways
- **Pattern Discovery**: Find communication patterns
- **Comparison**: Compare networks across conditions
- **Visualization**: Circle, chord, bubble, heatmap, river plots

## Common Use Cases

### Compare Conditions
```r
# Create objects for each condition
cellchat_list <- list(Control = cc_ctrl, Disease = cc_dis)
cellchat_merged <- mergeCellChat(cellchat_list)

# Compare
compareInteractions(cellchat_merged)
netVisual_diffInteraction(cellchat_merged)
```

### Specific Pathway Analysis
```r
# Visualize TGFb signaling
netVisual_chord_gene(cellchat, signaling = "TGFb")
plotGeneExpression(cellchat, signaling = "TGFb")
```

## Citation

```
Jin et al. (2021). Inference and analysis of cell-cell communication 
using CellChat. Nature Communications.
```

## Documentation

- Main skill: [SKILL.md](SKILL.md)
- Official docs: https://github.com/sqjin/CellChat
