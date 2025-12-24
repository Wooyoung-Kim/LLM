# NicheNet - Ligand-Receptor Analysis

Predict ligand-mediated cell-cell communication using prior knowledge of signaling networks.

## Overview

NicheNet identifies which ligands from sender cells most likely affect gene expression in receiver cells, providing mechanistic insights into intercellular communication.

## Quick Start

```r
library(nichenetr)
library(Seurat)
library(tidyverse)
source('scripts/nichenet_helpers.R')

# Load networks
networks <- load_nichenet_networks(species = "human", local_path = "data/nichenet/")

# Define populations
sender <- "Fibroblast"
receiver <- "Epithelial"

# Run analysis
results <- run_nichenet_analysis(
  seurat_obj = seurat_obj,
  sender_cells = sender,
  receiver_cells = receiver,
  geneset_oi = deg_genes,  # Your DEGs
  lr_network = networks$lr_network,
  ligand_target_matrix = networks$ligand_target_matrix
)

# Visualize
plot_ligand_activities(results$ligand_activities, top_n = 20)
plot_ligand_target_heatmap(results$ligand_target_links, geneset_oi = deg_genes)
```

## Key Features

- **Ligand Activity Prediction**: Rank ligands by regulatory potential
- **Target Gene Inference**: Identify which genes are regulated by each ligand  
- **Receptor Identification**: Find expressed receptors for active ligands
- **Visualization Tools**: Built-in plotting functions
- **Multi-condition Support**: Compare communication across conditions

## Common Use Cases

### Predict Communication in Disease
```r
# Get DEGs in diseased receiver cells
geneset_oi <- find_degs_in_receiver(seurat_obj, 
                                     receiver = "Epithelial",
                                     condition1 = "diseased",
                                     condition2 = "healthy")

# Run NicheNet
results <- run_nichenet_analysis(...)
```

### Identify Receptors
```r
# Get receptors for top ligands
receptors <- identify_receptors(
  ligands = results$best_ligands,
  lr_network = networks$lr_network,
  receiver_expressed = receiver_genes
)
```

## Citation

```
Browaeys et al. (2020). NicheNet: modeling intercellular communication 
by linking ligands to target genes. Nature Methods.
```

## Documentation

- Main skill: [SKILL.md](SKILL.md)
- Official docs: https://github.com/saeyslab/nichenetr
