# CellChat Workflow Examples

## Example 1: Basic CellChat Analysis

### Complete Workflow from Seurat to Communication Network

```r
library(CellChat)
library(Seurat)
library(patchwork)
library(dplyr)

# Load Seurat object
seurat_obj <- readRDS("seurat_annotated.rds")

# Create CellChat object
cellchat <- createCellChat(object = seurat_obj,
                           meta = seurat_obj@meta.data,
                           group.by = "celltype",
                           assay = "RNA")

# Load database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

# Use secreted signaling only
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

# Preprocessing
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Infer communication
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Pathway-level
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Compute centrality
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Visualize overall network
groupSize <- as.numeric(table(cellchat@idents))

par(mfrow = c(1,2))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Interaction strength")

# Save
saveRDS(cellchat, "cellchat_analysis.rds")
```

## Example 2: Focused Sender-Receiver Analysis

### Analyze Communication from Fibroblasts to Epithelial Cells

```r
library(CellChat)

# Load processed CellChat object
cellchat <- readRDS("cellchat_analysis.rds")

# Define sender and receiver
sender_celltypes <- c("Fibroblast_1", "Fibroblast_2")
receiver_celltypes <- c("Epithelial_1", "Epithelial_2")

# Get cluster IDs
sender_ids <- which(levels(cellchat@idents) %in% sender_celltypes)
receiver_ids <- which(levels(cellchat@idents) %in% receiver_celltypes)

# Extract interactions
df.net <- subsetCommunication(cellchat,
                              sources.use = sender_ids,
                              targets.use = receiver_ids)

# View top interactions
df.net %>%
  arrange(desc(prob)) %>%
  head(20)

# Circle plot for specific direction
groupSize <- as.numeric(table(cellchat@idents))

netVisual_circle(cellchat@net$weight,
                 sources.use = sender_ids,
                 targets.use = receiver_ids,
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Fibroblast to Epithelial")

# Chord diagram
netVisual_chord_gene(cellchat,
                     sources.use = sender_ids,
                     targets.use = receiver_ids,
                     slot.name = "netP",
                     legend.pos.x = 10)

# Bubble plot of ligand-receptor pairs
netVisual_bubble(cellchat,
                 sources.use = sender_ids,
                 targets.use = receiver_ids,
                 remove.isolate = FALSE,
                 angle.x = 45)
```

## Example 3: Pathway-Specific Analysis

### Analyze TGFβ, EGF, and WNT Signaling

```r
# Pathways to analyze
pathways.show <- c("TGFb", "EGF", "WNT")

# Chord diagrams for each pathway
pdf("pathway_chord_plots.pdf", width = 8, height = 8)
for (pathway in pathways.show) {
  netVisual_chord_gene(cellchat,
                       signaling = pathway,
                       slot.name = "netP",
                       legend.pos.x = 10,
                       title.name = pathway)
}
dev.off()

# Heatmaps
pdf("pathway_heatmaps.pdf", width = 10, height = 8)
par(mfrow = c(2,2))
for (pathway in pathways.show) {
  netVisual_heatmap(cellchat,
                    signaling = pathway,
                    color.heatmap = "Reds")
}
dev.off()

# Gene expression for pathways
plotGeneExpression(cellchat,
                   signaling = pathways.show,
                   enriched.only = FALSE)

# Signaling role network
for (pathway in pathways.show) {
  netAnalysis_signalingRole_network(cellchat,
                                    signaling = pathway,
                                    width = 10, height = 3)
}

# Bubble plot comparing pathways
netVisual_bubble(cellchat,
                 signaling = pathways.show,
                 remove.isolate = FALSE,
                 angle.x = 45)
ggsave("pathways_bubble.pdf", width = 12, height = 8)
```

## Example 4: Pattern Analysis

### Identify Communication Patterns

```r
library(CellChat)
library(NMF)
library(ggalluvial)

cellchat <- readRDS("cellchat_analysis.rds")

# Determine optimal number of patterns
selectK(cellchat, pattern = "outgoing")
selectK(cellchat, pattern = "incoming")

# Outgoing patterns (cell types as senders)
npatterns_out <- 3
cellchat <- identifyCommunicationPatterns(cellchat,
                                          pattern = "outgoing",
                                          k = npatterns_out,
                                          height = 9,
                                          width = 5)

# River plot for outgoing
netAnalysis_river(cellchat, pattern = "outgoing")
ggsave("outgoing_river.pdf", width = 8, height = 6)

# Dot plot for outgoing
netAnalysis_dot(cellchat, pattern = "outgoing")
ggsave("outgoing_dot.pdf", width = 8, height = 6)

# Incoming patterns (cell types as receivers)
npatterns_in <- 3
cellchat <- identifyCommunicationPatterns(cellchat,
                                          pattern = "incoming",
                                          k = npatterns_in,
                                          height = 9,
                                          width = 5)

# River plot for incoming
netAnalysis_river(cellchat, pattern = "incoming")
ggsave("incoming_river.pdf", width = 8, height = 6)

# Dot plot for incoming
netAnalysis_dot(cellchat, pattern = "incoming")
ggsave("incoming_dot.pdf", width = 8, height = 6)
```

## Example 5: Comparing Conditions

### Disease vs. Healthy Communication Networks

```r
library(CellChat)
library(patchwork)

# Load Seurat objects for each condition
seurat_healthy <- subset(seurat_obj, subset = condition == "healthy")
seurat_disease <- subset(seurat_obj, subset = condition == "disease")

# Function to process one condition
process_cellchat <- function(seurat_subset, condition_name) {
  
  cat(sprintf("Processing %s...\n", condition_name))
  
  # Create CellChat
  cc <- createCellChat(object = seurat_subset,
                       meta = seurat_subset@meta.data,
                       group.by = "celltype",
                       assay = "RNA")
  
  # Load database
  CellChatDB <- CellChatDB.human
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cc@DB <- CellChatDB.use
  
  # Process
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  
  return(cc)
}

# Process both conditions
cellchat_healthy <- process_cellchat(seurat_healthy, "Healthy")
cellchat_disease <- process_cellchat(seurat_disease, "Disease")

# Merge
cellchat_list <- list(Healthy = cellchat_healthy, Disease = cellchat_disease)
cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

# Compare number of interactions
compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1,2))
ggsave("compare_interactions_count.pdf", width = 6, height = 4)

# Compare interaction strength
compareInteractions(cellchat_merged, show.legend = FALSE,
                    group = c(1,2), measure = "weight")
ggsave("compare_interactions_weight.pdf", width = 6, height = 4)

# Differential interaction strength
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE)
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight")

# Identify changed signaling
pos.dataset <- "Disease"
features.name <- netAnalysis_signalingChanges_scatter(cellchat_merged)
ggsave("signaling_changes.pdf", width = 8, height = 8)

# Compare specific pathways
pathways.show <- c("TGFb", "EGF", "WNT", "NOTCH")

# Heatmap comparison
pdf("pathway_comparison_heatmaps.pdf", width = 12, height = 6)
par(mfrow = c(1,2))
for (i in 1:2) {
  netVisual_heatmap(cellchat_list[[i]],
                    signaling = pathways.show,
                    color.heatmap = "Reds")
  title(names(cellchat_list)[i])
}
dev.off()

# Save merged object
saveRDS(cellchat_merged, "cellchat_merged_comparison.rds")
```

## Example 6: Integration with Seurat Visualization

### Overlay CellChat Results on UMAP

```r
library(CellChat)
library(Seurat)
library(ggplot2)

# Get top pathways
pathways <- cellchat@netP$pathways
pathway_strengths <- sapply(pathways, function(p) {
  sum(cellchat@net$weight)
})

top_pathways <- names(sort(pathway_strengths, decreasing = TRUE))[1:6]

# For each pathway, get top ligands and receptors
for (pathway in top_pathways) {
  
  # Get pathway info
  pathway_genes <- CellChatDB.use$interaction %>%
    filter(pathway_name == pathway)
  
  ligands <- unique(pathway_genes$ligand)
  receptors <- unique(pathway_genes$receptor)
  
  # Visualize on UMAP
  if (length(ligands) > 0) {
    p_lig <- FeaturePlot(seurat_obj, features = ligands[1:min(4, length(ligands))],
                         ncol = 2)
    ggsave(sprintf("%s_ligands_umap.pdf", pathway), p_lig,
           width = 8, height = 8)
  }
  
  if (length(receptors) > 0) {
    p_rec <- FeaturePlot(seurat_obj, features = receptors[1:min(4, length(receptors))],
                         ncol = 2)
    ggsave(sprintf("%s_receptors_umap.pdf", pathway), p_rec,
           width = 8, height = 8)
  }
}

# Cell type composition in communication
DimPlot(seurat_obj, reduction = "umap", group.by = "celltype", label = TRUE)
ggsave("celltype_umap.pdf", width = 8, height = 6)
```

## Example 7: Export Results for Publication

### Generate All Key Figures

```r
# Create output directory
dir.create("cellchat_figures", showWarnings = FALSE)

# 1. Overall network
pdf("cellchat_figures/01_overall_network.pdf", width = 12, height = 6)
par(mfrow = c(1,2))
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Interaction strength")
dev.off()

# 2. Top pathways
top_pathways <- names(sort(pathway_strengths, decreasing = TRUE))[1:10]

pdf("cellchat_figures/02_top_pathways.pdf", width = 10, height = 12)
par(mfrow = c(3,2))
for (pathway in top_pathways[1:6]) {
  netVisual_aggregate(cellchat, signaling = pathway,
                      vertex.weight = groupSize,
                      layout = "circle")
}
dev.off()

# 3. Pathway analysis
pdf("cellchat_figures/03_pathway_bubble.pdf", width = 12, height = 8)
netVisual_bubble(cellchat, signaling = top_pathways,
                 remove.isolate = FALSE, angle.x = 45)
dev.off()

# 4. Signaling roles
pdf("cellchat_figures/04_signaling_roles.pdf", width = 12, height = 10)
netAnalysis_signalingRole_scatter(cellchat)
dev.off()

# 5. Pattern analysis
pdf("cellchat_figures/05_patterns.pdf", width = 10, height = 12)
par(mfrow = c(2,1))
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()

# 6. Export data table
write.csv(df.net, "cellchat_figures/interactions_table.csv", row.names = FALSE)

cat("All figures exported to cellchat_figures/\n")
```
