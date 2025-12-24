# NicheNet Workflow Examples

## Example 1: Basic Ligand Activity Analysis

### Predict Ligands Affecting Epithelial Cell DEGs

```r
library(nichenetr)
library(Seurat)
library(dplyr)
library(ggplot2)

# Load Seurat object
seurat_obj <- readRDS("seurat_integrated.rds")

# Load NicheNet networks
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

lr_network <- lr_network %>% distinct(from, to)

# Define sender and receiver
sender_cells <- "Fibroblast"
receiver_cells <- "Epithelial"

# Get DEGs in receiver cells (disease vs healthy)
Idents(seurat_obj) <- "condition"
receiver_subset <- subset(seurat_obj, subset = celltype == "Epithelial")

de_results <- FindMarkers(receiver_subset,
                          ident.1 = "diseased",
                          ident.2 = "healthy",
                          min.pct = 0.10)

geneset_oi <- rownames(de_results %>% 
                        filter(p_val_adj < 0.05 & avg_log2FC > 0.5))

cat(sprintf("Gene set of interest: %d genes\n", length(geneset_oi)))

# Get expressed genes
get_expressed <- function(cell_type, seurat_obj, pct = 0.10) {
  cells <- WhichCells(seurat_obj, idents = cell_type)
  expr <- GetAssayData(seurat_obj, slot = "data")[, cells]
  rownames(expr)[rowSums(expr > 0) >= pct * ncol(expr)]
}

sender_expressed <- get_expressed(sender_cells, seurat_obj)
receiver_expressed <- get_expressed(receiver_cells, seurat_obj)
background_expressed <- rownames(seurat_obj)[rowSums(GetAssayData(seurat_obj, slot = "counts")) > 0]

# Get potential ligands
ligands <- lr_network %>% pull(from) %>% unique()
expressed_ligands <- intersect(ligands, sender_expressed)
potential_ligands <- expressed_ligands %>% intersect(colnames(ligand_target_matrix))

cat(sprintf("Potential ligands: %d\n", length(potential_ligands)))

# Predict ligand activities
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected)

# Get top 20 ligands
best_upstream_ligands <- ligand_activities %>%
  top_n(20, aupr_corrected) %>%
  arrange(-aupr_corrected) %>%
  pull(test_ligand)

print(best_upstream_ligands)

# Visualize
library(ggplot2)
top20 <- ligand_activities %>% top_n(20, aupr_corrected) %>% arrange(aupr_corrected)
top20$test_ligand <- factor(top20$test_ligand, levels = top20$test_ligand)

p <- ggplot(top20, aes(x = aupr_corrected, y = test_ligand)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Ligand Activity", y = "", title = "Top 20 Predicted Ligands") +
  theme_minimal()

ggsave("ligand_activities.pdf", p, width = 6, height = 8)
```

## Example 2: Ligand-Target Network Analysis

### Identify Target Genes of Top Ligands

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

# Filter for genes of interest
active_ligand_target_links_genes_oi <- active_ligand_target_links %>%
  filter(target %in% geneset_oi)

# Prepare for visualization
ligand_target_matrix_vis <- active_ligand_target_links_genes_oi %>%
  spread(target, weight, fill = 0) %>%
  column_to_rownames("ligand") %>%
  as.matrix()

# Heatmap
library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(c(0, max(ligand_target_matrix_vis)), c("white", "red"))

ht <- Heatmap(
  ligand_target_matrix_vis,
  name = "Regulatory\nPotential",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title = "Target Genes",
  row_title = "Ligands"
)

pdf("ligand_target_heatmap.pdf", width = 12, height = 10)
draw(ht)
dev.off()
```

## Example 3: Ligand-Receptor Pair Analysis

### Identify Active L-R Pairs

```r
# Get receptors
receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors <- intersect(receptors, receiver_expressed)

# Get L-R pairs for top ligands
lr_network_top <- lr_network %>%
  filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>%
  distinct(from, to) %>%
  rename(ligand = from, receptor = to)

# Add ligand activity scores
lr_network_top_df <- lr_network_top %>%
  inner_join(ligand_activities %>% select(test_ligand, aupr_corrected),
             by = c("ligand" = "test_ligand")) %>%
  arrange(-aupr_corrected)

head(lr_network_top_df, 20)

# Visualize expression of L-R pairs
library(Seurat)

# Top pairs
top_lr <- lr_network_top_df %>% head(10)

# Extract ligands and receptors
ligands_to_plot <- unique(top_lr$ligand)
receptors_to_plot <- unique(top_lr$receptor)

# Subset Seurat object
sender_subset <- subset(seurat_obj, idents = sender_cells)
receiver_subset <- subset(seurat_obj, idents = receiver_cells)

# Dot plots
p1 <- DotPlot(sender_subset, features = ligands_to_plot) +
  coord_flip() +
  ggtitle("Ligands (Sender)")

p2 <- DotPlot(receiver_subset, features = receptors_to_plot) +
  coord_flip() +
  ggtitle("Receptors (Receiver)")

library(patchwork)
combined <- p1 | p2
ggsave("ligand_receptor_expression.pdf", combined, width = 14, height = 8)
```

## Example 4: Multi-Sender Analysis

### Compare Multiple Sender Cell Types

```r
# Define multiple sender types
sender_types <- c("Fibroblast", "Macrophage", "Endothelial")

# Run analysis for each sender
results_list <- lapply(sender_types, function(sender) {
  
  cat(sprintf("\nAnalyzing %s as sender...\n", sender))
  
  # Get expressed genes
  sender_expressed <- get_expressed(sender, seurat_obj)
  
  # Get potential ligands
  ligands <- lr_network %>% pull(from) %>% unique()
  expressed_ligands <- intersect(ligands, sender_expressed)
  potential_ligands <- expressed_ligands %>% intersect(colnames(ligand_target_matrix))
  
  # Predict activities
  ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  )
  
  # Add sender info
  ligand_activities$sender <- sender
  
  return(ligand_activities)
})

# Combine results
all_ligand_activities <- bind_rows(results_list)

# Get top ligands per sender
top_ligands_per_sender <- all_ligand_activities %>%
  group_by(sender) %>%
  top_n(10, aupr_corrected) %>%
  arrange(sender, -aupr_corrected)

# Visualize comparison
library(ggplot2)

p <- ggplot(top_ligands_per_sender, 
            aes(x = aupr_corrected, y = test_ligand, fill = sender)) +
  geom_bar(stat = "identity") +
  facet_wrap(~sender, scales = "free_y", ncol = 1) +
  labs(x = "Ligand Activity", y = "Ligand") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("multi_sender_comparison.pdf", p, width = 8, height = 12)
```

## Example 5: Condition-Specific Analysis

### Compare Disease vs Healthy Communication

```r
# Subset by condition
seurat_diseased <- subset(seurat_obj, subset = condition == "diseased")
seurat_healthy <- subset(seurat_obj, subset = condition == "healthy")

# Function to run NicheNet for one condition
run_nichenet_condition <- function(seurat_subset, condition_name) {
  
  # Get DEGs (compare to overall average or other approach)
  # Here we'll use genes highly expressed in receiver cells
  receiver_subset <- subset(seurat_subset, idents = receiver_cells)
  
  # Get top variable genes as proxy for interesting genes
  receiver_subset <- FindVariableFeatures(receiver_subset, nfeatures = 500)
  geneset_oi <- VariableFeatures(receiver_subset)
  
  # Get expressed genes
  sender_expressed <- get_expressed(sender_cells, seurat_subset)
  receiver_expressed <- get_expressed(receiver_cells, seurat_subset)
  background <- rownames(seurat_subset)[rowSums(GetAssayData(seurat_subset)) > 0]
  
  # Potential ligands
  ligands <- lr_network %>% pull(from) %>% unique()
  expressed_ligands <- intersect(ligands, sender_expressed)
  potential_ligands <- expressed_ligands %>% intersect(colnames(ligand_target_matrix))
  
  # Predict
  ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  )
  
  ligand_activities$condition <- condition_name
  
  return(ligand_activities)
}

# Run for both conditions
diseased_results <- run_nichenet_condition(seurat_diseased, "diseased")
healthy_results <- run_nichenet_condition(seurat_healthy, "healthy")

# Compare
combined_results <- bind_rows(diseased_results, healthy_results)

# Find condition-specific ligands
diseased_specific <- diseased_results %>%
  top_n(20, aupr_corrected) %>%
  pull(test_ligand)

healthy_specific <- healthy_results %>%
  top_n(20, aupr_corrected) %>%
  pull(test_ligand)

# Unique to disease
disease_only <- setdiff(diseased_specific, healthy_specific)
cat("Disease-specific ligands:\n")
print(disease_only)

# Compare top ligands
top_both <- combined_results %>%
  filter(test_ligand %in% c(diseased_specific, healthy_specific)) %>%
  select(test_ligand, aupr_corrected, condition) %>%
  spread(condition, aupr_corrected, fill = 0)

# Scatter plot
p <- ggplot(top_both, aes(x = healthy, y = diseased)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text_repel(aes(label = test_ligand), size = 3, max.overlaps = 15) +
  labs(x = "Ligand Activity (Healthy)", 
       y = "Ligand Activity (Diseased)",
       title = "Condition-Specific Ligand Activities") +
  theme_minimal()

ggsave("condition_comparison.pdf", p, width = 8, height = 8)
```

## Example 6: Integration with Seurat Visualization

### Overlay NicheNet Results on UMAP

```r
# Get top ligands and their expression
top_ligands <- best_upstream_ligands[1:6]

# Feature plots for ligands (in sender cells)
FeaturePlot(seurat_obj, features = top_ligands, ncol = 3)
ggsave("top_ligands_umap.pdf", width = 12, height = 8)

# Get top receptors
top_receptors <- lr_network_top_df %>%
  filter(ligand %in% top_ligands) %>%
  pull(receptor) %>%
  unique() %>%
  head(6)

# Feature plots for receptors (in receiver cells)
FeaturePlot(seurat_obj, features = top_receptors, ncol = 3)
ggsave("top_receptors_umap.pdf", width = 12, height = 8)

# Violin plots split by condition
VlnPlot(seurat_obj, features = top_ligands[1:3], 
        split.by = "condition", ncol = 3)
ggsave("top_ligands_violin.pdf", width = 12, height = 4)
```
