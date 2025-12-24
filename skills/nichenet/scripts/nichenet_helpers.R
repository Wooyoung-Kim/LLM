# Helper functions for NicheNet analysis

library(nichenetr)
library(dplyr)
library(Seurat)

#' Load NicheNet networks
#'
#' @param species Species ("human" or "mouse")
#' @param local_path Optional local path to saved networks
#' @return List containing lr_network, ligand_target_matrix, weighted_networks
#' @export
load_nichenet_networks <- function(species = "human", local_path = NULL) {
  
  if (!is.null(local_path)) {
    # Load from local files
    lr_network <- readRDS(file.path(local_path, 
                          paste0("lr_network_", species, "_21122021.rds")))
    ligand_target_matrix <- readRDS(file.path(local_path,
                          paste0("ligand_target_matrix_nsga2r_final",
                                 ifelse(species == "mouse", "_mouse", ""), ".rds")))
    weighted_networks <- readRDS(file.path(local_path,
                          paste0("weighted_networks_nsga2r_final",
                                 ifelse(species == "mouse", "_mouse", ""), ".rds")))
  } else {
    # Load from Zenodo
    base_url <- "https://zenodo.org/record/7074291/files/"
    
    lr_network <- readRDS(url(paste0(base_url, "lr_network_", species, "_21122021.rds")))
    
    if (species == "mouse") {
      ligand_target_matrix <- readRDS(url(paste0(base_url, "ligand_target_matrix_nsga2r_final_mouse.rds")))
      weighted_networks <- readRDS(url(paste0(base_url, "weighted_networks_nsga2r_final_mouse.rds")))
    } else {
      ligand_target_matrix <- readRDS(url(paste0(base_url, "ligand_target_matrix_nsga2r_final.rds")))
      weighted_networks <- readRDS(url(paste0(base_url, "weighted_networks_nsga2r_final.rds")))
    }
  }
  
  # Clean lr_network
  lr_network <- lr_network %>% distinct(from, to)
  
  cat(sprintf("Loaded NicheNet networks for %s:\n", species))
  cat(sprintf("  L-R pairs: %d\n", nrow(lr_network)))
  cat(sprintf("  Ligands in matrix: %d\n", ncol(ligand_target_matrix)))
  cat(sprintf("  Target genes: %d\n", nrow(ligand_target_matrix)))
  
  return(list(
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    weighted_networks = weighted_networks
  ))
}

#' Get expressed genes in specific cell population
#'
#' @param cell_ids Cell type name or cell barcodes
#' @param seurat_obj Seurat object
#' @param pct Minimum percentage of cells expressing gene
#' @return Vector of expressed gene names
#' @export
get_expressed_genes <- function(cell_ids, seurat_obj, pct = 0.10) {
  
  if (length(cell_ids) == 1 && is.character(cell_ids)) {
    cells <- WhichCells(seurat_obj, idents = cell_ids)
  } else {
    cells <- cell_ids
  }
  
  expression_data <- GetAssayData(seurat_obj, slot = "data")[, cells]
  
  expressed_genes <- rownames(expression_data)[
    rowSums(expression_data > 0) >= pct * ncol(expression_data)
  ]
  
  return(expressed_genes)
}

#' Run complete NicheNet analysis
#'
#' @param seurat_obj Seurat object
#' @param sender_cells Sender cell type or barcodes
#' @param receiver_cells Receiver cell type or barcodes
#' @param geneset_oi Genes of interest (e.g., DEGs)
#' @param lr_network Ligand-receptor network
#' @param ligand_target_matrix Ligand-target matrix
#' @param top_n Number of top ligands to return
#' @return List with results
#' @export
run_nichenet_analysis <- function(seurat_obj,
                                   sender_cells,
                                   receiver_cells,
                                   geneset_oi,
                                   lr_network,
                                   ligand_target_matrix,
                                   top_n = 20) {
  
  # Get expressed genes
  sender_expressed <- get_expressed_genes(sender_cells, seurat_obj, pct = 0.10)
  receiver_expressed <- get_expressed_genes(receiver_cells, seurat_obj, pct = 0.10)
  
  # Background genes
  background_expressed_genes <- rownames(seurat_obj)[
    rowSums(GetAssayData(seurat_obj, slot = "counts")) > 0
  ]
  
  # Get potential ligands
  ligands <- lr_network %>% pull(from) %>% unique()
  expressed_ligands <- intersect(ligands, sender_expressed)
  potential_ligands <- expressed_ligands %>% intersect(colnames(ligand_target_matrix))
  
  # Get receptors
  receptors <- lr_network %>% pull(to) %>% unique()
  expressed_receptors <- intersect(receptors, receiver_expressed)
  
  cat(sprintf("\nNicheNet Analysis:\n"))
  cat(sprintf("  Sender cells: %d\n", length(sender_cells)))
  cat(sprintf("  Receiver cells: %d\n", length(receiver_cells)))
  cat(sprintf("  Genes of interest: %d\n", length(geneset_oi)))
  cat(sprintf("  Potential ligands: %d\n", length(potential_ligands)))
  cat(sprintf("  Expressed receptors: %d\n\n", length(expressed_receptors)))
  
  # Predict ligand activities
  ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  )
  
  ligand_activities <- ligand_activities %>% arrange(-aupr_corrected)
  
  # Get top ligands
  best_upstream_ligands <- ligand_activities %>%
    top_n(top_n, aupr_corrected) %>%
    arrange(-aupr_corrected) %>%
    pull(test_ligand)
  
  # Get ligand-target links
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
  
  # Get ligand-receptor pairs
  lr_network_top <- lr_network %>%
    filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>%
    distinct(from, to) %>%
    rename(ligand = from, receptor = to) %>%
    inner_join(ligand_activities %>% select(test_ligand, aupr_corrected),
               by = c("ligand" = "test_ligand"))
  
  cat("Top 10 ligands:\n")
  print(head(ligand_activities, 10))
  
  return(list(
    ligand_activities = ligand_activities,
    best_ligands = best_upstream_ligands,
    ligand_target_links = active_ligand_target_links,
    lr_pairs = lr_network_top,
    expressed_ligands = expressed_ligands,
    expressed_receptors = expressed_receptors
  ))
}

#' Find DEGs in receiver cells between conditions
#'
#' @param seurat_obj Seurat object
#' @param receiver Receiver cell type
#' @param condition1 First condition
#' @param condition2 Second condition
#' @param condition_column Metadata column for condition
#' @return Vector of DEG names
#' @export
find_degs_in_receiver <- function(seurat_obj,
                                   receiver,
                                   condition1,
                                   condition2,
                                   condition_column = "condition") {
  
  # Subset to receiver cells
  receiver_subset <- subset(seurat_obj, subset = celltype == receiver)
  
  # Set identity to condition
  Idents(receiver_subset) <- condition_column
  
  # Find DEGs
  de_results <- FindMarkers(receiver_subset,
                            ident.1 = condition1,
                            ident.2 = condition2,
                            min.pct = 0.10,
                            logfc.threshold = 0.25)
  
  # Filter significant genes
  geneset_oi <- rownames(de_results %>% 
                          filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25))
  
  cat(sprintf("Found %d DEGs in %s (%s vs %s)\n", 
              length(geneset_oi), receiver, condition1, condition2))
  
  return(geneset_oi)
}

#' Plot ligand activities
#'
#' @param ligand_activities Data frame from predict_ligand_activities
#' @param top_n Number of top ligands to show
#' @export
plot_ligand_activities <- function(ligand_activities, top_n = 20) {
  library(ggplot2)
  
  top_ligands <- ligand_activities %>%
    top_n(top_n, aupr_corrected) %>%
    arrange(aupr_corrected)
  
  top_ligands$test_ligand <- factor(top_ligands$test_ligand,
                                     levels = top_ligands$test_ligand)
  
  p <- ggplot(top_ligands, aes(x = aupr_corrected, y = test_ligand)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(x = "Ligand Activity (AUPR)", y = "Ligand", 
         title = "Top Predicted Ligands") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10))
  
  return(p)
}
