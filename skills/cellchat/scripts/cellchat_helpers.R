# Helper functions for CellChat analysis

library(CellChat)
library(Seurat)
library(dplyr)

#' Create and preprocess CellChat object from Seurat
#'
#' @param seurat_obj Seurat object
#' @param group_by Column name for cell grouping
#' @param assay Assay to use (default: "RNA")
#' @param species Species ("human" or "mouse")
#' @param database_type Type of interactions ("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact", or "all")
#' @return Preprocessed CellChat object
#' @export
create_cellchat <- function(seurat_obj,
                             group_by = "celltype",
                             assay = "RNA",
                             species = "human",
                             database_type = "Secreted Signaling") {
  
  # Create CellChat object
  cellchat <- createCellChat(object = seurat_obj,
                             meta = seurat_obj@meta.data,
                             group.by = group_by,
                             assay = assay)
  
  # Load database
  if (species == "human") {
    CellChatDB <- CellChatDB.human
  } else {
    CellChatDB <- CellChatDB.mouse
  }
  
  # Subset database if needed
  if (database_type != "all") {
    CellChatDB.use <- subsetDB(CellChatDB, search = database_type)
  } else {
    CellChatDB.use <- CellChatDB
  }
  
  cellchat@DB <- CellChatDB.use
  
  cat(sprintf("Created CellChat object:\n"))
  cat(sprintf("  Cell groups: %d\n", length(unique(cellchat@idents))))
  cat(sprintf("  Database: %s (%s)\n", species, database_type))
  cat(sprintf("  Interactions: %d\n", nrow(CellChatDB.use$interaction)))
  
  return(cellchat)
}

#' Run complete CellChat workflow
#'
#' @param cellchat CellChat object
#' @param min_cells Minimum cells for filtering
#' @return Processed CellChat object
#' @export
run_cellchat_workflow <- function(cellchat, min_cells = 10) {
  
  cat("Running CellChat workflow...\n")
  
  # Preprocessing
  cat("  1. Preprocessing...\n")
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Compute communication probability
  cat("  2. Computing communication probability...\n")
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  
  # Pathway-level communication
  cat("  3. Computing pathway-level communication...\n")
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Aggregate network
  cat("  4. Aggregating network...\n")
  cellchat <- aggregateNet(cellchat)
  
  # Compute centrality
  cat("  5. Computing centrality...\n")
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  cat("CellChat workflow complete!\n")
  cat(sprintf("  Number of interactions: %d\n", nrow(cellchat@net$count)))
  cat(sprintf("  Total interaction weight: %.2f\n", sum(cellchat@net$weight)))
  
  return(cellchat)
}

#' Visualize overall communication network
#'
#' @param cellchat CellChat object
#' @param plot_type Type of plot ("circle", "chord", "heatmap")
#' @param measure Measure to visualize ("count" or "weight")
#' @export
visualize_overall_network <- function(cellchat, 
                                       plot_type = "circle",
                                       measure = "count") {
  
  groupSize <- as.numeric(table(cellchat@idents))
  
  if (plot_type == "circle") {
    par(mfrow = c(1,1))
    netVisual_circle(cellchat@net[[measure]],
                     vertex.weight = groupSize,
                     weight.scale = TRUE,
                     label.edge = FALSE,
                     title.name = paste("Interaction", measure))
  } else if (plot_type == "heatmap") {
    netVisual_heatmap(cellchat, measure = measure)
  }
}

#' Extract top signaling pathways
#'
#' @param cellchat CellChat object
#' @param top_n Number of top pathways
#' @return Data frame of top pathways
#' @export
get_top_pathways <- function(cellchat, top_n = 10) {
  
  pathways <- cellchat@netP$pathways
  pathway_strengths <- sapply(pathways, function(pathway) {
    sum(cellchat@net$weight)
  })
  
  top_pathways <- names(sort(pathway_strengths, decreasing = TRUE))[1:top_n]
  
  return(data.frame(
    pathway = top_pathways,
    strength = pathway_strengths[top_pathways]
  ))
}

#' Compare CellChat objects from different conditions
#'
#' @param cellchat_list List of CellChat objects
#' @param condition_names Names of conditions
#' @return Merged CellChat object
#' @export
compare_cellchat_conditions <- function(cellchat_list, condition_names) {
  
  names(cellchat_list) <- condition_names
  
  # Merge
  cellchat_merged <- mergeCellChat(cellchat_list, add.names = condition_names)
  
  cat(sprintf("Merged %d CellChat objects\n", length(cellchat_list)))
  
  # Compare interactions
  compareInteractions(cellchat_merged, show.legend = FALSE)
  
  return(cellchat_merged)
}

#' Extract specific sender-receiver interactions
#'
#' @param cellchat CellChat object
#' @param sources Source cell types
#' @param targets Target cell types
#' @return Data frame of filtered interactions
#' @export
get_sender_receiver_interactions <- function(cellchat, sources, targets) {
  
  df.net <- subsetCommunication(cellchat,
                                sources.use = sources,
                                targets.use = targets)
  
  cat(sprintf("Found %d interactions from %s to %s\n",
              nrow(df.net),
              paste(sources, collapse = ", "),
              paste(targets, collapse = ", ")))
  
  return(df.net)
}

#' Visualize specific pathways
#'
#' @param cellchat CellChat object
#' @param pathways Vector of pathway names
#' @param plot_type Type of plot ("chord", "bubble", "heatmap")
#' @export
visualize_pathways <- function(cellchat, pathways, plot_type = "chord") {
  
  for (pathway in pathways) {
    cat(sprintf("Visualizing %s...\n", pathway))
    
    if (plot_type == "chord") {
      netVisual_chord_gene(cellchat,
                           signaling = pathway,
                           slot.name = "netP",
                           legend.pos.x = 10,
                           title.name = pathway)
    } else if (plot_type == "bubble") {
      netVisual_bubble(cellchat,
                       signaling = pathway,
                       remove.isolate = FALSE)
    } else if (plot_type == "heatmap") {
      netVisual_heatmap(cellchat,
                        signaling = pathway,
                        color.heatmap = "Reds")
    }
  }
}

#' Run pattern analysis
#'
#' @param cellchat CellChat object
#' @param pattern Pattern type ("outgoing" or "incoming")
#' @param k Number of patterns
#' @return CellChat object with patterns
#' @export
identify_patterns <- function(cellchat, pattern = "outgoing", k = NULL) {
  
  # Suggest k if not provided
  if (is.null(k)) {
    selectK(cellchat, pattern = pattern)
    k <- readline(prompt = "Enter number of patterns (k): ")
    k <- as.numeric(k)
  }
  
  cat(sprintf("Identifying %d %s patterns...\n", k, pattern))
  
  cellchat <- identifyCommunicationPatterns(cellchat,
                                            pattern = pattern,
                                            k = k,
                                            height = 9)
  
  # Visualize
  netAnalysis_river(cellchat, pattern = pattern)
  netAnalysis_dot(cellchat, pattern = pattern)
  
  return(cellchat)
}
