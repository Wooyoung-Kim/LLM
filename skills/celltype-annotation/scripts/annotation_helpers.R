# Helper functions for cell type annotation

library(Seurat)
library(dplyr)
library(tibble)

#' Load ScType functions from GitHub
#'
#' @export
load_sctype <- function() {
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  message("ScType functions loaded successfully")
}

#' Run ScType annotation on Seurat object
#'
#' @param seurat_obj Seurat object with clustered data
#' @param tissue Tissue type for annotation
#' @param cluster_col Column name for clusters (default: "seurat_clusters")
#' @param assay Assay to use (default: "RNA")
#' @param db_path Path to ScType database (default: GitHub URL)
#' @param confidence_threshold Score threshold as fraction of ncells (default: 0.25)
#' @return Seurat object with sctype annotations
#' @export
run_sctype_annotation <- function(seurat_obj,
                                    tissue = "Immune system",
                                    cluster_col = "seurat_clusters",
                                    assay = "RNA",
                                    db_path = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",
                                    confidence_threshold = 0.25) {
  
  # Load ScType functions
  load_sctype()
  
  # Prepare gene sets
  gs_list <- gene_sets_prepare(db_path, tissue)
  
  # Check Seurat version
  seurat_package_v5 <- isFALSE('counts' %in% names(attributes(seurat_obj[[assay]])))
  
  # Extract scaled data
  scRNAseqData_scaled <- if (seurat_package_v5) {
    as.matrix(seurat_obj[[assay]]$scale.data)
  } else {
    as.matrix(seurat_obj[[assay]]@scale.data)
  }
  
  # Run ScType
  es.max <- sctype_score(
    scRNAseqData = scRNAseqData_scaled,
    scaled = TRUE,
    gs = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
  )
  
  # Aggregate by cluster
  cL_results <- do.call("rbind", lapply(unique(seurat_obj@meta.data[[cluster_col]]), function(cl) {
    cells_in_cluster <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[cluster_col]] == cl, ])
    es.max.cl <- sort(rowSums(es.max[, cells_in_cluster]), decreasing = TRUE)
    head(data.frame(
      cluster = cl,
      type = names(es.max.cl),
      scores = es.max.cl,
      ncells = length(cells_in_cluster)
    ), 10)
  }))
  
  # Get top prediction per cluster
  sctype_scores <- cL_results %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = scores)
  
  # Set low-confidence to Unknown
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells * confidence_threshold] <- "Unknown"
  
  # Add to Seurat
  seurat_obj$sctype <- ""
  for(j in unique(sctype_scores$cluster)) {
    cl_type <- sctype_scores[sctype_scores$cluster == j, ]
    seurat_obj$sctype[seurat_obj@meta.data[[cluster_col]] == j] <- as.character(cl_type$type[1])
  }
  
  # Add confidence scores
  seurat_obj$sctype_score <- NA
  for(j in unique(sctype_scores$cluster)) {
    cl_score <- sctype_scores[sctype_scores$cluster == j, ]
    seurat_obj$sctype_score[seurat_obj@meta.data[[cluster_col]] == j] <- as.numeric(as.character(cl_score$scores[1]))
  }
  
  message(sprintf("ScType annotation complete. %d clusters annotated.", nrow(sctype_scores)))
  
  return(list(
    seurat_obj = seurat_obj,
    scores = sctype_scores
  ))
}

#' Run Scibet annotation on Seurat object
#'
#' @param seurat_obj Seurat object
#' @param reference_path Path to Scibet reference model CSV
#' @param cluster_col Column name for clusters
#' @param n_genes Number of genes for feature selection
#' @return Seurat object with scibet annotations
#' @export
run_scibet_annotation <- function(seurat_obj,
                                    reference_path,
                                    cluster_col = "seurat_clusters",
                                    n_genes = 50) {
  
  library(scibet)
  
  # Load reference
  model <- readr::read_csv(reference_path)
  model <- pro.core(model)
  
  # Join layers (Seurat v5)
  seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
  
  # Prepare query
  query <- seurat_obj@assays$RNA$data %>%
    t() %>%
    as.data.frame()
  
  query_meta <- seurat_obj@meta.data %>%
    select(all_of(cluster_col))
  colnames(query_meta) <- "label"
  query <- cbind(query, query_meta)
  
  # Select genes
  etest_gene <- SelectGene(query, k = n_genes)
  
  # Predict
  query_labels <- query$label
  query <- query[, -ncol(query)]
  
  prd <- LoadModel(model)
  predicted_labels <- prd(query)
  
  # Add to Seurat
  seurat_obj$scibet <- predicted_labels
  
  # Get consensus per cluster
  cluster_celltype_map <- seurat_obj@meta.data %>%
    group_by(across(all_of(cluster_col)), scibet) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(across(all_of(cluster_col)), desc(count)) %>%
    group_by(across(all_of(cluster_col))) %>%
    slice(1) %>%
    ungroup()
  
  names(cluster_celltype_map)[names(cluster_celltype_map) == "scibet"] <- "scibet_consensus"
  
  seurat_meta <- seurat_obj@meta.data %>%
    rownames_to_column("cells") %>%
    left_join(cluster_celltype_map, by = cluster_col) %>%
    column_to_rownames("cells")
  
  seurat_obj@meta.data <- seurat_meta
  
  message(sprintf("Scibet annotation complete using %d genes.", n_genes))
  
  return(seurat_obj)
}

#' Run scMayoMap annotation
#'
#' @param markers Data frame with marker genes (from FindAllMarkers)
#' @param tissue Tissue type
#' @param seurat_obj Seurat object (optional, for adding annotations)
#' @param cluster_col Cluster column name
#' @return List with annotations and optionally updated Seurat object
#' @export
run_scmayomap_annotation <- function(markers,
                                      tissue,
                                      seurat_obj = NULL,
                                      cluster_col = "seurat_clusters") {
  
  library(scMayoMap)
  
  # Run scMayoMap
  obj <- scMayoMap(data = markers, tissue = tissue)
  mayo_markers <- obj$markers
  
  # Get top prediction per cluster
  top_markers <- mayo_markers %>%
    group_by(cluster) %>%
    slice_max(order_by = score, n = 1) %>%
    ungroup() %>%
    select(cluster, celltype, genes)
  
  colnames(top_markers) <- c(cluster_col, "mayo", "mayo_gene")
  
  # Add to Seurat if provided
  if (!is.null(seurat_obj)) {
    seurat_meta <- seurat_obj@meta.data %>%
      rownames_to_column("cells") %>%
      left_join(top_markers, by = cluster_col) %>%
      column_to_rownames("cells")
    
    seurat_obj@meta.data <- seurat_meta
    
    message(sprintf("scMayoMap annotation complete for tissue: %s", tissue))
    
    return(list(
      seurat_obj = seurat_obj,
      annotations = top_markers
    ))
  }
  
  return(list(annotations = top_markers))
}

#' Compare multiple annotation methods
#'
#' @param seurat_obj Seurat object with multiple annotation columns
#' @param methods Vector of column names to compare
#' @param cluster_col Cluster column
#' @return Comparison data frame
#' @export
compare_annotations <- function(seurat_obj,
                                  methods = c("sctype", "scibet_consensus", "mayo"),
                                  cluster_col = "seurat_clusters") {
  
  comparison <- seurat_obj@meta.data %>%
    group_by(across(all_of(cluster_col))) %>%
    summarise(
      across(all_of(methods), ~names(which.max(table(.)))),
      n_cells = n(),
      .groups = "drop"
    )
  
  # Calculate agreement
  if (length(methods) >= 2) {
    comparison$agreement_score <- apply(comparison[, methods], 1, function(x) {
      sum(x == x[1]) / length(x)
    })
  }
  
  return(comparison)
}

#' Create consensus annotation
#'
#' @param comparison_df Output from compare_annotations
#' @param methods Methods to consider
#' @param min_agreement Minimum agreement score (0-1)
#' @return Data frame with consensus annotations
#' @export
create_consensus <- function(comparison_df,
                               methods = c("sctype", "scibet_consensus", "mayo"),
                               min_agreement = 0.67) {
  
  consensus <- comparison_df %>%
    mutate(
      consensus_celltype = case_when(
        agreement_score >= min_agreement ~ .data[[methods[1]]],
        TRUE ~ "Manual_Review_Needed"
      ),
      consensus_confidence = case_when(
        agreement_score == 1.0 ~ "High",
        agreement_score >= 0.67 ~ "Medium",
        TRUE ~ "Low"
      )
    )
  
  return(consensus)
}

#' Visualize annotation comparison
#'
#' @param seurat_obj Seurat object
#' @param methods Methods to compare
#' @param reduction Reduction to use (default: "umap")
#' @export
plot_annotation_comparison <- function(seurat_obj,
                                        methods = c("sctype", "scibet_consensus", "mayo"),
                                        reduction = "umap") {
  
  library(patchwork)
  
  plots <- lapply(methods, function(method) {
    DimPlot(seurat_obj, reduction = reduction, label = TRUE, group.by = method) +
      ggtitle(method) +
      NoLegend()
  })
  
  wrap_plots(plots, ncol = length(methods))
}
