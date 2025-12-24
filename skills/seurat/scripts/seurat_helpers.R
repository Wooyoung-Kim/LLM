# Quick QC and filtering helper functions for Seurat

library(Seurat)
library(ggplot2)

#' Comprehensive QC metrics calculation
#'
#' @param seurat_obj Seurat object
#' @param species Species ("human" or "mouse")
#' @return Seurat object with QC metrics added
#' @export
add_qc_metrics <- function(seurat_obj, species = "human") {
  # Mitochondrial genes
  mt_pattern <- ifelse(species == "human", "^MT-", "^mt-")
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
  
  # Ribosomal genes
  ribo_pattern <- "^RP[SL]"
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = ribo_pattern)
  
  # Hemoglobin genes (for blood samples)
  if (species == "human") {
    seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^HB[^(P)]")
  }
  
  # Complexity (genes per UMI)
  seurat_obj[["log10GenesPerUMI"]] <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
  
  return(seurat_obj)
}

#' Plot comprehensive QC metrics
#'
#' @param seurat_obj Seurat object with QC metrics
#' @export
plot_qc_metrics <- function(seurat_obj) {
  library(patchwork)
  
  # Violin plots
  p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, pt.size = 0.1)
  
  # Scatter plots
  p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  p4 <- FeatureScatter(seurat_obj, feature1 = "nFeature_RNA", feature2 = "percent.mt")
  
  # Combine
  combined <- p1 / (p2 | p3 | p4)
  
  print(combined)
  
  # Print summary statistics
  cat("\nQC Metrics Summary:\n")
  cat("nFeature_RNA (genes):\n")
  print(summary(seurat_obj$nFeature_RNA))
  cat("\nnCount_RNA (UMIs):\n")
  print(summary(seurat_obj$nCount_RNA))
  cat("\npercent.mt:\n")
  print(summary(seurat_obj$percent.mt))
}

#' Suggest QC thresholds based on data distribution
#'
#' @param seurat_obj Seurat object with QC metrics
#' @export
suggest_qc_thresholds <- function(seurat_obj) {
  # Calculate MAD-based thresholds
  median_features <- median(seurat_obj$nFeature_RNA)
  mad_features <- mad(seurat_obj$nFeature_RNA)
  
  median_counts <- median(seurat_obj$nCount_RNA)
  mad_counts <- mad(seurat_obj$nCount_RNA)
  
  # Suggest thresholds
  cat("\nSuggested QC Thresholds:\n")
  cat(sprintf("nFeature_RNA: %d - %d\n", 
              max(200, median_features - 3*mad_features),
              median_features + 3 * mad_features))
  cat(sprintf("nCount_RNA: %d - %d\n",
              max(500, median_counts - 3*mad_counts),
              median_counts + 3*mad_counts))
  cat(sprintf("percent.mt: < %.1f\n", quantile(seurat_obj$percent.mt, 0.95)))
  
  cat("\nExample filter command:\n")
  cat(sprintf("seurat_obj <- subset(seurat_obj, subset = \n"))
  cat(sprintf("  nFeature_RNA > %d & nFeature_RNA < %d & \n",
              max(200, median_features - 3*mad_features),
              median_features + 3*mad_features))
  cat(sprintf("  percent.mt < %.1f)\n", quantile(seurat_obj$percent.mt, 0.95)))
}

#' Standard Seurat workflow
#'
#' @param seurat_obj Seurat object after QC filtering
#' @param n_pcs Number of PCs to use (default: 30)
#' @param resolution Clustering resolution (default: 0.5)
#' @param use_sct Use SCTransform instead of standard normalization (default: FALSE)
#' @return Processed Seurat object
#' @export
run_standard_workflow <- function(seurat_obj, n_pcs = 30, resolution = 0.5, use_sct = FALSE) {
  library(Seurat)
  
  if (use_sct) {
    # SCTransform workflow
    seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
  } else {
    # Standard workflow
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
    seurat_obj <- ScaleData(seurat_obj)
  }
  
  # PCA
  seurat_obj <- RunPCA(seurat_obj, npcs = n_pcs, verbose = FALSE)
  
  # Clustering
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_pcs)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  
  # UMAP
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_pcs)
  
  return(seurat_obj)
}

#' Export top markers to CSV
#'
#' @param markers_df Data frame from FindAllMarkers
#' @param top_n Number of top markers per cluster
#' @param output_file Output file path
#' @export
export_top_markers <- function(markers_df, top_n = 10, output_file = "top_markers.csv") {
  library(dplyr)
  
  top_markers <- markers_df %>%
    group_by(cluster) %>%
    top_n(n = top_n, wt = avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC))
  
  write.csv(top_markers, output_file, row.names = FALSE)
  cat(sprintf("Top %d markers per cluster exported to %s\n", top_n, output_file))
  
  return(top_markers)
}
