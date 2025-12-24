# Helper functions for pseudobulk DESeq2 analysis

library(Seurat)
library(DESeq2)
library(dplyr)
library(tibble)

#' Aggregate Seurat object to pseudobulk
#'
#' @param seurat_obj Seurat object
#' @param group_by Vector of metadata columns to group by
#' @param assay Assay to use
#' @return List with counts matrix and metadata
#' @export
aggregate_to_pseudobulk <- function(seurat_obj,
                                      group_by = c("condition", "replicate", "celltype"),
                                      assay = "RNA") {
  
  # Aggregate
  pseudobulk <- AggregateExpression(
    seurat_obj,
    assays = assay,
    slot = "counts",
    return.seurat = TRUE,
    group.by = group_by,
    verbose = TRUE
  )
  
  # Extract counts
  counts <- GetAssayData(pseudobulk, layer = "counts")
  
  # Verify integer counts
  if (!all(counts@x %% 1 == 0)) {
    warning("Counts are not integers! Check that you used slot='counts'")
  }
  
  # Create metadata
  samples <- colnames(counts)
  
  metadata <- tibble(sample = samples) %>%
    tidyr::separate(sample, into = group_by, sep = "_", remove = FALSE)
  
  cat(sprintf("Created %d pseudobulk samples\n", ncol(counts)))
  cat(sprintf("Groups: %s\n", paste(group_by, collapse = ", ")))
  
  return(list(
    counts = counts,
    metadata = metadata,
    pseudobulk = pseudobulk
  ))
}

#' Run DESeq2 for one cell type
#'
#' @param counts Count matrix
#' @param metadata Sample metadata
#' @param celltype Cell type to analyze
#' @param design_formula Formula for design (e.g., "~ condition")
#' @param contrast Contrast to test (e.g., c("condition", "treated", "control"))
#' @param shrink Use LFC shrinkage
#' @param shrink_type Type of shrinkage ("ashr", "apeglm", "normal")
#' @return Data frame with results
#' @export
run_deseq2_celltype <- function(counts,
                                 metadata,
                                 celltype,
                                 design_formula = "~ condition",
                                 contrast = c("condition", "treated", "control"),
                                 shrink = TRUE,
                                 shrink_type = "ashr") {
  
  # Filter to cell type
  keep_samples <- metadata$celltype == celltype
  if (sum(keep_samples) < 3) {
    warning(sprintf("Not enough samples for %s (need >= 3)", celltype))
    return(NULL)
  }
  
  counts_sub <- counts[, keep_samples, drop = FALSE]
  meta_sub <- metadata[keep_samples, , drop = FALSE]
  
  # Filter low-count genes
  keep_genes <- rowSums(counts_sub >= 1) >= 2
  counts_sub <- counts_sub[keep_genes, , drop = FALSE]
  
  if (nrow(counts_sub) < 10) {
    warning(sprintf("Too few genes for %s", celltype))
    return(NULL)
  }
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData = meta_sub,
    design = as.formula(design_formula)
  )
  
  # Run DESeq2
  dds <- DESeq(dds, fitType = "local")
  
  # Get results
  res <- results(dds, contrast = contrast)
  
  # Shrink LFC
  if (shrink) {
    if (shrink_type == "ashr" && requireNamespace("ashr", quietly = TRUE)) {
      res <- lfcShrink(dds, contrast = contrast, res = res, type = "ashr")
    } else if (shrink_type == "apeglm" && requireNamespace("apeglm", quietly = TRUE)) {
      # Need coefficient name for apeglm
      coef_name <- grep(paste0("^", contrast[1]), resultsNames(dds), value = TRUE)[1]
      res <- lfcShrink(dds, coef = coef_name, type = "apeglm")
    } else {
      res <- lfcShrink(dds, contrast = contrast, res = res, type = "normal")
    }
  }
  
  # Convert to data frame
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(
      celltype = celltype,
      contrast = paste(contrast[2], "vs", contrast[3])
    ) %>%
    arrange(padj)
  
  return(res_df)
}

#' Run DESeq2 for all cell types
#'
#' @param counts Count matrix
#' @param metadata Sample metadata
#' @param design_formula Design formula
#' @param contrast Contrast to test
#' @param parallel Run in parallel
#' @return Combined data frame with all results
#' @export
run_deseq2_all_celltypes <- function(counts,
                                      metadata,
                                      design_formula = "~ condition",
                                      contrast = c("condition", "treated", "control"),
                                      parallel = FALSE) {
  
  celltypes <- unique(metadata$celltype)
  
  cat(sprintf("Running DESeq2 for %d cell types...\n", length(celltypes)))
  
  if (parallel) {
    library(parallel)
    library(BiocParallel)
    register(MulticoreParam(workers = 4))
    
    results_list <- bplapply(celltypes, function(ct) {
      run_deseq2_celltype(counts, metadata, ct, design_formula, contrast)
    })
  } else {
    results_list <- lapply(celltypes, function(ct) {
      cat(sprintf("  Processing %s...\n", ct))
      run_deseq2_celltype(counts, metadata, ct, design_formula, contrast)
    })
  }
  
  # Remove NULL results and combine
  results_list <- results_list[!sapply(results_list, is.null)]
  combined_results <- bind_rows(results_list)
  
  cat(sprintf("Complete! Found %d total DEGs (padj < 0.05)\n",
              sum(combined_results$padj < 0.05, na.rm = TRUE)))
  
  return(combined_results)
}

#' Get significant DEGs
#'
#' @param deseq_results DESeq2 results data frame
#' @param padj_cutoff Adjusted p-value cutoff
#' @param lfc_cutoff Log fold change cutoff (absolute value)
#' @return Filtered data frame
#' @export
get_sig_degs <- function(deseq_results,
                          padj_cutoff = 0.05,
                          lfc_cutoff = 0.5) {
  
  sig_degs <- deseq_results %>%
    filter(padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff)
  
  cat(sprintf("Significant DEGs: %d\n", nrow(sig_degs)))
  cat(sprintf("  Upregulated: %d\n", sum(sig_degs$log2FoldChange > 0)))
  cat(sprintf("  Downregulated: %d\n", sum(sig_degs$log2FoldChange < 0)))
  
  return(sig_degs)
}

#' Export DESeq2 results
#'
#' @param deseq_results DESeq2 results
#' @param output_prefix Output file prefix
#' @export
export_deseq2_results <- function(deseq_results, output_prefix = "deseq2") {
  
  # All results
  write.csv(deseq_results, paste0(output_prefix, "_all_results.csv"), 
            row.names = FALSE)
  
  # Significant only
  sig <- get_sig_degs(deseq_results)
  write.csv(sig, paste0(output_prefix, "_significant.csv"), 
            row.names = FALSE)
  
  # Summary by cell type
  summary_df <- deseq_results %>%
    group_by(celltype, contrast) %>%
    summarise(
      total_genes = n(),
      sig_genes = sum(padj < 0.05, na.rm = TRUE),
      up_genes = sum(padj < 0.05 & log2FoldChange > 0.5, na.rm = TRUE),
      down_genes = sum(padj < 0.05 & log2FoldChange < -0.5, na.rm = TRUE),
      .groups = "drop"
    )
  
  write.csv(summary_df, paste0(output_prefix, "_summary.csv"), 
            row.names = FALSE)
  
  cat(sprintf("Exported results to %s_*.csv\n", output_prefix))
}
