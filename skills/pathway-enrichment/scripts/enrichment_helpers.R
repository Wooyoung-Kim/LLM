# Helper functions for pathway enrichment analysis

library(fgsea)
library(dplyr)

#' Prepare ranked gene list from DEG results
#'
#' @param deg_df Data frame with DEG results
#' @param gene_col Column name for gene symbols
#' @param rank_col Column name for ranking (e.g., "avg_log2FC", "stat")
#' @param remove_na Remove NA values
#' @return Named numeric vector, sorted
#' @export
make_ranked_list <- function(deg_df,
                              gene_col = "gene",
                              rank_col = "avg_log2FC",
                              remove_na = TRUE) {
  
  # Extract values
  ranks <- deg_df[[rank_col]]
  names(ranks) <- deg_df[[gene_col]]
  
  # Remove NA
  if (remove_na) {
    ranks <- ranks[!is.na(ranks)]
  }
  
  # Sort decreasing
  ranks <- sort(ranks, decreasing = TRUE)
  
  return(ranks)
}

#' Run fgsea and format results
#'
#' @param gene_ranks Named numeric vector
#' @param pathways List of gene sets
#' @param minSize Minimum pathway size
#' @param maxSize Maximum pathway size
#' @return Data frame with fgsea results
#' @export
run_fgsea_analysis <- function(gene_ranks,
                                pathways,
                                minSize = 15,
                                maxSize = 500,
                                nperm = 10000) {
  
  fgsea_res <- fgsea(
    pathways = pathways,
    stats = gene_ranks,
    minSize = minSize,
    maxSize = maxSize,
    nperm = nperm
  )
  
  # Sort by significance
  fgsea_res <- fgsea_res %>%
    arrange(padj, desc(abs(NES)))
  
  return(fgsea_res)
}

#' Plot top fgsea pathways
#'
#' @param fgsea_results Results from fgsea
#' @param top_n Number of top pathways
#' @param clean_names Remove prefix and format names
#' @export
plot_fgsea_bar <- function(fgsea_results,
                            top_n = 20,
                            clean_names = TRUE) {
  
  library(ggplot2)
  
  top_pathways <- fgsea_results %>%
    filter(padj < 0.05) %>%
    top_n(top_n, wt = abs(NES))
  
  if (clean_names) {
    top_pathways <- top_pathways %>%
      mutate(
        pathway = gsub("^HALLMARK_", "", pathway),
        pathway = gsub("^GOBP_|^GOCC_|^GOMF_", "", pathway),
        pathway = gsub("_", " ", pathway),
        pathway = tolower(pathway),
        pathway = stringr::str_to_title(pathway)
      )
  }
  
  p <- ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES)) +
    geom_col(aes(fill = -log10(padj))) +
    coord_flip() +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    labs(
      x = "Pathway",
      y = "Normalized Enrichment Score",
      fill = "-log10(padj)"
    ) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  return(p)
}

#' Load MSigDB gene sets
#'
#' @param species Species ("Homo sapiens" or "Mus musculus")
#' @param category MSigDB category ("H", "C2", "C5", etc.)
#' @param subcategory Subcategory (e.g., "CP:KEGG", "GO:BP")
#' @return List of gene sets
#' @export
load_msigdb_pathways <- function(species = "Homo sapiens",
                                  category = "H",
                                  subcategory = NULL) {
  
  library(msigdbr)
  
  if (is.null(subcategory)) {
    gene_sets <- msigdbr(species = species, category = category)
  } else {
    gene_sets <- msigdbr(species = species, 
                          category = category,
                          subcategory = subcategory)
  }
  
  pathways <- split(gene_sets$gene_symbol, gene_sets$gs_name)
  
  cat(sprintf("Loaded %d pathways from %s", 
              length(pathways), category))
  if (!is.null(subcategory)) {
    cat(sprintf(" (%s)", subcategory))
  }
  cat("\n")
  
  return(pathways)
}

#' Run enrichment for up and down genes separately
#'
#' @param deg_df DEG results
#' @param fc_col Log fold change column
#' @param pval_col Adjusted p-value column
#' @param fc_cutoff Fold change cutoff
#' @param pval_cutoff P-value cutoff
#' @param pathways Gene set list
#' @return List with up and down enrichment results
#' @export
run_directional_enrichment <- function(deg_df,
                                        fc_col = "avg_log2FC",
                                        pval_col = "p_val_adj",
                                        fc_cutoff = 0.5,
                                        pval_cutoff = 0.05,
                                        pathways) {
  
  # Upregulated genes
  up_genes <- deg_df %>%
    filter(.data[[fc_col]] > fc_cutoff & 
           .data[[pval_col]] < pval_cutoff) %>%
    pull(gene)
  
  # Downregulated genes
  down_genes <- deg_df %>%
    filter(.data[[fc_col]] < -fc_cutoff & 
           .data[[pval_col]] < pval_cutoff) %>%
    pull(gene)
  
  cat(sprintf("Up: %d genes, Down: %d genes\n", 
              length(up_genes), length(down_genes)))
  
  # Run enrichment
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  up_enrich <- enricher(
    gene = up_genes,
    TERM2GENE = stack(pathways),
    pvalueCutoff = 0.05
  )
  
  down_enrich <- enricher(
    gene = down_genes,
    TERM2GENE = stack(pathways),
    pvalueCutoff = 0.05
  )
  
  return(list(
    up = up_enrich,
    down = down_enrich,
    up_genes = up_genes,
    down_genes = down_genes
  ))
}
