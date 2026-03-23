library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  counts <- read.delim(counts_csv, row.names = 1)
  meta   <- read.csv(metafile_csv)
  meta_sub   <- meta[meta$timepoint %in% selected_times, ]
  counts_sub <- counts[, meta_sub$samplename]
  meta_sub$timepoint <- factor(meta_sub$timepoint, levels = c("vP0", setdiff(selected_times, "vP0")))
  se <- SummarizedExperiment(
    assays  = list(counts = as.matrix(counts_sub)),
    colData = meta_sub
  )
  
  return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  dds <- DESeqDataSet(se, design)
  dds <- DESeq(dds)
  res <- results(dds)
  res_df <- as.data.frame(res)
  
  return(list(
    dds = dds,
    res = res_df
  ))
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  result <- deseq2_res %>%
    rownames_to_column(var = "genes") %>%
    as_tibble() %>%
    mutate(volc_plot_status = case_when(
      padj < padj_threshold & log2FoldChange > 0 ~ "UP",
      padj < padj_threshold & log2FoldChange < 0 ~ "DOWN",
      TRUE ~ "NS"
    ))
  
  return(result)
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  plot <- ggplot(labeled_results, aes(x = pvalue)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    labs(
      title = "Histogram of Unadjusted P-values",
      x = "P-value",
      y = "Count"
    ) +
    theme_minimal()
  
  return(plot)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  sig_results <- labeled_results %>%
    filter(padj < padj_threshold)
  
  plot <- ggplot(sig_results, aes(x = log2FoldChange)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    labs(
      title = "Histogram of Log2 Fold Change (Significant Genes)",
      x = "Log2 Fold Change",
      y = "Count"
    ) +
    theme_minimal()
  
  return(plot)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  top_genes <- labeled_results %>%
    arrange(padj) %>%
    head(num_genes) %>%
    pull(genes)
  norm_counts <- counts(dds_obj, normalized = TRUE)
  
  plot_data <- norm_counts[top_genes, ] %>%
    as.data.frame() %>%
    rownames_to_column(var = "genes") %>%
    pivot_longer(-genes, names_to = "sample", values_to = "norm_counts")
  
  plot <- ggplot(plot_data, aes(x = genes, y = norm_counts, color = sample)) +
    geom_point(size = 2) +
    labs(
      title = paste("Normalized Counts for Top", num_genes, "Genes by Adjusted P-value"),
      x = "Gene",
      y = "Normalized Counts"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(plot)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  plot <- ggplot(labeled_results, aes(x = log2FoldChange, y = -log10(padj), color = volc_plot_status)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NS" = "grey")) +
    labs(
      title = "Volcano Plot",
      x = "Log2 Fold Change",
      y = "-log10(Adjusted P-value)",
      color = "Status"
    ) +
    theme_minimal()
  
  return(plot)
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  id2gene <- read.delim(id2gene_path, header = FALSE, col.names = c("ensembl", "symbol"))
  
  rnk_list <- labeled_results %>%
    filter(!is.na(log2FoldChange)) %>%
    left_join(id2gene, by = c("genes" = "ensembl")) %>%
    mutate(symbol = ifelse(is.na(symbol), genes, symbol)) %>%
    group_by(symbol) %>%
    slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(log2FoldChange)) %>%
    { setNames(.$log2FoldChange, .$symbol) }
  
  return(rnk_list)
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  pathways <- gmtPathways(gmt_file_path)
  
  fgsea_res <- fgsea(
    pathways = pathways,
    stats    = rnk_list,
    minSize  = min_size,
    maxSize  = max_size
  )
  
  return(as_tibble(fgsea_res))
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  top_pos <- fgsea_results %>%
    filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    head(num_paths)
  
  top_neg <- fgsea_results %>%
    filter(NES < 0) %>%
    arrange(NES) %>%
    head(num_paths)
  
  plot_data <- bind_rows(top_pos, top_neg) %>%
    arrange(NES) %>%
    mutate(pathway = factor(pathway, levels = pathway))
  
  plot <- ggplot(plot_data, aes(x = pathway, y = NES, fill = NES > 0)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                      labels = c("TRUE" = "Positive", "FALSE" = "Negative"),
                      name = "Direction") +
    coord_flip() +
    labs(
      title = paste("Top", num_paths, "Positive and Negative NES Pathways"),
      x = "Pathway",
      y = "NES"
    ) +
    theme_minimal()
  
  return(plot)
}

