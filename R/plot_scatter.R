
#' @name plot_scatter
#' @title plot_scatter
#' @description plot scatter plot for pairwise comparison
#'
#' @param expr data.frame of expressing values,
#'        rownames are gene names and colnames are sample/group names, all columns are expressing values,
#'        for better transform estimation, expr should be include all sample instead of only inculde column x and y, if possible
#' @param x x-axis name, corelated to one column name of expr
#' @param y y-axis name, corelated to one column name of expr
#' @param dge_genes character vector of differential changed genes, corelated to rownames of expr
#' @param transform_method matrix transform method, one of vst, rlog or none
#' @param filter_low_values numeric, filter rowMeans low value rows
#' @param size point size
#' @param alpha point transparency
#' @param title plot title
#' @param xlim two element numeric vector, restrict x axis, default: NULL
#' @param ylim two element numeric vector, restrict y axis, default: NULL
#' @param ns_resampling numeric, downsampling points
#' @param color character vector, map point color
#' @param dge_genes_color character vector, map dge_genes point color
#'
#' @importFrom dplyr %>% mutate filter if_else
#' @importFrom tibble has_rownames rownames_to_column as_tibble
#' @importFrom DESeq2 vst rlog
#' @importFrom rlang sym !! :=
#' @import ggplot2
#'
#' @return ggplot2 object
#'
#' @export
#'
plot_scatter <- function(expr,
                         x,
                         y,
                         dge_genes=NULL,
                         transform_method="vst", # vst, rlog, none
                         filter_low_values=1,
                         size=1.5,
                         alpha=0.7,
                         title="Pairwsie scatter point",
                         xlim=NULL,
                         ylim=NULL,
                         ns_resampling=1000,
                         color="grey50",
                         dge_genes_color="red2"
                         ) {

  if(!has_rownames(expr)) stop("expr has no valid rownames")

  if(!all(x %in% colnames(expr), y %in% colnames(expr))) stop("x or y was not in any of expr column")

  expr <- as.matrix(expr)

  if(filter_low_values > 0) remove_genes <- rownames(expr)[which(rowMeans(expr) <= filter_low_values)]

  # round for transformation
  expr <- round(expr)

  # transform
  switch(transform_method,
         vst={
           print("vst transformation")
           expr <- suppressMessages(vst(expr, blind = TRUE))
         },
         rlog={
           print("rlog transformation")
           names <- rownames(expr)
           expr <- suppressMessages(rlog(expr, blind = TRUE))
           rownames(expr) <- names
         },
         expr
  )

  # resample
  if(is.null(dge_genes)) {
    # remove low filter values and keep all dge genes
    keep_genes <- setdiff(rownames(expr), remove_genes) %>%
      sample(ns_resampling) %>%
      unique()
    # convert to df
    df <- expr[, c(x, y)] %>%
      as_tibble(rownames="name") %>%
      filter(name %in% keep_genes) %>%
      mutate(sig="ns")
  } else {
    # remove low filter values and keep all dge genes
    keep_genes <- setdiff(rownames(expr), remove_genes) %>%
      setdiff(dge_genes) %>%
      sample(ns_resampling) %>%
      c(dge_genes) %>%
      unique()
    # convert to df
    df <- expr[, c(x, y)] %>%
      as_tibble(rownames="name") %>%
      filter(name %in% keep_genes) %>%
      mutate(sig="ns") %>%   # add dge info
      mutate(sig=if_else(name %in% dge_genes, "dge", sig))
  }

  # set lims
  df$shape <- "inrange"
  if(!is.null(xlim)) {
    df <- df %>%
      mutate(shape=if_else(!!sym(x) > xlim[2] | !!sym(x) < xlim[1], "outrange", "inrange")) %>%
      mutate(!!sym(x):=if_else(!!sym(x) > xlim[2], xlim[2],
                              if_else(!!sym(x) < xlim[1], xlim[1], !!sym(x))))
  }
  if(!is.null(ylim)) {
    df <- df %>%
      mutate(shape=if_else(!!sym(y) > ylim[2] | !!sym(y) < ylim[1], "outrange", shape)) %>%
      mutate(!!sym(y):=if_else(!!sym(y) > ylim[2], ylim[2],
                              if_else(!!sym(y) < ylim[1], ylim[1], !!sym(y))))
  }

  # plot
  plot <- ggplot(df, aes(!!sym(x), !!sym(y))) +
    geom_point(aes(shape=shape, fill=sig), size=size, alpha=alpha, stroke=0) +
    scale_shape_manual(values=c(inrange=21, outrange=24)) +
    scale_fill_manual(values=c(ns=color, dge=dge_genes_color)) +
    ggtitle(title) +
    theme_linedraw()


  return(plot)
}




# not run
if(F) {
  mat <- readxl::read_xlsx("/Users/hh/Desktop/2020-01-13_smartSeq_lung_EC_PC_HH/gene_expression.xlsx") %>%
    .[, c(2, 20:31)] %>%
    filter(!is.na(external_gene_name)) %>%
    distinct(external_gene_name, .keep_all = T) %>%
    column_to_rownames("external_gene_name")
  mat_mean <- data.frame(pbs=rowMeans(mat[,1:6]), lps=rowMeans(mat[, 7:12]))
  rownames(mat_mean) <- rownames(mat)

  res <- readxl::read_xlsx("/Users/hh/Desktop/2020-01-13_smartSeq_lung_EC_PC_HH/statistics/2020-01-15_wald_TEST_lps.pc_pbs.pc.xlsx")
  genes <- filter(res, abs(log2FoldChange) > 2, padj < 0.01) %>%
    pull(external_gene_name)

  plot_scatter(mat, x="pbs_pc#3", y="lps_pc#3", dge_genes=genes, transform_method = "rlog", filter_low_values = 10, size=2)
  plot_scatter(mat_mean, x="pbs", y="lps", dge_genes=genes, transform_method = "vst",
               dge_genes_color = "darkblue", size=2, xlim=c(0, 10), ylim=c(0, 10))
}
