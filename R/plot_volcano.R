

#' @name plot_volcano
#' @title plot volcano
#' @description plot volcano plot for DESeq2 pariwise compare result
#'
#' @param res a data.drame from DESeq2 result
#' @param title plot title
#' @param baseMean_thred baseMean cutoff
#' @param fc_thred foldchange cutoff
#' @param padj_thred p.adjust cutoff
#' @param size point size
#' @param alpha point transparency
#' @param xlim two element numeric vector, restrict x axis
#' @param ylim two element numeric vector, restrict y axis
#' @param ns_resampling numeric, downsampling NS points
#' @param color three element character vector, map point color, ordered by down-regulated, ns, up-regulated
#'
#' @importFrom dplyr filter mutate if_else sample_n bind_rows %>%
#' @import ggplot2
#'
#' @return ggplot2 object
#'
#' @export
#'
plot_volcano <- function(res,
                         title="Volcano of DEGs",
                         baseMean_thred=50,
                         fc_thred=2,
                         padj_thred=0.05,
                         size=1.5,
                         alpha=0.7,
                         xlim=NULL, # c(-5, 5)
                         ylim=NULL, # c(0, 20)
                         ns_resampling=1000,
                         color=c("#0571b0", "#bababa", "#ca0020")
                         ) {
  df <- filter(res, baseMean > baseMean_thred, !is.na(entrezgene_id), !is.na(padj)) %>%
    mutate(`-log10(padj)`=-log10(padj))
  # take sig genes
  df.sig <- filter(df, padj < padj_thred, abs(log2FoldChange) > log2(fc_thred)) %>%
    mutate(direction=if_else(log2FoldChange > 0, "up", "down"))
  # sample un-sig genes,  fo reducing rendering points
  df.unsig <- filter(df, padj >= padj_thred) %>%
    sample_n(size=ns_resampling) %>%
    mutate(direction="ns")
  # merge sig and un-sig
  df.merge <- bind_rows(df.sig, df.unsig)

  # set lims
  df.merge$shape <- "inrange"
  if(!is.null(xlim)) {
    df.merge <- df.merge %>%
      mutate(shape=if_else(log2FoldChange > xlim[2] | log2FoldChange < xlim[1], "outrange", "inrange")) %>%
      mutate(log2FoldChange=if_else(log2FoldChange > xlim[2], xlim[2],
                                    if_else(log2FoldChange < xlim[1], xlim[1], log2FoldChange)))
    }
  if(!is.null(ylim)) {
    df.merge <- df.merge %>%
      mutate(shape=if_else(`-log10(padj)` > ylim[2] | `-log10(padj)` < ylim[1], "outrange", shape)) %>%
      mutate(`-log10(padj)`=if_else(`-log10(padj)` > ylim[2], ylim[2],
                                    if_else(`-log10(padj)` < ylim[1], ylim[1], `-log10(padj)`)))
  }

  # plot
  plot <- ggplot(df.merge, aes(log2FoldChange, `-log10(padj)`)) +
    geom_point(aes(shape=shape, fill=direction), size=size, alpha=alpha, stroke=0) +
    scale_shape_manual(values=c(inrange=21, outrange=24)) +
    scale_fill_manual(values=c(down=color[1], ns=color[2], up=color[3])) +
    geom_vline(xintercept=c(-log2(fc_thred), log2(fc_thred)), linetype="dashed") +
    geom_hline(yintercept=-log10(padj_thred), linetype="dashed") +
    labs(x="log2 (FoldChange)", y="-log10 (p.adj)", title=title, subtitle=str_glue(
           "padj_threshold: {padj_thred}, FoldChange_threshold: {fc_thred}, baseMean_threshold: {baseMean_thred}")) +
    theme_linedraw()

  return(plot)
}



# not run
if(F) {
  res <- readxl::read_xlsx("/Users/hh/Desktop/2020-01-13_smartSeq_lung_EC_PC_HH/statistics/2020-01-15_wald_TEST_lps.pc_pbs.pc.xlsx")

  p <- plot_volcano(res)
  cowplot::save_plot("~/Desktop/test.pdf", p, base_height = 3)

  plot_volcano(res, xlim = c(-5, 5), ylim = c(0, 10))

}
