
#' @name plot_enrich_ora
#' @title plot_enrich_ora
#' @description dotplot for enrich result from enrich_ORA
#'
#' @param enrichment data.frame, from ORA enrichment,
#'        result can be merged from multi enrichment, grouping info can be viewed by facet
#' @param show_category integer, show enriched category number
#' @param axis_x character, plot axis x by  which column,
#'        one of -log10(p.adjust), enrichFactor, Count
#' @param plot_type character, one of dot, bar
#' @param filter_by character, filter show_category by which column,
#'        one of -log10(p.adjust), enrichFactor, Count
#' @param color_by character, map point color, one of -log10(p.adjust), enrichFactor, Count
#' @param color_range two element character vector, mapping color range, from low to high
#' @param size_by character, map point size, one of -log10(p.adjust), enrichFactor, Count
#'        only works if plot_type is dot
#' @param size_range two element numeric vector, mapping point size, from low to high
#'        only works if plot_type is dot
#' @param text_len_limit integer, wrap y-axis text length
#' @param facet_by character, facet plot by which column in enrichment
#' @param facet_scales character, facet scales, one of free_y, free_x, "free, fixed,
#'        only works if facet_by is not NULL
#'
#' @import ggplot2
#' @importFrom dplyr arrange mutate slice vars group_by
#' @importFrom stringr str_wrap
#' @importFrom rlang sym !!
#'
#' @return ggplot2 object
#'
#' @export
#'
plot_enrich_ora <- function(enrichment,
                            show_category=10,
                            plot_type="dot",
                            axis_x="enrichFactor",
                            filter_by="-log10(p.adjust)",
                            color_by="-log10(p.adjust)",
                            color_range=c("darkblue", "red"),
                            size_by="Count",
                            size_range=c(2, 8),
                            text_len_limit=50,
                            facet_by=NULL,
                            facet_scales="free_y"
) {
  if(!is.data.frame(enrichment))
    stop("enrichment must be a tibble/data.frame")
  
  if(!(plot_type %in% c("dot", "bar")))
    stop("plot_type must be one of dor, bar")
  
  if(!(axis_x %in% c("-log10(p.adjust)", "enrichFactor", "Count")))
    stop("axis_x must be one of -log10(p.adjust), enrichFactor, Count")
  
  if(!(filter_by %in% c("-log10(p.adjust)", "enrichFactor", "Count")))
    stop("filter_by must be one of -log10(p.adjust), enrichFactor, Count")
  
  if(!(size_by %in% c("-log10(p.adjust)", "enrichFactor", "Count")))
    stop("size_by must be one of -log10(p.adjust), enrichFactor, Count")
  
  if(!(color_by %in% c("-log10(p.adjust)", "enrichFactor", "Count")))
    stop("color_by must be one of -log10(p.adjust), enrichFactor, Count")
  
  if(!is.null(facet_by)) {
    if(!(facet_scales %in% c("free", "free_x", "free_y", "fixed")))
      stop("facet_scales must be one of free, free_x, free_y, fixed")
  }
  
  # group facet
  if(!is.null(facet_by)) {
    data <- group_by(enrichment, !!sym(facet_by)) %>%
      mutate(`-log10(p.adjust)`=-log10(p.adjust)) %>%
      arrange(desc(!!sym(filter_by))) %>%
      dplyr::slice(seq_len(show_category)) %>%
      arrange(!!sym(axis_x)) %>%
      mutate(Description=factor(Description, levels=unique(.$Description)))
  } else {
    data <- mutate(enrichment, `-log10(p.adjust)`=-log10(p.adjust)) %>%
      arrange(desc(!!sym(filter_by))) %>%
      dplyr::slice(seq_len(show_category)) %>%
      arrange(!!sym(axis_x)) %>%
      mutate(Description=factor(Description, levels=unique(.$Description)))
  }
  
  if(plot_type=="dot") {
    plot <- ggplot(data, aes(!!sym(axis_x), Description)) +
      geom_segment(aes(yend=Description), xend=0, size=0.7, color="grey50", alpha=0.7) +
      geom_point(aes(size=!!sym(size_by), color=!!sym(color_by))) +
      scale_size_continuous(range=size_range) +
      scale_color_gradient(low=color_range[1], high=color_range[2]) +
      scale_y_discrete(labels=function(x) str_wrap(x, width=text_len_limit)) +
      theme_cowplot() +
      ylab("")
  } else {
    plot <- ggplot(data, aes(!!sym(axis_x), Description)) +
      geom_bar(aes(fill=!!sym(color_by)), stat="identity", color=NA, width=0.75) +
      scale_size_continuous(range=size_range) +
      scale_fill_gradient(low=color_range[1], high=color_range[2]) +
      scale_y_discrete(labels=function(x) str_wrap(x, width=text_len_limit)) +
      theme_cowplot() +
      ylab("")
  }
  
  # facet plot
  if(!is.null(facet_by)) {
    plot <- plot + facet_wrap(vars(!!sym(facet_by)), scales=facet_scales)
  }
  
  # return
  return(plot)
}



# for test
if(F) {
  res <- readxl::read_xlsx("/Users/hh/Desktop/2020-01-13_smartSeq_lung_EC_PC_HH/statistics/2020-01-15_wald_TEST_lps.pc_pbs.pc.xlsx", sheet=2)
  res <- arrange(res, desc(stat)) %>%
    filter(!is.na(entrezgene_id)) %>%
    distinct(entrezgene_id, .keep_all=T)
  
  ora <- enrich_ORA(res$external_gene_name[1:400], "mm", "symbol", db="KEGG")
  
  plot_enrich_ORA(ora@result, axis_x="enrichFactor", filter_by="enrichFactor",
                  color_by="-log10(p.adjust)", plot_type="dot", text_len_limit=20)
}

