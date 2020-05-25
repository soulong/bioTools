
#' @name plot_enrich_gsea
#' @title plot_enrich_gsea
#' @description dotplot for enrich result from enrich_GSEA
#'
#' @param enrichment data.frame, from GSEA enrichment,
#'        result can be merged from multi enrichment, grouping info can be viewed by facet
#' @param show_category integer, show enriched category number
#' @param plot_type character, one of dot, bar
#' @param filter_by character, filter show_category by which column,
#'        one of -log10(p.adjust), NES
#' @param color_by character, map point color, one of -log10(p.adjust), enrichFactor, Count
#' @param color_range two element character vector, mapping color range, from low to high
#' @param point_size point size, only works if plot_type is dot
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
plot_enrich_gsea <- function(enrichment,
                             show_category=10,
                             plot_type="dot",
                             filter_by="-log10(p.adjust)",
                             color_by="NES",
                             color_range=c("darkblue", "red"),
                             point_size=5,
                             text_len_limit=50,
                             facet_by=NULL,
                             facet_scales="free_y"
) {
  
  if(!is.data.frame(enrichment))
    stop("enrichment must be a tibble/data.frame")
  
  if(!(plot_type %in% c("dot", "bar")))
    stop("plot_type must be one of dor, bar")
  
  if(!(filter_by %in% c("-log10(p.adjust)", "NES")))
    stop("filter_by must be one of -log10(p.adjust), NES")
  
  if(!(color_by %in% c("-log10(p.adjust)", "NES")))
    stop("color_by must be one of -log10(p.adjust), NES")
  
  if(!is.null(facet_by)) {
    if(!(facet_scales %in% c("free", "free_x", "free_y", "fixed")))
      stop("facet_scales must be one of free, free_x, free_y, fixed")
  }
  
  # group facet
  if(!is.null(facet_by)) {
    data <- mutate(enrichment,
                   `-log10(p.adjust)`=-log10(p.adjust),
                   direction=if_else(NES>=0, "Up-regulated", "Down-regulated")) %>%
      group_by(direction, !!sym(facet_by)) %>%
      arrange(desc(abs(!!sym(filter_by)))) %>%
      dplyr::slice(seq_len(show_category)) %>%
      arrange(NES) %>%
      mutate(Description=factor(Description, levels=unique(.$Description)))
  } else {
    data <- mutate(enrichment,
                   `-log10(p.adjust)`=-log10(p.adjust),
                   direction=if_else(NES>=0, "Up-regulated", "Down-regulated")) %>%
      group_by(direction) %>%
      arrange(desc(abs(!!sym(filter_by)))) %>%
      dplyr::slice(seq_len(show_category)) %>%
      arrange(NES) %>%
      mutate(Description=factor(Description, levels=unique(.$Description)))
  }
  
  if(plot_type=="dot") {
    plot <- ggplot(data, aes(NES, Description)) +
      geom_segment(aes(yend=Description), xend=0, size=0.7, color="grey50", alpha=0.7) +
      geom_point(aes(size=!!sym(size_by), color=!!sym(color_by)), size=point_size) +
      scale_color_gradient2(low=color_range[1], mid="white", high=color_range[2], midpoint=0) +
      scale_y_discrete(labels=function(x) str_wrap(x, width=text_len_limit)) +
      theme_cowplot() +
      ylab("")
  } else {
    plot <- ggplot(data, aes(NES, Description)) +
      geom_bar(aes(fill=!!sym(color_by)), stat="identity", color=NA, width=0.75) +
      scale_fill_gradient2(low=color_range[1], mid="white", high=color_range[2], midpoint=0) +
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
  
  genelist <- res$stat
  names(genelist) <- res$entrezgene_id
  
  gsea <- enrich_GSEA(genelist, "mm", db="KEGG")
  gsea@result %>% as_tibble()
  gseaplot(gsea, gsea@result$ID[1])
  gseaplot2(gsea, gsea@result$ID[1:2])
  
  plot_enrich_GSEA(gsea$GO, filter_by="-log10(p.adjust)", facet_by="ontology", text_len_limit=30)
  plot_enrich_GSEA(gsea$KEGG, filter_by="NES", text_len_limit=30)
}

