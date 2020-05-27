
#' @name plot_enrich_gsea
#' @title plot_enrich_gsea
#' @description dotplot for enrich result from enrich_GSEA
#'
#' @param enrichment data.frame, from GSEA enrichment,
#'        result can be merged from multi enrichment, grouping info can be viewed by facet
#' @param show_category integer, show enriched category number
#' @param plot_type character, one of dot, bar
#' @param filter_by character, filter show_category by which column,
#'        one of -log10(p.adjust), NES, Count
#' @param color_by character, map point color, one of -log10(p.adjust), NES, Count
#' @param color_range two element character vector, mapping color range, from low to high,
#'        for gsea plot, the value of NSE = 0 will be white color
#' @param size_by character, map point size, currently, only Count avaible for gseaResult
#' @param size_range two element numeric vector, mapping point size, from low to high
#'        only works if plot_type is dot
#' @param text_len_limit integer, wrap y-axis text length
#' @param facet_by character, facet plot by which column in enrichment,
#'        default avaible value is direction, user use add aditional column in enrichment
#' @param facet_scales character, facet scales, one of free_y, free_x, "free, fixed,
#'        only works if facet_by is not NULL
#'
#' @import ggplot2
#' @importFrom dplyr arrange mutate slice vars group_by
#' @importFrom stringr str_wrap str_split
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
                             color_range=c("blue2", "red2"),
                             size_by="Count", # not changable by now
                             size_range=c(2, 6),
                             text_len_limit=40,
                             facet_by=NULL, # direction
                             facet_scales="free_y"
) {

  if(!(class(enrichment) == "data.frame"))
    stop("enrichment must be a data.frame (result slot from clusterProfiler)")

  if(!(plot_type %in% c("dot", "bar")))
    stop("plot_type must be one of dor, bar")

  if(!(filter_by %in% c("-log10(p.adjust)", "NES", "Count")))
    stop("filter_by must be one of -log10(p.adjust), NES, Count")

  if(!(color_by %in% c("-log10(p.adjust)", "NES", "Count")))
    stop("color_by must be one of -log10(p.adjust), NES, Count")

  if(!is.null(facet_by)) {
    if(!(facet_scales %in% c("free", "free_x", "free_y", "fixed")))
      stop("facet_scales must be one of free, free_x, free_y, fixed")
  }

  # count core enrichment genes
  enrichment$Count <-sapply(enrichment$core_enrichment, function(x) length(str_split(x, "/", simplify = T)))

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
    data <- mutate(enrichment, `-log10(p.adjust)`=-log10(p.adjust),
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
      geom_point(aes(size=Count, fill=!!sym(color_by)), shape=21, color="transparent") +
      scale_size_continuous(range=size_range) +
      scale_fill_gradient2(low=color_range[1], mid="white", high=color_range[2], midpoint=0) +
      scale_y_discrete(labels=function(x) str_wrap(x, width=text_len_limit)) +
      theme_classic() +
      ylab("")
  } else {
    plot <- ggplot(data, aes(NES, Description)) +
      geom_bar(aes(fill=!!sym(color_by)), stat="identity", color=NA, width=0.75) +
      scale_fill_gradient2(low=color_range[1], mid="white", high=color_range[2], midpoint=0) +
      scale_y_discrete(labels=function(x) str_wrap(x, width=text_len_limit)) +
      theme_classic() +
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

  gsea <- enrich_gsea(genelist, "mm", keytype = "entrezid", db="KEGG")
  gseaplot(gsea, gsea@result$ID[1])
  gseaplot2(gsea, gsea@result$ID[1:2])

  plot_enrich_gsea(gsea@result, filter_by="-log10(p.adjust)", color_by = "NES", text_len_limit=30, facet_by = "direction")
  plot_enrich_gsea(gsea@result, filter_by="NES", color_by = "-log10(p.adjust)", text_len_limit=30, plot_type = "dot")
}

