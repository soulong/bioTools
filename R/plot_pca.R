
#' @name plot_pca
#' @title plot pca
#' @description pca plot for a matrix
#'
#' @param data A data-matrix or data-frame containing numerical data only,
#'        variables are expected to be in the rows and samples in the columns by default
#' @param metadata data-frame containing metadata,
#'        strictly enforced that rownames(metadata) == colnames(mat). DEFAULT: NULL
#' @param color_by color point by, one of column in metadata
#' @param shape_by shape point by, one of column in metadata
#' @param center center the data before performing PCA, same as prcomp() 'center' parameter
#' @param scale scale the data before performing PCA, same as prcomp() 'scale' parameter
#'
#' @import PCAtools
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom ggplotify as.ggplot
#'
#' @return list of ggplot2 object, or list of pca object
#'
#' @export
#'
plot_pca <- function(data,
                     metadata=NULL,
                     color_by=NULL,
                     shape_by=NULL,
                     center=TRUE,
                     scale=FALSE,
                     return_data=FALSE
                     ) {
  # get pca object
  pca <- pca(data, metadata, center, scale, rank=10, removeVar)

  # return data
  if(return_data) return(pca)

  # screeplot
  scree <- screeplot(pca, sizeCumulativeSumLine=1, sizeCumulativeSumPoints=1.5, title=NULL) + theme_linedraw()

  n <- 5

  # eigencor plot
  eigencor <- suppressWarnings(
    eigencorplot(pca, metavars=colnames(metadata),
                           col=c("#008837","#f7f7f7","#7b3294"))
    ) %>% ggplotify::as.ggplot()
  # loading plpt
  loadings <- suppressMessages(plotloadings(pca, components = getComponents(pca, seq_len(n)),
                           shapeSizeRange=6, legendLabSize=6, legendIconSize=0.2,
                           col=c("#01665e","#f5f5f5","#8c510a"))) + theme_linedraw()
  # pairs plot
  if(!is.null(color_by)) {
    if(!is.null(shape_by)) {
      pairs <- suppressMessages(pairsplot(pca, components = getComponents(pca, seq_len(n)), triangle=T, trianglelabSize=12,
                         colby=color_by, shape=shape_by,
                         gridlines.major=F, gridlines.minor=F, axisLabSize=8))
    } else {
      pairs <- suppressMessages(pairsplot(pca, components = getComponents(pca, seq_len(n)), triangle=T, trianglelabSize=12,
                         colby=color_by,
                         gridlines.major=F, gridlines.minor=F, axisLabSize=8))
    }
  } else {
    if(!is.null(shape_by)) {
      pairs <- suppressMessages(pairsplot(pca, components = getComponents(pca, seq_len(n)), triangle=T, trianglelabSize=12,
                         shape=shape_by,
                         gridlines.major=F, gridlines.minor=F, axisLabSize=8))
    } else {
      pairs <- suppressMessages(pairsplot(pca, components = getComponents(pca, seq_len(n)), triangle=T, trianglelabSize=12,
                         gridlines.major=F, gridlines.minor=F, axisLabSize=8))
    }
  }

  # biplot
  if(!is.null(color_by)) {
    if(!is.null(shape_by)) {
      biplot <- biplot(pca, colby=color_by, shape=shape_by) + theme_linedraw()
    } else {
      biplot <- biplot(pca, colby=color_by) + theme_linedraw()
    }
  } else {
    if(!is.null(shape_by)) {
      biplot <- biplot(pca, shape=shape_by) + theme_linedraw()
    } else {
      biplot <- biplot(pca) + theme_linedraw()
    }
  }

  # return plot
  plot <- list(scree=scree, loadings=loadings, eigencor=eigencor, pairs=pairs, biplot=biplot)
  return(plot)
}



# not run
if(F) {
  data <- readxl::read_xlsx("/Users/hh/Desktop/2020-01-13_smartSeq_lung_EC_PC_HH/gene_expression.xlsx") %>%
    .[, c(2, 6:31)] %>%
    filter(!is.na(external_gene_name)) %>%
    distinct(external_gene_name, .keep_all = T) %>%
    column_to_rownames("external_gene_name")

  metadata <- data.frame(treat=c(rep("pbs", 7), rep("lps", 7), rep("pbs", 6), rep("lps", 6)),
                         cell=c(rep("ec", 14), rep("pc", 12)))
  rownames(metadata) <- colnames(data)

  pca <- plot_pca(data, metadata, color_by = "treat", shape_by = "cell")
  pca$scree
  pca$loadings
  pca$eigencor
  pca$pairs
  pca$biplot

  plot_pca(data, metadata, color_by = "treat", shape_by = "cell",return_data = T)
}

