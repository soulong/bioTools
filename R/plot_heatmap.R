
#' @name plot_heatmap
#' @title plot heatmap
#' @description plot heatmap using complexHeatmap
#'
#' @param matrix matrix, with unique rownames and colnames
#' @param group vecotr, used to anno column or cluster_withon_group, the order is same with matrix colname
#' @param scale_method character, scale method, one of "scale", "log1p", "none"
#' @param col_cluster_type character, cluster column method, one of "auto", "semi", "none"
#' @param split_row_num integer, divide row cluster into how many sub-clusters
#' @param color_map numeric vector, heatmap color mapping cut, length must be same with color_palette
#' @param color_palette character vector, color mapped to color_map, length must be same with color_map
#' @param row_names_size integer, rownames size
#' @param column_names_size nteger, colnames size
#' @param dist_method dist method, refer to dist
#' @param hclust_method hclust method, refer to hclust
#' @param ... additional parameters pass to ComplexHeatmap::Heatmap
#'
#' @importFrom ComplexHeatmap Heatmap cluster_within_group
#' @importFrom dendextend color_branches
#' @importFrom circlize colorRamp2
#' @importFrom tibble rownames_to_column as_tibble
#'
#' @return a list, contain heatmap and scaled ordered heatmap data
#'
#' @export
#'
plot_heatmap <- function(
  matrix,
  group=NULL,
  scale_method="scale",
  col_cluster_type="auto",
  split_row_num=1,
  color_map=c(-2, 0, 2),
  color_palette=c("#377eb8", "white", "#e41a1c"),
  row_names_size=2,
  column_names_size=8,
  dist_method="euclidean",
  hclust_method="complete",
  ...
) {

  # process matrix
  if(scale_method=="scale") {
    # remove much zero rows
    mat <- apply(matrix, 1, scale) %>%
      t() %>% na.omit()
    colnames(mat) <- colnames(matrix)

    na_rows <- attributes(mat)$na.action %>% names()
    if(length(na_rows)>=1) {
      print(paste0("some rows will be discarded, due to scaling:"))
      print(na_rows)
    }
  }

  if(scale_method=="log1p") {
    mat <- log2(matrix + 1)
  }

  if(scale_method=="none") {
    mat <- matrix
  }

  # column dendegram
  if(col_cluster_type=="semi") {
    if(is.null(group))
      stop("please provide valid group value!")
    set.seed(123)
    dend_col <- cluster_within_group(mat, factor=group)

  } else {
    if(col_cluster_type=="auto") {
      dend_col <- TRUE
    } else {
      dend_col <- FALSE
    }
  }

  # row dendegram
  set.seed(123)
  dend_row <- mat %>%
    dist(method=dist_method) %>%
    hclust(method=hclust_method) %>%
    as.dendrogram()

  # heatmap
  if(split_row_num==1) {
    ht <- Heatmap(mat, name="Scale",
                  col=colorRamp2(color_map, color_palette),
                  row_names_gp=gpar(fontsize=row_names_size),
                  column_names_gp=gpar(fontsize=column_names_size),
                  heatmap_legend_param=list(title_gp=gpar(fontsize=6, fontface="plain")),
                  column_names_rot=45,
                  cluster_columns=dend_col,
                  cluster_rows=dend_row,
                  row_dend_reorder=T,
                  ...)
  } else {
    # color cluster dendgrams
    dend_row <- color_branches(dend_row, k=split_row_num)
    ht <- Heatmap(mat, name="Scale",
                  col=colorRamp2(color_map, color_palette),
                  row_names_gp=gpar(fontsize=row_names_size),
                  column_names_gp=gpar(fontsize=column_names_size),
                  heatmap_legend_param=list(title_gp=gpar(fontsize=6, fontface="plain")),
                  column_names_rot=45,
                  cluster_columns=dend_col,
                  cluster_rows=dend_row,
                  row_split=split_row_num,
                  row_gap=unit(2, "mm"),
                  row_dend_reorder=T,
                  ...)
  }

  # export ordered matrix
  if(split_row_num==1) {
    ordered_row_export <- row_order(ht) %>% list()
  } else {
    ordered_row_export <- row_order(ht)
  }

  names(ordered_row_export) <- 1:split_row_num

  ordered_matlist <- lapply(ordered_row_export,
                            function(x) as_tibble(rownames_to_column(as.data.frame(mat[x, ]), "name")))
  # merge list to a df
  ordered_mat <- bind_rows(ordered_matlist, .id = "cluster")

  return(list(heatmap=ht, ordered_mat=ordered_mat))
}




# not run
if(F) {
  expr <- readxl::read_xlsx("/Users/hh/Desktop/2020-01-13_smartSeq_lung_EC_PC_HH/2020-01-15_gene_expression_trim.xlsx")
  vars <- matrixStats::rowVars(as.matrix(expr[, 6:17]))
  matrix <- expr[order(vars, decreasing = TRUE)[1:1000], c(2, 6:17)] %>%
    column_to_rownames("external_gene_name") %>%
    as.matrix()
  matrix <- matrix[rowMeans(matrix) > 1, ]
  dim(matrix)

  group <- c(rep("pbs", 6), rep("lps", 6))

  plot <- plot_heatmap(matrix, group, col_cluster_type = "semi", scale_method = "scale", split_row_num = 2)
  plot$heatmap
}
