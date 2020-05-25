
#' @name enrich_gsea
#' @title enrich_gsea
#' @description GSEA enrichment using annoHub db
#'
#' @param genelist character vector, must be gene symbols
#' @param keytype character, indicating genes type, one of entrezid, symbol
#' @param species character, must be one of hs, mm
#' @param db character, which internal db to enrich, one of GO, KEGG
#' @param ontology which GO ontology to use, one or more of BP, MF, CC,
#'        only works if db contains GO
#' @param pvalue_cutoff pvalueCutoff
#' @param padjust_method pvalue_cutoff method, default is "BH"
#' @param min_size min set size
#' @param max_size max set size
#' @param ... additional parameters pass to clusterProfiler::enricher
#'
#' @importFrom dplyr filter left_join distinct
#' @importFrom clusterProfiler GSEA
#'
#' @return enrichResult from enricher
#'
#' @export
#'
enrich_gsea <- function(genelist,
                        species="hs",
                        keytype="entrezid", # symbol
                        db="GO",
                        ontology=c("BP", "MF", "CC"),
                        pvalue_cutoff=0.1,
                        padjust_method="BH",
                        min_size=10,
                        max_size=400,
                        ...) {
  
  names <- names(genelist)
  if(is.null(names)) {
    stop("genelist must have names")
  } else {
    if(!length(unique(names))==length(names)) {
      stop("genelist names are not unique")
    }
  }
  
  if(!(species %in% c("hs", "mm")))
    stop("species must be one of hs, mm")
  
  if(!(keytype %in% c("entrezid", "symbol")))
    stop("keytype must be one of entrezid, symbol")
  
  if(!all(db %in% c("GO", "KEGG")))
    stop("db must be one of GO, KEGG")
  
  # rank genelist
  genelist <- genelist[order(genelist, decreasing=TRUE)]
  
  print(paste0("enrich on: ", db))
  
  if(db=="GO") {
    t2g <- filter(bioTools::geneset_go[[species]], ontology %in% ontology) %>%
      .[, c("goid", keytype)]
    t2n <- filter(bioTools::geneset_go[[species]], ontology %in% ontology) %>%
      .[, c("goid", "term")]
    
    enrich <- GSEA(genelist, seed=123,
                   TERM2GENE=t2g, TERM2NAME=t2n,
                   pvalueCutoff=pvalue_cutoff, pAdjustMethod=padjust_method,
                   minGSSize=min_size, maxGSSize=max_size, ...)
  }
  
  if(db=="KEGG") {
    t2g <- bioTools::geneset_kegg[[species]][, c("keggid", keytype)]
    t2n <- bioTools::geneset_kegg[[species]][, c("keggid", "term")]
    enrich <- GSEA(genelist, seed=123,
                   TERM2GENE=t2g, TERM2NAME=t2n,
                   pvalueCutoff=pvalue_cutoff, pAdjustMethod=padjust_method,
                   minGSSize=min_size, maxGSSize=max_size, ...)
  }
  
  # add metadata
  enrich@organism <- species
  enrich@setType <- db
  enrich@keytype <- keytype
  
  # convert to data.frame for gseaplot from enrichplot
  enrich@result <- as.data.frame(enrich@result)
  
  return(enrich)
}
