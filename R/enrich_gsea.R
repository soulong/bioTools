
#' @name enrich_gsea
#' @title enrich_gsea
#' @description GSEA enrichment using annoHub db
#'
#' @param genelist character vector, must be gene symbols
#' @param keytype character, indicating genes type, one of entrezid, symbol
#' @param species character, must be one of hs, mm
#' @param db character, which internal db to enrich, one of GO, KEGG
#' @param go_ontology which GO ontology to use, one or more of BP, MF, CC,
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
                        keytype="symbol", # symbol, entrezid
                        db="GO",
                        go_ontology="BP", # "BP", "MF", "CC"
                        pvalue_cutoff=0.1,
                        padjust_method="BH",
                        min_size=10,
                        max_size=400,
                        ...
                        ) {

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
    t2g <- filter(bioTools::geneset_go[[species]], ontology %in% go_ontology) %>%
      .[, c("goid", keytype)]
    t2n <- filter(bioTools::geneset_go[[species]], ontology %in% go_ontology) %>%
      .[, c("goid", "term")]
    print(t2n)
    enrich <- GSEA(genelist, seed=123,
                   TERM2GENE=t2g, TERM2NAME=t2n,
                   pvalueCutoff=pvalue_cutoff, pAdjustMethod=padjust_method,
                   minGSSize=min_size, maxGSSize=max_size, ...)
    # # add ontology, this will make enrichplot::gseaplot2 not work
    # ontology <- bioTools::geneset_go[[species]] %>%
    #   dplyr::select(goid, ontology) %>%
    #   distinct(goid, .keep_all = TRUE)
    # enrich@result <- left_join(enrich@result, ontology, by=c("ID"="goid"))
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



# not run
if(F) {
  # data
  dge <- readxl::read_xlsx("~/Desktop/2020-01-13_smartSeq_lung_EC_PC_HH/statistics/2020-01-15_wald_TEST_lps.pc_pbs.pc.xlsx", 2)
  dge <- filter(data, !is.na(external_gene_name), !is.na(entrezgene_id)) %>%
    arrange(desc(stat))
  genelist <- dge$stat
  names(genelist) <- dge$external_gene_name

  # enrich
  go <- enrich_gsea(genelist, "mm", keytype = "symbol", db = "GO", go_ontology = "CC")
  dotplot(go)
  enrichplot::gseaplot(go, go@result$ID[1], by="runningScore")
  enrichplot::gseaplot2(go, go@result$ID[1])
  enrichplot::gseaplot2(go, go@result$ID[1:2])
  gsea <- enrich_gsea(genelist, "mm", keytype = "symbol", db = "KEGG")
  enrichplot::gseaplot2(gsea, gsea@result$ID[1], title = gsea@result$Description[1])
}

