
#' @name enrich_ora
#' @title enrich_ora
#' @description ORA enrichment using internal db
#'
#' @param genes character vector, must be gene symbols
#' @param keytype character, indicating genes type, one of entrezid, symbol
#' @param species character, must be one of hs, mm
#' @param db character, which internal db to enrich, one of GO, KEGG
#' @param go_ontology which GO ontology to use, one or more of BP, MF, CC,
#'        only works if db contains GO
#' @param pvalue_cutoff pvalueCutoff
#' @param qvalue_cutoff qvalueCutoff
#' @param min_size min set size
#' @param max_size max set size
#' @param ... additional parameters pass to clusterProfiler::enricher
#'
#' @importFrom dplyr left_join filter distinct mutate
#' @importFrom clusterProfiler enricher
#'
#' @return enrichResult from enricher
#'
#' @export
#'
enrich_ora <- function(genes,
                       species="hs",
                       keytype="symbol", # symbol, entrezid
                       db="GO",
                       go_ontology="BP", # "BP", "MF", "CC"
                       pvalue_cutoff=0.05,
                       qvalue_cutoff=0.2,
                       min_size=10,
                       max_size=400,
                       ...
                       ) {

  if(!(species %in% c("hs", "mm")))
    stop("species must be one of hs, mm")

  if(!(keytype %in% c("entrezid", "symbol")))
    stop("keytype must be one of entrezid, symbol")

  if(!all(db %in% c("GO", "KEGG")))
    stop("db must be one of GO, KEGG")

  print(paste0("enrich on: ", db))

  if(db=="GO") {
    t2g <- filter(bioTools::geneset_go[[species]], ontology %in% go_ontology) %>%
      .[, c("goid", keytype)]
    t2n <- filter(bioTools::geneset_go[[species]], ontology %in% go_ontology) %>%
      .[, c("goid", "term")]

    enrich <- enricher(genes,
                       TERM2GENE=t2g, TERM2NAME=t2n,
                       pvalueCutoff=pvalue_cutoff, qvalueCutoff=qvalue_cutoff,
                       minGSSize=min_size, maxGSSize=max_size, ...)
    # add ontology
    ontology <- bioTools::geneset_go[[species]] %>%
      dplyr::select(goid, ontology) %>%
      distinct(goid, .keep_all = TRUE)
    enrich@result <- left_join(enrich@result, ontology, by=c("ID"="goid"))
  }

  if(db=="KEGG") {
    t2g <- bioTools::geneset_kegg[[species]][, c("keggid", keytype)]
    t2n <- bioTools::geneset_kegg[[species]][, c("keggid", "term")]
    enrich <- enricher(genes,
                       TERM2GENE=t2g, TERM2NAME=t2n,
                       pvalueCutoff=pvalue_cutoff, qvalueCutoff=qvalue_cutoff,
                       minGSSize=min_size, maxGSSize=max_size, ...)
  }

  # add enrichFactor
  enrich@result <- mutate(enrich@result,
                          enrichFactor=Count / as.numeric(sub("/\\d+", "", BgRatio)))

  # add metadata
  enrich@organism <- species
  enrich@keytype <- keytype

  # convert to data.frame for gseaplot from enrichplot
  enrich@result <- as.data.frame(enrich@result)

  return(enrich)
}



# not run
if(F) {
  # data
  dge <- readxl::read_xlsx("~/Desktop/2020-01-13_smartSeq_lung_EC_PC_HH/statistics/2020-01-15_wald_TEST_lps.pc_pbs.pc.xlsx", 2)
  dge <- filter(data, !is.na(external_gene_name), !is.na(entrezgene_id))
  colnames(dge)[1:5] <- c("ensembl","symbol","entrezid","biotype","description")

  # enrich
  ora <- enrich_ora(dge$symbol, "mm", keytype = "symbol", db = "GO", go_ontology = "BP")
  dotplot(ora)
}


