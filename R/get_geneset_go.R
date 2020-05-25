
#' @name get_geneset_go
#' @title get GO annotation
#' @description get GO annotation from bioconductor org.Xx.eg.db
#'
#' @param species character, must be one of hs, mm
#' @param ontology character, GO ontology to keep, one or more of "BP","MF","CC"
#' @param min_size integer, minimal set size
#' @param max_size integer, maxium set size
#' @param evidence_exclude character, excluded GO evidence,
#'        refer to http://geneontology.org/docs/guide-go-evidence-codes/
#'
#' @importFrom AnnotationDbi keys metadata
#' @importFrom dplyr distinct filter transmute %>%
#' @importFrom tibble as_tibble
#'
#' @return tibble
#'
#' @export
#'
get_geneset_go <- function(species="hs",
                           ontology=c("BP","MF","CC"),
                           min_size=5,
                           max_size=500,
                           evidence_exclude=c("NAS","NR","ND","IEA")
) {

  if(!species %in% c("hs", "mm")) stop("species must be one of hs, mm!")
  if(!all(ontology %in% c("BP","MF","CC"))) stop("ontology must be in BP, MF, CC!")

  if(species=="hs") {
    require(org.Hs.eg.db, quietly=T)
    species.go <- org.Hs.eg.db
  } else {
    require(org.Mm.eg.db, quietly=T)
    species.go <- org.Mm.eg.db
  }

  t2g <- suppressMessages(
    AnnotationDbi::select(species.go, keys(species.go),
                          c("GO", "ONTOLOGY", "SYMBOL", "EVIDENCE"), "ENTREZID")) %>%
    as_tibble() %>%
    filter(!(EVIDENCE %in% evidence_exclude), ONTOLOGY %in% ontology)

  require(GO.db, quietly=T)
  t2n <- suppressMessages(
    AnnotationDbi::select(GO.db, unique(t2g$GO),
                          c("GOID", "TERM"), "GOID")) %>%
    as_tibble()

  df <- left_join(t2g, t2n, by=c("GO"="GOID")) %>%
    transmute(goid=GO, entrezid=as.integer(ENTREZID), symbol=SYMBOL, term=TERM, ontology=ONTOLOGY) %>%
    distinct()

  # filter geneset size
  keep_go <- group_by(df, goid) %>%
    summarise(size=n()) %>%
    mutate(keep=if_else(size >= min_size & size <= max_size, "keep", "drop")) %>%
    filter(keep=="keep")
  df <- filter(df, goid %in% keep_go$goid)

  print(distinct(df, goid, .keep_all=T) %>% group_by(ontology) %>% summarise(set_size=n()))

  return(go=df)
}




# not run by user
if(F) {
  go_hs <- get_geneset_go("hs")
  go_mm <- get_geneset_go("mm")
  geneset_go <- list(hs=go_hs, mm=go_mm)
  # save to data/
  usethis::use_data(geneset_go, overwrite=T, version=3)
}




