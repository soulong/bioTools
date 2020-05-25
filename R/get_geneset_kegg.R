
#' @name get_geneset_kegg
#' @title get KEGG annotation
#' @description get KEGG annotation from KEGG API
#'
#' @param species a character, must be one of hs, mm
#'
#' @importFrom clusterProfiler download_KEGG bitr
#' @importFrom dplyr left_join select if_else transmute distinct
#' @importFrom tibble as_tibble
#'
#' @return tibble
#'
#' @export
#'
get_geneset_kegg <- function(species="hs"
) {
  
  if(!species %in% c("hs", "mm")) stop("species must be one of hs, mm!")
  
  print("build KEGG database using KEGG rest API")
  
  species_kegg <- if_else(species=="hs", "hsa", "mmu")
  kegg_list <- download_KEGG(species_kegg)
  
  if(species=="hs") {
    require(org.Hs.eg.db)
    orgdb <- "org.Hs.eg.db"
  } else {
    require(org.Mm.eg.db)
    orgdb <- "org.Mm.eg.db"
  }
  
  kegg2symbol <- suppressMessages(
    bitr(as.integer(kegg_list$KEGGPATHID2EXTID$to), "ENTREZID", "SYMBOL", orgdb))
  
  df <- left_join(kegg_list$KEGGPATHID2EXTID, kegg2symbol, by=c("to"="ENTREZID")) %>%
    left_join(kegg_list$KEGGPATHID2NAME, by=c("from"="from")) %>%
    as_tibble() %>%
    transmute(keggid=from, entrezid=as.integer(to.x), symbol=SYMBOL, term=to.y) %>%
    distinct()
  
  print(group_by(df, keggid) %>% summarise(size=n()))
  
  return(kegg=df)
}



# not run by user
if(F) {
  kegg_hs <- get_geneset_KEGG("hs")
  kegg_mm <- get_geneset_KEGG("mm")
  geneset_kegg <- list(hs=kegg_hs, mm=kegg_mm)
  # save to data/
  usethis::use_data(geneset_kegg, overwrite=T, version=3)
}
