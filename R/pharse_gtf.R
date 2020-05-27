
#' @name pharse_gtf
#' @title pharse genes from local gtf file
#' @description prase txdb, gene genes from local gtf file
#'
#' @param genecode_gtf character, file path to genecode GTF file
#' @param genecode_entrezid character, file path to genecode entrezid file
#'
#' @importFrom stringr str_split
#' @importFrom dplyr as_tibble filter distinct left_join rename arrange select
#' @importFrom readr read_table2
#' @importFrom rtracklayer import
#'
#' @return a list, contains tx2gene, genes
#'
#' @export
#'
pharse_gtf <- function(genecode_gtf,
                       genecode_entrezid
                       ) {
  print("read gtf file ...")
  gtf <- rtracklayer::import(genecode_gtf) %>%
    as_tibble() %>%
    filter(type %in% c("gene", "transcript")) %>%
    dplyr::select(type=type, length=width, chromosome=seqnames, strand=strand, source=source,
                  transcript_id, transcript_name, transcript_type, transcript_support_level, tag,
                  gene_id, gene_name, gene_type)
  print("read entrezid ...")
  entrezid <- readr::read_table2(genecode_entrezid, col_names = FALSE, progress = FALSE) %>%
    rename(transcript_id="X1", entrezid="X2")

  # transcripts
  print("get transcripts ...")
  transcripts <- filter(gtf, type=="transcript") %>%
    dplyr::select(transcript_id, gene_id, gene_name, gene_type) %>%
    distinct(transcript_id, .keep_all=TRUE) %>%
    left_join(entrezid) %>%
    arrange(gene_id)

  # remove version info
  print("remove version info ...")
  transcripts$transcript_id <- sapply(transcripts$transcript_id, function(x) str_split(x, "[.]", simplify = T)[1])
  transcripts$gene_id <- sapply(transcripts$gene_id, function(x) str_split(x, "[.]", simplify = T)[1])

  # tx2gene
  print("get tx2gene ...")
  tx2gene <- transcripts %>%
    dplyr::select(transcript_id, gene_id)

  # genes
  print("get genes ...")
  genes <- transcripts %>%
    dplyr::select(ensembl=gene_id, symbol=gene_name, entrezid, gene_type) %>%
    distinct(ensembl, .keep_all=TRUE)

  print("done!")
  return(list(tx2gene=tx2gene, genes=genes))
}



## not run
if(FALSE) {
  # human
  hs <- pharse_gtf(genecode_gtf="/Users/hh/Desktop/hs/annotation.gtf",
                   genecode_entrezid="/Users/hh/Desktop/hs/entrezid.txt")
  # mouse
  mm <- pharse_gtf(genecode_gtf="/Users/hh/Desktop/mm/annotation.gtf",
                   genecode_entrezid="/Users/hh/Desktop/mm/entrezid.txt")

  tx2gene <- list(hs=hs$tx2gene, mm=mm$tx2gene)
  genes <- list(hs=hs$genes, mm=mm$genes)

  # save data
  usethis::use_data(tx2gene, genes, overwrite=TRUE, version=3)
}


