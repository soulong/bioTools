
#' @title  make TF to target dataset from online data
#'
#' @description
#'
#' @param category must be c("trrust", "dorothea"), for now
#' @param dorothea_confidence dorothea filter level, "A" means validated at assay level
#'
#' @importFrom readr read_tsv
#' @importFrom magrittr set_colnames %>%
#' @importFrom rlist list.append
#' @importFrom dplyr select filter mutate if_else distinct bind_rows arrange
#'
#' @return a list containing two tibble: tf to target (hs and mm)
#'
#' @examples
#'
#'
make_TFdb <- function(category=c("trrust", "dorothea"),
                      dorothea_confidence=c("A") # A-E
                      ) {

  ## trrust
  if("trrust" %in% category) {
    # address
    species <- c(
      hs = "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
      mm = "https://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv"
    )

    # download TRRUST data
    trrust <- list()
    for(i in seq_along(species)) {

      print(paste0("process on: ", names(species)[i]))
      print("download TRRUST data from website")
      tmp_file <- tempfile()
      download.file(species[i], destfile = tmp_file)

      # read data
      trrust_data <- suppressMessages(read_tsv(tmp_file, col_names = FALSE, progress = FALSE)) %>%
        set_colnames(c("tf", "target", "regulation", "reference")) %>%
        mutate(reference=NULL, source="trrust") %>%
        distinct(tf, target, .keep_all = T)

      # merge species data to trrust list
      trrust <- list.append(trrust, trrust_data)
    }

    # set species names
    names(trrust) <- names(species)

  }


  ## dorothea
  if("dorothea" %in% category) {
    # install dorothea package
    if(!("dorothea" %in% installed.packages()[, 1, drop = TRUE])) {
      print("dorothea package not found, install it though Bioconductor")
      BiocManager::install("dorothea")

      # filter dorothea
      dorothea_hs <- filter(dorothea::dorothea_hs, confidence %in% dorothea_confidence) %>%
        select(tf, target, regulation=mor) %>%
        mutate(source = "dorothea") %>%
        distinct(tf, target, .keep_all = T) %>%
        mutate(regulation = if_else(regulation == 1, "Activation",
                                    if_else(regulation == -1, "Repression", "Unknown")))
      dorothea_mm <- filter(dorothea::dorothea_mm, confidence %in% dorothea_confidence) %>%
        select(tf, target, regulation=mor) %>%
        mutate(source = "dorothea") %>%
        distinct(tf, target, .keep_all = T) %>%
        mutate(regulation = if_else(regulation == 1, "Activation",
                                    if_else(regulation == -1, "Repression", "Unknown")))
      # merge to list
      dorothea <- list(hs=dorothea_hs, mm=dorothea_mm)
    }
  }


  ## merge datasets
  tf_hs <- bind_rows(trrust$hs, dorothea$hs) %>%
    arrange(regulation, desc(source)) %>% # optimal for Activation, Repression, then Unknown, ttrust dataset was also optimal for selection
    distinct(tf, target, .keep_all = T)
  tf_mm <- bind_rows(trrust$mm, dorothea$mm) %>%
    arrange(regulation, desc(source)) %>% # optimal for Activation, Repression, then Unknown, ttrust dataset was also optimal for selection
    distinct(tf, target, .keep_all = T)
  tf2target = list(hs = tf_hs, mm = tf_mm)

  return(tf2target)
}



# not run
if(F) {
  library(usethis)
  tf2target <- make_TFdb()
  use_data(tf2target, overwrite = TRUE, version = 3)
}

