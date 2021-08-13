#' Reads taxonomy file from mothur
#'
#' Reads taxonomy file from mothur and remove all the symbols such as quotation, comma, semicolon etc
#' @param tax.file A taxonomy file obtained from mothur
#' @return A table of formatted taxonomy file
#' @export

read.mothur.taxonomy <- function(tax.file) {
  tbl <- read.delim(tax.file, header = FALSE, row.names = 1) %>%
    .[2:nrow(.),]
  split <- strsplit(as.character(tbl$V3), ";", fixed = TRUE)
  kingdom <- sapply(split, "[", 1) %>% sub("\\([0-9.]+\\)", "", .)
  phylum <- sapply(split, "[", 2) %>% sub("\\([0-9.]+\\)", "", .)
  class <- sapply(split, "[", 3) %>% sub("\\([0-9.]+\\)", "", .)
  order <- sapply(split, "[", 4) %>% sub("\\([0-9.]+\\)", "", .)
  family <- sapply(split, "[", 5) %>% sub("\\([0-9.]+\\)", "", .)
  genus <- sapply(split, "[", 6) %>% sub("\\([0-9.]+\\)", "", .)
  species<- sapply(split, "[", 7) %>% sub("\\([0-9.]+\\)", "", .)
  tax_df <- data.frame(Count= tbl$V2, Kingdom=kingdom, Phylum = phylum, Class = class, Order = order, Family = family,
                       Genus = genus, Species = species) %>% `rownames<-`(rownames(tbl))
  return(tax_df)
}
