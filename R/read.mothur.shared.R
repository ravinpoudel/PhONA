#' Reads .shared file from mothur
#'
#' Reads in .shared file, a count matrix from mothur.
#' @param shared.file A .shared file
#' @return A table with rows as features(OTUs) and columns as samples
#' @export
read.mothur.shared <- function(shared.file){
  otu_data <-read.table(shared.file, header=T, row.names = 2, stringsAsFactors = FALSE) %>%
    .[,-c(1,2)]
  return(otu_data)
}
