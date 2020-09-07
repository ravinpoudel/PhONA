#' Run Lasso Models
#'
#' This function takes in phyloseq object, association matix, p value matrix
#' and create a combined OTU-OTU and OTU-Phenotype network. User can select model to
#' define OTU-Phontype assocaition.
#' @param x A phyloseq object which combined OTU count, taxonomy and metadata
#' @return A matrix of the infile
#' @export
read.mothur.shared <- function(shared.file){
  otu_data <-read.table(shared.file, header=T, row.names = 2, stringsAsFactors = FALSE) %>%
    .[,-c(1,2)]
  return(otu_data)
}
