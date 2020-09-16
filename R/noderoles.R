#' Run Node Roles usign rnetcarto
#'
#' This function takes in phyloseq object, association matix, p value matrix
#' and create a combined OTU-OTU and OTU-Phenotype network. User can select model to
#' define OTU-Phontype assocaition.
#' @param x A phyloseq object which combined OTU count, taxonomy and metadata
#' @return A matrix of the infile
#' @export


node.role <- function(graph_object){
  net_adj = get.adjacency(graph_object,sparse=FALSE) # here we don't consider weights- just take presence and absence, as negative wights get infinity value in the rnetcarto SA.
  return(netcarto(net_adj))
}
