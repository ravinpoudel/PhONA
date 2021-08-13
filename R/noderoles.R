#' Run Node Roles using rnetcarto
#'
#' This function takes in graph_object from igraph, and run rnetcarto function to assign roles based SA algorithm.
#' @param graph_object A graph object from igraph.
#' @return A plot with nodes classified into four groups based on within-module degree and among-module connectivity
#' @export


node.role <- function(graph_object){
  net_adj = get.adjacency(graph_object,sparse=FALSE) # here we don't consider weights- just take presence and absence, as negative wights get infinity value in the rnetcarto SA.
  return(netcarto(net_adj))
}
