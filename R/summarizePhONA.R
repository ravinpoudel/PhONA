#' Summarize graph object
#'
#' This function takes in phyloseq object, association matix, p value matrix
#' and create a combined OTU-OTU and OTU-Phenotype network. User can select model to
#' define OTU-Phontype assocaition.
#' @param x A phyloseq object which combined OTU count, taxonomy and metadata
#' @return A matrix of the infile
#' @export

summarizePhONA <- function (graph.object, roledf){
  top3Hub <- function(graph.object){
    aa <- authority.score(graph.object)$vector
    bb <- V(graph.object)
    dd <- list(bb,round(aa,3))
    df <-  data.frame(matrix(unlist(dd), nrow=length(V(graph.object)), byrow=F))
    ordf <- df[with(df, order(-df[,2])), ]
    Hub <- ordf[1:3,1]
    grabHub <- paste(Hub,collapse=";")
    HUbV <- ordf[1:3,2]
    grabHubV <- paste(HUbV,collapse=";")
    grab <- c(grabHub,grabHubV)
    return(grab)
  }
  role2 <- roledf
  nModules.SA <- length(unique(role2$module))
  edge<-ecount(graph.object)
  node<-vcount(graph.object)
  nodeDegree<-round((edge/node),3)
  connectance<-round((edge/(node^2)),3)
  avgpath <- round(average.path.length(graph.object),3)
  trans <- round(transitivity(graph.object, type = c("average")),3)
  wtc<-max((walktrap.community(graph.object))$membership)
  mod <- round(modularity(walktrap.community(graph.object)),3)
  top3hub<-top3Hub(graph.object)[1]
  top3hubv<-top3Hub(graph.object)[2]
  n_postiveL <-as.vector(E(graph.object)$weight_2 > 0 ) %>% na.omit() %>% as.vector() %>% sum()
  n_negativeL <- as.vector(E(graph.object)$weight_2 < 0 ) %>% na.omit() %>% as.vector() %>% sum()
  tabular <- data.frame(node,edge,nodeDegree,avgpath, trans, mod,connectance,wtc,
                        nModules.SA,top3hub,top3hubv, n_postiveL, n_negativeL)
  return(tabular)
}
