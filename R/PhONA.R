#' Run PhONA
#'
#' This function takes in phyloseq object, association matix, p value matrix
#' and create a combined OTU-OTU and OTU-Phenotype network. User can select model to
#' define OTU-Phontype assocaition.
#'
#'
#' @param physeqobj A phyloseq object which combined OTU count, taxonomy and metadata
#' @param cordata A pairwise square matrix defining OTU-OTU association
#' @param pdata A pairwise square matrix defining p-value for OTU-OTU association
#' @param model Model to define association between OTUs and Phenotype. Option available are "lm", "lasso".
#' In lasso option, we are using lasso model to reduce the number of features/OTUs. OTUs important to phenotype prediction
#'  were ranked using `varImp`. Selected OTUs were then for reduced GLM model.
#' @param OTU_OTU_pvalue Pvalue for OTU-OTU association
#' @param OTU_OTU_rvalue Level of OTU-OTU association
#' @param OTU_Phenotype_pvalue Pvalue for OTU-Phenotype association
#' @param definePhenotype Phenotype to be used. It is a column header from phenotype data
#' @param defineTreatment Select the treatment. It is same as the treatment name that the phyloseq object represents
#' @param coloredby Select taxonomic group to be used for coloring node. Options: Kingdom,Phylum, Class, Order, Family, Genus, Species
#' @param PhenoNodecolor Select color for phenotype node
#' @param PhenoNodesize Select node size for phenotype node
#' @param nodesize Select size for nodes other than phenotype node
#' @param Pheno2OTUedgecolor Select color of edge from OTU to phenotype node
#' @param netlayout Select layout for the network graph. All the layout options from igraph can be used
#' @return A PhONA representing OTU-OTU association as well as OTU-Phenotype association
#' @export
#' @examples
#'\strong{Running linear model based PhONA,where OTU-Phenotype association is defined using linear model}
#'PhONA(
#'physeqobj = phyobj,
#'cordata = sparcc.cor,
#'pdata = sparcc.pval,
#'model = "lm",
#'OTU_OTU_pvalue = 0.001,
#'OTU_OTU_rvalue = 0.6,
#'OTU_Phenotype_pvalue = 0.6,
#'definePhenotype = "Marketable",
#'defineTreatment = "Maxifort",
#'coloredby = "Phylum",
#'PhenoNodecolor = "yellow",
#'PhenoNodesize = 20,
#'PhenoNodelabel = "Yield",
#'nodesize = 10,
#'Pheno2OTUedgecolor = "black",
#'netlayout = layout.fruchterman.reingold)
#'
#'
#'\strong{Running lasso model based PhONA,where OTU-Phenotypeassociation is defined using lasso model and GLM}
#'
#'PhONA(
#'physeqobj = phyobj,
#'cordata = sparcc.cor,
#'pdata = sparcc.pval,
#'model = "lasso",
#'OTU_OTU_pvalue = 0.001,
#'OTU_OTU_rvalue = 0.6,
#'OTU_Phenotype_pvalue = 0.6,
#'definePhenotype = "Marketable",
#'defineTreatment = "Maxifort",
#'coloredby = "Phylum",
#'PhenoNodecolor = "yellow",
#'PhenoNodesize = 20,
#'PhenoNodelabel = "Yield",
#'nodesize = 10,
#'Pheno2OTUedgecolor = "black",
#'netlayout = layout.fruchterman.reingold)

PhONA <- function(physeqobj = physeq,
                  cordata = sparcc.cor,
                  pdata = sparcc.pval,
                  model = c("lm","lasso"),
                  OTU_OTU_pvalue = 0.001,
                  OTU_OTU_rvalue = 0.6,
                  OTU_Phenotype_pvalue = 0.6,
                  definePhenotype="Marketable",
                  defineTreatment="Maxifort",
                  coloredby="Phylum",
                  PhenoNodecolor="yellow",
                  PhenoNodesize=20,
                  PhenoNodelabel="Yield",
                  nodesize=10,
                  Pheno2OTUedgecolor = "black",
                  netlayout=layout.fruchterman.reingold){

  ###
  odata = t(otu_table(physeqobj))
  mdata = sample_data(physeqobj)
 # tdata = taxdata
  tdata = data.frame(tax_table(physeqobj)) # is too slow

  ###
  ####### check the validity of parater, flag if not

  possible_tax <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

  if (coloredby %in% possible_tax) {
    # create a color palltet using coloredby
    # number of unique color needed

    mod.vec <- tdata %>%
      pull(coloredby) %>%
      as.character()
    n <- length(unique(mod.vec))
    mod.num <- as.numeric(as.factor(mod.vec))
    ucolors <- distinctColorPalette(n)
    tdata["color"] <- sapply(mod.num, function(j) ucolors[j])
  } else {
    cat(" Provided coloredby: ", coloredby,"\n")

    stop("Please provide correct coloredby paramter. Available options are: Kingdom,
         Phylum, Class, Order, Family, Genus, Species")

  }

  ######

  p.yes <- pdata < OTU_OTU_pvalue

  p.yes.r <- cordata * p.yes

  p.yes.r <- abs(p.yes.r) > OTU_OTU_rvalue

  p.yes.rr <- p.yes.r * sparcc.cor

  adjm <- as.matrix(p.yes.rr)

  adjm <- as.matrix(adjm)

  ## Select only the taxa for the filtered OTUs by using rownames of otu.pval
  sel.tax <- tdata[rownames(adjm), , drop = FALSE]

  all.equal(rownames(sel.tax), rownames(adjm))


  net.grph = graph.adjacency(adjm, mode = "undirected", weighted = TRUE, diag = FALSE)
  edgew <- E(net.grph)$weight

  ##
  V(net.grph)$Kingdom <- as.character(sel.tax$Kingdom)
  V(net.grph)$Phylum <- as.character(sel.tax$Phylum)
  V(net.grph)$Class <- as.character(sel.tax$Class)
  V(net.grph)$Order <- as.character(sel.tax$Order)
  V(net.grph)$Family <- as.character(sel.tax$Family)
  V(net.grph)$Genus <- as.character(sel.tax$Genus)
  V(net.grph)$Species <- as.character(sel.tax$Species)
  V(net.grph)$color <- tdata$color

##########

  x = mdata %>% pull(definePhenotype)  ## generalize to metadata

##########

  if (model == "lm"){
   #source("R/model.linear.R")
   bb = model.linear(x, odata)
   bb["Treatment"] <- defineTreatment
  }

############
  if (model == "lasso"){
    #source("R/model.linear.R")
    bb = model.lasso(x, odata)
    bb["Treatment"] <- defineTreatment
  }

  #bb$relation <- rnorm(dim(bb)[1], mean = 0.2, sd = 0.2)



#################

  from <- rep(defineTreatment,dim(bb)[1])  ### how to get sel.rootstock

  df <-data.frame(FROM = from, TO=bb$otu, relation=bb$relation, pv =bb$pvalue)
  df.g <- graph.data.frame(d = df, directed = FALSE)
  E(df.g)$weight <- bb$relation
  #plot(df.g)
  E(df.g)$color = Pheno2OTUedgecolor
  E(df.g)$lty <- ifelse(E(df.g)$weight < 0, 2 ,1) # if negative dotted


  ###
  net.grph
  E(net.grph)$color = ifelse(E(net.grph)$weight < 0,"red","blue") # if negative red
  E(net.grph)$lty <- 1

  ##c combined two graphs

  net.two <- graph.union(df.g, net.grph)

  # Edge color
  color.g1 <- E(net.two)$color_1 %>% .[!is.na(.)]
  color.g2 <- E(net.two)$color_2 %>% .[!is.na(.)]
  combined_color <- c(color.g2,color.g1)
  E(net.two)$color <- combined_color


  # Edge link type
  linktype.g1 <- E(net.two)$lty_1 %>% .[!is.na(.)]
  linktype.g2 <- E(net.two)$lty_2 %>% .[!is.na(.)]
  combined_linktype <- c(linktype.g2,linktype.g1)
  E(net.two)$lty<- combined_linktype


  bad.vs3 <- V(net.two)[degree(net.two) == 0 | degree(net.two) == 1 ]

  # # remove isolated nodes
  net.two <- delete.vertices(net.two, bad.vs3)

  ##

  V(net.two)$color[1] <- PhenoNodecolor
  V(net.two)$size <- rep(nodesize,vcount(net.two))
  V(net.two)$size[1] <- PhenoNodesize
  V(net.two)$nName <- V(net.two)$Genus ## can be passed as option
  V(net.two)$nName[1] <- defineTreatment
  V(net.two)$vertex.label.size = rep(0.5,vcount(net.two))
  V(net.two)$vertex.label.size[1] <- 1
  V(net.two)$vertex.label_type = rep(1,vcount(net.two))
  V(net.two)$vertex.label_type[1] <- 2
  V(net.two)$nName[1] <- PhenoNodelabel


  par(mar = c(0, 0, 0, 0))
  plot(net.two,
       vertex.frame.color="black",
       edge.curved=F,
       layout=netlayout,
       vertex.label=V(net.two)$nName,
       vertex.label.color="black",
       #vertex.label.family="Times New Roman",
       vertex.label.font=V(net.two)$vertex.label_type,
       vertex.label.cex=V(net.two)$vertex.label.size)

  # legend_nodes <- unique(V(net.two)$Phylum)
  legend_nodes <- unique(get.vertex.attribute(net.two)[coloredby][[1]])
  legend_nodes[1] <-  PhenoNodelabel
  legend(x=1.1, y=1.1,legend_nodes, pch=21, pt.bg=unique( V(net.two)$color), pt.cex=2, cex=0.8, bty="n", ncol=1)
  #legend("bottom",legend_nodes, pch=21, pt.bg=unique( V(net.two)$color), pt.cex=2, cex=.8, bty="n", ncol=4)

  net.two
}


