#' Function to assign color to the phyloseq object
#'
#' This function takes in phyloseq object, association matix, p value matrix
#' and create a combined OTU-OTU and OTU-Phenotype network. User can select model to
#' define OTU-Phontype assocaition.
#' @param x A phyloseq object which combined OTU count, taxonomy and metadata
#' @return A matrix of the infile
#' @export


taxacolor <- function(phyobj, coloredby="x"){
  
  possible_tax <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  
  if (coloredby %in% possible_tax) {
    # create a color palltet using coloredby
    # number of unique color needed
    
    tc_tdata = data.frame(as(tax_table(phyobj), "matrix"))
    
    mod.vec <- tc_tdata %>%
      pull(coloredby) %>%
      as.character()
    n <- length(unique(mod.vec))
    mod.num <- as.numeric(as.factor(mod.vec))
    ucolors <- distinctColorPalette(n)
    
    # we need to pass the variable by reference
    # https://stackoverflow.com/questions/29190455/how-to-pass-parameters-by-reference-in-r
    eval.parent(substitute(colmn <-coloredby))
    
    getcolmn <- paste0(colmn, "_color")
    
    tc_tdata[getcolmn] <- sapply(mod.num, function(j) ucolors[j])
    # keep back as phyobj and return

    tax.mat <- tax_table(as.matrix(tc_tdata))
    
    otu.mat = otu_table(phyobj)
    sample.mat <- sample_data(phyobj)
    
    phyobj2 = phyloseq(otu.mat, tax.mat, sample.mat)
    return(phyobj2)
  } else {
    cat(" Provided coloredby: ", coloredby,"\n")
    
    stop("Provide correct coloredby paramter. Available options are: Kingdom,
         Phylum, Class, Order, Family, Genus, Species")
    }
}
