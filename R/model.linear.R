#' Run Linear Models
#'
#' This function takes in phyloseq object, association matix, p value matrix
#' and create a combined OTU-OTU and OTU-Phenotype network. User can select model to
#' define OTU-Phontype assocaition.
#' @param x A phyloseq object which combined OTU count, taxonomy and metadata
#' @return A matrix of the infile
#' @export

model.linear <- function(n, x, odata){
  mat <- list()
  for (i in 1:ncol(odata)){
    y <- odata[, i]
    which_otu <- colnames(odata)[i]
    mod <- lm(x ~ y)
    sum_mod <- summary.lm(mod)
    pv <- round(sum_mod$coefficients[8], 4)
    r2 <- round(sum_mod$r.squared,4)
    mat[[i]] <- list(num = i, otu = which_otu, pvalue = pv, relation = r2)}

  dd <- data.frame(matrix(unlist(mat), nrow=length(mat), byrow=T))
  colnames(dd) = names(mat[[1]])
  return(dd)
}

