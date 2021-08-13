#' Run Linear Models
#'
#' This function takes in phyloseq object, association matix, p value matrix
#' and create a combined OTU-OTU and OTU-Phenotype network. User can select model to
#' define OTU-Phontype assocaition.
#' @param n Number of iterations to run
#' @param phenotype A variable representing phenotype- which is used a reponse variable for linear regression
#' @param odata Count data with rows as samples, and columns as features/OTUs/ASVs
#' @return Features from linear regression
#' @export

model.linear <- function(n, phenotype, odata){
  mat <- list()
  for (i in 1:ncol(odata)){
    y <- odata[, i]
    which_otu <- colnames(odata)[i]
    mod <- lm(phenotype ~ y)
    sum_mod <- summary.lm(mod)
    pv <- round(sum_mod$coefficients[8], 4)
    r2 <- round(sum_mod$r.squared,4)
    mat[[i]] <- list(num = i, otu = which_otu, pvalue = pv, relation = r2)}

  dd <- data.frame(matrix(unlist(mat), nrow=length(mat), byrow=T))
  colnames(dd) = names(mat[[1]])
  return(dd)
}

