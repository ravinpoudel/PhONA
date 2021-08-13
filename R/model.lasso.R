#' Run Lasso Models
#'
#' This function takes in phyloseq object, association matix, p value matrix
#' and create a combined OTU-OTU and OTU-Phenotype network. User can select model to
#' define OTU-Phontype assocaition.
#' @param n Number of iterations to run
#' @param phenotype A variable representing phenotype- which is used a reponse variable for lasso regression
#' @param odata Count data with rows as samples, and columns as features/OTUs/ASVs
#' @return Features from lasso regression
#' @export

model.lasso <- function(n, phenotype, odata){
  df_list <- list()
  for (i in 1:n){
    # message ("Running Iteration Number::", i)
    phenotype = phenotype
    odata1 = odata
    tryCatch({df_list[[i]] = lasso(phenotype, odata1)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  # flat list to df
  df <- do.call(rbind,df_list)
  # # filter combined df based on unique otus
  df_uniq = df %>% distinct(otus, .keep_all = TRUE)
  df_uniq
}


