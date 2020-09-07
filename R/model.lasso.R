#' Run Lasso Models
#'
#' This function takes in phyloseq object, association matix, p value matrix
#' and create a combined OTU-OTU and OTU-Phenotype network. User can select model to
#' define OTU-Phontype assocaition.
#' @param x A phyloseq object which combined OTU count, taxonomy and metadata
#' @return A matrix of the infile
#' @export

model.lasso <- function(n, x, odata){
  df_list <- list()
  for (i in 1:n){
    # message ("Running Iteration Number::", i)
    x1 =x
    odata1 = odata
    tryCatch({df_list[[i]] = lasso(x1, odata1)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  # flat list to df
  df <- do.call(rbind,df_list)
  # # filter combined df based on unique otus
  df_uniq = df %>% distinct(otus, .keep_all = TRUE)
  df_uniq
}


